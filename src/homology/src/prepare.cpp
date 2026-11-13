#include "homology/homology.hpp"
#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QSet>
#include <QTextStream>
#include <QUrl>
#include <algorithm>
namespace homology {
static void fail(const QString &s) {
    throw Error(s);
}
static QString capture(const QString &s, const QString &re) {
    return QRegularExpression(re).match(s).captured(1);
}
QVector<Protein> flybaseProteins(const QString &fasta) {
    QVector<Protein> out;
    for (auto s : readFasta(fasta, true)) {
        QString gene = capture(s.header, "(?:^|[;,= ])(FBgn[0-9]+)(?:[,; ]|$)"),
                tr = capture(s.header, "(FBtr[0-9]+)");
        if (gene.isEmpty())
            fail("Missing FlyBase gene in " + s.id);
        out.push_back({s.id, gene, tr, "dmel", s.sequence, capture(s.header, "loc=([^;]+)")});
    }
    return out;
}
QVector<Protein> ncbiProteins(const QString &fasta, const QString &gff, const QString &species) {
    QFile f(gff);
    if (!f.open(QIODevice::ReadOnly))
        fail("Cannot read GFF " + gff);
    QTextStream in(&f);
    QMap<QString, QString> geneByProtein, trByProtein;
    while (!in.atEnd()) {
        QString l = in.readLine();
        if (l.startsWith("##FASTA"))
            break;
        if (l.startsWith('#'))
            continue;
        auto c = l.split('\t');
        if (c.size() != 9 || c[2] != "CDS")
            continue;
        QMap<QString, QString> attrs;
        for (auto field : c[8].split(';')) {
            int i = field.indexOf('=');
            if (i > 0)
                attrs[field.left(i)] = QUrl::fromPercentEncoding(field.mid(i + 1).toUtf8());
        }
        QString pid = attrs.value("protein_id"),
                gene = capture(attrs.value("Dbxref"), "(?:^|,)GeneID:([0-9]+)(?:,|$)");
        if (pid.isEmpty() || gene.isEmpty())
            continue;
        gene = "GeneID:" + gene;
        if (geneByProtein.contains(pid) && geneByProtein[pid] != gene)
            fail("Protein maps to multiple genes: " + pid);
        geneByProtein[pid] = gene;
        trByProtein[pid] = attrs.value("Parent");
    }
    QVector<Protein> out;
    for (auto s : readFasta(fasta, true)) {
        if (!geneByProtein.contains(s.id))
            fail("No gene mapping in GFF for " + s.id);
        out.push_back({s.id, geneByProtein[s.id], trByProtein[s.id], species, s.sequence, {}});
    }
    return out;
}
void writeProteins(const QString &dir, const QVector<Protein> &proteins) {
    if (!QDir().mkpath(dir))
        fail("Cannot create " + dir);
    QString fasta, map = "protein_id\tgene_id\n",
                   anno = "protein_id\tgene_id\ttranscript_id\tspecies\tlocation\tlength\n";
    QSet<QString> ids;
    for (auto p : proteins) {
        if (ids.contains(p.proteinId))
            fail("Duplicate protein " + p.proteinId);
        ids.insert(p.proteinId);
        fasta += ">" + p.proteinId + '\n' + p.sequence + '\n';
        map += p.proteinId + '\t' + p.geneId + '\n';
        anno += QStringList{p.proteinId, p.geneId,   p.transcriptId,
                            p.species,   p.location, QString::number(p.sequence.size())}
                    .join('\t') +
                '\n';
    }
    writeText(dir + "/proteins.faa", fasta.toUtf8());
    writeText(dir + "/protein_map.tsv", map.toUtf8());
    writeText(dir + "/annotation.tsv", anno.toUtf8());
}
void extractQueries(const QString &dir, const QStringList &genes, const QString &outputDir) {
    if (genes.isEmpty())
        fail("No requested genes");
    auto seqs = readFasta(dir + "/proteins.faa", true);
    QMap<QString, Sequence> byId;
    for (auto s : seqs)
        byId[s.id] = s;
    QFile f(dir + "/annotation.tsv");
    if (!f.open(QIODevice::ReadOnly))
        fail("Cannot read annotation.tsv");
    QTextStream in(&f);
    in.readLine();
    QMap<QString, Protein> best;
    while (!in.atEnd()) {
        auto c = in.readLine().split('\t');
        if (c.size() != 6)
            fail("Malformed annotation row");
        if (!genes.contains(c[1]))
            continue;
        if (!byId.contains(c[0]))
            fail("Missing annotation protein " + c[0]);
        Protein p{c[0], c[1], c[2], c[3], byId[c[0]].sequence, c[4]};
        if (!best.contains(p.geneId) || p.sequence.size() > best[p.geneId].sequence.size() ||
            (p.sequence.size() == best[p.geneId].sequence.size() &&
             p.proteinId < best[p.geneId].proteinId))
            best[p.geneId] = p;
    }
    QVector<Protein> selected;
    QSet<QString> seen;
    for (auto g : genes) {
        if (seen.contains(g))
            fail("Duplicate requested gene " + g);
        seen.insert(g);
        if (!best.contains(g))
            fail("Gene not found in annotation " + g);
        selected.push_back(best[g]);
    }
    if (QFileInfo::exists(outputDir))
        fail("Query output directory already exists");
    writeProteins(outputDir, selected);
    writeText(
        outputDir + "/selection.json",
        QJsonDocument(
            QJsonObject{
                {"policy", "longest_annotated_protein_per_gene;lexicographic_protein_id_tiebreak"},
                {"annotation_sha256", sha256(dir + "/annotation.tsv")},
                {"proteome_sha256", sha256(dir + "/proteins.faa")},
                {"genes", QJsonArray::fromStringList(genes)}})
            .toJson());
}
QJsonObject verifyTrial(const QString &source, const QString &proteinDir,
                        const QString &outputDir) {
    QFile f(source);
    if (!f.open(QIODevice::ReadOnly))
        fail("Cannot read trial source");
    QTextStream in(&f);
    QString dna;
    bool collect = false;
    while (!in.atEnd()) {
        QString l = in.readLine().trimmed();
        if (l.startsWith('>')) {
            if (collect)
                break;
            collect = l.startsWith(">dmel|FBgn0262810|");
            continue;
        }
        if (collect && !l.isEmpty())
            dna += l;
    }
    const QString bases = "TCAG",
                  aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    QMap<QString, QChar> code;
    int k = 0;
    for (QChar x : bases)
        for (QChar y : bases)
            for (QChar z : bases) {
                QString c;
                c += x;
                c += y;
                c += z;
                code[c] = aa[k++];
            }
    if (dna.isEmpty() || dna.size() % 3 || !dna.startsWith("ATG"))
        fail("Invalid provisional trial CDS");
    QString protein;
    for (int i = 0; i < dna.size(); i += 3) {
        auto c = dna.mid(i, 3);
        if (!code.contains(c))
            fail("Ambiguous trial codon");
        protein += code[c];
    }
    if (!protein.endsWith('*'))
        fail("Trial has no terminal stop");
    protein.chop(1);
    if (protein.contains('*'))
        fail("Trial contains internal stop");
    auto refs = flybaseProteins(proteinDir + "/source_protein.faa");
    QStringList matches;
    for (auto p : refs)
        if (p.geneId == "FBgn0262810" && p.sequence == protein)
            matches << p.proteinId;
    QJsonObject result{{"gene_id", "FBgn0262810"},
                       {"dna_length", dna.size()},
                       {"protein_length", protein.size()},
                       {"exact_annotated_protein_matches", QJsonArray::fromStringList(matches)},
                       {"verified", !matches.isEmpty()},
                       {"annotation_sha256", sha256(proteinDir + "/source_protein.faa")},
                       {"source_sha256", sha256(source)}};
    QDir().mkpath(outputDir);
    writeText(outputDir + "/trial_verification.json", QJsonDocument(result).toJson());
    return result;
}
void prepare(const QString &manifestPath, const QString &outputDir, const Hooks &hooks) {
    QJsonObject m = readJson(manifestPath);
    if (m["schema_version"].toInt() != 1)
        fail("Invalid source manifest");
    QDir().mkpath(outputDir);
    QSet<QString> ids;
    QJsonArray preparationStates;
    QStringList failures;
    auto saveStates = [&] {
        writeText(QDir(outputDir).filePath("preparation_status.json"),
                  QJsonDocument(QJsonObject{{"schema_version", 1}, {"sources", preparationStates}})
                      .toJson());
    };
    for (auto value : m["sources"].toArray()) {
        auto s = value.toObject();
        QString id = s["id"].toString();
        if (!QRegularExpression("^[A-Za-z0-9_-]+$").match(id).hasMatch() || ids.contains(id))
            fail("Invalid source id");
        ids.insert(id);
        if (!s.value("enabled").toBool(true))
            continue;
        if ((s["provider"] != "ncbi" && s["provider"] != "flybase") ||
            s["species"].toString().isEmpty() || s["assembly"].toString().isEmpty())
            fail("Source requires valid provider, species, assembly");
        if (s["provider"] == "flybase" &&
            s["files"].toObject().value("source_protein.faa").toString().isEmpty())
            fail("FlyBase source requires protein URL");
        QString dir = QDir(outputDir).absoluteFilePath(id);
        QDir().mkpath(dir);
        QString manifestFile = dir + "/source.json";
        try {
            if (QFileInfo::exists(manifestFile)) {
                auto previous = readJson(manifestFile);
                auto cleanSpec = [](QJsonObject object) {
                    for (const auto &key : object.keys())
                        if (object.value(key).isNull())
                            object.remove(key);
                    return object;
                };
                if (cleanSpec(previous["specification"].toObject()) != cleanSpec(s))
                    fail("Cached source specification changed: " + id);
                for (auto v : previous["files"].toArray()) {
                    auto x = v.toObject();
                    if (sha256(dir + "/" + x["name"].toString()) != x["sha256"].toString())
                        fail("Cached file checksum mismatch " + id);
                }
                if (hooks.progress)
                    hooks.progress("Verified cached source " + id);
                preparationStates.append(QJsonObject{
                    {"id", id}, {"status", previous.value("status").toString("ready")}});
                saveStates();
                continue;
            }
            auto download = [&](QString url, QString dest) {
                QUrl u(url);
                if (u.scheme() != "https")
                    fail("Only HTTPS source URLs supported");
                if (hooks.progress)
                    hooks.progress("Downloading " + id + ": " + QFileInfo(dest).fileName());
                process("curl",
                        {"-f", "-s", "-S", "-L", "--retry", "2", "--connect-timeout", "30",
                         "--max-time", "1800", "-o", dest + ".part", url},
                        1900, hooks);
                if (QFileInfo::exists(dest) && !QFile::remove(dest))
                    fail("Cannot replace " + dest);
                if (!QFile::rename(dest + ".part", dest))
                    fail("Cannot finalize download");
            };
            QString folder = s.value("base_url").toString(), prefix = s.value("prefix").toString();
            if (s["provider"] == "ncbi") {
                QString acc = s["assembly"].toString();
                auto match = QRegularExpression("^(GCF|GCA)_([0-9]{9})\\.([0-9]+)$").match(acc);
                if (!match.hasMatch())
                    fail("Invalid assembly accession");
                QString n = match.captured(2);
                QString parent = "https://ftp.ncbi.nlm.nih.gov/genomes/all/" + match.captured(1) +
                                 "/" + n.mid(0, 3) + "/" + n.mid(3, 3) + "/" + n.mid(6, 3) + "/";
                download(parent, dir + "/listing.html");
                QFile lf(dir + "/listing.html");
                if (!lf.open(QIODevice::ReadOnly))
                    fail("Cannot read directory listing");
                QString listing = QString::fromUtf8(lf.readAll());
                auto re = QRegularExpression("href=\"(" + QRegularExpression::escape(acc) +
                                             "_[^\"/]+)/\"");
                auto it = re.globalMatch(listing);
                QSet<QString> names;
                while (it.hasNext())
                    names.insert(it.next().captured(1));
                if (names.size() != 1)
                    fail("Cannot resolve unique assembly directory " + acc);
                prefix = *names.begin();
                folder = parent + prefix + "/";
            }
            QJsonArray fileRecords;
            auto record = [&](QString name, QString url) {
                fileRecords.append(QJsonObject{
                    {"name", name}, {"source_url", url}, {"sha256", sha256(dir + "/" + name)}});
            };
            QString checksums;
            if (s["provider"] == "ncbi") {
                download(folder + "md5checksums.txt", dir + "/md5checksums.txt");
                QFile cf(dir + "/md5checksums.txt");
                if (!cf.open(QIODevice::ReadOnly))
                    fail("Cannot read checksums");
                checksums = QString::fromUtf8(cf.readAll());
                record("md5checksums.txt", folder + "md5checksums.txt");
            }
            QJsonObject urls = s.value("files").toObject();
            bool hasProtein = true;
            if (s["provider"] == "ncbi") {
                hasProtein = checksums.contains("./" + prefix + "_protein.faa.gz") &&
                             checksums.contains("./" + prefix + "_genomic.gff.gz");
                if (hasProtein) {
                    urls["source_protein.faa"] = folder + prefix + "_protein.faa.gz";
                    urls["source_genomic.gff"] = folder + prefix + "_genomic.gff.gz";
                } else if (!s.value("include_genome").toBool(true)) {
                    fail("No protein annotation or requested genome for " + id);
                }
                if (s.value("include_genome").toBool(true))
                    urls["genome.fna"] = folder + prefix + "_genomic.fna.gz";
            }
            for (auto key : urls.keys()) {
                if (QFileInfo(key).fileName() != key || key == "." || key == "..")
                    fail("Invalid source output name");
                QString url = urls[key].toString();
                QString packed = dir + "/" + key + ".gz";
                download(url, packed);
                if (!checksums.isEmpty()) {
                    QString filename = QUrl(url).fileName(), expected;
                    for (auto line : checksums.split('\n')) {
                        auto cols = line.simplified().split(' ');
                        if (cols.size() == 2 && cols[1] == "./" + filename)
                            expected = cols[0];
                    }
                    if (expected.isEmpty())
                        fail("Missing MD5 for " + filename);
                    QFile cf(packed);
                    if (!cf.open(QIODevice::ReadOnly))
                        fail("Cannot read downloaded archive");
                    QCryptographicHash h(QCryptographicHash::Md5);
                    if (!h.addData(&cf) || QString::fromLatin1(h.result().toHex()) != expected)
                        fail("Downloaded MD5 mismatch " + filename);
                }
                process("gzip", {"-dc", packed}, 600, hooks, dir + "/" + key);
                record(key, url);
                record(key + ".gz", url);
            }
            int proteinCount = 0;
            if (hasProtein) {
                auto proteins =
                    s["provider"] == "flybase"
                        ? flybaseProteins(dir + "/source_protein.faa")
                        : ncbiProteins(dir + "/source_protein.faa", dir + "/source_genomic.gff",
                                       s["species"].toString());
                proteinCount = proteins.size();
                writeProteins(dir, proteins);
                for (auto name : {"proteins.faa", "protein_map.tsv", "annotation.tsv"})
                    record(name, "derived");
            }
            QString status = hasProtein ? "ready" : "genome_only";
            QJsonObject metadata{
                {"specification", s},
                {"retrieved_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)},
                {"status", status},
                {"protein_count", proteinCount},
                {"limitation",
                 hasProtein ? "" : "No protein/GFF annotation in publisher checksum manifest"},
                {"files", fileRecords}};
            writeText(manifestFile, QJsonDocument(metadata).toJson());
            preparationStates.append(QJsonObject{{"id", id}, {"status", status}});
            saveStates();
        } catch (const Error &error) {
            if (QString::fromUtf8(error.what()) == "cancelled")
                throw;
            failures.append(id + ": " + QString::fromUtf8(error.what()));
            preparationStates.append(
                QJsonObject{{"id", id}, {"status", "error"}, {"message", error.what()}});
            saveStates();
            if (hooks.progress)
                hooks.progress("Source error " + id + ": " + error.what());
        }
    }
    if (!failures.isEmpty())
        fail("Some sources could not be prepared:\n" + failures.join('\n'));
}
} // namespace homology
