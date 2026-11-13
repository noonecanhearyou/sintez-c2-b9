#include "homology/homology.hpp"
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDateTime>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QProcess>
#include <QProcessEnvironment>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSet>
#include <QTextStream>
#include <algorithm>
#include <cmath>
#include <map>

namespace homology {
static void fail(const QString &s) {
    throw Error(s);
}
void writeText(const QString &path, const QByteArray &bytes) {
    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly) || f.write(bytes) != bytes.size() || !f.commit())
        fail("Cannot write " + path);
}
QJsonObject readJson(const QString &path) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        fail("Cannot read " + path);
    QJsonParseError e;
    auto d = QJsonDocument::fromJson(f.readAll(), &e);
    if (e.error != QJsonParseError::NoError || !d.isObject())
        fail("Invalid JSON " + path + ": " + e.errorString());
    return d.object();
}
QString sha256(const QString &path) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        fail("Cannot hash " + path);
    QCryptographicHash h(QCryptographicHash::Sha256);
    if (!h.addData(&f))
        fail("Hash read error " + path);
    return h.result().toHex();
}
QVector<Sequence> readFasta(const QString &path, bool protein) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        fail("Cannot read FASTA " + path);
    QVector<Sequence> out;
    QSet<QString> ids;
    QTextStream in(&f);
    const QString alphabet = protein ? "ABCDEFGHIKLMNPQRSTVWXYZUOJ" : "ACGTRYSWKMBDHVN";
    auto finish = [&] {
        if (out.isEmpty())
            return;
        auto &s = out.last();
        s.sequence = s.sequence.toUpper();
        if (protein && s.sequence.endsWith('*'))
            s.sequence.chop(1);
        if (s.sequence.isEmpty())
            fail("Empty sequence " + s.id);
        for (QChar c : s.sequence)
            if (!alphabet.contains(c))
                fail("Invalid residue in " + s.id + ": " + c);
    };
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (line.isEmpty())
            continue;
        if (line.startsWith('>')) {
            finish();
            QString h = line.mid(1).trimmed();
            QString id = h.section(QRegularExpression("\\s+"), 0, 0);
            if (id.isEmpty() || ids.contains(id))
                fail("Empty/duplicate FASTA ID " + id);
            ids.insert(id);
            out.push_back({id, h, {}});
        } else {
            if (out.isEmpty())
                fail("Sequence before FASTA header");
            line.remove(QRegularExpression("\\s"));
            out.last().sequence += line;
        }
    }
    finish();
    if (out.isEmpty())
        fail("Empty FASTA " + path);
    return out;
}
QByteArray process(const QString &exe, const QStringList &args, int timeout, const Hooks &hooks,
                   const QString &stdoutFile) {
    if (timeout <= 0)
        fail("Process timeout must be positive");
    if (hooks.cancelled && hooks.cancelled->load())
        fail("cancelled");
    QProcess p;
    // Disable BLAST usage reporting: local sequences and queries remain local.
    auto env = QProcessEnvironment::systemEnvironment();
    env.insert("BLAST_USAGE_REPORT", "false");
    p.setProcessEnvironment(env);
    if (!stdoutFile.isEmpty())
        p.setStandardOutputFile(stdoutFile);
    p.start(exe, args);
    if (!p.waitForStarted(10000))
        fail("Cannot start " + exe + ": " + p.errorString());
    QElapsedTimer timer;
    timer.start();
    QByteArray output, errors;
    while (true) {
        bool done = p.waitForFinished(100);
        output += p.readAllStandardOutput();
        errors += p.readAllStandardError();
        if (errors.size() > 65536)
            errors = errors.right(65536);
        if (done || p.state() == QProcess::NotRunning)
            break;
        if ((hooks.cancelled && hooks.cancelled->load()) ||
            timer.elapsed() > qint64(timeout) * 1000) {
            p.kill();
            p.waitForFinished();
            fail(hooks.cancelled && hooks.cancelled->load() ? "cancelled"
                                                            : "Process timeout: " + exe);
        }
    }
    if (p.exitStatus() != QProcess::NormalExit || p.exitCode() != 0)
        fail(exe + " failed: " + QString::fromUtf8(errors));
    return output;
}
qint64 coveredLength(QVector<Interval> v) {
    if (v.isEmpty())
        return 0;
    for (const auto &x : v)
        if (x.start < 1 || x.end < x.start)
            fail("Invalid interval");
    std::sort(v.begin(), v.end(), [](auto a, auto b) { return a.start < b.start; });
    qint64 n = 0, l = v[0].start, r = v[0].end;
    for (int i = 1; i < v.size(); ++i) {
        if (v[i].start <= r + 1)
            r = std::max(r, v[i].end);
        else {
            n += r - l + 1;
            l = v[i].start;
            r = v[i].end;
        }
    }
    return n + r - l + 1;
}
static QString pathAt(const QDir &d, const QString &p) {
    return QFileInfo(p).isAbsolute() ? p : d.absoluteFilePath(p);
}
static QString cell(QString s) {
    s.replace('\t', ' ');
    s.replace('\n', ' ');
    s.replace('\r', ' ');
    return s;
}
static QString number(double n) {
    return QString::number(n, 'g', 10);
}
struct Hit {
    QString query, target, targetGene, species, type, assembly;
    double identity, evalue, bitscore;
    qint64 length, qlen, slen, qs, qe, ss, se;
    int frame;
    QString exclusion, quality;
    double qcov = 0, scov = 0;
    bool pass = false;
};
QJsonObject run(const QString &configPath, const Hooks &hooks) {
    const QJsonObject cfg = readJson(configPath);
    QDir base(QFileInfo(configPath).absolutePath());
    if (cfg.value("schema_version").toInt() != 1)
        fail("Unsupported schema_version");
    auto required = [&](QString k) {
        auto s = cfg.value(k).toString();
        if (s.isEmpty())
            fail("Missing " + k);
        return pathAt(base, s);
    };
    QString queryFile = required("query_fasta"), queryMapFile = required("query_map"),
            outDir = required("output_dir");
    auto queries = readFasta(queryFile, true);
    QMap<QString, QString> queryGenes;
    QFile qm(queryMapFile);
    if (!qm.open(QIODevice::ReadOnly))
        fail("Cannot read query map");
    QTextStream qms(&qm);
    if (qms.readLine() != "protein_id\tgene_id")
        fail("query_map header must be protein_id TAB gene_id");
    while (!qms.atEnd()) {
        auto c = qms.readLine().split('\t');
        if (c.size() != 2 || c[0].isEmpty() || c[1].isEmpty() || queryGenes.contains(c[0]))
            fail("Invalid/duplicate query map row");
        queryGenes[c[0]] = c[1];
    }
    QMap<QString, qint64> queryLengths;
    for (auto q : queries) {
        if (!queryGenes.contains(q.id))
            fail("Missing gene mapping " + q.id);
        queryLengths[q.id] = q.sequence.size();
    }
    if (queryGenes.size() != queries.size())
        fail("Query map has extra entries");
    QString querySpecies = cfg.value("query_species").toString();
    if (querySpecies.isEmpty())
        fail("Missing query_species");
    const auto filters = cfg.value("filters").toObject();
    double searchE = filters.value("search_evalue").toDouble(1e-3),
           acceptE = filters.value("accept_evalue").toDouble(1e-5);
    double minQ = filters.value("min_query_coverage").toDouble(50),
           minId = filters.value("min_identity").toDouble(0);
    int threads = cfg.value("threads").toInt(2), timeout = cfg.value("timeout_seconds").toInt(3600);
    int mtMode = cfg.value("mt_mode").toInt(0);
    if (mtMode < 0 || mtMode > 2)
        fail("mt_mode must be 0, 1 or 2");
    if (!std::isfinite(searchE) || !std::isfinite(acceptE) || searchE <= 0 || acceptE <= 0 ||
        acceptE > searchE || minQ < 0 || minQ > 100 || minId < 0 || minId > 100 || threads < 1 ||
        timeout < 1)
        fail("Invalid thresholds/threads/timeout");
    QString bin = cfg.value("blast_bin").toString();
    if (!bin.isEmpty())
        bin = pathAt(base, bin);
    auto executable = [&](const QString &n) { return bin.isEmpty() ? n : QDir(bin).filePath(n); };
    auto targets = cfg.value("targets").toArray();
    if (targets.isEmpty())
        fail("No targets");
    QSet<QString> targetIds;
    QJsonArray normalized;
    for (auto v : targets) {
        auto t = v.toObject();
        QString id = t.value("id").toString(), sp = t.value("species").toString(),
                kind = t.value("kind").toString();
        if (!QRegularExpression("^[A-Za-z0-9_-]+$").match(id).hasMatch() || targetIds.contains(id))
            fail("Invalid/duplicate target ID " + id);
        targetIds.insert(id);
        if (sp.isEmpty() || t.value("assembly").toString().isEmpty() ||
            (kind != "protein" && kind != "genome"))
            fail("Target requires species, assembly, kind");
        if (t.value("enabled").toBool(true)) {
            QString f = t.value("fasta").toString();
            if (f.isEmpty())
                fail("Missing target fasta");
            f = pathAt(base, f);
            if (!QFileInfo::exists(f))
                fail("Missing target " + f);
            t["fasta"] = f;
            QString m = t.value("protein_map").toString();
            if (kind == "protein" && sp == querySpecies && m.isEmpty())
                fail("Within-species search requires protein_map to exclude same-gene isoforms");
            if (!m.isEmpty())
                t["protein_map"] = pathAt(base, m);
        }
        normalized.append(t);
    }
    if (QFileInfo::exists(outDir))
        fail("Output already exists; choose a new output_dir: " + outDir);
    if (!QDir().mkpath(outDir + "/raw") || !QDir().mkpath(outDir + "/db"))
        fail("Cannot create output directory");
    QJsonObject manifest{
        {"schema_version", 1},
        {"module_version", "0.2.0"},
        {"qt_version", qVersion()},
        {"configuration_sha256", sha256(configPath)},
        {"blast_usage_reporting", false},
        {"started_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)},
        {"status", "running"},
        {"config", cfg},
        {"query_sha256", sha256(queryFile)},
        {"query_map_sha256", sha256(queryMapFile)},
        {"coordinate_system",
         "1-based inclusive; genome subject coordinates are on complete input contigs"}};
    auto saveManifest = [&] { writeText(outDir + "/run.json", QJsonDocument(manifest).toJson()); };
    saveManifest();
    QJsonArray provenance;
    QJsonObject versions;
    QStringList statusRows;
    bool anyErrors = false;
    int completed = 0, errors = 0, skipped = 0;
    QMap<QString, QSet<QString>> paralogs, proteinSpecies, genomeSpecies, partialGenes,
        partialProteinSpecies, partialGenomeSpecies;
    QMap<QString, Hit> bestProtein, bestGenome;
    QMap<QString, QString> bestProteinTarget, bestGenomeTarget;
    auto logCommand = [&](const QString &exe, const QStringList &args) {
        provenance.append(
            QJsonObject{{"executable", exe}, {"arguments", QJsonArray::fromStringList(args)}});
    };
    const QString statusHeader = "gene_id\tquery_id\ttarget_id\ttarget_species\tsearch_"
                                 "type\tstatus\traw_hsp_count\taccepted_hsp_count\tmessage\n";
    const QString hitHeader =
        "gene_id\tquery_id\ttarget_id\ttarget_species\tassembly\tsearch_type\ttarget_sequence_"
        "id\ttarget_gene_id\tidentity\talignment_length\tquery_length\ttarget_length\tquery_"
        "start\tquery_end\ttarget_start\ttarget_end\tstrand\tsubject_frame\tevalue\tbitscore\thsp_"
        "query_coverage\tquery_coverage\ttarget_coverage\tpasses_filters\texclusion_"
        "reason\tquality\n";
    QString hitText = hitHeader,
            locusText =
                "gene_id\tquery_id\ttarget_id\ttarget_species\tassembly\tcontig\tstart_1based\tend_"
                "1based\tstart_0based\tend_0based_exclusive\tstrand\tevalue\tbitscore\thsp_query_"
                "coverage\tquery_coverage\tquality\n";
    try {
        for (QString program : {"blastp", "tblastn", "makeblastdb"})
            versions[program] =
                QString::fromUtf8(process(executable(program), {"-version"}, 30, hooks)).trimmed();
        manifest["versions"] = versions;
        for (auto tv : normalized) {
            auto t = tv.toObject();
            QString id = t["id"].toString(), sp = t["species"].toString(),
                    kind = t["kind"].toString(), assembly = t["assembly"].toString();
            QString type = kind == "genome" ? "tblastn"
                                            : (sp == querySpecies ? "within_species_blastp"
                                                                  : "cross_species_blastp");
            auto addStatuses = [&](QString state, QString message, const QVector<Hit> &hits) {
                for (auto q : queries) {
                    int raw = 0, accepted = 0;
                    for (auto h : hits)
                        if (h.query == q.id) {
                            ++raw;
                            if (h.pass)
                                ++accepted;
                        }
                    QString st = state;
                    if (st == "completed")
                        st = accepted ? "completed_with_hits" : "completed_no_accepted_hits";
                    statusRows << QStringList{queryGenes[q.id],
                                              q.id,
                                              id,
                                              sp,
                                              type,
                                              st,
                                              QString::number(raw),
                                              QString::number(accepted),
                                              cell(message)}
                                      .join('\t');
                }
            };
            if (!t.value("enabled").toBool(true)) {
                ++skipped;
                addStatuses("not_run",
                            t.value("disabled_reason").toString("disabled in configuration"), {});
                continue;
            }
            if (hooks.progress)
                hooks.progress("Searching " + id + " (" + type + ")");
            QVector<Hit> hits;
            try {
                QString fasta = t["fasta"].toString(), db = outDir + "/db/" + id;
                auto info = t;
                info["fasta_sha256"] = sha256(fasta);
                provenance.append(QJsonObject{{"input", info}});
                QMap<QString, QString> targetGenes;
                if (kind == "protein") {
                    auto seqs = readFasta(fasta, true);
                    QSet<QString> ids;
                    for (auto s : seqs)
                        ids.insert(s.id);
                    QString mapPath = t["protein_map"].toString();
                    if (!mapPath.isEmpty()) {
                        QFile f(mapPath);
                        if (!f.open(QIODevice::ReadOnly))
                            fail("Cannot read protein_map");
                        QTextStream s(&f);
                        if (s.readLine() != "protein_id\tgene_id")
                            fail("Invalid protein_map header");
                        while (!s.atEnd()) {
                            auto c = s.readLine().split('\t');
                            if (c.size() != 2 || c[0].isEmpty() || c[1].isEmpty() ||
                                targetGenes.contains(c[0]))
                                fail("Invalid protein_map row");
                            targetGenes[c[0]] = c[1];
                        }
                        for (auto i : ids)
                            if (!targetGenes.contains(i))
                                fail("Unmapped target protein " + i);
                        provenance.append(
                            QJsonObject{{"protein_map", mapPath}, {"sha256", sha256(mapPath)}});
                    }
                }
                QStringList mk{"-in",  fasta, "-dbtype", kind == "protein" ? "prot" : "nucl",
                               "-out", db,    "-title",  id};
                logCommand(executable("makeblastdb"), mk);
                process(executable("makeblastdb"), mk, timeout, hooks);
                QString program = kind == "protein" ? "blastp" : "tblastn",
                        raw = outDir + "/raw/" + id + ".tsv";
                QStringList args{"-query",
                                 queryFile,
                                 "-db",
                                 db,
                                 "-evalue",
                                 number(searchE),
                                 "-num_threads",
                                 QString::number(threads),
                                 "-max_target_seqs",
                                 "1000000",
                                 "-outfmt",
                                 "6 qseqid sseqid pident length qlen slen qstart qend sstart send "
                                 "evalue bitscore sframe",
                                 "-out",
                                 raw};
                if (cfg.contains("mt_mode"))
                    args << "-mt_mode" << QString::number(mtMode);
                logCommand(executable(program), args);
                process(executable(program), args, timeout, hooks);
                QFile f(raw);
                if (!f.open(QIODevice::ReadOnly))
                    fail("Missing BLAST output");
                QTextStream s(&f);
                while (!s.atEnd()) {
                    auto c = s.readLine().split('\t');
                    if (c.size() != 13)
                        fail("Malformed BLAST output");
                    auto integer = [&](int i) {
                        bool ok = false;
                        qint64 x = c[i].toLongLong(&ok);
                        if (!ok)
                            fail("Invalid integer in BLAST output");
                        return x;
                    };
                    auto real = [&](int i) {
                        bool ok = false;
                        double x = c[i].toDouble(&ok);
                        if (!ok || !std::isfinite(x))
                            fail("Invalid number in BLAST output");
                        return x;
                    };
                    Hit h;
                    h.query = c[0];
                    h.target = c[1];
                    h.targetGene = targetGenes.value(h.target);
                    h.species = sp;
                    h.type = type;
                    h.assembly = assembly;
                    h.identity = real(2);
                    h.length = integer(3);
                    h.qlen = integer(4);
                    h.slen = integer(5);
                    h.qs = integer(6);
                    h.qe = integer(7);
                    h.ss = integer(8);
                    h.se = integer(9);
                    h.evalue = real(10);
                    h.bitscore = real(11);
                    h.frame = int(integer(12));
                    if (!queryLengths.contains(h.query) || h.qlen != queryLengths[h.query] ||
                        h.qs < 1 || h.qe < h.qs || h.qe > h.qlen || h.ss < 1 || h.se < 1 ||
                        h.ss > h.slen || h.se > h.slen || h.length < 1 || h.identity < 0 ||
                        h.identity > 100 || h.evalue < 0)
                        fail("Inconsistent BLAST coordinates/lengths");
                    if (type == "within_species_blastp" && queryGenes[h.query] == h.targetGene)
                        h.exclusion = "same_gene_or_isoform";
                    hits.push_back(h);
                }
                // Union coverage uses only statistically accepted HSPs passing identity.
                // Genomic HSPs are never joined across introns or distant loci automatically.
                QMap<QString, QVector<Interval>> qintervals, sintervals;
                for (auto h : hits)
                    if (h.exclusion.isEmpty() && h.evalue <= acceptE && h.identity >= minId) {
                        QString key = h.query + "\t" + h.target;
                        qintervals[key].push_back({h.qs, h.qe});
                        sintervals[key].push_back({std::min(h.ss, h.se), std::max(h.ss, h.se)});
                    }
                for (auto &h : hits) {
                    QString key = h.query + "\t" + h.target;
                    h.qcov = kind == "protein" ? 100.0 * coveredLength(qintervals[key]) / h.qlen
                                               : 100.0 * (h.qe - h.qs + 1) / h.qlen;
                    h.scov = kind == "protein" ? 100.0 * coveredLength(sintervals[key]) / h.slen
                                               : 100.0 * (std::abs(h.se - h.ss) + 1) / h.slen;
                    if (h.exclusion.isEmpty()) {
                        if (h.evalue > acceptE)
                            h.exclusion = "evalue";
                        else if (h.identity < minId)
                            h.exclusion = "identity";
                        else if (h.qcov < minQ)
                            h.exclusion = "query_coverage";
                    }
                    h.pass = h.exclusion.isEmpty();
                    h.quality = h.qcov >= 80 ? "broad_coverage_candidate" : "partial_candidate";
                    if (!h.pass)
                        h.quality = "filtered";
                    if (h.exclusion == "query_coverage") {
                        h.quality = "significant_partial_similarity";
                        QString gene = queryGenes[h.query];
                        if (type == "within_species_blastp")
                            partialGenes[gene].insert(h.targetGene);
                        else if (kind == "protein")
                            partialProteinSpecies[gene].insert(sp);
                        else
                            partialGenomeSpecies[gene].insert(sp);
                    }
                    QString strand = kind == "genome" ? (h.frame < 0 ? "-" : "+") : ".";
                    hitText += QStringList{queryGenes[h.query],
                                           h.query,
                                           id,
                                           sp,
                                           assembly,
                                           type,
                                           h.target,
                                           h.targetGene,
                                           number(h.identity),
                                           QString::number(h.length),
                                           QString::number(h.qlen),
                                           QString::number(h.slen),
                                           QString::number(h.qs),
                                           QString::number(h.qe),
                                           QString::number(std::min(h.ss, h.se)),
                                           QString::number(std::max(h.ss, h.se)),
                                           strand,
                                           QString::number(h.frame),
                                           number(h.evalue),
                                           number(h.bitscore),
                                           number(100.0 * (h.qe - h.qs + 1) / h.qlen),
                                           number(h.qcov),
                                           number(h.scov),
                                           h.pass ? "true" : "false",
                                           h.exclusion,
                                           h.quality}
                                   .join('\t') +
                               '\n';
                    if (h.pass) {
                        QString gene = queryGenes[h.query];
                        auto &best = kind == "protein" ? bestProtein : bestGenome;
                        auto &bestTarget = kind == "protein" ? bestProteinTarget : bestGenomeTarget;
                        if (!best.contains(gene) || h.bitscore > best[gene].bitscore) {
                            best[gene] = h;
                            bestTarget[gene] = id;
                        }
                        if (type == "within_species_blastp")
                            paralogs[gene].insert(h.targetGene);
                        else if (kind == "protein")
                            proteinSpecies[gene].insert(sp);
                        else {
                            genomeSpecies[gene].insert(sp);
                            locusText += QStringList{gene,
                                                     h.query,
                                                     id,
                                                     sp,
                                                     assembly,
                                                     h.target,
                                                     QString::number(std::min(h.ss, h.se)),
                                                     QString::number(std::max(h.ss, h.se)),
                                                     QString::number(std::min(h.ss, h.se) - 1),
                                                     QString::number(std::max(h.ss, h.se)),
                                                     strand,
                                                     number(h.evalue),
                                                     number(h.bitscore),
                                                     number(100.0 * (h.qe - h.qs + 1) / h.qlen),
                                                     number(h.qcov),
                                                     h.quality}
                                             .join('\t') +
                                         '\n';
                        }
                    }
                }
                ++completed;
                addStatuses("completed", "", hits);
            } catch (const Error &e) {
                if (QString::fromUtf8(e.what()) == "cancelled")
                    throw;
                anyErrors = true;
                ++errors;
                addStatuses("error", QString::fromUtf8(e.what()), {});
                if (hooks.progress)
                    hooks.progress("Error in " + id + ": " + e.what());
            }
            manifest["commands_and_inputs"] = provenance;
            saveManifest();
        }
        auto sorted = [](const QSet<QString> &s) {
            auto l = s.values();
            l.sort();
            return l.join(',');
        };
        QString summary =
            "gene_id\tparalog_candidate_count\tparalog_gene_ids\tprotein_hit_species\tgenome_hit_"
            "species\tsearch_status\tcompleted_targets\terror_targets\tnot_run_"
            "targets\tlimitations\tpartial_similarity_gene_ids\tpartial_protein_species\tpartial_"
            "genome_species\tbest_protein_target\tbest_protein_sequence_id\tbest_protein_"
            "evalue\tbest_genome_target\tbest_genome_contig\tbest_genome_evalue\n";
        QSet<QString> genes;
        for (auto g : queryGenes)
            genes.insert(g);
        auto gs = genes.values();
        gs.sort();
        QString state = anyErrors ? (completed ? "partial_error" : "error")
                                  : (skipped ? "partial_not_run" : "completed");
        for (auto g : gs) {
            QString limitations = "homology_candidates_only;orthology_unresolved;no_origin_"
                                  "classification;genomic_HSPs_not_joined;thresholds_provisional";
            for (auto q : queries)
                if (queryGenes[q.id] == g && q.sequence.size() < 30)
                    limitations += ";short_query_default_BLAST_settings";
            summary +=
                QStringList{g,
                            QString::number(paralogs[g].size()),
                            sorted(paralogs[g]),
                            sorted(proteinSpecies[g]),
                            sorted(genomeSpecies[g]),
                            state,
                            QString::number(completed),
                            QString::number(errors),
                            QString::number(skipped),
                            limitations,
                            sorted(partialGenes[g]),
                            sorted(partialProteinSpecies[g]),
                            sorted(partialGenomeSpecies[g]),
                            bestProteinTarget.value(g),
                            bestProtein.contains(g) ? bestProtein[g].target : QString(),
                            bestProtein.contains(g) ? number(bestProtein[g].evalue) : QString(),
                            bestGenomeTarget.value(g),
                            bestGenome.contains(g) ? bestGenome[g].target : QString(),
                            bestGenome.contains(g) ? number(bestGenome[g].evalue) : QString()}
                    .join('\t') +
                '\n';
        }
        writeText(outDir + "/homology_hits.tsv", hitText.toUtf8());
        writeText(outDir + "/homology_summary.tsv", summary.toUtf8());
        writeText(outDir + "/search_status.tsv",
                  (statusHeader + statusRows.join('\n') + '\n').toUtf8());
        writeText(outDir + "/candidate_loci.tsv", locusText.toUtf8());
        QJsonObject outputHashes;
        for (QString filename : {"homology_hits.tsv", "homology_summary.tsv", "candidate_loci.tsv",
                                 "search_status.tsv"})
            outputHashes[filename] = sha256(outDir + "/" + filename);
        manifest["output_sha256"] = outputHashes;
        manifest["query_count"] = queries.size();
        manifest["target_count"] = normalized.size();
        manifest["status"] = state;
        manifest["finished_utc"] = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);
        manifest["commands_and_inputs"] = provenance;
        saveManifest();
    } catch (const std::exception &e) {
        manifest["status"] = QString::fromUtf8(e.what()) == "cancelled" ? "cancelled" : "error";
        manifest["error"] = e.what();
        manifest["commands_and_inputs"] = provenance;
        saveManifest();
        throw;
    }
    return manifest;
}
} // namespace homology
