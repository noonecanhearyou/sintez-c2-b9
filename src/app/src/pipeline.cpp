#include "core.hpp"
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QTextStream>
#include <QUrl>
#include <QUuid>
#include <algorithm>
#include <orf/orf.hpp>
#include <synteny/synteny.hpp>
namespace go {
namespace {
QString abs(const QString &p) {
    return QFileInfo(p).absoluteFilePath();
}
QJsonObject json(const orf::J &value) {
    return QJsonDocument::fromJson(QByteArray::fromStdString(value.dump())).object();
}
void saveJson(QString path, QJsonObject v) {
    homology::writeText(path, QJsonDocument(v).toJson());
}
QJsonArray array(const QSet<QString> &s) {
    auto l = s.values();
    l.sort();
    return QJsonArray::fromStringList(l);
}
QString joins(const QSet<QString> &s) {
    auto l = s.values();
    l.sort();
    return l.join(',');
}
QSet<QString> fromList(const QString &s) {
    auto list = s.split(',', Qt::SkipEmptyParts);
    return QSet<QString>(list.begin(), list.end());
}
struct Queries {
    QString dir;
    QSet<QString> valid;
    QMap<QString, QString> native, missing;
};
QString location(QString contig, QString strand, QVector<QPair<qint64, qint64>> blocks) {
    std::sort(blocks.begin(), blocks.end());
    blocks.erase(std::unique(blocks.begin(), blocks.end()), blocks.end());
    QStringList b;
    for (auto x : blocks)
        b << QString::number(x.first) + ".." + QString::number(x.second);
    QString text = b.size() > 1 ? "join(" + b.join(',') + ")" : b.value(0);
    return contig + ":" + (strand == "-" ? "complement(" + text + ")" : text);
}
Queries queries(const Request &r, const QVector<Record> &input, const QJsonObject &source,
                const QString &out, std::function<void(QString)> warn) {
    Queries q;
    q.dir = out + "/queries";
    QString sourceDir = r.workspace + "/homology/data/" + source.value("id").toString();
    auto proteins = homology::readFasta(sourceDir + "/proteins.faa", true);
    QMap<QString, homology::Sequence> byId;
    QMap<QString, QStringList> bySequence;
    for (auto p : proteins) {
        byId[p.id] = p;
        bySequence[p.sequence] << p.id;
    }
    QMap<QString, Row> anno;
    for (auto row : table(sourceDir + "/annotation.tsv"))
        anno[row["protein_id"]] = row;
    QVector<homology::Protein> picked;
    QMap<QString, QString> aliases;
    QMap<QString, QString> proteinForGene;
    if (r.reference == "dmel") {
        QMap<QString, Row> best;
        for (auto row : anno) {
            QString gene = row["gene_id"];
            if (!r.genes.contains(gene))
                continue;
            if (!best.contains(gene) || row["length"].toInt() > best[gene]["length"].toInt() ||
                (row["length"] == best[gene]["length"] &&
                 row["protein_id"] < best[gene]["protein_id"]))
                best[gene] = row;
        }
        for (auto gene : r.genes) {
            if (!best.contains(gene)) {
                q.missing[gene] = "Ген отсутствует в аннотации референса";
                continue;
            }
            auto row = best[gene];
            picked << homology::Protein{row["protein_id"],
                                        gene,
                                        row["transcript_id"],
                                        r.reference,
                                        byId[row["protein_id"]].sequence,
                                        row["location"]};
            q.native[gene] = gene;
            aliases[gene] = gene;
        }
    } else {
        for (auto gene : r.genes) {
            QVector<Record> candidates;
            for (auto record : input)
                if (record.gene == gene && record.species == r.reference &&
                    !translate(record.dna).isEmpty())
                    candidates << record;
            bool raw = false;
            for (auto x : candidates)
                if (x.section == "raw" && !x.locus.contains(':'))
                    raw = true;
            QSet<QString> nativeGenes;
            QStringList matches;
            for (auto record : candidates) {
                if (raw && (record.section != "raw" || record.locus.contains(':')))
                    continue;
                for (auto id : bySequence.value(translate(record.dna))) {
                    if (!anno.contains(id))
                        continue;
                    nativeGenes.insert(anno[id]["gene_id"]);
                    matches << id;
                }
            }
            if (nativeGenes.size() != 1) {
                q.missing[gene] = nativeGenes.isEmpty()
                                      ? "Нет полного совпадения CDS с аннотированным белком "
                                        "выбранного референса"
                                      : "Несколько генов референса имеют совпадающие белки";
                continue;
            }
            QString native = *nativeGenes.begin();
            if (aliases.contains(native)) {
                q.missing[gene] = "Один ген референса соответствует нескольким входным кандидатам";
                continue;
            }
            std::sort(matches.begin(), matches.end(), [&](auto a, auto b) {
                if (byId[a].sequence.size() != byId[b].sequence.size())
                    return byId[a].sequence.size() > byId[b].sequence.size();
                return a < b;
            });
            QString pid = matches.first();
            auto row = anno[pid];
            q.native[gene] = native;
            aliases[native] = gene;
            proteinForGene[gene] = pid;
            picked << homology::Protein{pid,
                                        gene,
                                        row["transcript_id"].isEmpty() ? pid : row["transcript_id"],
                                        r.reference,
                                        byId[pid].sequence,
                                        {}};
        }
        QMap<QString, QVector<QPair<qint64, qint64>>> blocks;
        QMap<QString, QString> contigs, strands;
        QSet<QString> wanted;
        for (auto p : picked)
            wanted.insert(p.proteinId);
        QFile f(sourceDir + "/source_genomic.gff");
        if (!f.open(QIODevice::ReadOnly))
            throw homology::Error("Для смены референса требуется GFF-аннотация: " + sourceDir);
        QTextStream in(&f);
        while (!in.atEnd()) {
            auto line = in.readLine();
            if (line.startsWith("##FASTA"))
                break;
            if (line.startsWith('#'))
                continue;
            auto c = line.split('\t');
            if (c.size() != 9 || c[2] != "CDS")
                continue;
            auto m = QRegularExpression("(?:^|;)protein_id=([^;]+)").match(c[8]);
            QString pid = QUrl::fromPercentEncoding(m.captured(1).toUtf8());
            if (!wanted.contains(pid))
                continue;
            if (contigs.contains(pid) && (contigs[pid] != c[0] || strands[pid] != c[6]))
                throw homology::Error("CDS белка находится на нескольких контигах: " + pid);
            contigs[pid] = c[0];
            strands[pid] = c[6];
            bool a, b;
            auto l = c[3].toLongLong(&a), u = c[4].toLongLong(&b);
            if (!a || !b || l < 1 || u < l)
                throw homology::Error("Неверные координаты CDS в GFF");
            blocks[pid].append({l, u});
        }
        QVector<homology::Protein> valid;
        for (auto p : picked) {
            if (blocks[p.proteinId].isEmpty()) {
                q.missing[p.geneId] = "Для белка не найдены CDS-блоки в GFF";
                continue;
            }
            p.location = location(contigs[p.proteinId], strands[p.proteinId], blocks[p.proteinId]);
            valid << p;
        }
        picked = valid;
    }
    if (picked.isEmpty())
        throw homology::Error("Для выбранных генов не найдены однозначные аннотированные белки у " +
                              displaySpecies(r.reference) +
                              ". Выберите другие гены или организм. Неподтверждённые "
                              "последовательности не подставляются.");
    homology::writeProteins(q.dir, picked);
    QVector<Row> nativeRows, ownMap;
    for (auto p : picked) {
        q.valid.insert(p.geneId);
        nativeRows << Row{{"gene_id", p.geneId},
                          {"reference_gene_id", q.native[p.geneId]},
                          {"reference_protein_id", p.proteinId},
                          {"reference_species", r.reference}};
    }
    for (auto row : table(sourceDir + "/protein_map.tsv")) {
        QString native = row["gene_id"];
        if (aliases.contains(native) && q.valid.contains(aliases[native]))
            row["gene_id"] = aliases[native];
        ownMap << row;
    }
    writeTable(q.dir + "/primary_target_map.tsv", {"protein_id", "gene_id"}, ownMap);
    writeTable(q.dir + "/reference_gene_map.tsv",
               {"gene_id", "reference_gene_id", "reference_protein_id", "reference_species"},
               nativeRows);
    saveJson(q.dir + "/selection.json",
             {{"genes", array(q.valid)},
              {"policy", "canonical annotated peptide; exact CDS/protein matching for "
                         "changed reference"},
              {"reference_species", r.reference},
              {"annotation_sha256", homology::sha256(sourceDir + "/annotation.tsv")},
              {"proteome_sha256", homology::sha256(sourceDir + "/proteins.faa")}});
    for (auto gene : r.genes)
        if (!q.valid.contains(gene)) {
            if (!q.missing.contains(gene))
                q.missing[gene] = "Нет подтверждённой последовательности";
            warn(gene + ": " + q.missing[gene]);
        }
    return q;
}
void recomputeSummary(QString dir, const QVector<Row> &prototype) {
    auto hits = table(dir + "/homology_hits.tsv"), states = table(dir + "/search_status.tsv");
    QMap<QString, QSet<QString>> paralogs, protein, genome, partialP, partialG, partialOwn;
    QMap<QString, Row> bestP, bestG;
    QMap<QString, QMap<QString, int>> counts;
    for (auto s : states)
        counts[s["gene_id"]][s["status"].startsWith("completed_") ? "completed_targets"
                             : s["status"] == "not_run"           ? "not_run_targets"
                                                                  : "error_targets"]++;
    for (auto h : hits) {
        QString g = h["gene_id"], sp = h["target_species"], type = h["search_type"];
        if (h["passes_filters"] == "true") {
            if (type == "within_species_blastp")
                paralogs[g].insert(h["target_gene_id"]);
            else if (type == "tblastn")
                genome[g].insert(sp);
            else
                protein[g].insert(sp);
            auto &best = type == "tblastn" ? bestG : bestP;
            if (!best.contains(g) || h["bitscore"].toDouble() > best[g]["bitscore"].toDouble())
                best[g] = h;
        }
        if (h["quality"] == "significant_partial_similarity") {
            if (type == "within_species_blastp")
                partialOwn[g].insert(h["target_gene_id"]);
            else if (type == "tblastn")
                partialG[g].insert(sp);
            else
                partialP[g].insert(sp);
        }
    }
    QVector<Row> rows;
    for (auto row : prototype) {
        auto g = row["gene_id"];
        row["paralog_candidate_count"] = QString::number(paralogs[g].size());
        row["paralog_gene_ids"] = joins(paralogs[g]);
        row["protein_hit_species"] = joins(protein[g]);
        row["genome_hit_species"] = joins(genome[g]);
        row["partial_similarity_gene_ids"] = joins(partialOwn[g]);
        row["partial_protein_species"] = joins(partialP[g]);
        row["partial_genome_species"] = joins(partialG[g]);
        for (auto key : {"completed_targets", "error_targets", "not_run_targets"})
            row[key] = QString::number(counts[g][key]);
        row["search_status"] = counts[g]["error_targets"] ? "partial_error" : "completed";
        for (auto kind : {QString("protein"), QString("genome")}) {
            auto b = kind == "protein" ? bestP[g] : bestG[g];
            row["best_" + kind + "_target"] = b["target_id"];
            row["best_" + kind + (kind == "protein" ? "_sequence_id" : "_contig")] =
                b["target_sequence_id"];
            row["best_" + kind + "_evalue"] = b["evalue"];
        }
        rows << row;
    }
    if (rows.isEmpty())
        throw homology::Error("Пустая сводка гомологии");
    QFile f(dir + "/homology_summary.tsv");
    if (!f.open(QIODevice::ReadOnly))
        throw homology::Error("Не удалось прочитать таблицу результата");
    auto columns = QString::fromUtf8(f.readLine()).trimmed().split('\t');
    f.close();
    writeTable(dir + "/homology_summary.tsv", columns, rows);
}
} // namespace
Result execute(const Request &r, const homology::Hooks &hooks,
               std::function<void(int, QString)> progress) {
    auto check = [&] {
        if (hooks.cancelled && hooks.cancelled->load())
            throw homology::Error("cancelled");
    };
    if (r.genes.isEmpty())
        throw homology::Error("Выберите хотя бы один ген.");
    if (r.outputRoot.isEmpty())
        throw homology::Error("Выберите каталог сохранения.");
    for (auto p : {r.hal, r.alignment, r.assemblies})
        if (!QFileInfo(p).isFile())
            throw homology::Error("Исходный файл не найден: " + p);
    QString out = QDir(r.outputRoot)
                      .absoluteFilePath("LibreGenesis_" +
                                        QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss") +
                                        "_" + QUuid::createUuid().toString(QUuid::Id128).left(5));
    if (!QDir().mkpath(out))
        throw homology::Error("Нет прав на запись в каталог результатов.");
    QFile log(out + "/events.jsonl");
    if (!log.open(QIODevice::WriteOnly))
        throw homology::Error("Не удалось создать журнал событий");
    auto event = [&](QString severity, QString stage, QString code, QString text) {
        QJsonObject line{{"time_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs)},
                         {"severity", severity},
                         {"stage", stage},
                         {"code", code},
                         {"message", text}};
        log.write(QJsonDocument(line).toJson(QJsonDocument::Compact) + '\n');
        log.flush();
        if (hooks.progress)
            hooks.progress(severity + " | " + stage + " | " + text);
    };
    auto update = [&](int n, QString text) {
        if (progress)
            progress(n, text);
        check();
    };
    QJsonObject manifest{{"application", "LibreGenesis"},
                         {"version", "1.0.0"},
                         {"status", "running"},
                         {"started_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)},
                         {"reference_species", r.reference},
                         {"requested_genes", QJsonArray::fromStringList(r.genes)},
                         {"output_directory", out},
                         {"support_score_kind", "uncalibrated heuristic; not probability"}};
    auto save = [&] { saveJson(out + "/run.json", manifest); };
    save();
    event("INFO", "Запуск", "run.started", "Результаты сохраняются в " + out);
    try {
        update(2, "Проверка входных файлов");
        QFile hal(r.hal);
        if (!hal.open(QIODevice::ReadOnly))
            throw homology::Error("Нет доступа к HAL-файлу");
        if (hal.read(8) != QByteArray::fromHex("894844460d0a1a0a"))
            throw homology::Error("Выбранный HAL не является HDF5-выравниванием.");
        auto cat = catalog(r.alignment);
        QSet<QString> selection(r.genes.begin(), r.genes.end());
        for (auto g : r.genes)
            if (!cat.genes.contains(g))
                throw homology::Error("Ген отсутствует во входном файле: " + g);
        if (!cat.species.contains(r.reference))
            throw homology::Error("В последовательностях нет основного организма: " + r.reference);
        auto assemblyRows = table(r.assemblies, ',');
        QMap<QString, QString> csvAssemblies;
        for (auto a : assemblyRows) {
            if (a["Species"].isEmpty() || a["Assembly"].isEmpty())
                throw homology::Error("CSV должен содержать столбцы Species и Assembly.");
            csvAssemblies[a["Species"]] = a["Assembly"];
        }
        auto allSources = homology::readJson(r.workspace + "/homology/examples/sources.all.json");
        QJsonArray sources;
        QJsonObject primarySource;
        QSet<QString> scope;
        for (auto v : allSources["sources"].toArray()) {
            auto s = v.toObject();
            auto sp = s["species"].toString();
            if (sp == "aaeg" || (!r.comparisonSpecies.isEmpty() &&
                                 !r.comparisonSpecies.contains(sp) && sp != r.reference))
                continue;
            QString name = displaySpecies(sp);
            if (!csvAssemblies.contains(name))
                continue;
            QString acc = s["assembly"].toString();
            if (sp == "dmel") {
                if (!csvAssemblies[name].contains("r6.36"))
                    throw homology::Error("Для dmel требуется согласованная аннотация r6.36.");
            } else if (csvAssemblies[name] != acc)
                throw homology::Error("Версия сборки " + name +
                                      " в CSV отличается от конфигурации. Подготовьте "
                                      "согласованный источник.");
            sources.append(s);
            scope.insert(sp);
            if (sp == r.reference)
                primarySource = s;
        }
        if (primarySource.isEmpty() || scope.size() < 2)
            throw homology::Error("Нужны аннотация референса и хотя бы один вид сравнения из CSV.");
        saveJson(out + "/sources.json", {{"schema_version", 1}, {"sources", sources}});
        manifest["scope_species"] = array(scope);
        manifest["excluded_species"] = QJsonArray{"aaeg"};
        manifest["input_hashes"] = QJsonObject{{"alignment", homology::sha256(r.alignment)},
                                               {"assemblies", homology::sha256(r.assemblies)}};
        manifest["hal_file"] = QJsonObject{
            {"path", abs(r.hal)},
            {"size", double(QFileInfo(r.hal).size())},
            {"mtime_utc", QFileInfo(r.hal).lastModified().toUTC().toString(Qt::ISODate)}};
        save();
        update(5, "Подготовка белковых запросов");
        auto selectedRecords = records(r.alignment, selection);
        bool needPreparation = false;
        for (auto v : sources) {
            auto s = v.toObject();
            if (!QFileInfo::exists(r.workspace + "/homology/data/" + s["id"].toString() +
                                   "/source.json"))
                needPreparation = true;
        }
        if (needPreparation) {
            event("INFO", "Источники", "sources.download",
                  "Загрузка недостающих протеомов, геномов и аннотаций");
            homology::prepare(out + "/sources.json", r.workspace + "/homology/data", hooks);
        }
        auto q = queries(r, selectedRecords, primarySource, out, [&](QString text) {
            event("WARNING", "Референс", "reference.missing", text);
        });
        QJsonObject missing;
        for (auto it = q.missing.begin(); it != q.missing.end(); ++it)
            missing[it.key()] = it.value();
        manifest["unavailable_reference_genes"] = missing;
        manifest["analyzed_reference_genes"] = array(q.valid);
        save();
        auto halTool = [&](QString name) {
            QString p = QDir(r.halBin).filePath(name);
            if (!QFileInfo::exists(p))
                throw homology::Error("Не найден " + name +
                                      ". Укажите папку HAL-инструментов в настройках.");
            return p;
        };
        auto tree = QString::fromUtf8(
            homology::process(halTool("halStats"), {r.hal, "--tree"}, 120, hooks));
        parseTree(tree);
        homology::writeText(out + "/species_tree.nwk", tree.toUtf8());
        QString hDir = out + "/homology";
        QJsonObject hRun;
        bool reuse = false;
        QString cache = r.workspace + "/homology/results/without_aaeg";
        if (r.useVerified && r.reference == "dmel" && QFileInfo::exists(cache + "/run.json")) {
            try {
                auto original = homology::readJson(cache + "/run.json");
                if (original["status"] != "completed" ||
                    original["config"].toObject()["query_species"] != "dmel")
                    throw homology::Error("Референс кэша не совпадает");
                auto proof = homology::readJson(
                    r.workspace + "/homology/results/verification/trial_verification.json");
                if (proof["source_sha256"] != manifest["input_hashes"].toObject()["alignment"])
                    throw homology::Error("Входные последовательности изменились");
                for (auto filename : {"homology_hits.tsv", "homology_summary.tsv",
                                      "candidate_loci.tsv", "search_status.tsv"})
                    if (homology::sha256(cache + "/" + filename) !=
                        original["output_sha256"].toObject()[filename])
                        throw homology::Error("Нарушена целостность сохранённого результата");
                QMap<QString, QString> originalQueries;
                for (auto s : homology::readFasta(
                         r.workspace + "/homology/data/candidate_queries/proteins.faa", true))
                    originalQueries[s.id] = s.sequence;
                for (auto s : homology::readFasta(q.dir + "/proteins.faa", true))
                    if (originalQueries.value(s.id) != s.sequence)
                        throw homology::Error("Белковый запрос отличается от кэша");
                QDir().mkpath(hDir);
                auto cfg = original["config"].toObject();
                QJsonArray targets;
                QSet<QString> ids;
                for (auto t : cfg["targets"].toArray()) {
                    auto v = t.toObject();
                    if (scope.contains(v["species"].toString())) {
                        targets.append(v);
                        ids.insert(v["id"].toString());
                    }
                }
                for (auto filename : {"homology_hits.tsv", "homology_summary.tsv",
                                      "candidate_loci.tsv", "search_status.tsv"}) {
                    auto all = table(cache + "/" + filename);
                    QVector<Row> filtered;
                    for (auto row : all)
                        if (q.valid.contains(row["gene_id"]) &&
                            (!row.contains("target_id") || ids.contains(row["target_id"])))
                            filtered << row;
                    QFile f(cache + "/" + filename);
                    if (!f.open(QIODevice::ReadOnly))
                        throw homology::Error("Не удалось прочитать таблицу результата");
                    auto header = QString::fromUtf8(f.readLine()).trimmed().split('\t');
                    writeTable(hDir + "/" + filename, header, filtered);
                }
                recomputeSummary(hDir, table(hDir + "/homology_summary.tsv"));
                cfg["targets"] = targets;
                cfg["query_fasta"] = q.dir + "/proteins.faa";
                cfg["query_map"] = q.dir + "/protein_map.tsv";
                hRun = original;
                hRun["config"] = cfg;
                hRun["query_count"] = q.valid.size();
                hRun["query_sha256"] = homology::sha256(q.dir + "/proteins.faa");
                hRun["query_map_sha256"] = homology::sha256(q.dir + "/protein_map.tsv");
                hRun["reused_from"] =
                    QJsonObject{{"path", cache},
                                {"manifest_sha256", homology::sha256(cache + "/run.json")},
                                {"verification", "output hashes, query peptides, original TXT hash "
                                                 "and assembly versions"}};
                QJsonObject hashes;
                for (auto filename : {"homology_hits.tsv", "homology_summary.tsv",
                                      "candidate_loci.tsv", "search_status.tsv"})
                    hashes[filename] = homology::sha256(hDir + "/" + filename);
                hRun["output_sha256"] = hashes;
                saveJson(hDir + "/run.json", hRun);
                reuse = true;
                event("INFO", "Гомология", "cache.verified",
                      "Проверенные результаты гомологии отфильтрованы по выбранным "
                      "генам и видам");
            } catch (const std::exception &e) {
                event("WARNING", "Гомология", "cache.invalid",
                      QString::fromUtf8(e.what()) + ". Будет выполнен новый поиск.");
                if (QDir(hDir).exists()) {
                    hDir = out + "/homology_fresh";
                }
            }
        }
        update(12, "Поиск гомологов");
        if (!reuse) {
            QJsonArray targets;
            for (auto value : sources) {
                auto s = value.toObject();
                QString sp = s["species"].toString(), id = s["id"].toString(),
                        dir = r.workspace + "/homology/data/" + id;
                targets.append(QJsonObject{{"id", id + "_protein"},
                                           {"species", sp},
                                           {"assembly", s["assembly"]},
                                           {"kind", "protein"},
                                           {"fasta", dir + "/proteins.faa"},
                                           {"protein_map", sp == r.reference
                                                               ? q.dir + "/primary_target_map.tsv"
                                                               : dir + "/protein_map.tsv"}});
                if (sp == r.reference)
                    continue;
                QString genome = dir + "/genome.fna";
                if (!QFileInfo::exists(genome) && sp == "dmel") {
                    genome = out + "/dmel_from_hal.fna";
                    event("INFO", "Источники", "hal.export",
                          "Извлечение генома dmel из предоставленного HAL");
                    homology::process(halTool("hal2fasta"), {r.hal, "Drosophila_melanogaster"}, 600,
                                      hooks, genome);
                }
                if (!QFileInfo::exists(genome))
                    throw homology::Error("Геном не подготовлен: " + displaySpecies(sp));
                targets.append(QJsonObject{{"id", id + "_genome"},
                                           {"species", sp},
                                           {"assembly", s["assembly"]},
                                           {"kind", "genome"},
                                           {"fasta", genome}});
            }
            QJsonObject cfg{{"schema_version", 1},
                            {"query_species", r.reference},
                            {"query_fasta", q.dir + "/proteins.faa"},
                            {"query_map", q.dir + "/protein_map.tsv"},
                            {"output_dir", hDir},
                            {"blast_bin", r.blastBin},
                            {"threads", r.threads},
                            {"parallel_targets", 2},
                            {"timeout_seconds", 3600},
                            {"filters", QJsonObject{{"search_evalue", 1e-3},
                                                    {"accept_evalue", 1e-5},
                                                    {"min_query_coverage", 50},
                                                    {"min_identity", 0}}},
                            {"targets", targets}};
            saveJson(out + "/homology_config.json", cfg);
            std::atomic_int count{0};
            homology::Hooks h = hooks;
            h.progress = [&](QString message) {
                if (message.startsWith("Searching")) {
                    int n = ++count;
                    update(12 + qMin(30, n * 30 / qMax(1, targets.size())),
                           "Гомология: " + message);
                }
                event(message.startsWith("Error") ? "ERROR" : "INFO", "Гомология",
                      "homology.progress", message);
            };
            // Serial callback ownership avoids concurrent writes to the event log.
            cfg.remove("parallel_targets");
            saveJson(out + "/homology_config.json", cfg);
            hRun = homology::run(out + "/homology_config.json", h);
        }
        if (hRun["status"] != "completed")
            event("WARNING", "Гомология", "homology.incomplete",
                  "Некоторые поиски завершились с ошибками; уверенная классификация "
                  "запрещена");
        update(45, "Перенос геномных участков");
        QString synDir = out + "/synteny";
        QJsonObject synConfig{
            {"hal_file", abs(r.hal)},
            {"annotation_tsv", q.dir + "/annotation.tsv"},
            {"query_fasta", q.dir + "/proteins.faa"},
            {"query_map", q.dir + "/protein_map.tsv"},
            {"homology_dir", hDir},
            {"sources_manifest", out + "/sources.json"},
            {"assembly_csv", abs(r.assemblies)},
            {"hal_bin", r.halBin},
            {"output_dir", synDir},
            {"reference_species", r.reference},
            {"source_metadata",
             r.workspace + "/homology/data/" + primarySource["id"].toString() + "/source.json"},
            {"primary_flank", 2000},
            {"fallback_flank", 10000},
            {"timeout_seconds", 600}};
        saveJson(out + "/synteny_config.json", {{"schema_version", 1}, {"synteny", synConfig}});
        QJsonObject synRun;
        int synCount = 0;
        homology::Hooks sh = hooks;
        sh.progress = [&](QString m) {
            if (m.startsWith("SYNTENY"))
                update(45 + qMin(20, (++synCount) * 20 / qMax(1, scope.size() - 1)),
                       "Синтения: " + m.section(' ', 1));
            event("INFO", "Синтения", "synteny.progress", m);
        };
        try {
            synRun = synteny::run(out + "/synteny_config.json", sh);
        } catch (const std::exception &e) {
            check();
            synRun = {{"status", "failed"}, {"error", e.what()}};
            event("ERROR", "Синтения", "synteny.failed", QString::fromUtf8(e.what()));
        }
        update(67, "Проверка кодирующего потенциала");
        QString selected = out + "/selected_alignments.txt";
        subset(r.alignment, r.genes, selected, scope);
        orf::Settings settings;
        settings.input = selected.toStdString();
        settings.outputRoot = (out + "/orf").toStdString();
        settings.finalOut = (out + "/orf/gene_profiles.jsonl").toStdString();
        settings.referenceSpecies = r.reference.toStdString();
        settings.maxEvents = 200;
        settings.progressEvery = 500;
        int orfStage = 0, recordTotal = 0, comparisonTotal = 0;
        for (const auto &record : selectedRecords)
            if (scope.contains(record.species)) {
                ++recordTotal;
                if (record.section != "raw" && record.species != r.reference)
                    ++comparisonTotal;
            }
        QRegularExpression recordsProgress("Обработано записей: ([0-9]+)");
        QJsonObject orfRun;
        if (comparisonTotal == 0) {
            QDir().mkpath(out + "/orf");
            orfRun = {{"status", "unavailable"},
                      {"reason", "Во входном TXT нет сравнительных CDS для выбранных видов."}};
            saveJson(out + "/orf/pipeline_run.json", orfRun);
            event("WARNING", "ORF", "orf.unavailable", orfRun.value("reason").toString());
        } else
            try {
                orfRun = json(orf::runPipeline(
                    settings,
                    [&](const std::string &message) {
                        QString m = QString::fromStdString(message);
                        if (m.contains(": вычисление")) {
                            ++orfStage;
                            update(67 + (orfStage - 1) * 22 / 5,
                                   "ORF: этап " + QString::number(orfStage) + " из 5");
                        }
                        auto match = recordsProgress.match(m);
                        if (match.hasMatch() && orfStage > 0) {
                            int done = match.captured(1).toInt(),
                                total = orfStage <= 2 ? recordTotal : comparisonTotal;
                            int value =
                                67 + (orfStage - 1) * 22 / 5 + qMin(4, done * 4 / qMax(1, total));
                            update(qMin(89, value), QString("ORF: этап %1 из 5 · %2 / %3 записей")
                                                        .arg(orfStage)
                                                        .arg(done)
                                                        .arg(total));
                        }
                        event("INFO", "ORF", "orf.progress", m);
                    },
                    hooks.cancelled));
            } catch (const std::exception &error) {
                check();
                orfRun = {{"status", "failed"}, {"error", error.what()}};
                event("ERROR", "ORF", "orf.failed", QString::fromUtf8(error.what()));
            }
        update(90, "Объединение свидетельств и классификация");
        QMap<QString, QJsonObject> profiles;
        QFile pf(out + "/orf/gene_profiles.jsonl");
        if (orfRun.value("status") == "completed" && !pf.open(QIODevice::ReadOnly))
            throw homology::Error("Не удалось прочитать профили ORF");
        while (pf.isOpen() && !pf.atEnd()) {
            auto o = QJsonDocument::fromJson(pf.readLine()).object();
            if (!o.isEmpty())
                profiles[o["gene_id"].toString()] = o;
        }
        auto hits = table(hDir + "/homology_hits.tsv"), states = table(hDir + "/search_status.tsv");
        QVector<Row> synRows, evidenceRows;
        if (QFileInfo::exists(synDir + "/synteny_results.tsv"))
            synRows = table(synDir + "/synteny_results.tsv");
        if (QFileInfo::exists(synDir + "/synteny_evidence.tsv"))
            evidenceRows = table(synDir + "/synteny_evidence.tsv");
        QJsonArray results;
        QMap<QString, int> queryLengths;
        for (auto row : table(q.dir + "/annotation.tsv"))
            queryLengths[row["gene_id"]] = row["length"].toInt();
        QVector<Row> summary;
        QString primaryGenus = displaySpecies(r.reference).section(' ', 0, 0);
        QSet<QString> basal{"dper", "dpse", "dmir", "dwil", "sleb", "bdor", "blat",
                            "mdom", "agam", "bmor", "mpha", "amel", "aflo", "dmag"};
        for (auto gene : r.genes) {
            check();
            QSet<QString> cross, distant, strongOwn, preserved, broken, contexts, complex;
            bool complete = q.valid.contains(gene),
                 shortQuery = queryLengths.value(gene) > 0 && queryLengths.value(gene) < 30;
            QJsonObject speciesEvidence;
            auto speciesItem = [&](QString sp) {
                auto o = speciesEvidence.value(sp).toObject();
                o["species"] = sp;
                o["name"] = displaySpecies(sp);
                return o;
            };
            for (auto sp : scope) {
                auto o = speciesItem(sp);
                o["homology"] = "нет принятых совпадений";
                o["synteny"] = "не оценено";
                o["orf"] = "не оценено";
                o["primary"] = sp == r.reference;
                speciesEvidence[sp] = o;
            }
            for (auto state : states)
                if (state["gene_id"] == gene && !state["status"].startsWith("completed_")) {
                    complete = false;
                    auto o = speciesItem(state["target_species"]);
                    o["homology"] = "ошибка или поиск не выполнен";
                    speciesEvidence[state["target_species"]] = o;
                }
            int significantPartial = 0;
            for (auto h : hits)
                if (h["gene_id"] == gene) {
                    QString sp = h["target_species"];
                    auto o = speciesItem(sp);
                    if (h["quality"] == "significant_partial_similarity") {
                        ++significantPartial;
                        if (o["homology"] == "нет принятых совпадений")
                            o["homology"] = "частичное сходство";
                    }
                    if (h["query_length"].toInt() < 30)
                        shortQuery = true;
                    if (h["passes_filters"] == "true") {
                        o["homology"] = "принятое совпадение";
                        if (o.value("target_label").toString().isEmpty() ||
                            h["search_type"] == "cross_species_blastp" &&
                                h["bitscore"].toDouble() >
                                    o.value("best_protein_bitscore").toDouble()) {
                            o["target_label"] = h["target_gene_id"].isEmpty()
                                                    ? h["target_sequence_id"]
                                                    : h["target_gene_id"];
                            if (h["search_type"] == "cross_species_blastp")
                                o["best_protein_bitscore"] = h["bitscore"].toDouble();
                        }
                        o["hsp_count"] = o.value("hsp_count").toInt() + 1;
                        if (h["search_type"] != "within_species_blastp")
                            cross.insert(sp);
                        bool strong = h["evalue"].toDouble() <= 1e-10 &&
                                      h["query_coverage"].toDouble() >= 80 &&
                                      h["alignment_length"].toInt() >= 50 &&
                                      h["identity"].toDouble() >= 30;
                        if (h["search_type"] == "within_species_blastp" && strong &&
                            h["identity"].toDouble() >= 40)
                            strongOwn.insert(h["target_gene_id"]);
                        if (h["search_type"] == "cross_species_blastp" && strong &&
                            (displaySpecies(sp).section(' ', 0, 0) != primaryGenus ||
                             (r.reference == "dmel" || r.reference == "dsim" ||
                              r.reference == "dsec") &&
                                 basal.contains(sp)))
                            distant.insert(sp);
                    }
                    speciesEvidence[sp] = o;
                }
            bool synAvailable = false, multi = false;
            for (auto row : synRows)
                if (row["gene_id"] == gene) {
                    auto sp = row["species"];
                    auto o = speciesItem(sp);
                    o["synteny"] = row["status"];
                    if (row["status"] != "unavailable_in_hal" && row["status"] != "error")
                        synAvailable = true;
                    if (row["context_groups"].toInt() > 0)
                        contexts.insert(sp);
                    if (row["mapping_groups"].toInt() > 1 || row["status"].contains("complex")) {
                        complex.insert(sp);
                        multi = true;
                    }
                    o["cds_coverage"] =
                        row["cds_bases"].toDouble() > 0
                            ? row["cds_covered"].toDouble() / row["cds_bases"].toDouble()
                            : 0;
                    speciesEvidence[sp] = o;
                }
            auto profile = profiles.value(gene);
            for (auto v : profile.value("species_profiles").toArray()) {
                auto p = v.toObject();
                QString sp = p["species"].toString(), status = p["status"].toString();
                if (!scope.contains(sp))
                    continue;
                auto o = speciesItem(sp);
                o["orf"] = status;
                o["loci"] = p["loci"];
                if (status.startsWith("preserved_both"))
                    preserved.insert(sp);
                if (status == "predicted_interruption_both")
                    broken.insert(sp);
                if (status == "multiple_loci_uncertain" || status == "multiple_predictions")
                    multi = true;
                speciesEvidence[sp] = o;
            }
            QJsonObject ev{{"primary_valid", q.valid.contains(gene)},
                           {"homology_complete", complete},
                           {"synteny_available", synAvailable},
                           {"orf_available", !preserved.isEmpty() || !broken.isEmpty()},
                           {"strong_paralogs", strongOwn.size()},
                           {"strong_paralog_ids", array(strongOwn)},
                           {"accepted_cross_species", array(cross)},
                           {"strong_distant_species", array(distant)},
                           {"preserved_species", array(preserved)},
                           {"interrupted_species", array(broken)},
                           {"context_species", array(contexts)},
                           {"complex_species", array(complex)},
                           {"short_query", shortQuery},
                           {"query_length_aa", queryLengths.value(gene)},
                           {"multiple_loci", multi},
                           {"partial_hsp_count", significantPartial},
                           {"species", speciesEvidence},
                           {"reference_species", r.reference},
                           {"reference_gene_id", q.native.value(gene)},
                           {"unavailable_reason", q.missing.value(gene)}};
            auto verdict = classify(gene, ev);
            verdict["orf_profile"] = profile;
            results.append(verdict);
            summary << Row{{"gene_id", gene},
                           {"reference_gene_id", q.native.value(gene)},
                           {"reference_species", r.reference},
                           {"class", verdict["class"].toString()},
                           {"class_label", verdict["class_label"].toString()},
                           {"confidence_level", verdict["confidence_level"].toString()},
                           {"support_score", QString::number(verdict["support_score"].toInt())},
                           {"reason", verdict["reason"].toString()},
                           {"score_kind", "heuristic_not_probability"}};
            event("INFO", "Классификация", "classifier.assigned",
                  gene + " → " + verdict["class_label"].toString() + "; поддержка " +
                      QString::number(verdict["support_score"].toInt()));
        }
        update(96, "Сохранение результатов");
        homology::writeText(out + "/summary.json", QJsonDocument(results).toJson());
        writeTable(out + "/summary.tsv",
                   {"gene_id", "reference_gene_id", "reference_species", "class", "class_label",
                    "confidence_level", "support_score", "reason", "score_kind"},
                   summary);
        QString readme = "# LibreGenesis — результат анализа\n\nОсновной организм: " +
                         displaySpecies(r.reference) +
                         ".\n\nsummary.tsv — ген, гипотеза класса, уровень "
                         "доверия, эвристическая поддержка и обоснование. Баллы не "
                         "являются вероятностью происхождения.\n\nГомология: " +
                         QFileInfo(hDir).fileName() +
                         "/; синтения: synteny/; ORF: orf/. Полный журнал: "
                         "events.jsonl.\n\nAedes aegypti исключён из области "
                         "анализа. Отсутствующие в HAL виды не трактуются как "
                         "потеря гена. ORF-предсказания и геномные участки не "
                         "объявляются одним локусом без подтверждения координат.\n";
        homology::writeText(out + "/README.md", readme.toUtf8());
        manifest["status"] = "completed";
        manifest["finished_utc"] = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);
        manifest["modules"] = QJsonObject{{"homology", hRun["status"]},
                                          {"homology_reused", reuse},
                                          {"synteny", synRun["status"]},
                                          {"orf", orfRun["status"]}};
        QJsonObject outputHashes;
        for (QString filename : {"summary.json", "summary.tsv", "species_tree.nwk"})
            outputHashes[filename] = homology::sha256(out + "/" + filename);
        manifest["output_sha256"] = outputHashes;
        manifest["completion_kind"] = synRun.value("status") == "completed" &&
                                              orfRun.value("status") == "completed" &&
                                              hRun.value("status") == "completed"
                                          ? "complete"
                                          : "with_limitations";
        manifest["result_sha256"] = homology::sha256(out + "/summary.json");
        manifest["gene_count"] = results.size();
        save();
        event("INFO", "Готово", "run.completed", "Анализ завершён. Результаты: " + out);
        update(100, "Анализ завершён");
        return {out, tree, results, r.reference, manifest};
    } catch (const std::exception &e) {
        bool cancelled = hooks.cancelled && hooks.cancelled->load();
        manifest["status"] = cancelled ? "cancelled" : "failed";
        manifest["error"] = e.what();
        manifest["finished_utc"] = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);
        save();
        event(cancelled ? "WARNING" : "ERROR", "Остановка",
              cancelled ? "run.cancelled" : "run.failed", QString::fromUtf8(e.what()));
        throw homology::Error(
            (cancelled ? QString("Анализ отменён. ") : QString::fromUtf8(e.what()) + "\n") +
            "Журнал и промежуточные результаты: " + out);
    }
}
} // namespace go
