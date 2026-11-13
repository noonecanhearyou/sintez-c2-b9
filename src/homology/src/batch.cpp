#include "homology/homology.hpp"
#include <QDateTime>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QSet>
#include <QTextStream>
#include <QThread>
#include <algorithm>
#include <future>
#include <mutex>
#include <vector>
namespace homology {
namespace {
using Row = QMap<QString, QString>;
struct Table {
    QString header;
    QVector<Row> rows;
};
Table readTable(const QString &path) {
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly))
        throw Error("Cannot read " + path);
    QTextStream in(&f);
    Table t;
    t.header = in.readLine();
    auto names = t.header.split('\t');
    while (!in.atEnd()) {
        auto line = in.readLine();
        if (line.isEmpty())
            continue;
        auto values = line.split('\t');
        if (values.size() != names.size())
            throw Error("Malformed table " + path);
        Row row;
        for (int i = 0; i < names.size(); ++i)
            row[names[i]] = values[i];
        t.rows.append(row);
    }
    return t;
}
QString absolute(const QDir &base, const QString &path) {
    return QFileInfo(path).isAbsolute() ? path : base.absoluteFilePath(path);
}
QString sortedList(const QSet<QString> &values) {
    auto list = values.values();
    list.sort();
    return list.join(',');
}
void addList(QSet<QString> &dest, const QString &value, QChar separator = ',') {
    for (auto x : value.split(separator, Qt::SkipEmptyParts))
        dest.insert(x);
}
} // namespace
QJsonObject mergeRuns(const QString &configPath, const QStringList &runs,
                      const QString &outputDir) {
    auto config = readJson(configPath);
    QDir base(QFileInfo(configPath).absolutePath());
    auto qhash = sha256(absolute(base, config.value("query_fasta").toString()));
    auto mhash = sha256(absolute(base, config.value("query_map").toString()));
    if (QFileInfo::exists(outputDir + "/run.json"))
        throw Error("Merged output already exists");
    QDir().mkpath(outputDir);
    QMap<QString, QString> headers, texts;
    QMap<QString, Row> summary;
    QMap<QString, QSet<QString>> lists[7];
    QStringList listNames{"paralog_gene_ids",
                          "protein_hit_species",
                          "genome_hit_species",
                          "partial_similarity_gene_ids",
                          "partial_protein_species",
                          "partial_genome_species",
                          "limitations"};
    QMap<QString, double> proteinScore, genomeScore;
    QSet<QString> targetIds, statusKeys;
    QJsonArray children, commands;
    int completed = 0, errors = 0, notRun = 0;
    QJsonObject versions;
    QString started;
    for (auto dir : runs) {
        auto manifest = readJson(dir + "/run.json");
        auto state = manifest.value("status").toString();
        QString time = manifest.value("started_utc").toString();
        if (time.isEmpty())
            for (auto nested : manifest.value("children").toArray()) {
                auto nestedTime = readJson(nested.toObject().value("path").toString())
                                      .value("started_utc")
                                      .toString();
                if (!nestedTime.isEmpty() && (time.isEmpty() || nestedTime < time))
                    time = nestedTime;
            }
        if (!time.isEmpty() && (started.isEmpty() || time < started))
            started = time;
        if (state != "completed" && state != "partial_error" && state != "partial_not_run" &&
            state != "error")
            throw Error("Incomplete child run " + dir);
        if (manifest.value("query_sha256").toString() != qhash ||
            manifest.value("query_map_sha256").toString() != mhash ||
            manifest.value("config").toObject().value("filters") != config.value("filters"))
            throw Error("Incompatible child inputs or filters");
        auto childVersions = manifest.value("versions").toObject();
        if (versions.isEmpty())
            versions = childVersions;
        else if (versions != childVersions)
            throw Error("BLAST version mismatch");
        for (auto v : manifest.value("config").toObject().value("targets").toArray()) {
            auto id = v.toObject().value("id").toString();
            if (targetIds.contains(id))
                throw Error("Duplicate target across child runs " + id);
            targetIds.insert(id);
        }
        children.append(QJsonObject{{"path", QFileInfo(dir + "/run.json").absoluteFilePath()},
                                    {"sha256", sha256(dir + "/run.json")}});
        for (auto v : manifest.value("commands_and_inputs").toArray())
            commands.append(v);
        for (QString name : {"homology_hits.tsv", "candidate_loci.tsv", "search_status.tsv"}) {
            auto t = readTable(dir + "/" + name);
            if (!headers.contains(name)) {
                headers[name] = t.header;
                texts[name] = t.header + '\n';
            } else if (headers[name] != t.header)
                throw Error("Table schema mismatch");
            auto columns = t.header.split('\t');
            for (auto row : t.rows) {
                QStringList cells;
                for (auto column : columns)
                    cells << row[column];
                texts[name] += cells.join('\t') + '\n';
                if (name == "search_status.tsv") {
                    QString key = row["query_id"] + '\t' + row["target_id"];
                    if (statusKeys.contains(key))
                        throw Error("Duplicate search status");
                    statusKeys.insert(key);
                }
                if (name == "homology_hits.tsv" && row["passes_filters"] == "true") {
                    auto gene = row["gene_id"];
                    bool genome = row["search_type"] == "tblastn";
                    auto &scores = genome ? genomeScore : proteinScore;
                    double score = row["bitscore"].toDouble();
                    if (!scores.contains(gene) || score > scores[gene]) {
                        scores[gene] = score;
                        auto &dest = summary[gene];
                        QString p = genome ? "best_genome_" : "best_protein_";
                        dest[p + "target"] = row["target_id"];
                        dest[p + (genome ? "contig" : "sequence_id")] = row["target_sequence_id"];
                        dest[p + "evalue"] = row["evalue"];
                    }
                }
            }
        }
        auto t = readTable(dir + "/homology_summary.tsv");
        if (!headers.contains("homology_summary.tsv"))
            headers["homology_summary.tsv"] = t.header;
        else if (headers["homology_summary.tsv"] != t.header)
            throw Error("Summary schema mismatch");
        bool first = true;
        for (auto row : t.rows) {
            QString gene = row["gene_id"];
            auto &dest = summary[gene];
            dest["gene_id"] = gene;
            for (int i = 0; i < listNames.size(); ++i)
                addList(lists[i][gene], row[listNames[i]], i == 6 ? ';' : ',');
            for (QString field : {"completed_targets", "error_targets", "not_run_targets"})
                dest[field] = QString::number(dest[field].toInt() + row[field].toInt());
            if (first) {
                completed += row["completed_targets"].toInt();
                errors += row["error_targets"].toInt();
                notRun += row["not_run_targets"].toInt();
                first = false;
            }
        }
    }
    QSet<QString> expectedTargets;
    for (auto v : config.value("targets").toArray())
        expectedTargets.insert(v.toObject().value("id").toString());
    if (targetIds != expectedTargets)
        throw Error("Merged target scope does not match configuration");
    auto queryCount =
        readFasta(absolute(base, config.value("query_fasta").toString()), true).size();
    if (statusKeys.size() != queryCount * targetIds.size())
        throw Error("Missing query/target states");
    QString state = errors ? (completed ? "partial_error" : "error")
                           : (notRun ? "partial_not_run" : "completed");
    auto columns = headers["homology_summary.tsv"].split('\t');
    QString mergedSummary = headers["homology_summary.tsv"] + '\n';
    for (auto it = summary.begin(); it != summary.end(); ++it) {
        auto &row = it.value();
        QString gene = it.key();
        for (int i = 0; i < listNames.size(); ++i) {
            row[listNames[i]] = sortedList(lists[i][gene]);
            if (i == 6)
                row[listNames[i]].replace(',', ';');
        }
        row["paralog_candidate_count"] = QString::number(lists[0][gene].size());
        row["search_status"] = state;
        QStringList cells;
        for (auto col : columns)
            cells << row[col];
        mergedSummary += cells.join('\t') + '\n';
    }
    for (auto it = texts.begin(); it != texts.end(); ++it)
        writeText(outputDir + "/" + it.key(), it.value().toUtf8());
    writeText(outputDir + "/homology_summary.tsv", mergedSummary.toUtf8());
    QJsonObject outputHashes;
    for (QString name :
         {"homology_hits.tsv", "homology_summary.tsv", "candidate_loci.tsv", "search_status.tsv"})
        outputHashes[name] = sha256(outputDir + "/" + name);
    QJsonObject manifest{
        {"schema_version", 1},
        {"started_utc", started},
        {"configuration_sha256", sha256(configPath)},
        {"output_sha256", outputHashes},
        {"module_version", "0.2.0"},
        {"status", state},
        {"config", config},
        {"query_sha256", qhash},
        {"query_map_sha256", mhash},
        {"versions", versions},
        {"children", children},
        {"commands_and_inputs", commands},
        {"finished_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)},
        {"target_count", targetIds.size()},
        {"query_count", queryCount},
        {"coordinate_system", "1-based inclusive; genome coordinates on full input contigs"}};
    writeText(outputDir + "/run.json", QJsonDocument(manifest).toJson());
    return manifest;
}
QJsonObject runBatch(const QString &configPath, const Hooks &hooks) {
    auto config = readJson(configPath);
    QDir base(QFileInfo(configPath).absolutePath());
    QString out = absolute(base, config.value("output_dir").toString());
    if (config.value("output_dir").toString().isEmpty())
        throw Error("Missing output_dir");
    bool resume = config.value("resume_batch").toBool(false);
    if (QFileInfo::exists(out) && !resume)
        throw Error("Batch output already exists");
    if (QFileInfo::exists(out + "/run.json"))
        throw Error("Batch already merged; choose a new output_dir");
    QDir().mkpath(out + "/configs");
    QDir().mkpath(out + "/targets");
    auto targets = config.value("targets").toArray();
    int jobs = config.value("parallel_targets").toInt(3);
    if (jobs < 1 || jobs > 8 || targets.isEmpty())
        throw Error("Invalid parallel_targets/targets");
    QSet<QString> ids;
    QStringList paths;
    for (auto v : targets) {
        auto id = v.toObject().value("id").toString();
        if (!QRegularExpression("^[A-Za-z0-9_-]+$").match(id).hasMatch() || ids.contains(id))
            throw Error("Invalid target ID");
        ids.insert(id);
        paths << out + "/targets/" + id;
    }
    std::atomic_int next{0};
    std::mutex mutex;
    QStringList errors;
    QJsonArray statuses;
    auto saveStatus = [&] {
        writeText(out + "/batch_status.json",
                  QJsonDocument(QJsonObject{{"targets", statuses}}).toJson());
    };
    Hooks safe = hooks;
    safe.progress = [&](const QString &message) {
        std::lock_guard<std::mutex> lock(mutex);
        if (hooks.progress)
            hooks.progress(message);
    };
    auto worker = [&] {
        for (;;) {
            int index = next.fetch_add(1);
            if (index >= targets.size())
                return;
            auto target = targets[index].toObject();
            QString id = target.value("id").toString(), dir = paths[index];
            try {
                if (hooks.cancelled && hooks.cancelled->load())
                    throw Error("cancelled");
                QString readyFile = target.value("source_ready_file").toString();
                if (!readyFile.isEmpty()) {
                    readyFile = absolute(base, readyFile);
                    QElapsedTimer timer;
                    timer.start();
                    int waitSeconds = config.value("wait_for_sources_seconds").toInt(0);
                    if (!QFileInfo::exists(readyFile))
                        safe.progress("Waiting for prepared source: " + id);
                    while (!QFileInfo::exists(readyFile)) {
                        if (hooks.cancelled && hooks.cancelled->load())
                            throw Error("cancelled");
                        if (timer.elapsed() >= qint64(waitSeconds) * 1000)
                            throw Error("Source not prepared: " + id);
                        QThread::msleep(200);
                    }
                }
                auto child = config;
                child.remove("parallel_targets");
                child.remove("resume_batch");
                for (QString field : {"query_fasta", "query_map", "blast_bin"})
                    if (!child.value(field).toString().isEmpty())
                        child[field] = absolute(base, child.value(field).toString());
                for (QString field : {"fasta", "protein_map"})
                    if (!target.value(field).toString().isEmpty())
                        target[field] = absolute(base, target.value(field).toString());
                child["targets"] = QJsonArray{target};
                child["output_dir"] = dir;
                QString childPath = out + "/configs/" + id + ".json";
                writeText(childPath, QJsonDocument(child).toJson());
                if (resume && QFileInfo::exists(dir + "/run.json")) {
                    auto old = readJson(dir + "/run.json");
                    bool compatible = (old.value("status") == "completed" ||
                                       old.value("status") == "partial_not_run") &&
                                      old.value("config").toObject() == child &&
                                      old.value("query_sha256").toString() ==
                                          sha256(child.value("query_fasta").toString()) &&
                                      old.value("query_map_sha256").toString() ==
                                          sha256(child.value("query_map").toString());
                    bool targetHash = !target.value("enabled").toBool(true);
                    bool mapHash = !target.value("enabled").toBool(true) ||
                                   target.value("protein_map").toString().isEmpty();
                    if (compatible)
                        for (auto v : old.value("commands_and_inputs").toArray()) {
                            auto record = v.toObject();
                            if (record.value("protein_map").toString() ==
                                    target.value("protein_map").toString() &&
                                !target.value("protein_map").toString().isEmpty() &&
                                record.value("sha256").toString() ==
                                    sha256(target.value("protein_map").toString()))
                                mapHash = true;
                            auto info = record.value("input").toObject();
                            if (info.value("id").toString() == id &&
                                info.value("fasta_sha256").toString() ==
                                    sha256(target.value("fasta").toString()))
                                targetHash = true;
                        }
                    if (!compatible || !targetHash || !mapHash)
                        throw Error("Existing child cannot be safely resumed: " + id);
                    safe.progress("Reusing completed target " + id);
                } else {
                    run(childPath, safe);
                }
                std::lock_guard<std::mutex> lock(mutex);
                statuses.append(
                    QJsonObject{{"id", id},
                                {"status", "finished"},
                                {"run_status", readJson(dir + "/run.json").value("status")}});
                saveStatus();
            } catch (const std::exception &error) {
                std::lock_guard<std::mutex> lock(mutex);
                errors << id + ": " + QString::fromUtf8(error.what());
                statuses.append(
                    QJsonObject{{"id", id}, {"status", "error"}, {"message", error.what()}});
                saveStatus();
            }
        }
    };
    std::vector<std::future<void>> workers;
    for (int i = 0; i < std::min(jobs, int(targets.size())); ++i)
        workers.push_back(std::async(std::launch::async, worker));
    for (auto &worker : workers)
        worker.get();
    if (!errors.isEmpty())
        throw Error(errors.join('\n'));
    return mergeRuns(configPath, paths, out);
}
} // namespace homology
