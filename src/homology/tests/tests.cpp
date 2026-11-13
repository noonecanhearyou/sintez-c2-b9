#include "homology/homology.hpp"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <QTextStream>
using namespace homology;
static void check(bool condition, const char *msg) {
    if (!condition)
        throw Error(msg);
}
static QString contents(QString path) {
    QFile f(path);
    check(f.open(QIODevice::ReadOnly), "cannot read test output");
    return QString::fromUtf8(f.readAll());
}
static QJsonObject target(QString id, QString species, QString kind, QString fasta,
                          QString map = {}) {
    return {{"id", id},     {"species", species}, {"assembly", "fixture_v1"},
            {"kind", kind}, {"fasta", fasta},     {"protein_map", map}};
}
static QString reverseComplement(QString dna) {
    QString r;
    for (auto it = dna.crbegin(); it != dna.crend(); ++it) {
        int i = QString("ACGT").indexOf(*it);
        check(i >= 0, "bad fixture DNA");
        r += QString("TGCA")[i];
    }
    return r;
}
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    QTextStream out(stdout);
    try {
        check(coveredLength({{1, 10}, {5, 15}, {20, 22}}) == 18,
              "union coverage double counts overlaps");
        check(coveredLength({{1, 3}, {4, 5}}) == 5, "adjacent interval union");
        QTemporaryDir tmp;
        check(tmp.isValid(), "temporary directory");
        QString d = tmp.path();
        writeText(d + "/bad.faa", ">same\nMPEPTIDE\n>same\nMPEPTIDE\n");
        bool rejected = false;
        try {
            readFasta(d + "/bad.faa", true);
        } catch (const Error &) {
            rejected = true;
        }
        check(rejected, "duplicate FASTA IDs accepted");
        writeText(d + "/bad.faa", ">id\nMP*AAA\n");
        rejected = false;
        try {
            readFasta(d + "/bad.faa", true);
        } catch (const Error &) {
            rejected = true;
        }
        check(rejected, "internal stop accepted");
        QVector<Protein> proteins{{"z", "g1", "t2", "sp", "MPEPTIDE", {}},
                                  {"a", "g1", "t1", "sp", "MPEPTIDE", {}},
                                  {"b", "g2", "t3", "sp", "MPEPTIDELONG", {}}};
        writeProteins(d + "/proteome", proteins);
        extractQueries(d + "/proteome", {"g1"}, d + "/selected");
        check(readFasta(d + "/selected/proteins.faa", true)[0].id == "a",
              "isoform selection not deterministic");
        std::atomic_bool cancelled{true};
        Hooks h;
        h.cancelled = &cancelled;
        rejected = false;
        try {
            process("unused", {}, 1, h);
        } catch (const Error &e) {
            rejected = QString(e.what()) == "cancelled";
        }
        check(rejected, "cancellation ignored");
        QString bin = qEnvironmentVariable("HOMOLOGY_BLAST_BIN");
        if (bin.isEmpty())
            bin =
                QDir(QCoreApplication::applicationDirPath()).absoluteFilePath("../tools/blast/bin");
        check(QFile::exists(bin + "/blastp"),
              "BLAST required for integration tests; set HOMOLOGY_BLAST_BIN");
        QString alphabet = "ACDEFGHIKLMNPQRSTVWY", query, negative, dna;
        unsigned state = 1234567;
        QStringList codons{"GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG",
                           "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"};
        for (int i = 0; i < 160; ++i) {
            state = state * 1664525u + 1013904223u;
            int j = int((state >> 8) % 20);
            query += alphabet[j];
            dna += codons[j];
            state = state * 1664525u + 1013904223u;
            negative += alphabet[int((state >> 8) % 20)];
        }
        writeText(d + "/queries.faa",
                  (">q1\n" + query + "\n>qnegative\n" + negative + '\n').toUtf8());
        writeText(d + "/query_map.tsv", "protein_id\tgene_id\nq1\tg1\nqnegative\tgnegative\n");
        writeText(d + "/within.faa",
                  (">q1iso\n" + query + "\n>paralogA\n" + query + "\n>paralogB\n" + query +
                   "\n>partial\n" + query.mid(20, 60) + '\n')
                      .toUtf8());
        writeText(d + "/within.tsv",
                  "protein_id\tgene_id\nq1iso\tg1\nparalogA\tg2\nparalogB\tg2\npartial\tg3\n");
        writeText(d + "/cross.faa", (">homolog\n" + query + '\n').toUtf8());
        writeText(d + "/genome.fna", (">minus_contig\n" + QString(60, 'A') +
                                      reverseComplement(dna) + QString(60, 'A') + '\n')
                                         .toUtf8());
        QJsonObject cfg{
            {"schema_version", 1},
            {"query_species", "sp1"},
            {"query_fasta", "queries.faa"},
            {"query_map", "query_map.tsv"},
            {"blast_bin", bin},
            {"output_dir", "out"},
            {"threads", 1},
            {"filters", QJsonObject{{"search_evalue", 1e-3},
                                    {"accept_evalue", 1e-5},
                                    {"min_query_coverage", 50}}},
            {"targets", QJsonArray{target("within", "sp1", "protein", "within.faa", "within.tsv"),
                                   target("cross", "sp2", "protein", "cross.faa"),
                                   target("genome", "sp2", "genome", "genome.fna")}}};
        auto writeConfig = [&] { writeText(d + "/config.json", QJsonDocument(cfg).toJson()); };
        writeConfig();
        auto result = run(d + "/config.json");
        check(result["status"] == "completed", "integration run incomplete");
        auto hits = contents(d + "/out/homology_hits.tsv");
        check(hits.contains("same_gene_or_isoform"), "same-gene isoform not excluded");
        check(hits.contains("significant_partial_similarity"), "partial significant evidence lost");
        auto summary = contents(d + "/out/homology_summary.tsv");
        check(summary.contains("g1\t1\tg2\tsp2\tsp2\tcompleted"),
              "gene-level paralogs not deduplicated");
        check(summary.contains("gnegative\t0\t\t\t\tcompleted"), "negative control summary wrong");
        auto statuses = contents(d + "/out/search_status.tsv");
        check(statuses.contains("completed_no_accepted_hits"), "no-hit status missing");
        auto loci = contents(d + "/out/candidate_loci.tsv");
        check(loci.contains("minus_contig\t61\t540\t60\t540\t-"),
              "reverse-strand or BED coordinates incorrect");
        rejected = false;
        try {
            run(d + "/config.json");
        } catch (const Error &) {
            rejected = true;
        }
        check(rejected, "existing result overwritten");
        auto targets = cfg["targets"].toArray();
        auto disabled = target("disabled", "sp3", "protein", "");
        disabled["enabled"] = false;
        targets.append(disabled);
        cfg["targets"] = targets;
        cfg["output_dir"] = "disabled_out";
        writeConfig();
        check(run(d + "/config.json")["status"] == "partial_not_run", "disabled target hidden");
        check(contents(d + "/disabled_out/search_status.tsv").contains("not_run"),
              "not_run status missing");
        auto broken = target("broken", "sp4", "genome", "genome.fna");
        writeText(d + "/broken.fna", "not a fasta file\n");
        broken["fasta"] = "broken.fna";
        targets.append(broken);
        cfg["targets"] = targets;
        cfg["output_dir"] = "partial_out";
        writeConfig();
        check(run(d + "/config.json")["status"] == "partial_error",
              "target failure not propagated");
        check(contents(d + "/partial_out/search_status.tsv").contains("\terror\t"),
              "error status missing");
        // Cached specifications from the old downloader contain optional null keys.
        QDir().mkpath(d + "/cache/cached");
        writeText(d + "/cache/cached/file.txt", "cached biological data");
        QJsonObject source{{"id", "cached"},
                           {"provider", "ncbi"},
                           {"species", "sp"},
                           {"assembly", "GCF_000000001.1"}};
        auto oldSource = source;
        oldSource["base_url"] = QJsonValue::Null;
        writeText(
            d + "/cache/cached/source.json",
            QJsonDocument(
                QJsonObject{{"specification", oldSource},
                            {"files", QJsonArray{QJsonObject{
                                          {"name", "file.txt"},
                                          {"sha256", sha256(d + "/cache/cached/file.txt")}}}}})
                .toJson());
        writeText(d + "/sources.json",
                  QJsonDocument(QJsonObject{{"schema_version", 1}, {"sources", QJsonArray{source}}})
                      .toJson());
        prepare(d + "/sources.json", d + "/cache");
        check(contents(d + "/cache/preparation_status.json").contains("ready"),
              "old cache with null optional fields rejected");
        auto batchTargets = cfg["targets"].toArray();
        while (batchTargets.size() > 3)
            batchTargets.removeLast();
        cfg["targets"] = batchTargets;
        cfg["output_dir"] = "batch_out";
        cfg["parallel_targets"] = 2;
        cfg["mt_mode"] = 1;
        writeConfig();
        check(runBatch(d + "/config.json")["status"] == "completed", "parallel batch failed");
        check(contents(d + "/batch_out/homology_summary.tsv")
                  .contains("g1\t1\tg2\tsp2\tsp2\tcompleted\t3\t0\t0"),
              "batch summary differs from serial search");
        auto manifest = readJson(d + "/batch_out/run.json");
        check(manifest["query_count"].toInt() == 2 && manifest["target_count"].toInt() == 3,
              "batch lost queries or targets");
        rejected = false;
        try {
            mergeRuns(d + "/config.json",
                      {d + "/batch_out/targets/within", d + "/batch_out/targets/within"},
                      d + "/bad_merge");
        } catch (const Error &) {
            rejected = true;
        }
        check(rejected, "merge accepts duplicate targets");
        auto waitingTarget = batchTargets[0].toObject();
        waitingTarget["source_ready_file"] = "missing_source.json";
        cfg["targets"] = QJsonArray{waitingTarget};
        cfg["parallel_targets"] = 1;
        cfg["wait_for_sources_seconds"] = 0;
        cfg["output_dir"] = "waiting_out";
        writeConfig();
        rejected = false;
        try {
            runBatch(d + "/config.json");
        } catch (const Error &) {
            rejected = true;
        }
        check(rejected &&
                  contents(d + "/waiting_out/batch_status.json").contains("Source not prepared"),
              "unprepared source not reported as an error");
        QDir().mkpath(d + "/cache/bad_cached");
        auto badSource = source;
        badSource["id"] = "bad_cached";
        writeText(d + "/cache/bad_cached/file.txt", "changed data");
        writeText(
            d + "/cache/bad_cached/source.json",
            QJsonDocument(QJsonObject{{"specification", badSource},
                                      {"files", QJsonArray{QJsonObject{{"name", "file.txt"},
                                                                       {"sha256", "incorrect"}}}}})
                .toJson());
        writeText(d + "/sources.json",
                  QJsonDocument(QJsonObject{{"schema_version", 1},
                                            {"sources", QJsonArray{badSource, source}}})
                      .toJson());
        rejected = false;
        try {
            prepare(d + "/sources.json", d + "/cache");
        } catch (const Error &) {
            rejected = true;
        }
        auto states = readJson(d + "/cache/preparation_status.json")["sources"].toArray();
        check(rejected && states.size() == 2 && states[0].toObject()["status"] == "error" &&
                  states[1].toObject()["status"] == "ready",
              "source failure prevents other sources from being checked");
        out << "PASS: interval union, FASTA validation, isoform selection, cancellation, "
               "BLASTP/TBLASTN, gene deduplication, negative control, reverse coordinates, output "
               "protection, disabled targets, partial errors\n";
        return 0;
    } catch (const std::exception &e) {
        QTextStream(stderr) << e.what() << Qt::endl;
        return 1;
    }
}
