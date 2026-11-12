#include "homology/homology.hpp"
#include <QCoreApplication>
#include <QFile>
#include <QJsonDocument>
#include <QMap>
#include <QTextStream>
using namespace homology;
using Row = QMap<QString, QString>;
static void check(bool ok, const QString &message) {
    if (!ok)
        throw Error(message);
}
static QVector<Row> table(const QString &path) {
    QFile f(path);
    check(f.open(QIODevice::ReadOnly), "Cannot read " + path);
    QTextStream in(&f);
    auto columns = in.readLine().split('\t');
    QVector<Row> rows;
    while (!in.atEnd()) {
        QString line = in.readLine();
        if (line.isEmpty())
            continue;
        auto cells = line.split('\t');
        check(cells.size() == columns.size(), "Column mismatch " + path);
        Row row;
        for (int i = 0; i < columns.size(); ++i)
            row[columns[i]] = cells[i];
        rows.push_back(row);
    }
    return rows;
}
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    auto a = app.arguments();
    try {
        check(a.size() == 4,
              "Usage: homology-validate VALIDATED_DIR CANDIDATES_DIR VERIFICATION.json");
        auto validation = readJson(a[3]);
        check(validation["verified"].toBool(), "Trial protein not verified");
        for (auto dir : {a[1], a[2]})
            check(readJson(dir + "/run.json")["status"] == "completed",
                  "Incomplete real run " + dir);
        auto controls = table(a[1] + "/homology_summary.tsv");
        check(controls.size() == 6, "Expected six control queries");
        QMap<QString, Row> control;
        for (auto r : controls)
            control[r["gene_id"]] = r;
        check(control["FBgn0000042"]["paralog_gene_ids"].split(',').contains("FBgn0000043"),
              "Act5C → Act42A positive control missed");
        check(control["FBgn0000577"]["partial_similarity_gene_ids"].split(',').contains(
                  "FBgn0001269"),
              "en/inv partial similarity missed");
        for (QString g : {"FBgn0262810", "FBgn0025616", "FBgn0029535"}) {
            check(control[g]["protein_hit_species"].split(',').contains("dsim"),
                  "Missing expected dsim protein hit " + g);
            check(control[g]["genome_hit_species"].split(',').contains("dsim"),
                  "Missing expected dsim genome hit " + g);
        }
        auto summary = table(a[2] + "/homology_summary.tsv");
        check(summary.size() == 780, "Expected 780 candidates");
        int paralogGenes = 0, proteinGenes = 0, genomeGenes = 0, shortGenes = 0, partialGenes = 0;
        for (auto r : summary) {
            check(r["search_status"] == "completed", "Candidate incomplete");
            check(r["completed_targets"] == "3" && r["error_targets"] == "0",
                  "Candidate missing target");
            if (r["paralog_candidate_count"].toInt() > 0)
                ++paralogGenes;
            if (!r["protein_hit_species"].isEmpty())
                ++proteinGenes;
            if (!r["genome_hit_species"].isEmpty())
                ++genomeGenes;
            if (r["limitations"].contains("short_query"))
                ++shortGenes;
            if (!r["partial_similarity_gene_ids"].isEmpty() ||
                !r["partial_protein_species"].isEmpty() || !r["partial_genome_species"].isEmpty())
                ++partialGenes;
        }
        auto statuses = table(a[2] + "/search_status.tsv");
        check(statuses.size() == 780 * 3, "Expected one status per query and target");
        for (auto dir : {a[1], a[2]}) {
            for (auto r : table(dir + "/homology_hits.tsv"))
                if (r["search_type"] == "within_species_blastp" && r["passes_filters"] == "true")
                    check(r["gene_id"] != r["target_gene_id"], "Self-gene passed filters");
            for (auto r : table(dir + "/candidate_loci.tsv")) {
                auto l = r["start_1based"].toLongLong(), u = r["end_1based"].toLongLong();
                check(l >= 1 && u >= l && r["start_0based"].toLongLong() == l - 1 &&
                          r["end_0based_exclusive"].toLongLong() == u,
                      "Coordinate conversion failed");
            }
        }
        QJsonObject report{{"status", "passed"},
                           {"candidate_genes", summary.size()},
                           {"candidate_target_searches", statuses.size()},
                           {"genes_with_accepted_paralog_candidates", paralogGenes},
                           {"genes_with_accepted_dsim_protein_hits", proteinGenes},
                           {"genes_with_accepted_dsim_genome_hits", genomeGenes},
                           {"genes_with_significant_partial_similarity", partialGenes},
                           {"short_queries_under_30_aa", shortGenes},
                           {"control_queries", controls.size()},
                           {"Act5C_to_Act42A", "passed"},
                           {"en_inv_partial_similarity", "passed"},
                           {"self_gene_exclusion", "passed"},
                           {"coordinate_conversion", "passed"},
                           {"trial_annotation_match", "FBpp0297196"},
                           {"scope", "dmel_r6.36 + dsim_GCF_000754195.2 only"}};
        QTextStream(stdout) << QJsonDocument(report).toJson();
        return 0;
    } catch (const std::exception &e) {
        QTextStream(stderr) << e.what() << Qt::endl;
        return 1;
    }
}
