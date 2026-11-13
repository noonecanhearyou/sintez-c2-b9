#include "homology/homology.hpp"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QSet>
#include <QTextStream>
using namespace homology;
using Row = QMap<QString, QString>;
static void check(bool ok, const QString &message) {
    if (!ok)
        throw Error(message);
}
static QVector<Row> table(QString path) {
    QFile f(path);
    check(f.open(QIODevice::ReadOnly), "Cannot read " + path);
    QTextStream in(&f);
    auto names = in.readLine().split('\t');
    QVector<Row> rows;
    while (!in.atEnd()) {
        auto l = in.readLine();
        if (l.isEmpty())
            continue;
        auto values = l.split('\t');
        check(values.size() == names.size(), "Table column mismatch " + path);
        Row row;
        for (int i = 0; i < names.size(); ++i)
            row[names[i]] = values[i];
        rows << row;
    }
    return rows;
}
static QSet<QString> list(QString value) {
    QSet<QString> r;
    for (auto x : value.split(',', Qt::SkipEmptyParts))
        r.insert(x);
    return r;
}
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    auto args = app.arguments();
    try {
        check(args.size() == 4, "Usage: homology-audit SOURCES.json DATA_DIR RESULT_DIR");
        QString dir = args[3];
        auto manifest = readJson(dir + "/run.json");
        auto cfg = manifest["config"].toObject();
        auto targets = cfg["targets"].toArray();
        auto sources = readJson(args[1])["sources"].toArray();
        check(manifest["status"] == "completed" || manifest["status"] == "partial_not_run",
              "Unresolved search errors");
        QMap<QString, Row> targetMap;
        QSet<QString> species, sourceSpecies;
        int enabled = 0, proteinTargets = 0, genomeTargets = 0;
        QJsonArray inventory;
        for (auto v : sources) {
            auto s = v.toObject();
            QString id = s["id"].toString();
            sourceSpecies.insert(s["species"].toString());
            auto provenance = readJson(args[2] + "/" + id + "/source.json");
            check(provenance["specification"].toObject()["assembly"] == s["assembly"],
                  "Source assembly mismatch");
            inventory.append(QJsonObject{{"species", s["species"]},
                                         {"assembly", s["assembly"]},
                                         {"status", provenance.value("status").toString("ready")}});
        }
        for (auto v : targets) {
            auto t = v.toObject();
            QString id = t["id"].toString();
            check(!targetMap.contains(id), "Duplicate target ID");
            targetMap[id] = {{"species", t["species"].toString()},
                             {"assembly", t["assembly"].toString()},
                             {"kind", t["kind"].toString()}};
            species.insert(t["species"].toString());
            if (t.value("enabled").toBool(true))
                ++enabled;
            if (t["kind"] == "protein")
                ++proteinTargets;
            else
                ++genomeTargets;
        }
        check(species == sourceSpecies && species.size() == sources.size(),
              "Source species scope mismatch");
        check(targets.size() == 2 * sources.size() - 1 && proteinTargets == sources.size() &&
                  genomeTargets == sources.size() - 1,
              "Unexpected target scope");
        auto summary = table(dir + "/homology_summary.tsv");
        check(summary.size() == 780, "Expected 780 genes");
        QSet<QString> genes;
        for (auto row : summary) {
            check(!genes.contains(row["gene_id"]), "Duplicate summary gene");
            genes.insert(row["gene_id"]);
            check(row["error_targets"] == "0" && row["completed_targets"].toInt() == enabled,
                  "Incomplete gene summary");
        }
        auto statuses = table(dir + "/search_status.tsv");
        check(statuses.size() == 780 * targets.size(), "Missing query/target states");
        QSet<QString> statusKeys;
        int errors = 0, skipped = 0;
        for (auto row : statuses) {
            QString key = row["query_id"] + '\t' + row["target_id"];
            check(!statusKeys.contains(key), "Duplicate status");
            statusKeys.insert(key);
            check(genes.contains(row["gene_id"]) && targetMap.contains(row["target_id"]),
                  "Unknown gene/target in statuses");
            if (row["status"] == "error")
                ++errors;
            if (row["status"] == "not_run")
                ++skipped;
        }
        check(errors == 0, "Search error states present");
        QMap<QString, QSet<QString>> paralogs, proteinSpecies, genomeSpecies, partialGenes,
            partialProtein, partialGenome, proteinBySpecies, genomeBySpecies, partialBySpecies;
        auto hits = table(dir + "/homology_hits.tsv");
        for (auto row : hits) {
            QString g = row["gene_id"], sp = row["target_species"], type = row["search_type"],
                    tid = row["target_id"];
            check(genes.contains(g) && targetMap.contains(tid), "Unknown hit gene/target");
            check(targetMap[tid]["species"] == sp && targetMap[tid]["assembly"] == row["assembly"],
                  "Hit assembly mismatch");
            if (row["passes_filters"] == "true") {
                check(row["evalue"].toDouble() <=
                          cfg["filters"].toObject()["accept_evalue"].toDouble(),
                      "Accepted E-value exceeds threshold");
                check(row["query_coverage"].toDouble() + 1e-7 >=
                          cfg["filters"].toObject()["min_query_coverage"].toDouble(),
                      "Accepted coverage below threshold");
                if (type == "within_species_blastp") {
                    check(g != row["target_gene_id"], "Own gene accepted");
                    paralogs[g].insert(row["target_gene_id"]);
                    proteinBySpecies[sp].insert(g);
                } else if (type == "cross_species_blastp") {
                    proteinSpecies[g].insert(sp);
                    proteinBySpecies[sp].insert(g);
                } else {
                    genomeSpecies[g].insert(sp);
                    genomeBySpecies[sp].insert(g);
                }
            }
            if (row["quality"] == "significant_partial_similarity") {
                partialBySpecies[sp].insert(g);
                if (type == "within_species_blastp")
                    partialGenes[g].insert(row["target_gene_id"]);
                else if (type == "cross_species_blastp")
                    partialProtein[g].insert(sp);
                else
                    partialGenome[g].insert(sp);
            }
        }
        int withParalogs = 0, withProtein = 0, withGenome = 0, withPartial = 0, shortQueries = 0,
            noAcceptedCross = 0;
        for (auto row : summary) {
            auto g = row["gene_id"];
            check(list(row["paralog_gene_ids"]) == paralogs[g] &&
                      row["paralog_candidate_count"].toInt() == paralogs[g].size(),
                  "Paralog summary mismatch");
            check(list(row["protein_hit_species"]) == proteinSpecies[g] &&
                      list(row["genome_hit_species"]) == genomeSpecies[g],
                  "Cross-species summary mismatch");
            check(list(row["partial_similarity_gene_ids"]) == partialGenes[g] &&
                      list(row["partial_protein_species"]) == partialProtein[g] &&
                      list(row["partial_genome_species"]) == partialGenome[g],
                  "Partial evidence lost during merge");
            if (!paralogs[g].isEmpty())
                ++withParalogs;
            if (!proteinSpecies[g].isEmpty())
                ++withProtein;
            if (!genomeSpecies[g].isEmpty())
                ++withGenome;
            if (!partialGenes[g].isEmpty() || !partialProtein[g].isEmpty() ||
                !partialGenome[g].isEmpty())
                ++withPartial;
            if (row["limitations"].contains("short_query"))
                ++shortQueries;
            if (proteinSpecies[g].isEmpty() && genomeSpecies[g].isEmpty())
                ++noAcceptedCross;
        }
        auto loci = table(dir + "/candidate_loci.tsv");
        for (auto row : loci) {
            qint64 l = row["start_1based"].toLongLong(), r = row["end_1based"].toLongLong();
            check(l >= 1 && r >= l && row["start_0based"].toLongLong() == l - 1 &&
                      row["end_0based_exclusive"].toLongLong() == r,
                  "Invalid locus/BED coordinates");
            check(row["strand"] == "+" || row["strand"] == "-", "Missing locus strand");
        }
        QString speciesTable = "species\tgenes_with_accepted_protein_hits\tgenes_with_accepted_"
                               "genome_hits\tgenes_with_partial_similarity\n";
        auto spList = species.values();
        spList.sort();
        for (auto sp : spList)
            speciesTable += QStringList{sp, QString::number(proteinBySpecies[sp].size()),
                                        QString::number(genomeBySpecies[sp].size()),
                                        QString::number(partialBySpecies[sp].size())}
                                .join('\t') +
                            '\n';
        writeText(dir + "/species_summary.tsv", speciesTable.toUtf8());
        QJsonObject report{{"status", "passed"},
                           {"species_count", species.size()},
                           {"target_count", targets.size()},
                           {"gene_count", genes.size()},
                           {"query_target_states", statuses.size()},
                           {"error_states", errors},
                           {"not_run_states", skipped},
                           {"raw_hsps", hits.size()},
                           {"accepted_genomic_hsps", loci.size()},
                           {"genes_with_paralog_candidates", withParalogs},
                           {"genes_with_cross_species_protein_hits", withProtein},
                           {"genes_with_cross_species_genome_hits", withGenome},
                           {"genes_with_partial_similarity", withPartial},
                           {"genes_without_accepted_cross_species_hits", noAcceptedCross},
                           {"short_queries_under_30aa", shortQueries},
                           {"excluded_species", cfg.value("excluded_species")},
                           {"source_inventory", inventory}};
        writeText(dir + "/validation.json", QJsonDocument(report).toJson());
        QTextStream(stdout) << QJsonDocument(report).toJson();
        return 0;
    } catch (const std::exception &e) {
        QTextStream(stderr) << e.what() << Qt::endl;
        return 1;
    }
}
