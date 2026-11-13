#include "core.hpp"
#include "denovo.hpp"
#include <QFile>
#include <QTemporaryDir>
#include <QtTest>
#include <orf/orf.hpp>
#include <synteny/synteny.hpp>
class Tests : public QObject {
    Q_OBJECT
  private slots:
    void translation() {
        QCOMPARE(go::translate("ATGAAATAA"), QString("MK"));
        QVERIFY(go::translate("ATGTAATAA").isEmpty());
        QVERIFY(go::translate("ATGNNNTAA").isEmpty());
    }
    void tree() {
        auto t = go::parseTree("((A:1,B:2)x:3,C:4)root;");
        QCOMPARE(go::leafNames(t), QStringList({"A", "B", "C"}));
        QVERIFY_EXCEPTION_THROWN(go::parseTree("(A,B)"), homology::Error);
    }
    void rules() {
        QJsonObject ev{{"primary_valid", true},
                       {"homology_complete", true},
                       {"orf_available", true},
                       {"synteny_available", true},
                       {"context_species", QJsonArray{"dsim", "dsec"}},
                       {"interrupted_species", QJsonArray{"dsim", "dsec"}}};
        QCOMPARE(go::classify("g", ev)["class"].toString(), QString("ambiguous"));
        ev["ancestral_noncoding_verified"] = true;
        QCOMPARE(go::classify("g", ev)["class"].toString(), QString("de_novo"));
        ev["ancestral_noncoding_verified"] = false;
        ev["strong_paralogs"] = 1;
        QCOMPARE(go::classify("g", ev)["class"].toString(), QString("duplication"));
        ev["strong_distant_species"] = QJsonArray{"dwil", "dpse"};
        QCOMPARE(go::classify("g", ev)["class"].toString(), QString("ambiguous"));
        ev["strong_paralogs"] = 0;
        QCOMPARE(go::classify("g", ev)["class"].toString(), QString("ancient_diverged"));
        ev["homology_complete"] = false;
        QCOMPARE(go::classify("g", ev)["class"].toString(), QString("ambiguous"));
        ev["homology_complete"] = true;
        ev["other_mechanism_verified"] = true;
        QCOMPARE(go::classify("g", ev)["class"].toString(), QString("other"));
        ev["primary_valid"] = false;
        QCOMPARE(go::classify("g", ev)["support_score"].toInt(), 0);
    }
    void branchScreen() {
        const QString gene = "FBgn0000001";
        const QString ref = "ATG" + QString("GCT").repeated(127) + "TAA";
        const QString outer = "contig-:101:180";
        QVector<go::Record> cohort{
            {gene, "dmel", "raw", ">dmel|FBgn0000001|score_1E10|shift_No|1", gene, ref},
            {gene, "dwil", "raw", ">dwil|annotated|score_1E10|shift_No|1", "annotated", ref},
            {gene, "sleb", "raw", ">sleb|" + outer, outer, "ATG"},
            {gene, "sleb", "spaln", ">sleb|" + outer + "|score_60|shift_No|1", outer, "ATG"},
            {gene, "sleb", "genewise", ">sleb|" + outer + "|score_1.6|shift_No|1", outer, "ATG"}};
        QSet<QString> scope{"dmel", "dwil", "sleb"};
        QVector<go::Row> loci{{{"gene_id", gene}, {"species", "sleb"}, {"contig", "contig"},
                              {"strand", "-"}, {"start_0based", "100"},
                              {"end_0based_exclusive", "180"}, {"provenance", "context_only"}}};
        QVector<go::Row> states;
        for (const auto &species : scope)
            states.append({{"gene_id", gene}, {"target_species", species},
                           {"search_type", "cross_species_blastp"},
                           {"status", "completed_no_accepted_hits"}});
        QVector<go::Row> hits{{{"gene_id", gene}, {"target_species", "dwil"},
                             {"search_type", "cross_species_blastp"}, {"passes_filters", "true"}}};
        auto screen = [&] { return go::denovoBranches(cohort, {{gene, 128}}, scope, loci, hits, states).value(gene); };
        auto evidence = screen();
        QCOMPARE(evidence["status"].toString(), QString("preliminary_candidate"));
        QCOMPARE(evidence["branch"].toString(), QString("Br8"));
        QJsonObject ev{{"primary_valid", true}, {"homology_complete", true},
                       {"accepted_cross_species", QJsonArray{"dwil"}},
                       {"de_novo_branch_screen", evidence}};
        QCOMPARE(go::classify(gene, ev)["class"].toString(), QString("de_novo"));
        hits[0]["target_species"] = "sleb";
        QCOMPARE(screen()["status"].toString(), QString("contradictory_local_homology"));
        loci[0]["strand"] = "+";
        QCOMPARE(screen()["status"].toString(), QString("outgroup_coordinates_unverified"));
    }
    void catalogAndSubset() {
        QTemporaryDir d;
        QString txt;
        for (auto section : {QString("raw"), QString("spaln"), QString("genewise")}) {
            txt += "Section " +
                   (section == "raw"     ? QString("I. raw")
                    : section == "spaln" ? QString("II. spaln")
                                         : QString("III. genewise")) +
                   "\n";
            for (auto gene : {QString("FBgn1"), QString("FBgn2")})
                txt += (section == "raw" ? QString("raw sequence") : section) + " for dmel_" +
                       gene + "\n>dmel|" + gene +
                       "|score_1E10|shift_No|1\nATGAAATAA\n>dsim|other|score_1E10|"
                       "shift_No|1\nATGAAATAA\n";
        }
        homology::writeText(d.path() + "/in.txt", txt.toUtf8());
        auto c = go::catalog(d.path() + "/in.txt");
        QCOMPARE(c.genes.size(), 2);
        QCOMPARE(c.species.size(), 2);
        go::subset(d.path() + "/in.txt", {"FBgn2"}, d.path() + "/out.txt");
        QCOMPARE(go::catalog(d.path() + "/out.txt").genes, QStringList{"FBgn2"});
        QCOMPARE(go::records(d.path() + "/out.txt").size(), 6);
    }
    void coords() {
        auto l = synteny::parseLocation("NC_1:complement(join(10..20,30..40))");
        QCOMPARE(l.strand, QChar('-'));
        QCOMPARE(l.blocks.size(), 2);
        QCOMPARE(synteny::coveredLength(l.blocks), qint64(22));
    }
};
QTEST_GUILESS_MAIN(Tests)
#include "tests.moc"
