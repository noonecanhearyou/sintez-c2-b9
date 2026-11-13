#include "window.hpp"
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QGraphicsEllipseItem>
#include <QJsonDocument>
#include <QMimeData>
#include <QTemporaryDir>
#include <QtTest>
class GuiTests : public QObject {
    Q_OBJECT
  private slots:
    void picker() {
        GenePicker p;
        p.setGenes({"FBgn1", "FBgn2", "FBgn3", "FBgn4"});
        p.show();
        p.click();
        QVERIFY(p.findChild<QFrame *>("popup")->isVisible());
        p.selectGenes({"FBgn1", "FBgn3"});
        QCOMPARE(p.selected().size(), 2);
        auto *search = p.findChild<QLineEdit *>();
        QVERIFY(search);
        search->setText("2");
        QCOMPARE(p.selected(), QStringList({"FBgn1", "FBgn3"}));
        auto *list = p.findChild<QListWidget *>();
        QVERIFY(list->item(0)->isHidden());
        QVERIFY(!list->item(1)->isHidden());
        list->item(1)->setCheckState(Qt::Checked);
        QCOMPARE(p.selected().size(), 3);
    }
    void treeSelectionAndExport() {
        EvidenceTree tree;
        tree.resize(760, 450);
        QJsonObject species{{"dmel", QJsonObject{{"homology", "принятое совпадение"},
                                                 {"orf", "preserved_both_identical"}}},
                            {"dsim", QJsonObject{{"homology", "принятое совпадение"}}}};
        tree.show();
        tree.setResult("(Drosophila_melanogaster:1,Drosophila_simulans:1)root;",
                       {{"species", species}, {"reference_gene_id", "FBgn1"}}, "FBgn1", "dmel");
        QTest::qWait(100);
        QSignalSpy spy(&tree, &EvidenceTree::speciesSelected);
        bool clicked = false;
        for (auto *item : tree.scene()->items())
            if (item->data(0) == "dmel" && dynamic_cast<QGraphicsEllipseItem *>(item)) {
                auto p = tree.mapFromScene(item->sceneBoundingRect().center());
                QTest::mouseClick(tree.viewport(), Qt::LeftButton, Qt::NoModifier, p);
                clicked = true;
                break;
            }
        QVERIFY(clicked);
        QCOMPARE(spy.count(), 1);
        tree.fitTree();
        const double minimum = tree.transform().m11();
        for (int i = 0; i < 80; ++i) tree.zoomOut();
        QVERIFY(qAbs(tree.transform().m11() - minimum) < .001);
        for (int i = 0; i < 80; ++i) tree.zoomIn();
        QVERIFY(tree.transform().m11() <= qMax(4.0, minimum) + .001);
        auto visible = tree.mapToScene(tree.viewport()->rect()).boundingRect();
        QVERIFY(visible.intersects(tree.sceneRect()));
        tree.fitTree();
        QTemporaryDir d;
        QVERIFY(tree.savePng(d.path() + "/tree.png"));
        QVERIFY(QFileInfo(d.path() + "/tree.png").size() > 1000);
    }
    void savedResult() {
        QTemporaryDir d;
        auto gene = go::classify(
            "FBgn1",
            {{"primary_valid", true},
             {"homology_complete", true},
             {"accepted_cross_species", QJsonArray{"dsim"}},
             {"species", QJsonObject{{"dmel", QJsonObject{{"primary", true}}},
                                     {"dsim", QJsonObject{{"homology", "принятое совпадение"}}}}}});
        homology::writeText(d.path() + "/summary.json", QJsonDocument(QJsonArray{gene}).toJson());
        homology::writeText(d.path() + "/species_tree.nwk",
                            "(Drosophila_melanogaster:1,Drosophila_simulans:1)root;");
        homology::writeText(
            d.path() + "/run.json",
            QJsonDocument(
                QJsonObject{{"application", "GeneOrigin"},
                            {"status", "completed"},
                            {"reference_species", "dmel"},
                            {"result_sha256", homology::sha256(d.path() + "/summary.json")}})
                .toJson());
        MainWindow w;
        w.loadSaved(d.path());
        QCOMPARE(w.resultDirectory(), d.path());
        auto *tree = w.findChild<EvidenceTree *>();
        QVERIFY(tree && !tree->scene()->items().isEmpty());
        homology::writeText(d.path() + "/summary.json", "[]");
        QVERIFY_EXCEPTION_THROWN(go::loadResult(d.path()), homology::Error);
    }
    void fileDropAndCancel() {
        QString workspace = QString::fromUtf8(GENEORIGIN_WORKSPACE);
        if (!QFileInfo::exists(workspace + "/Acalyptratae20.hal"))
            QSKIP("Real data not installed");
        MainWindow w;
        w.show();
        QMimeData data;
        data.setUrls(
            {QUrl::fromLocalFile(workspace + "/Acalyptratae20.hal"),
             QUrl::fromLocalFile(workspace + "/supplementary_file_unannotated_alignments.txt"),
             QUrl::fromLocalFile(workspace + "/supplementary_file_genome_assembly.csv")});
        QDragEnterEvent enter(QPoint(50, 50), Qt::CopyAction, &data, Qt::LeftButton,
                              Qt::NoModifier);
        QApplication::sendEvent(&w, &enter);
        QVERIFY(enter.isAccepted());
        QDropEvent drop(QPointF(50, 50), Qt::CopyAction, &data, Qt::LeftButton, Qt::NoModifier);
        QApplication::sendEvent(&w, &drop);
        auto *picker = w.findChild<GenePicker *>();
        QVERIFY(picker);
        QTRY_COMPARE_WITH_TIMEOUT(picker->all().size(), 780, 10000);
        w.selectGenes({"FBgn0262810"});
        auto *ref = w.findChild<QComboBox *>("reference");
        QVERIFY(ref);
        ref->setCurrentIndex(ref->findData("dsim"));
        QCOMPARE(ref->currentData().toString(), QString("dsim"));
        ref->setCurrentIndex(ref->findData("dmel"));
        QTemporaryDir output;
        w.findChild<QLineEdit *>("outputRoot")->setText(output.path());
        w.startAnalysis();
        QVERIFY(w.running());
        QPushButton *cancel = nullptr;
        for (auto *b : w.findChildren<QPushButton *>())
            if (b->text() == "Отменить анализ")
                cancel = b;
        QVERIFY(cancel);
        cancel->click();
        QTRY_VERIFY_WITH_TIMEOUT(!w.running(), 10000);
        auto dirs = QDir(output.path()).entryList(QDir::Dirs | QDir::NoDotAndDotDot);
        QCOMPARE(dirs.size(), 1);
        QCOMPARE(homology::readJson(output.path() + "/" + dirs.first() + "/run.json")["status"]
                     .toString(),
                 QString("cancelled"));
    }
};
QTEST_MAIN(GuiTests)
#include "gui_tests.moc"
