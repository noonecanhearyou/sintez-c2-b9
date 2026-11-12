#include "homology/homology.hpp"
#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QRegularExpression>
#include <QSet>
#include <QTextStream>
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    auto a = app.arguments();
    QTextStream out(stdout), err(stderr);
    try {
        homology::Hooks hooks;
        hooks.progress = [&](const QString &s) { err << s << Qt::endl; };
        if (a.size() == 3 && (a[1] == "run" || a[1] == "run-batch")) {
            auto r =
                a[1] == "run-batch" ? homology::runBatch(a[2], hooks) : homology::run(a[2], hooks);
            out << QJsonDocument(r).toJson();
            return r["status"].toString() == "error" || r["status"].toString() == "partial_error"
                       ? 2
                       : 0;
        }
        if (a.size() >= 4 && a[1] == "merge") {
            auto cfg = homology::readJson(a[2]);
            QString dest = cfg.value("output_dir").toString();
            if (dest.isEmpty())
                throw homology::Error("Missing output_dir");
            if (!QFileInfo(dest).isAbsolute())
                dest = QDir(QFileInfo(a[2]).absolutePath()).absoluteFilePath(dest);
            auto result = homology::mergeRuns(a[2], a.mid(3), dest);
            out << QJsonDocument(result).toJson();
            return result.value("status").toString() == "completed" ? 0 : 2;
        }
        if (a.size() == 4 && a[1] == "prepare") {
            homology::prepare(a[2], a[3], hooks);
            return 0;
        }
        if (a.size() == 5 && a[1] == "extract-candidates") {
            QFile f(a[3]);
            if (!f.open(QIODevice::ReadOnly))
                throw homology::Error("Cannot open candidate source");
            QTextStream in(&f);
            QStringList genes;
            QSet<QString> seen;
            while (!in.atEnd()) {
                auto l = in.readLine();
                if (l.startsWith("Section II"))
                    break;
                auto m = QRegularExpression("^raw sequence for dmel_(FBgn[0-9]+)\\s*$").match(l);
                if (m.hasMatch() && !seen.contains(m.captured(1))) {
                    seen.insert(m.captured(1));
                    genes << m.captured(1);
                }
            }
            homology::extractQueries(a[2], genes, a[4]);
            out << "Extracted " << genes.size() << " genes from Section I" << Qt::endl;
            return 0;
        }
        if (a.size() == 5 && a[1] == "extract") {
            homology::extractQueries(a[2], a[3].split(',', Qt::SkipEmptyParts), a[4]);
            return 0;
        }
        if (a.size() == 5 && a[1] == "verify-trial") {
            auto r = homology::verifyTrial(a[2], a[3], a[4]);
            out << QJsonDocument(r).toJson();
            return r["verified"].toBool() ? 0 : 2;
        }
        err << "Usage:\n  homology-cli prepare SOURCES.json DATA_DIR\n  homology-cli extract "
               "PROTEIN_DIR gene1,gene2 QUERY_DIR\n  homology-cli extract-candidates PROTEIN_DIR "
               "SOURCE.txt QUERY_DIR\n  homology-cli verify-trial SOURCE.txt FLYBASE_DIR "
               "OUTPUT_DIR\n  homology-cli run CONFIG.json\n  homology-cli run-batch CONFIG.json\n "
               " homology-cli merge CONFIG.json RUN_DIR...\n";
        return 1;
    } catch (const std::exception &e) {
        err << e.what() << Qt::endl;
        return 2;
    }
}
