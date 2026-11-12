#pragma once
#include <QJsonArray>
#include <QMap>
#include <QSet>
#include <homology/homology.hpp>
#include <memory>
namespace go {
using Row = QMap<QString, QString>;
struct Record {
    QString gene, species, section, header, locus, dna;
};
struct Catalog {
    QStringList genes, species;
    QMap<QString, int> counts;
};
struct Request {
    QString hal, alignment, assemblies, reference = "dmel", outputRoot, workspace, blastBin, halBin;
    QStringList genes, comparisonSpecies;
    bool useVerified = true;
    int threads = 2;
};
struct Result {
    QString directory, tree;
    QJsonArray genes;
    QString reference;
    QJsonObject manifest;
};
struct Node {
    QString name;
    double length = 0;
    QVector<std::shared_ptr<Node>> children;
};
Catalog catalog(const QString &path);
QVector<Record> records(const QString &path, const QSet<QString> &genes = {});
void subset(const QString &source, const QStringList &genes, const QString &output,
            const QSet<QString> &species = {});
QVector<Row> table(const QString &path, QChar delimiter = '\t');
void writeTable(const QString &path, const QStringList &columns, const QVector<Row> &rows);
QMap<QString, QString> speciesNames();
QString displaySpecies(const QString &code);
QString translate(const QString &dna);
std::shared_ptr<Node> parseTree(const QString &text);
QStringList leafNames(const std::shared_ptr<Node> &node);
QJsonObject classify(const QString &gene, const QJsonObject &evidence);
QString className(const QString &code);
QString findHalBin();
QString findBlastBin(const QString &workspace);
Result execute(const Request &, const homology::Hooks &hooks = {},
               std::function<void(int, QString)> progress = {});
Result loadResult(const QString &directory);
} // namespace go
