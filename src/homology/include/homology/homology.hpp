#pragma once
#include <QJsonObject>
#include <QString>
#include <QStringList>
#include <QVector>
#include <atomic>
#include <functional>
#include <stdexcept>

namespace homology {
struct Error : std::runtime_error {
    explicit Error(const QString &s) : std::runtime_error(s.toStdString()) {}
};
struct Sequence {
    QString id, header, sequence;
};
struct Protein {
    QString proteinId, geneId, transcriptId, species, sequence, location;
};
struct Hooks {
    // run() blocks: call from a Qt worker thread. Queue UI updates to the GUI thread.
    std::atomic_bool *cancelled = nullptr;
    std::function<void(const QString &)> progress;
};
QVector<Sequence> readFasta(const QString &path, bool protein);
void writeText(const QString &path, const QByteArray &bytes);
QJsonObject readJson(const QString &path);
QString sha256(const QString &path);
QByteArray process(const QString &executable, const QStringList &args, int timeoutSeconds,
                   const Hooks &hooks = {}, const QString &stdoutFile = {});
QVector<Protein> flybaseProteins(const QString &fasta);
QVector<Protein> ncbiProteins(const QString &fasta, const QString &gff, const QString &species);
void writeProteins(const QString &dir, const QVector<Protein> &proteins);
// Select longest protein per gene; ties resolved by protein ID, deterministically.
void extractQueries(const QString &proteinDir, const QStringList &genes, const QString &outputDir);
QJsonObject verifyTrial(const QString &source, const QString &proteinDir, const QString &outputDir);
// Downloads only explicit sources from the manifest. Cached files must match recorded hashes.
void prepare(const QString &manifestPath, const QString &outputDir, const Hooks &hooks = {});
// Config paths are relative to the config file. Existing output directories are rejected.
QJsonObject run(const QString &configPath, const Hooks &hooks = {});
QJsonObject runBatch(const QString &configPath, const Hooks &hooks = {});
QJsonObject mergeRuns(const QString &configPath, const QStringList &runs, const QString &outputDir);
struct Interval {
    qint64 start, end;
};
qint64 coveredLength(QVector<Interval> intervals);
} // namespace homology
