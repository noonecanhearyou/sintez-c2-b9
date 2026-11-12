#pragma once

#include <homology/homology.hpp>
#include <QMap>

namespace synteny {

struct Interval {
    qint64 start, end;
};

struct Location {
    QString sequence;
    QChar strand = '+';
    QVector<Interval> blocks;
};

struct Block {
    QString name, sequence;
    QChar orientation;
    Interval source, target;
    QString sourceSequence;
};

Location parseLocation(const QString &text);
QVector<Block> parsePsl(const QString &line);
qint64 coveredLength(QVector<Interval> intervals);
bool flankOrder(Interval left, Interval right, QChar orientation);
bool compatible(const QString &assembly, const QString &otherAssembly,
                const QString &sequence, const QString &otherSequence,
                QChar strand, QChar otherStrand);
QJsonObject run(const QString &pipelineConfigPath, const homology::Hooks &hooks = {});

} // namespace synteny
