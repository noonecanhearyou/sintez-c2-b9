#pragma once
#include "core.hpp"

namespace go {
// Branch numbers are specific to the D. melanogaster analysis in Fig. 1a.
QMap<QString, QJsonObject> denovoBranches(const QVector<Record> &cohort,
                                        const QMap<QString, int> &queryLengths,
                                        const QSet<QString> &scope,
                                        const QVector<Row> &syntenicLoci,
                                        const QVector<Row> &hits,
                                        const QVector<Row> &searchStates);
}
