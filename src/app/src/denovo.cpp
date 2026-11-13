#include "denovo.hpp"
#include <QRegularExpression>
#include <algorithm>
#include <cmath>
#include <optional>

namespace go {
namespace {
constexpr double cutoff = 1e-6;
const QVector<QStringList> groups{
    {"dmel"}, {"dsim", "dsec"}, {"dyak", "dere"}, {"dfic"},
    {"drho", "dele", "deug", "dbia", "dtak"}, {"dkik"}, {"dana"},
    {"dmir", "dpse", "dper"}, {"dwil"}, {"sleb"}, {"bdor", "blat"}};

int branch(const QString &species) {
    for (int i = 0; i < groups.size(); ++i)
        if (groups[i].contains(species))
            return i;
    return -1;
}

QString branchName(int index) {
    return index == 10 ? "Bactrocera" : QString("Br%1").arg(index);
}

struct Interval {
    QString contig, strand;
    qint64 start, end;
};

std::optional<Interval> inputInterval(const QString &label) {
    static const QRegularExpression pattern(R"(^(.+)([+-]):(\d+):(\d+)$)");
    const auto match = pattern.match(label);
    if (!match.hasMatch())
        return {};
    bool startOk, endOk;
    const auto start = match.captured(3).toLongLong(&startOk);
    const auto end = match.captured(4).toLongLong(&endOk);
    if (!startOk || !endOk || start < 1 || end < start)
        return {};
    // extract_raw_hit.py slices source DNA at start - 1 through end.
    return Interval{match.captured(1), match.captured(2), start - 1, end};
}

bool overlaps(const Interval &a, const Interval &b) {
    return a.contig == b.contig && a.strand == b.strand &&
           std::max(a.start, b.start) < std::min(a.end, b.end);
}

struct Curve {
    double baseline, fastAmplitude, fastRate, slowAmplitude, slowRate;
    double at(int length) const {
        return baseline + fastAmplitude * std::exp(-fastRate * length) +
               slowAmplitude * std::exp(-slowRate * length);
    }
};

struct Model { Curve mean, sd; };
// Peng & Zhao, archived extract_alnscore_genewise.py, alnscore models.
constexpr Model models[]{
    {{-.039279974, 2.736299221, .015644801, 15.302980734, .088245361},
     {0, .128773380, .001565254, 1.049427652, .036506722}},
    {{.010964339, .148062602, .015192716, 10.840972709, .322977380},
     {.000000102, .190068618, .020887105, .020123856, .001636796}}
};

struct Locus {
    QString species, label;
    Interval interval{};
    QVector<double> scores[2];
    bool raw = false;
};

struct Gene {
    int length = 0;
    QSet<QString> annotated;
    QMap<QString, Locus> loci;
};

struct MethodResult {
    int coding = -1, outer = -1;
    double p = 1, adjusted = 1;
};

QString key(const QString &gene, const QString &species) {
    return gene + '\t' + species;
}

int referenceLength(const QString &dna) {
    auto sequence = dna;
    sequence.remove('-');
    if (sequence.isEmpty() || sequence.size() % 3 ||
        sequence.contains(QRegularExpression("[^ACGT]")))
        return 0;
    const auto end = sequence.right(3);
    return sequence.size() / 3 - (end == "TAA" || end == "TAG" || end == "TGA");
}

std::optional<double> probability(const Locus &locus, int method, int length) {
    if (!locus.raw || length <= 0 || locus.scores[method].size() != 1)
        return {};
    const auto &model = models[method];
    const double normalized = locus.scores[method].first() / length;
    return .5 * std::erfc((normalized - model.mean.at(length)) /
                          (std::sqrt(2.0) * model.sd.at(length)));
}

QMap<QString, MethodResult> assessMethod(const QMap<QString, Gene> &genes, int method) {
    QMap<QString, MethodResult> results;
    for (auto it = genes.cbegin(); it != genes.cend(); ++it) {
        const auto &gene = it.value();
        if (!gene.length || !gene.annotated.contains("dmel"))
            continue;
        QMap<int, double> groupsP;
        for (const auto &species : gene.annotated)
            groupsP[branch(species)] = 0;
        for (const auto &locus : gene.loci) {
            const auto p = probability(locus, method, gene.length);
            if (p)
                groupsP[branch(locus.species)] = std::min(groupsP.value(branch(locus.species), 1.), *p);
        }
        MethodResult result;
        for (auto group = groupsP.cbegin(); group != groupsP.cend(); ++group)
            if (group.value() < cutoff)
                result.coding = group.key();
        auto outer = groupsP.upperBound(result.coding);
        if (outer != groupsP.end()) {
            result.outer = outer.key();
            result.p = outer.value();
        } else if (result.coding >= 0) {
            result.p = groupsP.value(result.coding);
        }
        results[it.key()] = result;
    }
    auto order = results.keys();
    std::stable_sort(order.begin(), order.end(), [&](const QString &a, const QString &b) {
        return results.value(a).p < results.value(b).p;
    });
    double adjusted = 1;
    for (qsizetype i = order.size(); i > 0; --i) {
        auto &result = results[order[i - 1]];
        adjusted = std::min(adjusted, result.p * order.size() / i);
        result.adjusted = adjusted;
    }
    return results;
}

QJsonObject methodJson(const MethodResult &result) {
    return {{"furthest_coding_branch", result.coding}, {"first_outer_branch", result.outer},
            {"p_raw", result.p}, {"p_bh", result.adjusted}};
}
}

QMap<QString, QJsonObject> denovoBranches(const QVector<Record> &cohort,
                                        const QMap<QString, int> &queryLengths,
                                        const QSet<QString> &scope,
                                        const QVector<Row> &syntenicLoci,
                                        const QVector<Row> &hits,
                                        const QVector<Row> &searchStates) {
    QMap<QString, Gene> genes;
    for (const auto &record : cohort) {
        if (branch(record.species) < 0)
            continue;
        auto &gene = genes[record.gene];
        const auto interval = inputInterval(record.locus);
        if (!interval) {
            if (record.section == "raw" && !record.dna.isEmpty() &&
                record.header.section('|', 2, 2) == "score_1E10") {
                gene.annotated.insert(record.species);
                if (record.species == "dmel" && record.locus == record.gene)
                    gene.length = referenceLength(record.dna);
            }
            continue;
        }
        auto &locus = gene.loci[key(record.species, record.locus)];
        locus.species = record.species;
        locus.label = record.locus;
        locus.interval = *interval;
        if (record.section == "raw") {
            locus.raw = !record.dna.isEmpty();
        } else if (record.section == "spaln" || record.section == "genewise") {
            const auto field = record.header.section('|', 2, 2);
            bool ok;
            const double score = field.mid(6).toDouble(&ok);
            if (field.startsWith("score_") && ok && std::isfinite(score) && std::abs(score) < 1e10)
                locus.scores[record.section == "spaln" ? 0 : 1].append(score);
        }
    }
    const auto spaln = assessMethod(genes, 0), genewise = assessMethod(genes, 1);

    QMap<QString, QVector<Interval>> mapped;
    for (const auto &row : syntenicLoci)
        if (row["provenance"] == "direct_CDS" || row["provenance"] == "context_only")
            mapped[key(row["gene_id"], row["species"])].append(
                {row["contig"], row["strand"], row["start_0based"].toLongLong(),
                 row["end_0based_exclusive"].toLongLong()});
    QMap<QString, QVector<Row>> accepted;
    for (const auto &row : hits)
        if (row["passes_filters"] == "true")
            accepted[row["gene_id"]].append(row);
    QMap<QString, QSet<QString>> completedProteinSearches;
    for (const auto &row : searchStates)
        if (row["status"].startsWith("completed_") && row["search_type"].endsWith("blastp"))
            completedProteinSearches[row["gene_id"]].insert(row["target_species"]);

    QMap<QString, QJsonObject> results;
    for (auto it = queryLengths.cbegin(); it != queryLengths.cend(); ++it) {
        const auto &id = it.key();
        const auto &gene = genes.value(id);
        const auto a = spaln.value(id), b = genewise.value(id);
        const int coding = a.coding;
        QJsonObject evidence{{"status", "insufficient_branch_evidence"},
                             {"spaln", methodJson(a)}, {"genewise", methodJson(b)},
                             {"spaln_bh_cohort_size", spaln.size()}, {"genewise_bh_cohort_size", genewise.size()},
                             {"threshold", cutoff}, {"external_filters", "not_run"}};
        results[id] = evidence;
        if (it.value() != gene.length || coding < 1 || coding > 9 || coding != b.coding ||
            a.outer < 0 || b.outer < 0 || a.adjusted < cutoff || b.adjusted < cutoff)
            continue;
        evidence["branch"] = branchName(coding);
        QJsonArray outer;
        QMap<QString, QVector<Interval>> weakIntervals;
        for (const auto &locus : gene.loci) {
            if (branch(locus.species) <= coding || !scope.contains(locus.species))
                continue;
            const auto p = probability(locus, 0, gene.length);
            const auto q = probability(locus, 1, gene.length);
            if (!p || !q || *p < cutoff || *q < cutoff)
                continue;
            const auto regions = mapped.value(key(id, locus.species));
            if (std::none_of(regions.cbegin(), regions.cend(), [&](const auto &region) {
                    return overlaps(locus.interval, region);
                }))
                continue;
            weakIntervals[locus.species].append(locus.interval);
            outer.append(QJsonObject{{"species", locus.species}, {"locus", locus.label},
                                     {"spaln_p_raw", *p}, {"genewise_p_raw", *q}});
        }
        evidence["outgroup_loci"] = outer;
        if (outer.isEmpty()) {
            evidence["status"] = "outgroup_coordinates_unverified";
            results[id] = evidence;
            continue;
        }
        QSet<QString> blockers;
        for (const auto &hit : accepted.value(id)) {
            const auto sp = hit["target_species"];
            const auto type = hit["search_type"];
            if (type == "within_species_blastp") {
                blockers.insert(sp);
                continue;
            }
            if (branch(sp) >= 0 && branch(sp) <= coding)
                continue;
            if (type == "tblastn") {
                const auto start = hit["target_start"].toLongLong(), end = hit["target_end"].toLongLong();
                Interval aligned{hit["target_sequence_id"], hit["strand"],
                                 std::min(start, end) - 1, std::max(start, end)};
                const auto intervals = weakIntervals.value(sp);
                if (std::any_of(intervals.cbegin(), intervals.cend(), [&](const auto &locus) {
                        return aligned.contig == locus.contig && aligned.strand == locus.strand &&
                               aligned.start >= locus.start && aligned.end <= locus.end;
                    }))
                    continue;
            }
            blockers.insert(sp);
        }
        auto blocked = blockers.values();
        blocked.sort();
        evidence["contradictory_hit_species"] = QJsonArray::fromStringList(blocked);
        auto missing = (scope - completedProteinSearches.value(id)).values();
        missing.sort();
        evidence["missing_local_protein_searches"] = QJsonArray::fromStringList(missing);
        evidence["status"] = !blockers.isEmpty() ? "contradictory_local_homology" :
                             !missing.isEmpty() ? "local_search_incomplete" : "preliminary_candidate";
        results[id] = evidence;
    }
    return results;
}
}
