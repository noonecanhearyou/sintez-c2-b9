#include "common.hpp"
#include <array>
namespace orf {
J compareCds(std::string r, std::string q, const Settings &s) {
    r = upper(r);
    q = upper(q);
    auto skipped = [](const char *status, const char *reason) {
        return J{{"status", status}, {"reason", reason}};
    };
    if (r.empty())
        return skipped("skipped_empty_reference", "Пустая референсная CDS");
    if (q.empty())
        return skipped("skipped_empty_query", "Пустая сравниваемая CDS");
    if (r.find_first_not_of("ACGT") != r.npos || q.find_first_not_of("ACGT") != q.npos)
        return skipped("skipped_non_acgt",
                       "Последовательность содержит неоднозначные или недопустимые нуклеотиды");
    if ((int)r.size() > s.maxAlignLength || (int)q.size() > s.maxAlignLength)
        return skipped("skipped_size_limit", "Превышена максимальная длина последовательности");
    if ((long long)r.size() * q.size() > (unsigned long long)s.maxAlignCells)
        return skipped("skipped_size_limit",
                       "Превышен лимит произведения длин последовательностей");
    int n = r.size(), m = q.size();
    using Scores = std::array<int, 3>;
    using Trace = std::array<uint8_t, 3>;
    constexpr int neg = -1000000000;
    std::vector<Scores> prev(m + 1), curr(m + 1);
    std::vector<Trace> trace(size_t(n + 1) * (m + 1));
    auto choose = [](Scores a) {
        uint8_t k = 0;
        if (a[1] > a[k])
            k = 1;
        if (a[2] > a[k])
            k = 2;
        return k;
    };
    prev[0] = {0, neg, neg};
    for (int j = 1; j <= m; ++j) {
        prev[j] = {neg, neg, -10 - (j - 1)};
        trace[j][2] = j == 1 ? 0 : 2;
    }
    for (int i = 1; i <= n; ++i) {
        curr[0] = {neg, -10 - (i - 1), neg};
        trace[size_t(i) * (m + 1)][1] = i == 1 ? 0 : 1;
        for (int j = 1; j <= m; ++j) {
            auto &t = trace[size_t(i) * (m + 1) + j];
            Scores a = prev[j - 1];
            t[0] = choose(a);
            curr[j][0] = a[t[0]] + (r[i - 1] == q[j - 1] ? 4 : -2);
            a = {prev[j][0] - 10, prev[j][1] - 1, prev[j][2] - 10};
            t[1] = choose(a);
            curr[j][1] = a[t[1]];
            a = {curr[j - 1][0] - 10, curr[j - 1][1] - 10, curr[j - 1][2] - 1};
            t[2] = choose(a);
            curr[j][2] = a[t[2]];
        }
        prev.swap(curr);
    }
    int state = choose(prev[m]), score = prev[m][state], i = n, j = m;
    std::string ops;
    while (i || j) {
        auto next = trace[size_t(i) * (m + 1) + j][state];
        ops += char('0' + state);
        if (state == 0) {
            --i;
            --j;
        } else if (state == 1)
            --i;
        else
            --j;
        state = next;
    }
    std::reverse(ops.begin(), ops.end());
    J mismatches = J::array(), gaps = J::array();
    int matched = 0, mismatch = 0, pairs = 0, gq = 0, gr = 0, events = 0, nontriplet = 0,
        triplet = 0, terminal = 0;
    std::string ar, aq;
    i = j = 0;
    for (size_t p = 0; p < ops.size();) {
        size_t e = p + 1;
        while (e < ops.size() && ops[e] == ops[p])
            ++e;
        int len = e - p;
        char op = ops[p];
        if (op == '0') {
            for (int k = 0; k < len; ++k) {
                ++pairs;
                if (r[i + k] == q[j + k])
                    ++matched;
                else {
                    ++mismatch;
                    if ((int)mismatches.size() < s.maxEvents)
                        mismatches.push_back({{"reference_position", i + k},
                                              {"query_position", j + k},
                                              {"reference_base", std::string(1, r[i + k])},
                                              {"query_base", std::string(1, q[j + k])}});
                }
            }
            ar += r.substr(i, len);
            aq += q.substr(j, len);
            i += len;
            j += len;
        } else {
            bool internal;
            J event;
            if (op == '1') {
                gq += len;
                internal = i > 0 && i + len < n && j > 0 && j < m;
                event = {{"kind", "gap_in_query"},    {"reference_start", i},
                         {"reference_end", i + len},  {"query_position", j},
                         {"length_nt", len},          {"internal", internal},
                         {"nontriplet", len % 3 != 0}};
                ar += r.substr(i, len);
                aq += std::string(len, '-');
                i += len;
            } else {
                gr += len;
                internal = j > 0 && j + len < m && i > 0 && i < n;
                event = {{"kind", "gap_in_reference"}, {"reference_position", i},
                         {"query_start", j},           {"query_end", j + len},
                         {"length_nt", len},           {"internal", internal},
                         {"nontriplet", len % 3 != 0}};
                ar += std::string(len, '-');
                aq += q.substr(j, len);
                j += len;
            }
            ++events;
            if (internal) {
                if (len % 3)
                    ++nontriplet;
                else
                    ++triplet;
            } else
                ++terminal;
            if ((int)gaps.size() < s.maxEvents)
                gaps.push_back(event);
        }
        p = e;
    }
    double identity = pairs ? double(matched) / pairs : 0, rc = double(pairs) / n,
           qc = double(pairs) / m;
    bool quality = identity >= s.minIdentity && rc >= s.minCoverage && qc >= s.minCoverage;
    J preview;
    if (s.previewColumns)
        preview = {{"reference", ar.substr(0, s.previewColumns)},
                   {"query", aq.substr(0, s.previewColumns)},
                   {"truncated", (int)ops.size() > s.previewColumns}};
    return {{"alignment_score", score / 2.0},
            {"alignment_columns", ops.size()},
            {"matched_nt", matched},
            {"mismatch_nt", mismatch},
            {"aligned_pairs", pairs},
            {"identity", round6(identity)},
            {"reference_coverage", round6(rc)},
            {"query_coverage", round6(qc)},
            {"gap_in_query_nt", gq},
            {"gap_in_reference_nt", gr},
            {"gap_event_count", events},
            {"internal_nontriplet_gap_events", nontriplet},
            {"internal_triplet_gap_events", triplet},
            {"terminal_gap_events", terminal},
            {"mismatch_examples", mismatches},
            {"gap_examples", gaps},
            {"mismatch_examples_truncated", mismatch > (int)mismatches.size()},
            {"gap_examples_truncated", events > (int)gaps.size()},
            {"alignment_preview", preview},
            {"quality_ok", quality},
            {"status", quality ? "aligned" : "low_quality_alignment"},
            {"reason", quality ? J() : J("Identity или покрытие ниже заданного порога")},
            {"interpretation",
             {{"internal_nontriplet_gaps_are", "potential_frame_disruption_candidates"},
              {"evolutionary_direction_inferred", false},
              {"codon_aware_alignment", false},
              {"quality_passed", quality}}}};
}
void comparisonStage(Context &c) {
    auto d = c.root / "4_results_third", input = c.root / "3_results_second/orf_records.jsonl";
    std::map<Key, std::vector<J>> refs, profiles;
    int refCount = 0;
    jsonl(input, [&](J r) {
        c.check();
        if (r["section"] == "raw" || r["species"] != c.s.referenceSpecies)
            return;
        auto &v = refs[key(r, {"gene_id", "section"})];
        if (v.size() < 2)
            v.push_back({{"sequence", r["sequence"]},
                         {"sequence_kind", r["sequence_kind"]},
                         {"cds_status", r["cds_validation"]["status"]},
                         {"locus", r["locus"]}});
        ++refCount;
    });
    int refPairs = refs.size();
    auto out = output(d / "comparisons.jsonl");
    Csv csv(d / "comparison_summary.csv", schema["comparative_orf"]["SUMMARY_FIELDS"]);
    J statuses = J::object(), methods = J::object(), agreements = J::object();
    std::set<std::string> genes;
    int count = 0, nontriplet = 0;
    jsonl(input, [&](J r) {
        if (r["section"] == "raw" || r["species"] == c.s.referenceSpecies)
            return;
        c.tick(++count);
        auto method = str(r, "section"), seq = str(r, "sequence");
        auto pk = key(r, {"gene_id", "species", "locus", "section"});
        profiles[pk].push_back({{"length_nt", seq.size()},
                                {"cds_status", r["cds_validation"]["status"]},
                                {"shift_label", r["source_metadata"]["shift_label"]},
                                {"sha256", hash(seq)}});
        auto &v = refs[key(r, {"gene_id", "section"})];
        J ref, comparison;
        auto skip = [&](const char *s, const char *t) {
            comparison = {{"status", s}, {"reason", t}};
        };
        if (v.empty())
            skip("skipped_missing_reference", "Референсная CDS выбранного вида отсутствует");
        else if (v.size() > 1)
            skip("skipped_multiple_references",
                 "Найдено несколько референсных последовательностей");
        else {
            ref = v[0];
            if ((ref["sequence_kind"] != "annotated_cds" && !(c.s.referenceSpecies != "dmel" && ref["sequence_kind"] == "predicted_cds")) ||
                ref["cds_status"] != "complete_candidate")
                skip("skipped_reference_noncanonical",
                     "Референсная CDS не прошла проверку структуры");
            else if (r["sequence_kind"] != "annotated_cds" && r["sequence_kind"] != "predicted_cds")
                skip("skipped_non_cds_type", "Последовательность не является подготовленной CDS");
            else
                comparison = compareCds(str(ref, "sequence"), seq, c.s);
        }
        J result = copyFields(r, {"gene_id", "species", "locus"});
        result["method"] = method;
        result["reference"] = {{"species", c.s.referenceSpecies},
                               {"length_nt", ref.is_null() ? J() : J(str(ref, "sequence").size())},
                               {"cds_status", ref.is_null() ? J() : ref["cds_status"]}};
        result["target"] = {{"length_nt", seq.size()},
                            {"sequence_kind", r["sequence_kind"]},
                            {"cds_status", r["cds_validation"]["status"]},
                            {"reported_shift", r["source_metadata"]["shift_label"]},
                            {"reported_score", r["source_metadata"]["score_label"]}};
        result["comparison"] = comparison;
        out << result.dump() << '\n';
        genes.insert(str(r, "gene_id"));
        inc(methods, method);
        inc(statuses, str(comparison, "status"));
        if (comparison["status"] == "aligned")
            nontriplet += num(comparison, "internal_nontriplet_gap_events");
        J row = copyFields(result, {"gene_id", "species", "method", "locus"});
        row["reference_length_nt"] = result["reference"]["length_nt"];
        row["target_length_nt"] = seq.size();
        row["reference_cds_status"] = result["reference"]["cds_status"];
        row["target_cds_status"] = result["target"]["cds_status"];
        row["reported_shift"] = result["target"]["reported_shift"];
        for (auto &f : schema["comparative_orf"]["SUMMARY_FIELDS"]) {
            auto k = f.get<std::string>();
            if (comparison.contains(k))
                row[k] = comparison[k];
        }
        row["alignment_status"] = comparison["status"];
        row["quality_ok"] = comparison.value("quality_ok", J(false));
        csv.write(row);
    });
    if (!count)
        throw std::runtime_error("Нет сравнений");
    std::set<Key> keys;
    for (auto &[k, v] : profiles)
        keys.insert({k[0], k[1], k[2]});
    Csv acsv(d / "tool_agreement.csv", schema["comparative_orf"]["AGREEMENT_FIELDS"]);
    for (auto &k : keys) {
        auto &a = profiles[{k[0], k[1], k[2], "spaln"}],
             &b = profiles[{k[0], k[1], k[2], "genewise"}];
        J row = {{"gene_id", k[0]},
                 {"species", k[1]},
                 {"locus", k[2]},
                 {"spaln_count", a.size()},
                 {"genewise_count", b.size()}};
        std::string type;
        if (a.size() > 1 || b.size() > 1)
            type = "ambiguous_multiple_hits";
        else if (a.empty())
            type = "missing_spaln";
        else if (b.empty())
            type = "missing_genewise";
        else {
            bool same = a[0]["sha256"] == b[0]["sha256"] && a[0]["length_nt"] == b[0]["length_nt"],
                 status = a[0]["cds_status"] == b[0]["cds_status"];
            type = same     ? "exact_sequence"
                   : status ? "same_cds_status_only"
                            : "different_cds_status";
            for (auto m : {"spaln", "genewise"}) {
                auto &v = std::string(m) == "spaln" ? a[0] : b[0];
                for (auto f : {"length_nt", "cds_status", "shift_label"})
                    row[std::string(m) + "_" + f] = v[f];
            }
            row["same_sequence"] = same;
            row["same_cds_status"] = status;
            row["same_shift_label"] = a[0]["shift_label"] == b[0]["shift_label"];
        }
        row["agreement_type"] = type;
        inc(agreements, type);
        acsv.write(row);
    }
    save(d / "comparison_report.json", {{"stage", "Comparative ORF Analysis Stage 3"},
                                        {"input_file", input.string()},
                                        {"gene_filter", c.s.gene.empty() ? J() : J(c.s.gene)},
                                        {"reference_species", c.s.referenceSpecies},
                                        {"methods", {"spaln", "genewise"}},
                                        {"settings",
                                         {{"max_align_length", c.s.maxAlignLength},
                                          {"max_align_cells", c.s.maxAlignCells},
                                          {"min_identity", c.s.minIdentity},
                                          {"min_coverage", c.s.minCoverage},
                                          {"max_events", c.s.maxEvents},
                                          {"preview_columns", c.s.previewColumns},
                                          {"alignment_mode", "global"},
                                          {"match_score", 2.0},
                                          {"mismatch_score", -1.0},
                                          {"open_gap_score", -5.0},
                                          {"extend_gap_score", -0.5}}},
                                        {"reference_records", refCount},
                                        {"reference_gene_method_pairs", refPairs},
                                        {"processed_comparisons", count},
                                        {"unique_genes", genes.size()},
                                        {"comparisons_per_method", methods},
                                        {"alignment_status_counts", statuses},
                                        {"tool_agreement_counts", agreements},
                                        {"accepted_internal_nontriplet_gap_events", nontriplet},
                                        {"coordinate_system", "0-based_half-open"},
                                        {"limitations", schema["comparative_orf"]["limitations"]}});
}
} // namespace orf
