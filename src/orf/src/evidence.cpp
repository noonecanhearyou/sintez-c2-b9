#include "common.hpp"
namespace orf {
J analyzeGapPath(const J &c, int maxDetails) {
    J r = {{"status", "not_available"},
           {"events_complete", false},
           {"phase_final_mod3", nullptr},
           {"phase_restorations", nullptr},
           {"phase_excursions", nullptr},
           {"codon_boundary_gap_events", nullptr},
           {"off_boundary_gap_events", nullptr},
           {"internal_nontriplet_observed", nullptr},
           {"details", J::array()}};
    auto finish = [&](const char *s) {
        r["status"] = s;
        return r;
    };
    if (str(c, "status") != "aligned" && str(c, "status") != "low_quality_alignment")
        return r;
    if (!c.contains("gap_examples") || !c["gap_examples"].is_array() ||
        !c.contains("gap_event_count") || !c["gap_event_count"].is_number_integer() ||
        num(c, "gap_event_count") < 0)
        return finish("missing_gap_data");
    auto events = c["gap_examples"];
    int total = num(c, "gap_event_count");
    if (c.value("gap_examples_truncated", J()) == true || (int)events.size() != total)
        return finish("incomplete_gap_examples");
    int phase = 0, rest = 0, exc = 0, bound = 0, off = 0, non = 0, trip = 0, terminal = 0;
    J details = J::array();
    for (auto &e : events) {
        if (!e.is_object())
            return finish("invalid_gap_data");
        auto kind = str(e, "kind");
        if ((kind != "gap_in_query" && kind != "gap_in_reference") || !e.contains("length_nt") ||
            !e["length_nt"].is_number_integer() || num(e, "length_nt") <= 0 ||
            !e.contains("internal") || !e["internal"].is_boolean())
            return finish("invalid_gap_data");
        int len = num(e, "length_nt");
        if (e["internal"] == false) {
            ++terminal;
            continue;
        }
        auto rk = kind == "gap_in_query" ? "reference_start" : "reference_position",
             qk = kind == "gap_in_query" ? "query_position" : "query_start";
        if (!e.contains(rk) || !e[rk].is_number_integer() || !e.contains(qk) ||
            !e[qk].is_number_integer() || num(e, rk) < 0 || num(e, qk) < 0)
            return finish("invalid_gap_data");
        int rp = num(e, rk), qp = num(e, qk);
        bool at = rp % 3 == 0 && qp % 3 == 0;
        if (at)
            ++bound;
        else
            ++off;
        if (len % 3)
            ++non;
        else
            ++trip;
        int before = phase;
        phase = ((phase + (kind == "gap_in_query" ? -len : len)) % 3 + 3) % 3;
        if (!before && phase)
            ++exc;
        if (before && !phase)
            ++rest;
        if ((int)details.size() < maxDetails)
            details.push_back({{"kind", kind},
                               {"length_nt", len},
                               {"reference_position", rp},
                               {"query_position", qp},
                               {"at_codon_boundary", at},
                               {"phase_before", before},
                               {"phase_after", phase}});
    }
    if (c.value("internal_nontriplet_gap_events", J()) != non ||
        c.value("internal_triplet_gap_events", J()) != trip ||
        c.value("terminal_gap_events", J()) != terminal)
        return finish("inconsistent_gap_counts");
    r.update({{"status", !total   ? "no_gaps"
                         : !non   ? "in_frame_gaps_only"
                         : !phase ? "phase_restored"
                                  : "phase_unresolved"},
              {"events_complete", true},
              {"phase_final_mod3", phase},
              {"phase_restorations", rest},
              {"phase_excursions", exc},
              {"codon_boundary_gap_events", bound},
              {"off_boundary_gap_events", off},
              {"internal_nontriplet_observed", non},
              {"details", details}});
    return r;
}
static std::vector<std::string> parseCsvLine(const std::string &s) {
    std::vector<std::string> v;
    std::string t;
    bool quoted = false;
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '"') {
            if (quoted && i + 1 < s.size() && s[i + 1] == '"') {
                t += '"';
                ++i;
            } else
                quoted = !quoted;
        } else if (c == ',' && !quoted) {
            v.push_back(t);
            t.clear();
        } else if (c != '\r')
            t += c;
    }
    v.push_back(t);
    return v;
}
void evidenceStage(Context &c) {
    auto d = c.root / "5_results_fourth";
    std::map<Key, J> profiles;
    std::map<Key, int> counts;
    int profileCount = 0;
    jsonl(c.root / "3_results_second/orf_records.jsonl", [&](J r) {
        c.check();
        if (r["section"] == "raw")
            return;
        ++profileCount;
        auto k = key(r, {"gene_id", "species", "section", "locus"});
        ++counts[k];
        profiles[k] = {
            {"length_nt", str(r, "sequence").size()},
            {"cds_status", r["cds_validation"]["status"]},
            {"internal_stop_count", r["cds_validation"].value("internal_stop_count", J())},
            {"longest_complete_aa", r["orf_scan"]["longest_complete_aa"]},
            {"reported_shift", r["source_metadata"]["shift_label"]}};
    });
    std::map<Key, std::string> agreements;
    std::ifstream af(c.root / "4_results_third/tool_agreement.csv");
    std::string line;
    std::getline(af, line);
    auto header = parseCsvLine(line);
    while (std::getline(af, line)) {
        auto cells = parseCsvLine(line);
        J r;
        for (size_t i = 0; i < header.size() && i < cells.size(); ++i)
            r[header[i]] = cells[i];
        auto k = key(r, {"gene_id", "species", "locus"});
        agreements[k] = agreements.count(k) ? "duplicate_agreement_rows" : str(r, "agreement_type");
    }
    auto out = output(d / "evidence_records.jsonl");
    Csv csv(d / "evidence_summary.csv", schema["orf_evidence"]["SUMMARY_FIELDS"]);
    J evidenceCounts = J::object(), alignCounts = J::object(), gapCounts = J::object(),
      matchCounts = J::object();
    struct Gene {
        int comparisons = 0, non = 0;
        std::set<std::string> species;
        std::set<Key> loci, exact;
        J classes = J::object();
    };
    std::map<std::string, Gene> genes;
    int count = 0, non = 0;
    jsonl(c.root / "4_results_third/comparisons.jsonl", [&](J r) {
        c.tick(++count);
        auto k = key(r, {"gene_id", "species", "method", "locus"});
        J p = profiles.count(k) ? profiles[k] : J::object();
        auto t = r["target"], comp = r["comparison"];
        std::string match = !counts.count(k) ? "missing_stage2_profile"
                            : counts[k] != 1 ? "ambiguous_multiple_profiles"
                                             : "matched";
        if (match == "matched") {
            if (p["length_nt"] != t.value("length_nt", J()))
                match = "length_mismatch";
            else if (p["cds_status"] != t.value("cds_status", J()))
                match = "cds_status_mismatch";
            else if (p["reported_shift"] != t.value("reported_shift", J()))
                match = "shift_label_mismatch";
        }
        if (counts[k] != 1)
            p = J::object();
        auto ak = key(r, {"gene_id", "species", "locus"});
        std::string agreement = agreements.count(ak) ? agreements[ak] : "missing_agreement";
        auto gap = analyzeGapPath(comp, c.s.maxGapDetails);
        std::string cls = "insufficient_evidence", reason;
        auto status = str(t, "cds_status"), gs = str(gap, "status");
        if (match != "matched")
            reason = "stage2_stage3_data_mismatch";
        else if (comp["status"] != "aligned" || comp.value("quality_ok", J()) != true)
            reason = "alignment_did_not_pass_quality_thresholds";
        else if (status == "interrupted") {
            cls = "potential_cds_disruption";
            reason = "internal_stop_in_predicted_cds";
        } else if (status == "complete_candidate") {
            if (gs == "no_gaps" || gs == "in_frame_gaps_only") {
                cls = "preserved_cds_candidate";
                reason = "intact_cds_without_nontriplet_internal_gaps";
            } else if (gs == "phase_restored" || gs == "phase_unresolved") {
                cls = "conflicting_signals";
                reason = "intact_cds_but_nontriplet_alignment_gaps";
            } else
                reason = "gap_path_not_fully_evaluable";
        } else
            reason = "target_cds_status_" + (status.empty() ? "None" : status);
        J result = copyFields(r, {"gene_id", "species", "method", "locus"});
        result["stage2_profile_match"] = match;
        result["cds"] = {{"status", t.value("cds_status", J())},
                         {"internal_stop_count", p.value("internal_stop_count", J())},
                         {"longest_complete_aa", p.value("longest_complete_aa", J())},
                         {"reported_shift", t.value("reported_shift", J())}};
        result["tool_agreement"] = agreement;
        result["alignment"] =
            copyFields(comp, {"status", "identity", "reference_coverage", "query_coverage",
                              "internal_nontriplet_gap_events"});
        result["gap_path"] = gap;
        result["evidence"] = {{"class", cls},
                              {"reason", reason},
                              {"evolutionary_direction_inferred", false},
                              {"confirmed_frameshift", false}};
        out << result.dump() << '\n';
        J row = copyFields(r, {"gene_id", "species", "method", "locus"});
        row.update(result["alignment"]);
        row.update(copyFields(result["cds"],
                              {"internal_stop_count", "longest_complete_aa", "reported_shift"}));
        row.update(copyFields(gap, {"phase_final_mod3", "phase_restorations",
                                    "codon_boundary_gap_events", "off_boundary_gap_events"}));
        row.update({{"profile_match", match},
                    {"cds_status", t["cds_status"]},
                    {"tool_agreement", agreement},
                    {"alignment_status", comp["status"]},
                    {"gap_path_status", gs},
                    {"evidence_class", cls},
                    {"reason", reason}});
        csv.write(row);
        inc(evidenceCounts, cls);
        inc(alignCounts, str(comp, "status"));
        inc(gapCounts, gs);
        inc(matchCounts, match);
        int ng = comp["status"] == "aligned" ? num(comp, "internal_nontriplet_gap_events") : 0;
        non += ng;
        auto &g = genes[str(r, "gene_id")];
        ++g.comparisons;
        g.species.insert(str(r, "species"));
        g.loci.insert({str(r, "species"), str(r, "locus")});
        if (agreement == "exact_sequence")
            g.exact.insert({str(r, "species"), str(r, "locus")});
        inc(g.classes, cls);
        g.non += ng;
    });
    if (!count)
        throw std::runtime_error("Нет сравнений");
    Csv gcsv(d / "gene_summary.csv", schema["orf_evidence"]["GENE_FIELDS"]);
    for (auto &[id, g] : genes) {
        J row = {{"gene_id", id},
                 {"comparisons", g.comparisons},
                 {"unique_species", g.species.size()},
                 {"unique_loci", g.loci.size()},
                 {"exact_agreement_loci", g.exact.size()},
                 {"nontriplet_event_observations", g.non}};
        for (auto cl : {"preserved_cds_candidate", "potential_cds_disruption",
                        "conflicting_signals", "insufficient_evidence"})
            row[cl] = num(g.classes, cl);
        gcsv.write(row);
    }
    save(d / "evidence_report.json",
         {{"stage", "ORF Evidence Analysis Stage 4"},
          {"orf_input", (c.root / "3_results_second/orf_records.jsonl").string()},
          {"comparisons_input", (c.root / "4_results_third/comparisons.jsonl").string()},
          {"agreement_input", (c.root / "4_results_third/tool_agreement.csv").string()},
          {"gene_filter", c.s.gene.empty() ? J() : J(c.s.gene)},
          {"processed_comparisons", count},
          {"unique_genes", genes.size()},
          {"evidence_class_counts", evidenceCounts},
          {"alignment_status_counts", alignCounts},
          {"gap_path_status_counts", gapCounts},
          {"profile_match_counts", matchCounts},
          {"nontriplet_event_observations", non},
          {"settings", {{"max_gap_details", c.s.maxGapDetails}}},
          {"limitations", schema["orf_evidence"]["limitations"]}});
}
} // namespace orf
