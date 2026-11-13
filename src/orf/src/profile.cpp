#include "common.hpp"
namespace orf {
static bool valid(const J &r) {
    if (r["profile_match"] != "matched")
        return false;
    auto e = str(r, "evidence_class");
    if (e == "preserved_cds_candidate" || e == "conflicting_signals")
        return r["cds_status"] == "complete_candidate" && r["alignment_status"] == "aligned";
    if (e == "potential_cds_disruption")
        return r["cds_status"] == "interrupted" && r["alignment_status"] == "aligned";
    return true;
}
J classifyLocus(const J &records) {
    if (records.empty())
        throw std::invalid_argument("Empty locus");
    J row = copyFields(records[0], {"gene_id", "species", "locus"});
    std::vector<J> a, b;
    std::set<std::string> agreements;
    int gaps = 0, phase = 0;
    for (auto &r : records) {
        (r["method"] == "spaln" ? a : b).push_back(r);
        agreements.insert(r["tool_agreement"].dump());
        gaps += num(r, "nontriplet_gap_observations");
        phase += num(r, "phase_restored_observations");
    }
    J agreement = agreements.size() == 1 ? J::parse(*agreements.begin()) : J("inconsistent");
    row.update({{"status", nullptr},
                {"reason", nullptr},
                {"tool_agreement", agreement},
                {"spaln_records", a.size()},
                {"genewise_records", b.size()},
                {"spaln_evidence", a.size() == 1 ? a[0]["evidence_class"] : J()},
                {"genewise_evidence", b.size() == 1 ? b[0]["evidence_class"] : J()},
                {"spaln_cds_status", a.size() == 1 ? a[0]["cds_status"] : J()},
                {"genewise_cds_status", b.size() == 1 ? b[0]["cds_status"] : J()},
                {"comparison_records", records.size()},
                {"nontriplet_gap_observations", gaps},
                {"phase_restored_observations", phase}});
    auto finish = [&](const char *s, const char *r) {
        row["status"] = s;
        row["reason"] = r;
        return row;
    };
    if (a.size() > 1 || b.size() > 1)
        return finish("multiple_predictions", "More than one prediction for the same method");
    if (a.empty() || b.empty())
        return finish("missing_method", "Only one prediction method is represented");
    if (!valid(a[0]))
        return finish("inconsistent_metadata", "Inconsistent Spaln evidence");
    if (!valid(b[0]))
        return finish("inconsistent_metadata", "Inconsistent GeneWise evidence");
    if (agreement != "exact_sequence" && agreement != "same_cds_status_only" &&
        agreement != "different_cds_status")
        return finish("inconsistent_metadata", "Missing or inconsistent method agreement");
    bool equal = a[0]["cds_status"] == b[0]["cds_status"];
    if ((agreement == "exact_sequence" || agreement == "same_cds_status_only") && !equal)
        return finish("inconsistent_metadata", "Agreement contradicts CDS statuses");
    if (agreement == "different_cds_status" && equal)
        return finish("inconsistent_metadata", "Different-CDS agreement contradicts statuses");
    auto x = str(a[0], "evidence_class"), y = str(b[0], "evidence_class");
    if (x == "preserved_cds_candidate" && y == x) {
        if (agreement == "exact_sequence")
            return finish("preserved_both_identical", "Both methods predict identical intact CDS");
        if (agreement == "same_cds_status_only")
            return finish("preserved_both_different",
                          "Both predict intact CDS but sequences differ");
        return finish("inconsistent_metadata", "Unexpected agreement for preserved CDS");
    }
    if (x == "potential_cds_disruption" && y == x)
        return finish("predicted_interruption_both", "Internal stops predicted by both methods");
    if ((x == "preserved_cds_candidate" && y == "potential_cds_disruption") ||
        (y == "preserved_cds_candidate" && x == "potential_cds_disruption"))
        return finish("discordant_methods", "Predictions contradict each other");
    if (x == "conflicting_signals" || y == "conflicting_signals")
        return finish("conflicting_alignment", "Nontriplet gap signal conflicts with CDS evidence");
    if (x == "insufficient_evidence" || y == "insufficient_evidence") {
        if (x == "preserved_cds_candidate" || y == "preserved_cds_candidate" ||
            x == "potential_cds_disruption" || y == "potential_cds_disruption")
            return finish("partial_evidence", "Only one method provides usable evidence");
        return finish("insufficient_evidence", "Both methods lack sufficient evidence");
    }
    return finish("insufficient_evidence", "No supported aggregation rule");
}
void profileStage(Context &c) {
    auto d = c.root / "6_results_fifth";
    std::map<Key, J> groups;
    int count = 0;
    jsonl(c.root / "5_results_fourth/evidence_records.jsonl", [&](J r) {
        c.tick(++count);
        J compact = copyFields(r, {"gene_id", "species", "method", "locus", "tool_agreement"});
        compact["evidence_class"] = r["evidence"]["class"];
        compact["cds_status"] = r["cds"]["status"];
        compact["profile_match"] = r["stage2_profile_match"];
        compact["alignment_status"] = r["alignment"]["status"];
        compact["reported_shift"] = r["cds"]["reported_shift"];
        compact["nontriplet_gap_observations"] =
            r["alignment"]["status"] == "aligned"
                ? num(r["alignment"], "internal_nontriplet_gap_events")
                : 0;
        compact["phase_restored_observations"] = int(r["alignment"]["status"] == "aligned" &&
                                                     r["gap_path"]["status"] == "phase_restored");
        auto &g = groups[key(r, {"gene_id", "species", "locus"})];
        if (g.is_null())
            g = J::array();
        g.push_back(compact);
    });
    if (!count)
        throw std::runtime_error("Нет сравнений");
    Csv lcsv(d / "locus_profiles.csv", schema["orf_profile"]["LOCUS_FIELDS"]),
        scsv(d / "species_profiles.csv", schema["orf_profile"]["SPECIES_FIELDS"]),
        gcsv(d / "gene_summary.csv", schema["orf_profile"]["GENE_FIELDS"]);
    std::map<Key, std::vector<J>> species;
    J lc = J::object(), sc = J::object();
    int gaps = 0, phase = 0;
    for (auto &[k, v] : groups) {
        auto row = classifyLocus(v);
        lcsv.write(row);
        inc(lc, str(row, "status"));
        gaps += num(row, "nontriplet_gap_observations");
        phase += num(row, "phase_restored_observations");
        species[{k[0], k[1]}].push_back(row);
    }
    std::map<std::string, std::vector<J>> genes;
    for (auto &[k, rows] : species) {
        J r = {{"gene_id", k[0]},
               {"species", k[1]},
               {"status", rows.size() == 1 ? rows[0]["status"] : J("multiple_loci_uncertain")},
               {"locus_count", rows.size()}};
        std::vector<std::string> loci;
        for (auto &row : rows) {
            for (auto f : {"comparison_records", "nontriplet_gap_observations",
                           "phase_restored_observations"})
                inc(r, f, num(row, f));
            loci.push_back(str(row, "locus"));
        }
        std::sort(loci.begin(), loci.end());
        std::string joined;
        for (auto &l : loci) {
            if (!joined.empty())
                joined += ';';
            joined += l;
        }
        r["loci"] = joined;
        scsv.write(r);
        inc(sc, str(r, "status"));
        genes[k[0]].push_back(r);
    }
    auto out = output(d / "gene_profiles.jsonl");
    for (auto &[id, rows] : genes) {
        J p = {
            {"gene_id", id},
            {"species_assessed", rows.size()},
            {"loci_assessed", 0},
            {"comparison_records", 0},
            {"species_status_counts", J::object()},
            {"nontriplet_gap_observations", 0},
            {"phase_restored_observations", 0},
            {"species_profiles", J::array()},
            {"origin_classification", nullptr},
            {"limitations",
             {"Only observed loci are evaluated.", "Missing species do not imply gene absence.",
              "Multiple loci are not resolved.", "Spaln and GeneWise are correlated predictions.",
              "Gap observations are not verified mutations.",
              "Ancestral coding status is not inferred."}}};
        for (auto &r : rows) {
            inc(p, "loci_assessed", num(r, "locus_count"));
            for (auto f : {"comparison_records", "nontriplet_gap_observations",
                           "phase_restored_observations"})
                inc(p, f, num(r, f));
            inc(p["species_status_counts"], str(r, "status"));
            std::vector<std::string> loci;
            std::string text = str(r, "loci");
            size_t a = 0, b;
            while ((b = text.find(';', a)) != text.npos) {
                loci.push_back(text.substr(a, b - a));
                a = b + 1;
            }
            loci.push_back(text.substr(a));
            p["species_profiles"].push_back({{"species", r["species"]},
                                             {"status", r["status"]},
                                             {"locus_count", r["locus_count"]},
                                             {"loci", loci}});
        }
        out << p.dump() << '\n';
        J row = p;
        for (auto &f : schema["orf_profile"]["GENE_FIELDS"]) {
            auto k = f.get<std::string>();
            if (k.rfind("species_", 0) == 0 && k != "species_assessed")
                row[k] = num(p["species_status_counts"], k.substr(8));
        }
        gcsv.write(row);
    }
    save(d / "profile_report.json",
         {{"stage", "ORF Profile Aggregation Stage 5"},
          {"input_file", (c.root / "5_results_fourth/evidence_records.jsonl").string()},
          {"gene_filter", c.s.gene.empty() ? J() : J(c.s.gene)},
          {"comparison_records", count},
          {"unique_loci", groups.size()},
          {"gene_species_pairs", species.size()},
          {"unique_genes", genes.size()},
          {"locus_status_counts", lc},
          {"species_status_counts", sc},
          {"nontriplet_gap_observations", gaps},
          {"phase_restored_observations", phase},
          {"origin_classification_performed", false},
          {"limitations", schema["orf_profile"]["limitations"]}});
}
} // namespace orf
