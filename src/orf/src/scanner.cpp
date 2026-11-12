#include "common.hpp"
#include <queue>
#include <tuple>
namespace orf {
std::string reverseComplement(std::string s) {
    s = upper(s);
    std::string a = "ACGTRYSWKMBDHVN", b = "TGCAYRSWMKVHDBN";
    for (char &c : s) {
        auto p = a.find(c);
        if (p != a.npos)
            c = b[p];
    }
    std::reverse(s.begin(), s.end());
    return s;
}
static std::string translate(const std::string &s) {
    std::string aa;
    const std::string bases = "TCAG",
                      code = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    for (size_t p = 0; p + 3 <= s.size(); p += 3) {
        int i = 0;
        for (int k = 0; k < 3; ++k)
            i = i * 4 + bases.find(s[p + k]);
        aa += code[i];
    }
    return aa;
}
J scanOrfs(std::string s, int minAa, int topK, int preview) {
    if (minAa < 1 || topK < 0 || preview < 0)
        throw std::invalid_argument("Invalid scan settings");
    s = upper(s);
    int n = s.size(), serial = 0, longest = 0;
    J all = zeros(), passed = zeros(), frames = J::object();
    using Priority = std::tuple<int, int, int, int>;
    struct Candidate {
        Priority priority;
        int serial;
        J data;
    };
    auto cmp = [](const Candidate &a, const Candidate &b) {
        return std::tie(a.priority, a.serial) > std::tie(b.priority, b.serial);
    };
    std::priority_queue<Candidate, std::vector<Candidate>, decltype(cmp)> heap(cmp);
    for (char strand : {'+', '-'}) {
        auto seq = strand == '+' ? s : reverseComplement(s);
        for (int offset = 0; offset < 3; ++offset) {
            std::string frame = std::string(1, strand) + std::to_string(offset + 1);
            frames[frame] = zeros();
            std::vector<int> starts;
            int p = offset;
            auto reg = [&](int start, int end, const std::string &type, J stopCodon) {
                int aa = (end - start) / 3 - (type == "complete");
                if (aa < 0)
                    return;
                inc(all, type);
                if (type == "complete")
                    longest = std::max(longest, aa);
                if (aa < minAa)
                    return;
                inc(passed, type);
                inc(frames[frame], type);
                int a = strand == '+' ? start : n - end, b = strand == '+' ? end : n - start;
                J c = {{"frame", frame},
                       {"strand_relative_to_input", std::string(1, strand)},
                       {"type", type},
                       {"start_nt", a},
                       {"end_nt", b},
                       {"length_nt", end - start},
                       {"aa_length", aa},
                       {"start_codon", "ATG"},
                       {"stop_codon", stopCodon}};
                if (!topK)
                    return;
                Candidate item{{aa, type == "complete", -a, -b}, ++serial, c};
                if ((int)heap.size() < topK)
                    heap.push(item);
                else if (item.priority > heap.top().priority) {
                    heap.pop();
                    heap.push(item);
                }
            };
            for (; p + 3 <= n; p += 3) {
                auto codon = seq.substr(p, 3);
                if (codon.find_first_not_of("ACGT") != codon.npos) {
                    for (int a : starts)
                        reg(a, p, "ambiguous_boundary", nullptr);
                    starts.clear();
                } else if (stop(codon)) {
                    for (int a : starts)
                        reg(a, p + 3, "complete", codon);
                    starts.clear();
                } else if (codon == "ATG")
                    starts.push_back(p);
            }
            for (int a : starts)
                reg(a, p, "open_end", nullptr);
        }
    }
    std::vector<Candidate> sorted;
    while (!heap.empty()) {
        sorted.push_back(heap.top());
        heap.pop();
    }
    std::reverse(sorted.begin(), sorted.end());
    J top = J::array();
    for (auto &c : sorted) {
        auto &v = c.data;
        auto dna = s.substr(num(v, "start_nt"), num(v, "length_nt"));
        if (v["strand_relative_to_input"] == "-")
            dna = reverseComplement(dna);
        auto pep = translate(dna);
        if (v["type"] == "complete" && !pep.empty())
            pep.pop_back();
        v["peptide_preview"] = pep.substr(0, preview);
        v["peptide_truncated"] = (int)pep.size() > preview;
        top.push_back(v);
    }
    return {{"min_aa", minAa},
            {"counts_all", all},
            {"counts_passing_min_aa", passed},
            {"counts_by_frame_passing_min_aa", frames},
            {"longest_complete_aa", longest},
            {"top_orfs", top}};
}
J validateCds(std::string s, const std::string &kind) {
    if (kind != "annotated_cds" && kind != "predicted_cds")
        return {{"status", "not_applicable"}, {"reason_codes", J::array()}};
    s = upper(s);
    int n = s.size(), rem = n % 3, amb = 0;
    bool start = n >= 3 && s.substr(0, 3) == "ATG",
         terminal = rem == 0 && n >= 3 && stop(s.substr(n - 3));
    J positions = J::array(), reasons = J::array();
    for (int p = 0; p + 3 <= n; p += 3) {
        auto codon = s.substr(p, 3);
        if (codon.find_first_not_of("ACGT") != codon.npos)
            ++amb;
        if (stop(codon) && !(terminal && p == n - 3))
            positions.push_back(p);
    }
    int bad = invalid(s);
    if (!n)
        return {{"status", "empty"},
                {"reason_codes", {"empty_cds"}},
                {"start_atg", false},
                {"terminal_stop", false},
                {"length_mod3", 0},
                {"internal_stop_count", 0},
                {"internal_stop_positions", J::array()},
                {"ambiguous_codon_count", 0},
                {"invalid_nt", 0}};
    if (rem)
        reasons.push_back("length_not_multiple_of_3");
    if (!start)
        reasons.push_back("start_not_atg");
    if (!terminal)
        reasons.push_back("terminal_stop_missing");
    if (!positions.empty())
        reasons.push_back("internal_stop_codons");
    if (amb)
        reasons.push_back("ambiguous_codons");
    if (bad)
        reasons.push_back("invalid_characters");
    int count = positions.size();
    while (positions.size() > 50)
        positions.erase(positions.end() - 1);
    return {{"status", count                        ? "interrupted"
                       : bad || amb                 ? "uncertain"
                       : rem || !start || !terminal ? "partial_or_noncanonical"
                                                    : "complete_candidate"},
            {"reason_codes", reasons},
            {"start_atg", start},
            {"terminal_stop", terminal},
            {"length_mod3", rem},
            {"internal_stop_count", count},
            {"internal_stop_positions", positions},
            {"internal_stop_positions_truncated", count > 50},
            {"ambiguous_codon_count", amb},
            {"invalid_nt", bad}};
}
void scannerStage(Context &c) {
    auto d = c.root / "3_results_second";
    auto out = output(d / "orf_records.jsonl");
    Csv csv(d / "orf_summary.csv", schema["orf_scanner"]["SUMMARY_FIELDS"]);
    J sections = J::object(), statuses = J::object();
    std::set<std::string> genes;
    int count = 0, withComplete = 0, withAmb = 0, withBad = 0, total = 0, passing = 0;
    jsonl(c.root / "2_results_first/records.jsonl", [&](J r) {
        c.tick(++count);
        auto s = upper(str(r, "sequence"));
        r["sequence"] = s;
        auto scan = scanOrfs(s, c.s.minAa, c.s.topK, c.s.maxPeptideAa),
             cds = validateCds(s, str(r, "sequence_kind"));
        J result =
            copyFields(r, {"gene_id", "species", "section", "locus", "sequence_kind", "sequence"});
        result["source_metadata"] =
            copyFields(r, {"original_header", "header_line", "score_label", "shift_label",
                           "exon_count", "strand", "start", "end", "issues"});
        result["orf_scan"] = scan;
        result["cds_validation"] = cds;
        out << result.dump() << '\n';
        genes.insert(str(r, "gene_id"));
        inc(sections, str(r, "section"));
        inc(statuses, str(cds, "status"));
        int a = num(scan["counts_all"], "complete"),
            b = num(scan["counts_passing_min_aa"], "complete");
        total += a;
        passing += b;
        withComplete += b > 0;
        withAmb += ambiguous(s) > 0;
        withBad += invalid(s) > 0;
        J row = copyFields(r, {"gene_id", "species", "section", "locus", "sequence_kind"});
        row.update(
            {{"length_nt", s.size()},
             {"ambiguous_nt", ambiguous(s)},
             {"invalid_nt", invalid(s)},
             {"reported_shift", r["shift_label"]},
             {"reported_score", r["score_label"]},
             {"cds_status", cds["status"]},
             {"complete_orfs_all", a},
             {"complete_orfs_min_aa", b},
             {"open_end_orfs_min_aa", scan["counts_passing_min_aa"]["open_end"]},
             {"ambiguous_boundary_min_aa", scan["counts_passing_min_aa"]["ambiguous_boundary"]},
             {"longest_complete_aa", scan["longest_complete_aa"]}});
        for (auto p : {std::pair{"start_atg", "cds_start_atg"},
                       {"terminal_stop", "cds_terminal_stop"},
                       {"length_mod3", "cds_length_mod3"},
                       {"internal_stop_count", "cds_internal_stop_count"},
                       {"ambiguous_codon_count", "cds_ambiguous_codons"}})
            row[p.second] = cds.value(p.first, J());
        J t = scan["top_orfs"].empty() ? J::object() : scan["top_orfs"][0];
        for (auto p : {std::pair{"frame", "top_orf_frame"},
                       {"type", "top_orf_type"},
                       {"aa_length", "top_orf_aa"},
                       {"start_nt", "top_orf_start_nt"},
                       {"end_nt", "top_orf_end_nt"}})
            row[p.second] = t.value(p.first, J());
        csv.write(row);
    });
    if (!count)
        throw std::runtime_error("Нет записей");
    save(d / "orf_report.json",
         {{"stage", "ORF Analysis Stage 2"},
          {"input_file", (c.root / "2_results_first/records.jsonl").string()},
          {"gene_filter", c.s.gene.empty() ? J() : J(c.s.gene)},
          {"settings",
           {{"min_aa", c.s.minAa},
            {"top_k", c.s.topK},
            {"max_peptide_aa", c.s.maxPeptideAa},
            {"start_codon", "ATG"},
            {"stop_codons", {"TAA", "TAG", "TGA"}},
            {"genetic_code", 1}}},
          {"scanned_jsonl_records", count},
          {"processed_records", count},
          {"unique_genes", genes.size()},
          {"records_per_section", sections},
          {"cds_status_counts", statuses},
          {"records_with_complete_orf_min_aa", withComplete},
          {"total_complete_orfs_all_lengths", total},
          {"total_complete_orfs_min_aa", passing},
          {"records_with_ambiguous_bases", withAmb},
          {"records_with_invalid_bases", withBad},
          {"coordinate_system", "0-based_half-open"},
          {"strand_definition", "relative_to_input_sequence"}});
}
} // namespace orf
