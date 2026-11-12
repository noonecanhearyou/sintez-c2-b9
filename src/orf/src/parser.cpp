#include "common.hpp"
#include <regex>
namespace orf {
static J parseRecord(const std::string &section, const std::string &gene, const std::string &header,
                     int line, const std::string &sequence) {
    std::vector<std::string> f;
    size_t a = 1, b;
    while ((b = header.find('|', a)) != header.npos) {
        f.push_back(header.substr(a, b - a));
        a = b + 1;
    }
    f.push_back(header.substr(a));
    J issues = J::array();
    if (f.size() != 2 && f.size() != 5)
        issues.push_back("header_wrong_field_count");
    if (section != "raw" && f.size() != 5)
        issues.push_back("prediction_missing_metadata");
    std::string species = f[0], locus = f.size() > 1 ? f[1] : "";
    if (species.empty() || locus.empty())
        issues.push_back("header_missing_species_or_locus");
    std::smatch m;
    bool loc = std::regex_match(locus, m, std::regex(R"(^(.+)([+-]):(\d+):(\d+)$)"));
    J contig, strand, start, end, score, shift, exons;
    if (loc) {
        contig = m[1].str();
        strand = m[2].str();
        start = std::stoll(m[3].str());
        end = std::stoll(m[4].str());
        if (start.get<long long>() >= end.get<long long>())
            issues.push_back("genomic_locus_coordinates_invalid");
    }
    std::string kind = "unknown";
    if (section == "raw" && f.size() == 2) {
        kind = "genomic_region";
        if (!loc)
            issues.push_back("genomic_locus_format_invalid");
    } else if (f.size() == 5)
        kind = loc ? "predicted_cds" : "annotated_cds";
    if (f.size() == 5) {
        if (f[2].rfind("score_", 0) == 0) {
            score = f[2].substr(6);
            try {
                size_t n;
                std::stod(score.get<std::string>(), &n);
                if (n != score.get<std::string>().size())
                    throw std::invalid_argument("score");
            } catch (...) {
                issues.push_back("score_not_numeric");
            }
        } else
            issues.push_back("score_field_invalid");
        if (f[3].rfind("shift_", 0) == 0) {
            shift = f[3].substr(6);
            if (shift != "Yes" && shift != "No")
                issues.push_back("shift_field_invalid");
        } else
            issues.push_back("shift_field_invalid");
        try {
            size_t n;
            auto v = std::stoll(f[4], &n);
            if (n != f[4].size())
                throw std::invalid_argument("exons");
            exons = v;
            if (v < 1)
                issues.push_back("exon_count_invalid");
        } catch (...) {
            issues.push_back("exon_count_invalid");
        }
    }
    auto s = upper(sequence);
    if (s.empty())
        issues.push_back(section == "raw" ? "empty_sequence" : "empty_prediction");
    if (invalid(s))
        issues.push_back("invalid_dna_characters");
    return {{"gene_id", gene},
            {"section", section},
            {"species", species},
            {"locus", locus},
            {"sequence_kind", kind},
            {"sequence", s},
            {"length_nt", s.size()},
            {"ambiguous_nt", ambiguous(s)},
            {"invalid_nt", invalid(s)},
            {"score_label", score},
            {"shift_label", shift},
            {"exon_count", exons},
            {"contig", contig},
            {"strand", strand},
            {"start", start},
            {"end", end},
            {"original_header", header},
            {"header_line", line},
            {"issues", issues}};
}
void parserStage(Context &c) {
    auto d = c.root / "2_results_first";
    auto out = output(d / "records.jsonl");
    Csv csv(d / "summary.csv", schema["orf_parser"]["CSV_FIELDS"]);
    std::ifstream in(c.s.input);
    if (!in)
        throw std::runtime_error("Входной файл не найден");
    std::string section, gene, header, seq, line;
    int lineno = 0, headerLine = 0, count = 0, groups = 0, ambRecords = 0, ambBases = 0;
    J sections = J::object(), counts = J::object(), warnings = J::object(), errors = J::object(),
      we = J::array(), ee = J::array();
    std::map<std::string, std::map<std::string, int>> geneGroups;
    std::map<std::string, std::map<std::string, J>> hits;
    auto flush = [&] {
        if (header.empty())
            return;
        J r = parseRecord(section, gene, header, headerLine, seq);
        header.clear();
        seq.clear();
        if (!c.s.gene.empty() && gene != c.s.gene)
            return;
        c.tick(++count);
        inc(counts, section);
        if (!hits[gene].count(section))
            hits[gene][section] = J::object();
        inc(hits[gene][section], str(r, "species"));
        if (num(r, "ambiguous_nt")) {
            ++ambRecords;
            ambBases += num(r, "ambiguous_nt");
        }
        for (auto &i : r["issues"]) {
            bool warn = i == "empty_prediction";
            inc(warn ? warnings : errors, i.get<std::string>());
            auto &examples = warn ? we : ee;
            if (examples.size() < 15)
                examples.push_back({{"line", headerLine},
                                    {"gene", gene},
                                    {"section", section},
                                    {"species", r["species"]},
                                    {"issue", i}});
        }
        out << r.dump() << '\n';
        std::string joined;
        for (auto &i : r["issues"]) {
            if (!joined.empty())
                joined += ';';
            joined += i.get<std::string>();
        }
        r["issues"] = joined;
        csv.write(r);
    };
    std::regex sr(R"(^Section\s+(I|II|III)\.\s+)", std::regex::icase),
        gr(R"(^(raw sequence|spaln|genewise) for dmel_(\S+)$)", std::regex::icase);
    std::smatch m;
    while (std::getline(in, line)) {
        ++lineno;
        c.check();
        line = trim(line);
        if (line.empty())
            continue;
        if (line.rfind("Section ", 0) == 0) {
            flush();
            if (!std::regex_search(line, m, sr))
                throw std::runtime_error("Неизвестная секция, строка " + std::to_string(lineno));
            auto roman = upper(m[1].str());
            section = roman == "I" ? "raw" : roman == "II" ? "spaln" : "genewise";
            gene.clear();
            inc(sections, section);
        } else if (std::regex_match(line, m, gr)) {
            if (section.empty())
                throw std::runtime_error("Группа вне секции");
            flush();
            auto group = upper(m[1].str());
            if (group != (section == "raw" ? "RAW SEQUENCE" : upper(section)))
                throw std::runtime_error("Неверная группа");
            gene = m[2];
            ++groups;
            if (c.s.gene.empty() || gene == c.s.gene) {
                ++geneGroups[section][gene];
                hits[gene][section];
            }
        } else if (line[0] == '>') {
            if (section.empty() || gene.empty())
                throw std::runtime_error("FASTA вне группы");
            flush();
            header = line;
            headerLine = lineno;
        } else {
            if (header.empty())
                throw std::runtime_error("Неожиданный текст, строка " + std::to_string(lineno));
            seq += line;
        }
    }
    if (in.bad())
        throw std::runtime_error("Ошибка чтения");
    flush();
    auto addError = [&](const std::string &type, const std::string &g, const std::string &s) {
        inc(errors, type);
        if (ee.size() < 15)
            ee.push_back({{"gene", g}, {"section", s}, {"issue", type}});
    };
    for (auto &[g, by] : hits) {
        for (auto s : {"raw", "spaln", "genewise"})
            if (!by.count(s))
                addError("gene_missing_in_section", g, s);
        if (by.count("raw") && by.count("spaln") && by.count("genewise"))
            for (auto s : {"spaln", "genewise"})
                if (by["raw"] != by[s])
                    addError("species_hit_count_mismatch", g, s);
    }
    J genesPer = J::object();
    for (auto s : {"raw", "spaln", "genewise"}) {
        if (num(sections, s) != 1)
            inc(errors, "missing_or_repeated_section");
        sections[s] = num(sections, s);
        counts[s] = num(counts, s);
        genesPer[s] = geneGroups[s].size();
        for (auto &[g, n] : geneGroups[s])
            if (n > 1)
                inc(errors, "repeated_gene_group");
    }
    if (!c.s.gene.empty() && !hits.count(c.s.gene))
        inc(errors, "requested_gene_not_found");
    int ne = 0, nw = 0;
    for (auto &v : errors)
        ne += v.get<int>();
    for (auto &v : warnings)
        nw += v.get<int>();
    save(d / "validation.json", {{"input_file", c.s.input},
                                 {"gene_filter", c.s.gene.empty() ? J() : J(c.s.gene)},
                                 {"sections_seen", sections},
                                 {"genes_per_section", genesPer},
                                 {"unique_genes", hits.size()},
                                 {"records_per_section", counts},
                                 {"selected_records", count},
                                 {"groups_scanned", groups},
                                 {"records_with_iupac_ambiguity", ambRecords},
                                 {"total_ambiguous_bases", ambBases},
                                 {"warning_count", nw},
                                 {"warning_counts", warnings},
                                 {"warning_examples", we},
                                 {"error_count", ne},
                                 {"error_counts", errors},
                                 {"error_examples", ee}});
    if (ne)
        throw std::runtime_error("Парсер обнаружил ошибки: " + std::to_string(ne) +
                                 ". См. validation.json");
}
} // namespace orf
