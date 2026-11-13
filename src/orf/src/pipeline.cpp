#include "common.hpp"
#include <QCoreApplication>
#include <QLockFile>
#include <QSaveFile>
namespace orf {
static void atomicSave(const fs::path &p, const std::string &bytes) {
    fs::create_directories(p.parent_path());
    QSaveFile f(QString::fromStdString(p.string()));
    if (!f.open(QIODevice::WriteOnly) ||
        f.write(bytes.data(), bytes.size()) != (qint64)bytes.size() || !f.commit())
        throw std::runtime_error("Cannot publish " + p.string());
}
J runPipeline(const Settings &settings, Progress progress, const std::atomic_bool *cancel) {
    Settings s = settings;
    if (s.input.empty() || s.outputRoot.empty() || s.finalOut.empty())
        throw std::invalid_argument("Требуются вход, каталог результатов и итоговый файл");
    if (s.mode != "fresh" && s.mode != "resume")
        throw std::invalid_argument("Режим: fresh или resume");
    if (s.minAa < 1 || s.topK < 0 || s.maxPeptideAa < 0 || s.maxAlignLength < 1 ||
        s.maxAlignCells < 1 || s.maxEvents < 0 || s.previewColumns < 0 || s.maxGapDetails < 0 ||
        s.progressEvery < 0 || !std::isfinite(s.minIdentity) || !std::isfinite(s.minCoverage) ||
        s.minIdentity < 0 || s.minIdentity > 1 || s.minCoverage < 0 || s.minCoverage > 1)
        throw std::invalid_argument("Некорректные параметры");
    s.input = fs::weakly_canonical(s.input).string();
    s.outputRoot = fs::weakly_canonical(s.outputRoot).string();
    s.finalOut = fs::weakly_canonical(s.finalOut).string();
    if (!fs::is_regular_file(s.input))
        throw std::runtime_error("Входной файл не найден");
    if (fs::path(s.finalOut).extension() != ".jsonl" && fs::path(s.finalOut).extension() != ".csv")
        throw std::invalid_argument("Итоговый файл: .jsonl или .csv");
    fs::path root = s.outputRoot;
    std::vector<std::vector<std::string>> outputs = {
        {"2_results_first/records.jsonl", "2_results_first/summary.csv",
         "2_results_first/validation.json"},
        {"3_results_second/orf_records.jsonl", "3_results_second/orf_summary.csv",
         "3_results_second/orf_report.json"},
        {"4_results_third/comparisons.jsonl", "4_results_third/comparison_summary.csv",
         "4_results_third/tool_agreement.csv", "4_results_third/comparison_report.json"},
        {"5_results_fourth/evidence_records.jsonl", "5_results_fourth/evidence_summary.csv",
         "5_results_fourth/gene_summary.csv", "5_results_fourth/evidence_report.json"},
        {"6_results_fifth/locus_profiles.csv", "6_results_fifth/species_profiles.csv",
         "6_results_fifth/gene_profiles.jsonl", "6_results_fifth/gene_summary.csv",
         "6_results_fifth/profile_report.json"}};
    std::set<fs::path> reserved = {root / "pipeline_manifest.json", root / "pipeline_run.json",
                                   root / ".pipeline.lock"};
    for (auto &v : outputs)
        for (auto &f : v)
            reserved.insert(fs::weakly_canonical(root / f));
    if (s.input == s.finalOut || reserved.count(s.input) || reserved.count(s.finalOut))
        throw std::invalid_argument("Входной, итоговый и служебные файлы должны различаться");
    fs::create_directories(root);
    QLockFile lock(QString::fromStdString((root / ".pipeline.lock").string()));
    lock.setStaleLockTime(0);
    if (!lock.tryLock(0))
        throw std::runtime_error("Каталог уже используется другим запуском");
    Context c{s, progress, cancel, root};
    J config = {{"gene", s.gene}, {"reference_species", s.referenceSpecies},
                {"min_aa", s.minAa},
                {"top_k", s.topK},
                {"max_peptide_aa", s.maxPeptideAa},
                {"max_align_length", s.maxAlignLength},
                {"max_align_cells", s.maxAlignCells},
                {"min_identity", s.minIdentity},
                {"min_coverage", s.minCoverage},
                {"max_events", s.maxEvents},
                {"preview_columns", s.previewColumns},
                {"max_gap_details", s.maxGapDetails}};
    auto executable = QCoreApplication::applicationFilePath().toStdString();
    std::string signature = hash(config.dump() + fileHash(s.input) + fileHash(executable));
    J manifest = J::object();
    if (s.mode == "resume" && fs::exists(root / "pipeline_manifest.json")) {
        try {
            manifest = load(root / "pipeline_manifest.json");
            if (!manifest.is_object() || !manifest.contains("signature") ||
                !manifest["signature"].is_string() || !manifest.contains("stages") ||
                !manifest["stages"].is_object())
                manifest = J::object();
        } catch (...) {
            manifest = J::object();
        }
    }
    bool reuse = s.mode == "resume" && manifest.value("signature", "") == signature;
    J result = {{"status", "running"},       {"computed_stages", 0}, {"reused_stages", 0},
                {"stages", J::array()},      {"settings", config},   {"input", s.input},
                {"final_output", s.finalOut}};
    if (!reuse)
        manifest = {{"signature", signature}, {"stages", J::object()}};
    std::vector<std::string> names = {"parser", "scanner", "comparative", "evidence", "profile"};
    std::vector<std::function<void(Context &)>> functions = {
        parserStage, scannerStage, comparisonStage, evidenceStage, profileStage};
    try {
        for (size_t i = 0; i < names.size(); ++i) {
            c.check();
            auto name = names[i];
            bool valid =
                reuse && manifest["stages"].contains(name) && manifest["stages"][name].is_object();
            if (valid)
                for (auto &f : outputs[i])
                    if (!fs::is_regular_file(root / f) ||
                        manifest["stages"][name].value(f, "") != fileHash(root / f)) {
                        valid = false;
                        break;
                    }
            if (valid) {
                inc(result, "reused_stages");
                if (progress)
                    progress(name + ": сохранённый результат проверен");
            } else {
                reuse = false;
                if (progress)
                    progress(name + ": вычисление");
                functions[i](c);
                J hashes = J::object();
                for (auto &f : outputs[i])
                    hashes[f] = fileHash(root / f);
                manifest["stages"][name] = hashes;
                atomicSave(root / "pipeline_manifest.json", manifest.dump(2));
                inc(result, "computed_stages");
            }
            result["stages"].push_back(
                {{"stage", name}, {"action", valid ? "reused" : "computed"}});
        }
        c.check();
        auto source = root / (fs::path(s.finalOut).extension() == ".csv"
                                  ? "6_results_fifth/gene_summary.csv"
                                  : "6_results_fifth/gene_profiles.jsonl");
        std::ifstream f(source);
        std::string data((std::istreambuf_iterator<char>(f)), {});
        atomicSave(s.finalOut, data);
        result["status"] = "completed";
        result["final_sha256"] = fileHash(s.finalOut);
        atomicSave(root / "pipeline_run.json", result.dump(2));
        return result;
    } catch (const std::exception &e) {
        result["status"] = cancel && cancel->load() ? "cancelled" : "failed";
        result["error"] = e.what();
        atomicSave(root / "pipeline_run.json", result.dump(2));
        throw;
    }
}
} // namespace orf
