#pragma once
#include <atomic>
#include <functional>
#include <nlohmann/json.hpp>
#include <string>
namespace orf {
using J = nlohmann::json;
struct Settings {
    std::string input, outputRoot, finalOut, gene, mode = "fresh", referenceSpecies = "dmel";
    int minAa = 25, topK = 10, maxPeptideAa = 120, maxAlignLength = 4000;
    long long maxAlignCells = 4000000;
    double minIdentity = .6, minCoverage = .6;
    int maxEvents = 25, previewColumns = 90, maxGapDetails = 10, progressEvery = 500;
};
using Progress = std::function<void(const std::string &)>;
std::string reverseComplement(std::string s);
J scanOrfs(std::string s, int minAa = 25, int topK = 10, int preview = 120);
J validateCds(std::string s, const std::string &kind);
J compareCds(std::string reference, std::string query, const Settings &settings);
J analyzeGapPath(const J &, int maxDetails = 10);
J classifyLocus(const J &records);
J runPipeline(const Settings &, Progress = {}, const std::atomic_bool *cancel = nullptr);
} // namespace orf
