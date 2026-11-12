#include "homology/job.hpp"
namespace homology {
void HomologyJob::start() {
    Hooks hooks;
    hooks.cancelled = &cancelled_;
    hooks.progress = [this](const QString &s) { emit progress(s); };
    try {
        emit succeeded(readJson(config_).value("parallel_targets").toInt(1) > 1
                           ? runBatch(config_, hooks)
                           : run(config_, hooks));
    } catch (const std::exception &e) {
        emit failed(QString::fromUtf8(e.what()));
    }
    emit finished();
}
} // namespace homology
