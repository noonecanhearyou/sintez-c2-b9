#pragma once
#include "homology/homology.hpp"
#include <QObject>
namespace homology {
// Construct without parent, moveToThread(workerThread), then call start there.
// requestCancel() may be called directly from any thread; it only sets an atomic flag.
class HomologyJob final : public QObject {
    Q_OBJECT
  public:
    explicit HomologyJob(QString configPath, QObject *parent = nullptr)
        : QObject(parent), config_(std::move(configPath)) {}
    void requestCancel() noexcept {
        cancelled_.store(true);
    }
  public slots:
    void start();
  signals:
    void progress(QString message);
    void succeeded(QJsonObject runManifest);
    void failed(QString message);
    void finished();

  private:
    QString config_;
    std::atomic_bool cancelled_{false};
};
} // namespace homology
