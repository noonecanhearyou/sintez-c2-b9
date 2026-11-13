#pragma once
#include "core.hpp"
#include "tree.hpp"
#include <QCheckBox>
#include <QComboBox>
#include <QFutureWatcher>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMainWindow>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QStackedWidget>
#include <QTableWidget>
#include <QTextBrowser>
#include <QToolButton>
#include <atomic>
class GenePicker final : public QToolButton {
    Q_OBJECT
  public:
    explicit GenePicker(QWidget *parent = nullptr);
    void setGenes(QStringList genes);
    QStringList selected() const;
    void selectGenes(QStringList genes);
    QStringList all() const {
        return genes_;
    }
  signals:
    void selectionChanged();

  private:
    void updateLabel();
    QStringList genes_;
    QSet<QString> selected_;
    QFrame *popup_;
    QLineEdit *search_;
    QListWidget *list_;
};
class MainWindow final : public QMainWindow {
    Q_OBJECT
  public:
    explicit MainWindow(QString workspace = QString::fromUtf8(GENEORIGIN_WORKSPACE));
    void loadProjectInputs();
    void loadSaved(QString directory);
    void startAnalysis();
    void selectGenes(QStringList genes);
    bool running() const {
        return busy_;
    }
    QString resultDirectory() const {
        return result_.directory;
    }
    void saveScreenshot(QString path);

  protected:
    void dragEnterEvent(QDragEnterEvent *) override;
    void dropEvent(QDropEvent *) override;
    void closeEvent(QCloseEvent *) override;

  private:
    QWidget *dataPage();
    QWidget *progressPage();
    QWidget *resultsPage();
    QWidget *card(QString title, QString subtitle = {});
    void acceptPaths(QStringList paths);
    void indexGenes();
    void showGene(QString gene);
    void appendEvent(QString text);
    void setBusy(bool busy);
    void chooseSpecies();
    void setPage(int page);
    void error(QString message);
    void filterResults();
    QString workspace_;
    go::Catalog catalog_;
    go::Result result_;
    QStringList scope_;
    std::atomic_bool cancel_{false};
    bool busy_ = false, closing_ = false;
    int generation_ = 0;
    QStackedWidget *pages_;
    QPushButton *nav_[3];
    QLabel *status_, *selectionHint_, *resultHint_, *progressText_, *progressSub_, *verdictTitle_,
        *confidence_, *speciesDetail_;
    QLineEdit *paths_[3], *output_, *blast_, *hal_, *resultSearch_;
    GenePicker *genes_;
    QComboBox *reference_, *resultGene_, *classFilter_;
    QCheckBox *reuse_, *consoleToggle_;
    QPushButton *start_, *cancelButton_, *folderButton_;
    QProgressBar *progress_, *score_;
    QPlainTextEdit *console_;
    QLabel *steps_[4];
    EvidenceTree *tree_;
    QTextBrowser *explanation_;
    QTableWidget *speciesTable_, *allGenes_;
    QFutureWatcher<go::Result> *watcher_;
    QString sessionLog_;
};
