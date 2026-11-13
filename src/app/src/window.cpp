#include "window.hpp"
#include <QApplication>
#include <QCloseEvent>
#include <QCompleter>
#include <QDateTime>
#include <QDesktopServices>
#include <QDialogButtonBox>
#include <QDir>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QException>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QJsonDocument>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QPainter>
#include <QScrollArea>
#include <QSettings>
#include <QSplitter>
#include <QStyle>
#include <QTabWidget>
#include <QTimer>
#include <QUrl>
#include <QVBoxLayout>
#include <QtConcurrent>
namespace {
class GuiError final : public QException {
  public:
    explicit GuiError(QString text) : message_(text.toUtf8()) {}
    void raise() const override {
        throw *this;
    }
    GuiError *clone() const override {
        return new GuiError(*this);
    }
    const char *what() const noexcept override {
        return message_.constData();
    }

  private:
    QByteArray message_;
};
QLabel *label(QString text, QString name = {}, QWidget *parent = nullptr) {
    auto *l = new QLabel(text, parent);
    if (!name.isEmpty())
        l->setObjectName(name);
    l->setWordWrap(true);
    l->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);
    if (name == "badge" || name == "eyebrow" || name == "h1" || name == "h2")
        l->setWordWrap(false);
    return l;
}
QPushButton *button(QString text, QString name = {}) {
    auto *b = new QPushButton(text);
    if (!name.isEmpty())
        b->setObjectName(name);
    b->setCursor(Qt::PointingHandCursor);
    b->setMinimumHeight(36);
    return b;
}
QVBoxLayout *vertical(QWidget *w, int margin = 20) {
    auto *l = new QVBoxLayout(w);
    l->setContentsMargins(margin, margin, margin, margin);
    l->setSpacing(12);
    return l;
}
QString orfName(QString s) {
    return QMap<QString, QString>{
        {"preserved_both_identical", "CDS сохранена · оба метода"},
        {"preserved_both_different", "CDS сохранена · последовательности различаются"},
        {"predicted_interruption_both", "Предсказано нарушение CDS"},
        {"multiple_loci_uncertain", "Несколько локусов"},
        {"discordant_methods", "Методы противоречат друг другу"},
        {"conflicting_alignment", "Противоречие выравнивания"},
        {"partial_evidence", "Частичные данные"},
        {"insufficient_evidence", "Недостаточно данных"},
        {"missing_method", "Нет одного метода"}}
        .value(s, s);
}
QString synName(QString s) {
    return QMap<QString, QString>{
        {"direct_CDS_candidate", "Отображение CDS"},
        {"direct_CDS_complex", "Неоднозначное отображение CDS"},
        {"inferred_context", "Двусторонний контекст"},
        {"one_sided_or_discordant_context", "Неполный / противоречивый контекст"},
        {"no_alignment", "Отображение не найдено"},
        {"unavailable_in_hal", "Вида нет в HAL"}}
        .value(s, s);
}
} // namespace
GenePicker::GenePicker(QWidget *p) : QToolButton(p) {
    setObjectName("genePicker");
    setToolButtonStyle(Qt::ToolButtonTextOnly);
    setMinimumHeight(43);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    setCursor(Qt::PointingHandCursor);
    popup_ = new QFrame(this, Qt::Popup);
    popup_->setObjectName("popup");
    popup_->setMinimumSize(430, 390);
    auto *l = vertical(popup_, 15);
    search_ = new QLineEdit;
    search_->setPlaceholderText("Поиск по идентификатору гена…");
    l->addWidget(search_);
    auto *actions = new QHBoxLayout;
    auto *all = button("Выбрать видимые"), *none = button("Снять выбор");
    actions->addWidget(all);
    actions->addWidget(none);
    l->addLayout(actions);
    list_ = new QListWidget;
    list_->setAlternatingRowColors(true);
    l->addWidget(list_);
    l->addWidget(label("Выбор сохраняется при поиске. Для закрытия нажмите вне списка.", "muted"));
    connect(this, &QToolButton::clicked, this, [this] {
        popup_->resize(qMax(430, width()), 430);
        popup_->move(mapToGlobal(QPoint(0, height() + 5)));
        popup_->show();
        search_->setFocus();
    });
    connect(search_, &QLineEdit::textChanged, this, [this](QString s) {
        for (int i = 0; i < list_->count(); ++i)
            list_->item(i)->setHidden(!list_->item(i)->text().contains(s, Qt::CaseInsensitive));
    });
    connect(list_, &QListWidget::itemChanged, this, [this](QListWidgetItem *i) {
        if (i->checkState() == Qt::Checked)
            selected_.insert(i->text());
        else
            selected_.remove(i->text());
        updateLabel();
        emit selectionChanged();
    });
    connect(all, &QPushButton::clicked, this, [this] {
        for (int i = 0; i < list_->count(); ++i)
            if (!list_->item(i)->isHidden())
                list_->item(i)->setCheckState(Qt::Checked);
    });
    connect(none, &QPushButton::clicked, this, [this] {
        for (int i = 0; i < list_->count(); ++i)
            if (!list_->item(i)->isHidden())
                list_->item(i)->setCheckState(Qt::Unchecked);
    });
    updateLabel();
}
void GenePicker::setGenes(QStringList genes) {
    genes_ = genes;
    selected_.clear();
    list_->blockSignals(true);
    list_->clear();
    for (auto g : genes) {
        auto *item = new QListWidgetItem(g, list_);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Unchecked);
    }
    list_->blockSignals(false);
    selectGenes(genes.mid(0, 3));
}
void GenePicker::selectGenes(QStringList genes) {
    QSet<QString> s(genes.begin(), genes.end());
    list_->blockSignals(true);
    selected_.clear();
    for (int i = 0; i < list_->count(); ++i) {
        auto *item = list_->item(i);
        bool chosen = s.contains(item->text());
        item->setCheckState(chosen ? Qt::Checked : Qt::Unchecked);
        if (chosen)
            selected_.insert(item->text());
    }
    list_->blockSignals(false);
    updateLabel();
    emit selectionChanged();
}
QStringList GenePicker::selected() const {
    QStringList out;
    for (auto g : genes_)
        if (selected_.contains(g))
            out << g;
    return out;
}
void GenePicker::updateLabel() {
    setText(genes_.isEmpty()
                ? "Сначала загрузите последовательности    ⌄"
                : QString("Выбрано %1 из %2 генов    ⌄").arg(selected_.size()).arg(genes_.size()));
}
MainWindow::MainWindow(QString workspace) : workspace_(workspace) {
    setWindowTitle("LibreGenesis — происхождение генов");
    resize(1450, 940);
    setMinimumSize(1120, 740);
    setAcceptDrops(true);
    QDir().mkpath(workspace_ + "/app_logs");
    sessionLog_ = workspace_ + "/app_logs/session_" +
                  QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss") + ".jsonl";
    setStyleSheet(R"(
 QMainWindow,QWidget{background:#F3F6F2;color:#203F39;font-family:'Inter','Arial';font-size:13px;}QWidget#sidebar{background:#102F2E;color:#E8F1EA;}QLabel#brand{background:transparent;color:#F0F7F1;font-size:25px;font-weight:700;}QLabel#sideText{background:transparent;color:#91ADA3;font-size:12px;}QLabel#sideTitle{background:transparent;color:#6F9688;font-size:10px;font-weight:700;letter-spacing:2px;}QPushButton#nav{background:transparent;color:#9FBBB0;border:0;text-align:left;padding:14px 18px;border-radius:10px;font-size:14px;}QPushButton#nav[active="true"]{background:#244943;color:#F3F9F4;}QPushButton#nav:hover{background:#20443D;}QFrame#card{background:white;border:1px solid #DCE6DF;border-radius:14px;}QLabel{background:transparent;}QLabel#eyebrow{color:#4B8E7A;font-weight:700;letter-spacing:2px;font-size:10px;}QLabel#h1{font-size:29px;font-weight:700;color:#163F37;}QLabel#h2{font-size:17px;font-weight:700;}QLabel#muted{color:#71867D;font-size:12px;}QLabel#badge{color:#357963;background:#E5EEE6;border-radius:12px;padding:7px 12px;font-size:10px;font-weight:700;}QLabel#fileTag{color:#427B66;background:#EAF1E9;border-radius:8px;padding:9px;font-size:11px;font-weight:700;}QLineEdit,QComboBox,QToolButton#genePicker{background:white;border:1px solid #CFDDD3;border-radius:8px;padding:10px;color:#27483B;selection-background-color:#BCDCCB;}QLineEdit:focus,QComboBox:focus{border:1px solid #278B72;}QComboBox::drop-down{border:0;width:26px;}QPushButton{background:white;border:1px solid #D5E1D7;border-radius:8px;padding:7px 13px;color:#315747;font-weight:600;}QPushButton:hover{background:#EDF4EC;border-color:#87AE9D;}QPushButton:disabled{background:#EDF1EC;color:#A0B0A5;border-color:#E1E8E0;}QPushButton#primary{background:#1B755E;color:white;border:0;padding:12px 22px;border-radius:9px;font-size:14px;}QPushButton#primary:hover{background:#125C4B;}QPushButton#primary:disabled{background:#A8BFB3;color:#EAF2EB;}QPushButton#danger{color:#AD5D43;border-color:#DECCBC;}QProgressBar{border:0;background:#E4ECE3;border-radius:5px;text-align:center;color:#456250;min-height:11px;}QProgressBar::chunk{background:#25856A;border-radius:5px;}QPlainTextEdit{background:#152E29;color:#C0DBCA;border:0;border-radius:10px;padding:12px;font-family:'Menlo','Consolas';font-size:11px;}QTextBrowser{background:white;border:0;font-size:13px;}QListWidget,QTableWidget{background:white;border:0;alternate-background-color:#F5F8F3;selection-background-color:#DDECE1;selection-color:#244D3C;outline:0;}QHeaderView::section{background:#EDF3EB;color:#647E6B;border:0;border-bottom:1px solid #DCE5D9;padding:9px;font-size:11px;font-weight:600;}QTabWidget::pane{background:transparent;border:0;}QTabBar::tab{background:transparent;color:#7C9181;padding:11px 20px;border-bottom:2px solid transparent;}QTabBar::tab:selected{color:#236B50;border-bottom:2px solid #236B50;font-weight:700;}QScrollArea{border:0;background:transparent;}QFrame#popup{background:white;border:1px solid #C8D9CD;border-radius:12px;}QCheckBox{spacing:8px;color:#597462;}QSplitter::handle{background:#F3F6F2;width:8px;}QMenuBar,QMenu{background:white;}QToolTip{background:#173C30;color:#F2F7EF;border:0;padding:7px;}
 )");
    QPixmap icon(128, 128);
    icon.fill(Qt::transparent);
    {
        QPainter p(&icon);
        p.setRenderHint(QPainter::Antialiasing);
        p.setPen(QPen(QColor("#1B755E"), 9, Qt::SolidLine, Qt::RoundCap));
        p.drawLine(28, 64, 68, 64);
        p.drawLine(68, 32, 68, 96);
        p.drawLine(68, 32, 105, 32);
        p.drawLine(68, 96, 105, 96);
        for (auto c : {QPoint(28, 64), QPoint(105, 32), QPoint(105, 96)}) {
            p.setBrush(QColor("#D1B16D"));
            p.setPen(Qt::NoPen);
            p.drawEllipse(c, 10, 10);
        }
    }
    setWindowIcon(QIcon(icon));
    auto *root = new QWidget;
    auto *row = new QHBoxLayout(root);
    row->setContentsMargins(0, 0, 0, 0);
    row->setSpacing(0);
    auto *side = new QWidget;
    side->setObjectName("sidebar");
    side->setFixedWidth(240);
    auto *sl = vertical(side, 23);
    sl->setSpacing(10);
    auto *mark = label("⌁", "brand");
    mark->setStyleSheet("font-size:42px;color:#B6D2B0;");
    sl->addWidget(mark);
    sl->addWidget(label("LibreGenesis", "brand"));
    sl->addSpacing(42);
    sl->addWidget(label("РАБОЧАЯ ОБЛАСТЬ", "sideTitle"));
    QStringList navNames{"01    Данные", "02    Выполнение", "03    Результаты"};
    for (int i = 0; i < 3; ++i) {
        nav_[i] = button(navNames[i], "nav");
        sl->addWidget(nav_[i]);
        connect(nav_[i], &QPushButton::clicked, this, [this, i] {
            if (i == 2 && result_.genes.isEmpty())
                return;
            setPage(i);
        });
    }
    sl->addStretch();
    sl->addWidget(label("Бригада 9", "sideText"));
    auto *repository = label(R"(<a href="https://github.com/noonecanhearyou/sintez-c2-b9" style="color:#B6D2B0">Репозиторий на GitHub ↗</a>)", "sideText");
    repository->setOpenExternalLinks(true);
    sl->addWidget(repository);
    row->addWidget(side);
    auto *right = new QWidget;
    auto *rl = vertical(right, 28);
    auto *top = new QHBoxLayout;
    top->addStretch();
    status_ = label("ГОТОВО К РАБОТЕ", "badge");
    top->addWidget(status_);
    auto *open = button("Открыть результат");
    top->addWidget(open);
    connect(open, &QPushButton::clicked, this, [this] {
        if (busy_)
            return;
        auto p =
            QFileDialog::getExistingDirectory(this, "Открыть результат LibreGenesis", workspace_);
        if (!p.isEmpty())
            loadSaved(p);
    });
    auto *topBar = new QWidget;
    topBar->setObjectName("topBar");
    topBar->setLayout(top);
    top->setContentsMargins(0, 0, 0, 0);
    rl->addWidget(topBar);
    pages_ = new QStackedWidget;
    pages_->addWidget(dataPage());
    pages_->addWidget(progressPage());
    pages_->addWidget(resultsPage());
    rl->addWidget(pages_, 1);
    row->addWidget(right, 1);
    setCentralWidget(root);
    for (auto *field : findChildren<QLineEdit *>())
        field->setAcceptDrops(false);
    auto *file = menuBar()->addMenu("Файл");
    file->addAction(
        "Загрузить исходные файлы…", this,
        [this] {
            acceptPaths(QFileDialog::getOpenFileNames(this, "Выберите HAL, TXT и CSV", workspace_,
                                                      "Геномные данные (*.hal *.txt *.csv)"));
        },
        QKeySequence::Open);
    file->addAction("Открыть результат…", open, &QPushButton::click);
    file->addSeparator();
    file->addAction("Выход", this, &QWidget::close, QKeySequence::Quit);
    watcher_ = new QFutureWatcher<go::Result>(this);
    connect(watcher_, &QFutureWatcher<go::Result>::finished, this, [this] {
        bool cancelled = cancel_.load();
        setBusy(false);
        try {
            result_ = watcher_->result();
            resultGene_->clear();
            for (auto v : result_.genes)
                resultGene_->addItem(v.toObject()["gene_id"].toString());
            allGenes_->setRowCount(result_.genes.size());
            int i = 0;
            for (auto v : result_.genes) {
                auto g = v.toObject();
                QStringList cells{g["gene_id"].toString(), g["class_label"].toString(),
                                  g["confidence_level"].toString(),
                                  QString::number(g["support_score"].toInt()),
                                  g["reason"].toString()};
                for (int j = 0; j < cells.size(); ++j)
                    allGenes_->setItem(i, j, new QTableWidgetItem(cells[j]));
                allGenes_->item(i, 1)->setData(Qt::UserRole, g["class"].toString());
                ++i;
            }
            resultHint_->setText(QString("%1 генов · %2 · результаты сохранены")
                                     .arg(result_.genes.size())
                                     .arg(go::displaySpecies(result_.reference)));
            if (resultGene_->count())
                showGene(resultGene_->currentText());
            auto modules = result_.manifest.value("modules").toObject();
            status_->setText((modules.value("synteny") == "completed" &&
                              modules.value("orf") == "completed" &&
                              modules.value("homology") == "completed")
                                 ? "АНАЛИЗ ЗАВЕРШЁН"
                                 : "ЗАВЕРШЕНО С ОГРАНИЧЕНИЯМИ");
            setPage(2);
        } catch (const std::exception &e) {
            status_->setText(cancelled ? "ОТМЕНЕНО" : "НУЖНО ВНИМАНИЕ");
            if (cancelled) {
                progressText_->setText("Анализ отменён. Записанные данные сохранены.");
                appendEvent("WARNING | Отмена | " + QString::fromUtf8(e.what()));
            } else if (!closing_)
                error(QString::fromUtf8(e.what()));
        }
        if (closing_)
            QTimer::singleShot(0, this, &QWidget::close);
    });
    setPage(0);
    appendEvent("INFO | Приложение | Рабочая область готова");
}
QWidget *MainWindow::card(QString title, QString subtitle) {
    auto *frame = new QFrame;
    frame->setObjectName("card");
    auto *l = vertical(frame);
    l->setAlignment(Qt::AlignTop);
    if (!title.isEmpty())
        l->addWidget(label(title, "h2"));
    if (!subtitle.isEmpty())
        l->addWidget(label(subtitle, "muted"));
    return frame;
}
QWidget *MainWindow::dataPage() {
    auto *page = new QWidget;
    auto *l = vertical(page, 0);
    auto *grid = new QGridLayout;
    grid->setSpacing(18);
    auto *input =
        card("01  Исходные данные", "Перетащите три файла в окно или выберите их через меню ОС.");
    auto *il = qobject_cast<QVBoxLayout *>(input->layout());
    auto *dropHint = label("↓\nПеретащите HAL · TXT · CSV сюда", "muted");
    dropHint->setAlignment(Qt::AlignCenter);
    dropHint->setMinimumHeight(110);
    dropHint->setStyleSheet("QLabel {border:1px dashed "
                            "#B8D0BC;border-radius:12px;background:#F7FAF5;color:"
                            "#6B8C73;font-size:14px;}");
    il->setAlignment(Qt::Alignment{});
    il->addWidget(dropHint, 2);
    QStringList tags{"HAL", "TXT", "CSV"},
        names{"Полногеномное выравнивание", "Последовательности кандидатов",
              "Организмы и версии сборок"};
    for (int i = 0; i < 3; ++i) {
        auto *r = new QHBoxLayout;
        auto *tag = label(tags[i], "fileTag");
        tag->setFixedWidth(52);
        r->addWidget(tag);
        auto *c = new QVBoxLayout;
        c->setSpacing(4);
        c->addWidget(label(names[i]));
        paths_[i] = new QLineEdit;
        paths_[i]->setPlaceholderText("Файл не выбран");
        paths_[i]->setObjectName("input" + QString::number(i));
        c->addWidget(paths_[i]);
        r->addLayout(c, 1);
        auto *pick = button("Выбрать");
        r->addWidget(pick);
        il->addStretch(1);
        il->addLayout(r);
        connect(pick, &QPushButton::clicked, this, [this, i, tags] {
            QString ext = tags[i].toLower();
            auto p = QFileDialog::getOpenFileName(this, "Выбрать " + tags[i], workspace_,
                                                  "Файлы " + tags[i] + " (*." + ext + ")");
            if (!p.isEmpty()) {
                paths_[i]->setText(p);
                if (i == 1)
                    indexGenes();
            }
        });
        connect(paths_[i], &QLineEdit::editingFinished, this, [this, i] {
            if (i == 1 && !paths_[i]->text().isEmpty())
                indexGenes();
        });
    }
    auto *project = button("Использовать файлы проекта");
    il->addStretch(1);
    il->addWidget(project);
    connect(project, &QPushButton::clicked, this, &MainWindow::loadProjectInputs);
    grid->addWidget(input, 0, 0, 2, 1);
    auto *selection = card("02  Область исследования",
                           "Выберите гены и организм, относительно которого выполняется пересчёт.");
    auto *se = qobject_cast<QVBoxLayout *>(selection->layout());
    se->addWidget(label("Исследуемые гены", "muted"));
    genes_ = new GenePicker;
    genes_->setObjectName("genePicker");
    se->addWidget(genes_);
    se->addWidget(label("Основной организм", "muted"));
    reference_ = new QComboBox;
    reference_->setObjectName("reference");
    reference_->addItem("Сначала загрузите TXT");
    se->addWidget(reference_);
    selectionHint_ =
        label("Список генов появится после проверки файла последовательностей.", "muted");
    se->addWidget(selectionHint_);
    auto *species = button("Виды организмов…");
    se->addWidget(species);
    connect(species, &QPushButton::clicked, this, &MainWindow::chooseSpecies);
    connect(genes_, &GenePicker::selectionChanged, this, [this] {
        selectionHint_->setText(
            QString("Выбрано генов: %1. Смена референса пересчитывает все три блока.")
                .arg(genes_->selected().size()));
    });
    grid->addWidget(selection, 0, 1);
    auto *storage = card("03  Сохранение", "Для каждого запуска создаётся отдельная папка. "
                                           "Исходные файлы не перезаписываются.");
    auto *st = qobject_cast<QVBoxLayout *>(storage->layout());
    auto *rr = new QHBoxLayout;
    output_ = new QLineEdit;
    output_->setObjectName("outputRoot");
    output_->setPlaceholderText("Рядом с файлом последовательностей");
    rr->addWidget(output_, 1);
    auto *browse = button("Папка…");
    rr->addWidget(browse);
    st->addLayout(rr);
    connect(browse, &QPushButton::clicked, this, [this] {
        auto p = QFileDialog::getExistingDirectory(this, "Каталог результатов", workspace_);
        if (!p.isEmpty())
            output_->setText(p);
    });
    reuse_ = new QCheckBox("Повторно использовать результаты с проверенными входами");
    reuse_->setChecked(true);
    reuse_->setToolTip("Сохранённая гомология dmel используется только после "
                       "проверки хэшей, белковых запросов и версий сборок. Смена "
                       "референса требует нового поиска.");
    st->addWidget(reuse_);
    auto *options = new QWidget;
    auto *ol = vertical(options, 0);
    ol->addWidget(label("Папка BLAST+", "muted"));
    blast_ = new QLineEdit(go::findBlastBin(workspace_));
    ol->addWidget(blast_);
    ol->addWidget(label("Папка HAL-инструментов", "muted"));
    hal_ = new QLineEdit(go::findHalBin());
    ol->addWidget(hal_);
    st->addWidget(options);
    options->hide();
    grid->addWidget(storage, 1, 1);
    grid->setColumnStretch(0, 5);
    grid->setColumnStretch(1, 4);
    l->addLayout(grid, 1);
    auto *bottom = new QHBoxLayout;
    bottom->addStretch();
    start_ = button("Запустить анализ  →", "primary");
    start_->setObjectName("primary");
    bottom->addWidget(start_);
    connect(start_, &QPushButton::clicked, this, &MainWindow::startAnalysis);
    l->addLayout(bottom);
    return page;
}
QWidget *MainWindow::progressPage() {
    auto *page = new QWidget;
    auto *l = vertical(page, 0);
    progressSub_ = label("Анализ выполняется в фоне. Окно остаётся доступным.", "muted");
    l->addWidget(progressSub_);
    auto *c = card("Ход анализа");
    auto *cl = qobject_cast<QVBoxLayout *>(c->layout());
    progressText_ = label("Подготовка", "h2");
    cl->addWidget(progressText_);
    progress_ = new QProgressBar;
    progress_->setRange(0, 100);
    progress_->setTextVisible(false);
    progress_->setMinimumHeight(12);
    cl->addWidget(progress_);
    QStringList names{"Гомология · белки и геномы", "Синтения · соответствующие участки",
                      "ORF · кодирующий потенциал", "Классификация и сохранение"};
    for (int i = 0; i < 4; ++i) {
        steps_[i] = label("○   " + names[i], "muted");
        cl->addWidget(steps_[i]);
    }
    l->addWidget(c);
    consoleToggle_ = new QCheckBox("Показывать журнал этапов");
    consoleToggle_->setChecked(true);
    l->addWidget(consoleToggle_);
    console_ = new QPlainTextEdit;
    console_->setReadOnly(true);
    console_->setMaximumBlockCount(2500);
    l->addWidget(console_, 1);
    connect(consoleToggle_, &QCheckBox::toggled, console_, &QWidget::setVisible);
    auto *r = new QHBoxLayout;
    r->addStretch();
    cancelButton_ = button("Отменить анализ", "danger");
    r->addWidget(cancelButton_);
    connect(cancelButton_, &QPushButton::clicked, this, [this] {
        cancel_.store(true);
        cancelButton_->setEnabled(false);
        progressText_->setText("Отмена: завершаем текущую операцию…");
        appendEvent("WARNING | Отмена | Запрошена остановка. Уже записанные файлы "
                    "сохраняются.");
    });
    l->addLayout(r);
    return page;
}
QWidget *MainWindow::resultsPage() {
    auto *page = new QWidget;
    auto *l = vertical(page, 0);
    auto *head = new QHBoxLayout;
    head->addStretch();
    folderButton_ = button("Открыть папку результатов");
    head->addWidget(folderButton_);
    connect(folderButton_, &QPushButton::clicked, this, [this] {
        if (!result_.directory.isEmpty())
            QDesktopServices::openUrl(QUrl::fromLocalFile(result_.directory));
    });
    resultHint_ = label("Результаты появятся после выполнения анализа.", "muted");

    auto *controls = new QHBoxLayout;
    controls->addWidget(label("Ген"));
    resultGene_ = new QComboBox;
    resultGene_->setEditable(false);
    resultGene_->setObjectName("resultGene");
    resultGene_->setStyleSheet("QComboBox::down-arrow {image:url(:/geneorigin/chevron.png);width:12px;height:8px;}");
    resultGene_->setMaxVisibleItems(20);
    resultGene_->setInsertPolicy(QComboBox::NoInsert);
    resultGene_->setMinimumWidth(245);

    controls->addWidget(resultGene_);
    controls->addStretch();
    auto *legend = label("<span style='color:#168D77'>● CDS сохранена</span> "
                         "&nbsp; <span style='color:#4382B4'>● Сходство</span> "
                         "&nbsp; <span style='color:#D8963D'>● Конфликт</span>",
                         "muted");
    legend->setWordWrap(false);
    controls->addWidget(legend);
    l->addLayout(controls);
    connect(resultGene_, &QComboBox::currentTextChanged, this, &MainWindow::showGene);
    auto *tabs = new QTabWidget;
    l->addWidget(tabs, 1);
    auto *view = new QWidget;
    auto *vl = vertical(view, 0);
    auto *split = new QSplitter(Qt::Horizontal);
    auto *tc = card("");
    tc->setToolTip("Дерево видов из HAL. Цвет отражает данные выбранного гена, "
                   "а не реконструкцию предкового гена.");
    auto *tl = qobject_cast<QVBoxLayout *>(tc->layout());
    auto *tools = new QHBoxLayout;
    tools->addWidget(label("Дерево видов", "h2"));
    tools->addStretch();
    for (auto name : {QString("−"), QString("+"), QString("Сохранить PNG")}) {
        auto *b = button(name);
        tools->addWidget(b);
        connect(b, &QPushButton::clicked, this, [this, name] {
            if (name == "−")
                tree_->zoomOut();
            else if (name == "+")
                tree_->zoomIn();
            else {
                auto path = QFileDialog::getSaveFileName(
                    this, "Сохранить дерево", result_.directory + "/tree.png", "PNG (*.png)");
                if (!path.isEmpty() && !tree_->savePng(path))
                    error("Не удалось сохранить изображение дерева.");
            }
        });
    }
    tools->addStretch();
    tl->addLayout(tools);
    tree_ = new EvidenceTree;
    tl->addWidget(tree_, 1);
    speciesDetail_ =
        label("Колесо — масштаб · перетаскивание — перемещение · нажатие — данные вида", "muted");
    speciesDetail_->hide();
    split->addWidget(tc);
    auto *verdict = card("Рабочая гипотеза");
    auto *ve = qobject_cast<QVBoxLayout *>(verdict->layout());
    verdictTitle_ = label("Ambiguous", "h1");
    ve->addWidget(verdictTitle_);
    confidence_ = label("Недостаточно данных", "muted");
    ve->addWidget(confidence_);
    score_ = new QProgressBar;
    score_->setRange(0, 100);
    score_->setTextVisible(false);
    ve->addWidget(score_);
    ve->addWidget(
        label("Поддержка 0–100 — эвристический балл, не вероятность истинности.", "muted"));
    explanation_ = new QTextBrowser;
    explanation_->setOpenExternalLinks(false);
    ve->addWidget(explanation_, 1);
    verdict->setMinimumWidth(330);
    split->setHandleWidth(12);
    split->setChildrenCollapsible(false);
    vl->addWidget(split, 1);
    speciesTable_ = new QTableWidget;
    speciesTable_->setColumnCount(4);
    speciesTable_->setHorizontalHeaderLabels({"Организм", "Гомология", "Синтения", "ORF"});
    speciesTable_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    speciesTable_->setWordWrap(true);
    connect(speciesTable_->horizontalHeader(), &QHeaderView::sectionResized, this,
            [this] { speciesTable_->resizeRowsToContents(); });
    speciesTable_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    speciesTable_->setSelectionBehavior(QAbstractItemView::SelectRows);
    speciesTable_->setAlternatingRowColors(true);
    auto *speciesCard = card("Данные видов");
    speciesCard->layout()->addWidget(speciesTable_);
    split->addWidget(speciesCard);
    split->setStretchFactor(0, 45);
    split->setStretchFactor(1, 55);
    split->setSizes({450, 550});
    tabs->addTab(view, "Дерево");
    tabs->addTab(verdict, "Гипотеза");
    auto *about = new QWidget;
    auto *aboutLayout = vertical(about);
    aboutLayout->addWidget(label("Сведения об анализе", "h2"));
    aboutLayout->addWidget(resultHint_);
    aboutLayout->addLayout(head);
    aboutLayout->addStretch();
    tabs->addTab(about, "Об анализе");
    connect(tree_, &EvidenceTree::speciesSelected, this, [this](QString sp) {
        for (int i = 0; i < speciesTable_->rowCount(); ++i)
            if (speciesTable_->item(i, 0)->data(Qt::UserRole).toString() == sp) {
                speciesTable_->selectRow(i);
                speciesTable_->scrollToItem(speciesTable_->item(i, 0));
                speciesDetail_->setText(go::displaySpecies(sp) + " · " +
                                        speciesTable_->item(i, 1)->text() + " · " +
                                        speciesTable_->item(i, 3)->text());
                break;
            }
    });
    auto *tablePage = new QWidget;
    auto *tp = vertical(tablePage, 0);
    auto *filters = new QHBoxLayout;
    resultSearch_ = new QLineEdit;
    resultSearch_->setPlaceholderText("Поиск гена или фрагмента обоснования…");
    filters->addWidget(resultSearch_, 1);
    classFilter_ = new QComboBox;
    classFilter_->addItem("Все классы", "");
    for (auto c : {"de_novo", "duplication", "ancient_diverged", "other", "ambiguous"})
        classFilter_->addItem(go::className(c), c);
    filters->addWidget(classFilter_);
    tp->addLayout(filters);
    allGenes_ = new QTableWidget;
    allGenes_->setColumnCount(5);
    allGenes_->setHorizontalHeaderLabels(
        {"Ген", "Класс-гипотеза", "Доверие", "Поддержка", "Обоснование"});
    allGenes_->horizontalHeader()->setStretchLastSection(true);
    allGenes_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    allGenes_->setSelectionBehavior(QAbstractItemView::SelectRows);
    allGenes_->setAlternatingRowColors(true);
    tp->addWidget(allGenes_, 1);
    tabs->addTab(tablePage, "Все гены");
    connect(resultSearch_, &QLineEdit::textChanged, this, &MainWindow::filterResults);
    connect(classFilter_, &QComboBox::currentIndexChanged, this, &MainWindow::filterResults);
    connect(allGenes_, &QTableWidget::cellDoubleClicked, this, [this, tabs](int row, int) {
        if (allGenes_->item(row, 0)) {
            resultGene_->setCurrentText(allGenes_->item(row, 0)->text());
            tabs->setCurrentIndex(0);
        }
    });
    auto *files = new QWidget;
    auto *fl = vertical(files);
    fl->addWidget(label("Файлы для проверки и передачи команде", "h2"));
    fl->addWidget(label("summary.tsv / summary.json — ген, гипотеза класса, доверие, поддержка и "
                        "основания.\nhomology/ — совпадения и состояния поиска.\nsynteny/ — "
                        "переносы, контекст, координаты и ограничения.\norf/ — пять этапов "
                        "анализа и профили локусов.\nevents.jsonl — отдельный журнал INFO / "
                        "WARNING / ERROR.\nrun.json — входы, область анализа, версии и статусы "
                        "модулей.",
                        "muted"));
    auto *exportFile = button("Сохранить копию summary.tsv…");
    fl->addWidget(exportFile);
    connect(exportFile, &QPushButton::clicked, this, [this] {
        if (result_.directory.isEmpty())
            return;
        auto path = QFileDialog::getSaveFileName(this, "Экспорт итоговой таблицы",
                                                 QDir::homePath() + "/gene_origin_summary.tsv",
                                                 "TSV (*.tsv)");
        if (!path.isEmpty()) {
            QFile f(result_.directory + "/summary.tsv");
            if (!f.open(QIODevice::ReadOnly)) {
                error("Итоговая таблица недоступна.");
                return;
            }
            try {
                homology::writeText(path, f.readAll());
            } catch (const std::exception &e) {
                error(e.what());
            }
        }
    });
    fl->addWidget(label("Отсутствие совпадений не доказывает De novo. Неприменимые методы "
                        "и неоднозначные локусы сохранены явно.",
                        "muted"));
    fl->addStretch();
    tabs->addTab(files, "Файлы и журнал");
    return page;
}
void MainWindow::setPage(int page) {
    pages_->setCurrentIndex(page);
    if (auto *topBar = findChild<QWidget *>("topBar"))
        topBar->setVisible(page != 2);
    for (int i = 0; i < 3; ++i) {
        nav_[i]->setProperty("active", i == page);
        nav_[i]->style()->unpolish(nav_[i]);
        nav_[i]->style()->polish(nav_[i]);
    }
    if (page == 2)
        QTimer::singleShot(0, tree_, &EvidenceTree::fitTree);
}
void MainWindow::loadProjectInputs() {
    acceptPaths({workspace_ + "/Acalyptratae20.hal",
                 workspace_ + "/supplementary_file_unannotated_alignments.txt",
                 workspace_ + "/supplementary_file_genome_assembly.csv"});
}
void MainWindow::acceptPaths(QStringList paths) {
    if (busy_)
        return;
    bool txt = false;
    for (auto p : paths) {
        auto ext = QFileInfo(p).suffix().toLower();
        int i = ext == "hal" ? 0 : ext == "txt" ? 1 : ext == "csv" ? 2 : -1;
        if (i < 0 || !QFileInfo(p).isFile()) {
            error("Ожидаются три локальных файла: HAL, TXT и CSV. Не распознан: " + p);
            continue;
        }
        paths_[i]->setText(QFileInfo(p).absoluteFilePath());
        txt |= i == 1;
    }
    if (txt)
        indexGenes();
    setPage(0);
}
void MainWindow::indexGenes() {
    QString path = paths_[1]->text();
    if (path.isEmpty())
        return;
    int generation = ++generation_;
    start_->setEnabled(false);
    selectionHint_->setText("Читаем список генов и организмов…");
    auto *watcher = new QFutureWatcher<go::Catalog>(this);
    connect(watcher, &QFutureWatcher<go::Catalog>::finished, this, [this, watcher, generation] {
        if (generation != generation_) {
            watcher->deleteLater();
            return;
        }
        try {
            catalog_ = watcher->result();
            genes_->setGenes(catalog_.genes);
            reference_->clear();
            for (auto sp : catalog_.species)
                reference_->addItem(go::displaySpecies(sp) + "  ·  " + sp, sp);
            int index = reference_->findData("dmel");
            if (index >= 0)
                reference_->setCurrentIndex(index);
            start_->setEnabled(true);
            appendEvent(QString("INFO | Входы | Найдено %1 генов, %2 видов")
                            .arg(catalog_.genes.size())
                            .arg(catalog_.species.size()));
        } catch (const std::exception &e) {
            genes_->setGenes({});
            error(e.what());
        }
        watcher->deleteLater();
    });
    watcher->setFuture(QtConcurrent::run([path] {
        try {
            return go::catalog(path);
        } catch (const std::exception &e) {
            throw GuiError(QString::fromUtf8(e.what()));
        }
    }));
}
void MainWindow::setBusy(bool busy) {
    busy_ = busy;
    start_->setEnabled(!busy);
    cancelButton_->setEnabled(busy);
    genes_->setEnabled(!busy);
    reference_->setEnabled(!busy);
    for (auto p : paths_)
        p->setReadOnly(busy);
    nav_[0]->setEnabled(!busy);
    nav_[2]->setEnabled(!busy && !result_.genes.isEmpty());
}
void MainWindow::selectGenes(QStringList genes) {
    genes_->selectGenes(genes);
}
void MainWindow::startAnalysis() {
    if (busy_)
        return;
    if (genes_->selected().isEmpty()) {
        error("Выберите хотя бы один ген из списка.");
        return;
    }
    for (auto p : paths_)
        if (!QFileInfo(p->text()).isFile()) {
            error("Выберите три существующих файла HAL, TXT и CSV.");
            return;
        }
    if (reference_->currentData().toString().isEmpty()) {
        error("Выберите основной организм.");
        return;
    }
    go::Request request;
    request.hal = paths_[0]->text();
    request.alignment = paths_[1]->text();
    request.assemblies = paths_[2]->text();
    request.genes = genes_->selected();
    request.reference = reference_->currentData().toString();
    request.workspace = workspace_;
    request.outputRoot = output_->text().trimmed().isEmpty()
                             ? QFileInfo(request.alignment).absolutePath()
                             : output_->text().trimmed();
    request.blastBin = blast_->text().trimmed();
    request.halBin = hal_->text().trimmed();
    request.useVerified = reuse_->isChecked();
    request.comparisonSpecies = scope_;
    cancel_.store(false);
    console_->clear();
    progress_->setValue(0);
    setBusy(true);
    setPage(1);
    status_->setText("ВЫПОЛНЯЕТСЯ АНАЛИЗ");
    progressSub_->setText(QString("%1 генов · основной организм: %2")
                              .arg(request.genes.size())
                              .arg(go::displaySpecies(request.reference)));
    appendEvent("INFO | Запуск | Начат новый анализ");
    watcher_->setFuture(QtConcurrent::run([this, request] {
        homology::Hooks hooks;
        hooks.cancelled = &cancel_;
        hooks.progress = [this](QString text) {
            QMetaObject::invokeMethod(
                this, [this, text] { appendEvent(text); }, Qt::QueuedConnection);
        };
        auto progress = [this](int value, QString text) {
            QMetaObject::invokeMethod(
                this,
                [this, value, text] {
                    progress_->setValue(value);
                    progressText_->setText(text);
                    int current = value < 45 ? 0 : value < 67 ? 1 : value < 90 ? 2 : 3;
                    for (int i = 0; i < 4; ++i) {
                        QString t = steps_[i]->text();
                        t.remove(0, 1);
                        steps_[i]->setText((i < current ? "✓" : i == current ? "●" : "○") + t);
                        steps_[i]->setStyleSheet(i == current ? "color:#18755C;font-weight:700;"
                                                              : "color:#82968A;");
                    }
                },
                Qt::QueuedConnection);
        };
        try {
            return go::execute(request, hooks, progress);
        } catch (const std::exception &e) {
            throw GuiError(QString::fromUtf8(e.what()));
        }
    }));
}
void MainWindow::appendEvent(QString text) {
    QString rendered = QDateTime::currentDateTime().toString("HH:mm:ss") + "  " + text;
    if (console_)
        console_->appendPlainText(rendered);
    auto parts = text.split(" | ");
    QJsonObject row{{"time_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs)},
                    {"severity", parts.value(0, "INFO")},
                    {"stage", parts.value(1, "Приложение")},
                    {"message", parts.mid(2).join(" | ")}};
    QFile f(sessionLog_);
    if (f.open(QIODevice::WriteOnly | QIODevice::Append))
        f.write(QJsonDocument(row).toJson(QJsonDocument::Compact) + '\n');
}
void MainWindow::error(QString text) {
    appendEvent("ERROR | Ошибка | " + text);
    QMessageBox box(this);
    box.setIcon(QMessageBox::Warning);
    box.setWindowTitle("LibreGenesis — нужно внимание");
    QString friendly = text.section('\n', 0, 0);
    if (text.contains("Cannot start"))
        friendly = "Не удалось запустить внешний инструмент. Проверьте папки "
                   "BLAST+ и HAL.";
    else if (text.contains("Process timeout"))
        friendly = "Внешний инструмент превысил время ожидания. Уменьшите область "
                   "анализа или проверьте журнал.";
    else if (text.contains("hash mismatch", Qt::CaseInsensitive))
        friendly = "Контрольная сумма данных не совпадает. Результат изменён или "
                   "относится к другому запуску.";
    else if (text.contains("Invalid JSON") || text.contains("json.exception"))
        friendly = "Не удалось прочитать данные одного из блоков. Проверьте "
                   "целостность входов и результатов.";
    else if (text.contains("bad_alloc"))
        friendly = "Недостаточно памяти. Выберите меньше генов и повторите анализ.";
    box.setText(friendly);
    box.setInformativeText("Проверьте входные данные и настройки. Подробности "
                           "сохраняются в журнале.");
    box.setDetailedText(text + "\n\nЖурнал приложения: " + sessionLog_);
    box.exec();
}
void MainWindow::showGene(QString gene) {
    QJsonObject record;
    for (auto v : result_.genes)
        if (v.toObject()["gene_id"] == gene) {
            record = v.toObject();
            break;
        }
    if (record.isEmpty())
        return;
    auto ev = record["evidence"].toObject();
    tree_->setResult(result_.tree, ev, gene, result_.reference);
    verdictTitle_->setText(record["class_label"].toString());
    confidence_->setText(record["confidence_level"].toString() + " · поддержка " +
                         QString::number(record["support_score"].toInt()) + " / 100");
    score_->setValue(record["support_score"].toInt());
    QString html = "<p style='color:#294B3B;line-height:1.5'>" +
                   record["reason"].toString().toHtmlEscaped() +
                   "</p><p style='color:#658370;font-size:11px'><b>ПРОВЕРЕННЫЕ "
                   "ОСНОВАНИЯ</b></p>";
    for (auto v : record["rule_trace"].toArray()) {
        auto r = v.toObject();
        html +=
            "<p style='font-size:12px;color:" + QString(r["met"].toBool() ? "#267854" : "#86988A") +
            "'>" + (r["met"].toBool() ? QString("✓  ") : QString("○  ")) +
            r["text"].toString().toHtmlEscaped() + "</p>";
    }
    html += "<p style='color:#8E7954;font-size:11px'><b>ОГРАНИЧЕНИЯ</b></p>";
    for (auto v : record["limitations"].toArray())
        html +=
            "<p style='font-size:11px;color:#7C856D'>• " + v.toString().toHtmlEscaped() + "</p>";
    if (!ev["unavailable_reason"].toString().isEmpty())
        html += "<p>" + ev["unavailable_reason"].toString().toHtmlEscaped() + "</p>";
    explanation_->setHtml(html);
    auto species = ev["species"].toObject();
    speciesTable_->setRowCount(species.size());
    int i = 0;
    for (auto sp : species.keys()) {
        auto o = species[sp].toObject();
        QStringList cells{go::displaySpecies(sp), o["homology"].toString(),
                          synName(o["synteny"].toString()), orfName(o["orf"].toString())};
        for (int j = 0; j < 4; ++j)
            speciesTable_->setItem(i, j, new QTableWidgetItem(cells[j]));
        speciesTable_->item(i, 0)->setData(Qt::UserRole, sp);
        ++i;
    }
    speciesTable_->resizeRowsToContents();
}
void MainWindow::filterResults() {
    for (int i = 0; i < allGenes_->rowCount(); ++i) {
        auto *id = allGenes_->item(i, 0);
        if (!id)
            continue;
        bool text =
            id->text().contains(resultSearch_->text(), Qt::CaseInsensitive) ||
            allGenes_->item(i, 4)->text().contains(resultSearch_->text(), Qt::CaseInsensitive);
        bool type = classFilter_->currentIndex() == 0 ||
                    allGenes_->item(i, 1)->data(Qt::UserRole).toString() == classFilter_->currentData().toString();
        allGenes_->setRowHidden(i, !text || !type);
    }
}
void MainWindow::chooseSpecies() {
    QDialog dialog(this);
    dialog.setWindowTitle("Виды организмов");
    dialog.resize(470, 570);
    auto *l = vertical(&dialog);
    l->addWidget(label("Референс добавляется автоматически. Aedes aegypti "
                       "исключён из текущей области анализа.",
                       "muted"));
    auto *list = new QListWidget;
    l->addWidget(list, 1);
    auto names = go::speciesNames();
    for (auto sp : names.keys()) {
        if (sp == "aaeg")
            continue;
        auto *i = new QListWidgetItem(names[sp], list);
        i->setData(Qt::UserRole, sp);
        i->setFlags(i->flags() | Qt::ItemIsUserCheckable);
        i->setCheckState(scope_.isEmpty() || scope_.contains(sp) ? Qt::Checked : Qt::Unchecked);
    }
    auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    l->addWidget(buttons);
    connect(buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
    if (dialog.exec() == QDialog::Accepted) {
        scope_.clear();
        for (int i = 0; i < list->count(); ++i)
            if (list->item(i)->checkState() == Qt::Checked)
                scope_ << list->item(i)->data(Qt::UserRole).toString();
        if (scope_.isEmpty())
            scope_ << reference_->currentData().toString();
    }
}
void MainWindow::loadSaved(QString dir) {
    if (busy_)
        return;
    try {
        result_ = go::loadResult(dir);
        resultGene_->clear();
        for (auto v : result_.genes)
            resultGene_->addItem(v.toObject()["gene_id"].toString());
        allGenes_->setRowCount(result_.genes.size());
        int i = 0;
        for (auto v : result_.genes) {
            auto g = v.toObject();
            QStringList cells{g["gene_id"].toString(), g["class_label"].toString(),
                              g["confidence_level"].toString(),
                              QString::number(g["support_score"].toInt()), g["reason"].toString()};
            for (int j = 0; j < cells.size(); ++j)
                allGenes_->setItem(i, j, new QTableWidgetItem(cells[j]));
            allGenes_->item(i, 1)->setData(Qt::UserRole, g["class"].toString());
            ++i;
        }
        resultHint_->setText(QString("%1 генов · %2 · сохранённый результат")
                                 .arg(result_.genes.size())
                                 .arg(go::displaySpecies(result_.reference)));
        nav_[2]->setEnabled(true);
        auto modules = result_.manifest.value("modules").toObject();
        status_->setText(modules.isEmpty() ? "СОХРАНЁННЫЙ РЕЗУЛЬТАТ"
                         : (modules.value("synteny") == "completed" &&
                            modules.value("orf") == "completed" &&
                            modules.value("homology") == "completed")
                             ? "АНАЛИЗ ЗАВЕРШЁН"
                             : "ЗАВЕРШЕНО С ОГРАНИЧЕНИЯМИ");
        setPage(2);
        if (resultGene_->count())
            showGene(resultGene_->currentText());
    } catch (const std::exception &e) {
        error(e.what());
    }
}
void MainWindow::dragEnterEvent(QDragEnterEvent *e) {
    if (!busy_ && e->mimeData()->hasUrls())
        e->acceptProposedAction();
}
void MainWindow::dropEvent(QDropEvent *e) {
    QStringList paths;
    for (auto url : e->mimeData()->urls())
        if (url.isLocalFile())
            paths << url.toLocalFile();
    acceptPaths(paths);
    e->acceptProposedAction();
}
void MainWindow::closeEvent(QCloseEvent *e) {
    if (busy_) {
        closing_ = true;
        cancel_.store(true);
        cancelButton_->setEnabled(false);
        progressText_->setText("Завершение текущей операции перед закрытием…");
        e->ignore();
    } else
        e->accept();
}
void MainWindow::saveScreenshot(QString path) {
    grab().save(path);
}
