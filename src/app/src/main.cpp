#include "window.hpp"
#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QStandardPaths>
#include <QTextStream>
#include <QTimer>
#ifdef Q_OS_MACOS
#include <mach-o/dyld.h>
#include <vector>
#endif
int main(int argc, char **argv) {
    // Register the bundle's plug-ins before QApplication initializes Cocoa.
    // Finder launch and headless launch can resolve qt.conf prefixes differently.
#ifdef Q_OS_MACOS
    uint32_t executablePathSize = 0;
    _NSGetExecutablePath(nullptr, &executablePathSize);
    std::vector<char> executablePath(executablePathSize);
    if (_NSGetExecutablePath(executablePath.data(), &executablePathSize) == 0) {
        const QFileInfo executable(QString::fromUtf8(executablePath.data()));
        const QString plugins = QDir(executable.absolutePath()).absoluteFilePath("../PlugIns");
        if (QDir(plugins).exists())
            QCoreApplication::addLibraryPath(QDir(plugins).canonicalPath());
    }
#endif
    QApplication app(argc, argv);
    app.setStyle("Fusion");
    QPalette palette;
    palette.setColor(QPalette::Window, QColor("#F3F6F2"));
    palette.setColor(QPalette::Base, Qt::white);
    palette.setColor(QPalette::WindowText, QColor("#203F39"));
    palette.setColor(QPalette::Text, QColor("#203F39"));
    palette.setColor(QPalette::Button, Qt::white);
    palette.setColor(QPalette::ButtonText, QColor("#203F39"));
    palette.setColor(QPalette::Highlight, QColor("#23775C"));
    palette.setColor(QPalette::HighlightedText, Qt::white);
    app.setPalette(palette);
    app.setApplicationName("LibreGenesis");
    app.setOrganizationName("LibreGenesis");
    QCommandLineParser p;
    p.setApplicationDescription("LibreGenesis — объединённый анализ происхождения генов");
    p.addHelpOption();
    p.addVersionOption();
    app.setApplicationVersion("1.0.0");
    auto opt = [&](QString name, QString help, QString value = {}, QString initial = {}) {
        p.addOption(QCommandLineOption(QStringList{name}, help, value, initial));
    };
    opt("cli", "Выполнить анализ без окна");
    opt("hal", "HAL-выравнивание", "path");
    opt("input", "TXT последовательностей", "path");
    opt("assemblies", "CSV сборок", "path");
    opt("reference", "Основной организм", "code", "dmel");
    opt("genes", "Идентификаторы через запятую", "list");
    opt("scope", "Коды видов сравнения через запятую", "list");
    opt("output", "Каталог сохранения", "path");
    opt("workspace", "Каталог проекта и биологических баз", "path",
        QString::fromUtf8(GENEORIGIN_WORKSPACE));
    opt("fresh", "Выполнить новый BLAST");
    opt("open-result", "Открыть сохранённый результат", "directory");
    opt("demo-inputs", "Заполнить входы файлами проекта");
    opt("screenshot", "Сохранить снимок окна", "path");
    opt("smoke-test", "Проверить запуск окна");
    p.process(app);
    QString workspace = p.value("workspace");
    if (!QFileInfo::exists(workspace + "/homology/examples/sources.all.json")) {
        workspace =
            QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) + "/workspace";
        QDir().mkpath(workspace + "/homology/examples");
        QFile source(":/geneorigin/sources.all.json");
        if (source.open(QIODevice::ReadOnly))
            homology::writeText(workspace + "/homology/examples/sources.all.json",
                                source.readAll());
    }
    if (p.isSet("cli")) {
        try {
            go::Request r;
            r.workspace = workspace;
            r.hal = p.value("hal");
            r.alignment = p.value("input");
            r.assemblies = p.value("assemblies");
            r.reference = p.value("reference");
            r.genes = p.value("genes").split(',', Qt::SkipEmptyParts);
            if (p.value("genes") == "all")
                r.genes = go::catalog(r.alignment).genes;
            r.comparisonSpecies = p.value("scope").split(',', Qt::SkipEmptyParts);
            r.outputRoot = p.value("output");
            r.useVerified = !p.isSet("fresh");
            r.blastBin = go::findBlastBin(workspace);
            r.halBin = go::findHalBin();
            homology::Hooks hooks;
            hooks.progress = [](QString m) { QTextStream(stderr) << m << Qt::endl; };
            auto result = go::execute(r, hooks, [](int n, QString m) {
                QTextStream(stderr) << n << "% " << m << Qt::endl;
            });
            QTextStream(stdout) << result.directory << Qt::endl;
            return 0;
        } catch (const std::exception &e) {
            QTextStream(stderr) << e.what() << Qt::endl;
            return 2;
        }
    }
    MainWindow w(workspace);
    app.setWindowIcon(w.windowIcon());
    w.show();
    if (p.isSet("demo-inputs"))
        w.loadProjectInputs();
    if (p.isSet("open-result"))
        w.loadSaved(p.value("open-result"));
    if (p.isSet("screenshot")) {
        QString path = p.value("screenshot");
        QTimer::singleShot(2200, &w, [&w, path] { w.saveScreenshot(path); });
    }
    if (p.isSet("smoke-test"))
        QTimer::singleShot(2600, &app, &QCoreApplication::quit);
    return app.exec();
}
