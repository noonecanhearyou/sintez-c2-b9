# Подключение к Qt-приложению

Библиотека зависит от Qt Core и не создаёт окно. Для подключения из родительского CMake-проекта:

```cmake
add_subdirectory(homology)
target_link_libraries(your_qt_app PRIVATE homology)
```

Синхронный API `homology::run(configPath, hooks)` должен выполняться вне GUI-потока. Он возвращает QJsonObject манифеста, бросает `homology::Error` при фатальной ошибке/отмене и сохраняет частичные ошибки отдельных целей в результате.

Готовый адаптер для QThread:

```cpp
#include <homology/job.hpp>
#include <QThread>

// Этот блок выполняется в GUI-потоке; window — существующее окно.
auto* thread = new QThread(window);
auto* job = new homology::HomologyJob(configPath);
job->moveToThread(thread);
QObject::connect(thread, &QThread::started, job, &homology::HomologyJob::start);
QObject::connect(job, &homology::HomologyJob::progress,
                 window, &MainWindow::showHomologyProgress);
QObject::connect(job, &homology::HomologyJob::succeeded,
                 window, &MainWindow::loadHomologyResults);
QObject::connect(job, &homology::HomologyJob::failed,
                 window, &MainWindow::showHomologyError);
QObject::connect(job, &homology::HomologyJob::finished,
                 thread, &QThread::quit);
QObject::connect(job, &homology::HomologyJob::finished,
                 job, &QObject::deleteLater);
QObject::connect(thread, &QThread::finished,
                 thread, &QObject::deleteLater);
thread->start();
```

Указанные методы MainWindow — методы вашего приложения; это пример подключения, не отдельное готовое GUI.

`succeeded` означает, что функция вернула манифест: он может содержать `partial_error`. Окно обязано проверить `status` перед показом успеха. Для таблиц подходит собственная QAbstractTableModel; не загружайте многомиллионные HSP-таблицы целиком в виджеты.

Отмена: пока объект job гарантированно жив, вызовите `job->requestCancel()` напрямую. Метод меняет только atomic-флаг и допускается из другого потока. Не отправляйте отмену обычным queued-slot в занятый рабочий поток: тот не обрабатывает очередь во время поиска. Для безопасного жизненного цикла храните активную задачу в контроллере, отключайте кнопку отмены после завершения и перед уничтожением окна отмените задачу и дождитесь потока. Не уничтожайте работающий QThread. Проверка отмены во время внешнего процесса выполняется примерно каждые 100 мс; чтение/хэширование входов не прерывается посередине.

Для загрузки данных `prepare()` также следует использовать рабочий поток. Она поддерживает те же Hooks и не требует Qt Network: внешние curl/gzip запускаются через QProcess без оболочки. В будущем её можно заменить загрузчиком приложения, сохранив формат source.json и входные FASTA/TSV.
