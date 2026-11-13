# Подключение ORF к приложению из трёх модулей

Общему приложению достаточно добавить каталог через `add_subdirectory(orf_cpp)`
и связать свою цель с `orf_core`. Текущее окно — рабочая оболочка для ORF с двумя
местами подключения остальных модулей; алгоритмы гомологии и геномного контекста
должны предоставляться соответствующими участниками команды.

```cmake
add_subdirectory(orf_cpp)
target_link_libraries(your_application PRIVATE orf_core)
```

```cpp
#include <orf/orf.hpp>

orf::Settings settings;
settings.input = "/data/supplementary_file_unannotated_alignments.txt";
settings.outputRoot = "/results/orf";
settings.finalOut = "/results/orf_profiles.jsonl";
settings.gene = ""; // Все гены.
settings.maxEvents = 200;
std::atomic_bool cancel{false};
auto report = orf::runPipeline(settings, [](const std::string& text) {
    // Callback вызывается в потоке вычислений. Передавайте обновления GUI через queued connection.
}, &cancel);
```

Вызов синхронный: запускайте его в рабочем потоке. В процессе должен существовать
`QCoreApplication` или `QApplication`. Ошибки передаются исключениями
`std::exception`; отмена также останавливает вызов исключением и фиксируется как
`cancelled` в `pipeline_run.json`. Отмена не прерывает текущее выравнивание.
На один каталог результатов действует межпроцессная блокировка.

## Контракт данных

Идентификатор гена — поле `gene_id`. Итоговый JSONL содержит одну запись на ген:

- `species_assessed`, `loci_assessed`, `comparison_records` — объём наблюдений;
- `species_status_counts` — распределение статусов;
- `species_profiles` — виды, их статусы и список локусов;
- `nontriplet_gap_observations`, `phase_restored_observations` — наблюдения по выравниваниям;
- `origin_classification: null` — классификацию выполняет объединяющий модуль;
- `limitations` — ограничения интерпретации.

Наблюдения разрывов могут повторяться по двум методам. Не складывайте их как
число независимых мутаций. `multiple_loci_uncertain` сохраняет неоднозначность,
автоматический выбор ортолога не выполняется.

Для детализации используйте `6_results_fifth/locus_profiles.csv`, затем
`5_results_fourth/evidence_records.jsonl` и `4_results_third/comparisons.jsonl`.
Идентичность локуса: **gene_id + species + locus**. Для записи отдельного метода
добавляется **method** (`spaln`/`genewise`). Сопоставление только по виду некорректно.

Координаты ORF и разрывов — **0-based, half-open** относительно входной
последовательности, а не генома. Знак рамки относится к переданной последовательности.
Геномные координаты заголовка сохраняются как исходные метаданные.

Промежуточные пути и имена файлов повторяют Python. Служебный manifest содержит
собственный C++-формат, несовместимый с Python. Финальный файл публикуется атомарно.
Изменение версии исполняемого файла инвалидирует кэш. При статическом связывании
`orf_core` это также учитывает изменения библиотеки.
