# Проверка C++-модуля и сравнение с Python

Все команды ниже выполняются из корня `orf_cpp`.
Для запуска приложения Python не нужен; для проверки эталона используется
окружение `/home/stass/miniconda3/envs/hac_orf_stas/bin/python` с Biopython 1.88.

## 1. Автоматические проверки

```bash
export LD_LIBRARY_PATH="$PWD/runtime/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export QT_PLUGIN_PATH="$PWD/runtime/plugins"
ctest --test-dir build --output-on-failure
```

Ожидается **3/3 PASS**:

- `core_algorithms`: вложенные и обратные ORF, IUPAC, границы, CDS, выравнивание, лимиты;
- `pipeline_integration`: реальный запуск пяти этапов на искусственных данных,
  повторное использование кэша, повреждение файла, изменение параметров,
  неизвестный ген, ошибочный ввод, защита ранее опубликованного результата;
- `gui_startup`: запуск Qt GUI с offscreen-платформой.

Эта проверка запуска GUI не заменяет ручную проверку выбора файлов и отмены.

Проверка исходного Python-стенда:

```bash
PYTHONDONTWRITEBYTECODE=1 /home/stass/miniconda3/envs/hac_orf_stas/bin/python \
  -m pytest -q -p no:cacheprovider ../test_orf
```

На использованной версии стенда: **66 passed**.

Дифференциальные проверки функций:

```bash
PYTHONDONTWRITEBYTECODE=1 /home/stass/miniconda3/envs/hac_orf_stas/bin/python \
  tests/differential.py --binary "$PWD/build/orf_workbench" \
  --reference ../test_orf --out tests/differential_result.json
```

Ожидается **6382 passed, 0 failed**. Проверяются все поля результатов, а не
только итоговые статусы. В набор входят искусственные случаи, случайные
последовательности с фиксированным seed, варианты лимитов и превью,
перебор коротких пар с равными оценками выравниваний, исходные тестовые
сценарии агрегации локусов и регрессии округления.

## 2. Проверка одного гена

Создайте отдельный Python-эталон с теми же параметрами:

```bash
/home/stass/miniconda3/envs/hac_orf_stas/bin/python ../test_orf/run_orf_pipeline.py \
  --input ../test_orf/1_raw_data/supplementary_file_unannotated_alignments.txt \
  --output-root validation/python_single \
  --final-out validation/python_FBgn0262810.jsonl \
  --gene FBgn0262810 --mode fresh --max-events 200

./run.sh --cli \
  --input ../test_orf/1_raw_data/supplementary_file_unannotated_alignments.txt \
  --output-root validation/single \
  --final-out validation/FBgn0262810.jsonl \
  --gene FBgn0262810 --mode fresh --max-events 200

python3 tools/compare_results.py \
  --python-root validation/python_single --cpp-root validation/single \
  --out validation/comparison_single.json
```

Ожидается **19 passed, 0 failed**. В GUI откройте оба итоговых JSONL и сравните
поля `species_status_counts`, `species_profiles`, число локусов и наблюдений.
Для ручной проверки конкретного вида найдите его в `species_profiles`, затем
в `locus_profiles.csv`, затем соответствующие записи двух методов в
`evidence_records.jsonl` и `comparisons.jsonl`.

Для `FBgn0262810` при анализе Spaln ожидаются, в частности, CDS dsim/dsec
длиной 387 nt без разрывов и CDS dyak длиной 384 nt с внутренним разрывом 3 nt.
Итоговый статус вида учитывает и GeneWise: не выводите его только по Spaln.

## 3. Проверка всех генов

В `../test_orf/pipeline_results_all` уже сохранён полный Python-эталон.
Параметры сверяются по его `orf_report.json` и `comparison_report.json`:
min_aa=25, top_k=10, max_peptide_aa=120, max_align_length=4000,
max_align_cells=4000000, min_identity=0.6, min_coverage=0.6,
max_events=200, preview_columns=90, max_gap_details=10.

```bash
./run.sh --cli \
  --input ../test_orf/1_raw_data/supplementary_file_unannotated_alignments.txt \
  --output-root validation/full \
  --final-out validation/all_gene_profiles.jsonl \
  --mode fresh --max-events 200

python3 tools/compare_results.py \
  --python-root ../test_orf/pipeline_results_all --cpp-root validation/full \
  --out validation/comparison.json
```

Ожидается **19 passed, 0 failed**. Сравнение потоковое; весь набор JSONL
не загружается в память. Проверяются порядок и число записей, каждый объект,
каждая CSV-ячейка и все отчётные показатели. Из JSON-отчётов исключены только
поля путей `input_file`, `orf_input`, `comparisons_input`, `agreement_input`.
Для JSON допускается другая расстановка ключей и пробелов; для чисел CSV —
эквивалентная десятичная запись без численного допуска.

Контрольные объёмы:

| Показатель | Ожидается |
|---|---:|
| Гены | 780 |
| Исходные записи | 51 033 |
| Записи в каждой секции | 17 011 |
| Межвидовые сравнения | 32 354 |
| Локусы | 16 177 |
| Сочетания ген — вид | 12 954 |
| Локусы preserved_both_identical | 3 090 |
| Локусы conflicting_alignment | 3 616 |
| Локусы insufficient_evidence | 8 870 |

## 4. Resume и обработка ошибок

Повторите полную команду с `--mode resume`: ожидаются
`computed_stages=0`, `reused_stages=5` в `validation/full/pipeline_run.json`.
Не меняйте исполняемый файл между запусками: его SHA-256 входит в подпись кэша.

В отдельной тестовой копии каталога удалите или измените `comparisons.jsonl`.
При `resume` ожидается повторное использование первых двух этапов и пересчёт
последних трёх. Изменение параметров запуска пересчитывает все этапы.
Запуск с неизвестным геном завершается ошибкой, предыдущий итоговый файл
остаётся доступным. Два запуска в один каталог одновременно запрещены.

Ручная проверка GUI: запуск одного гена → просмотр профиля → полный запуск →
«Отменить» → дождаться остановки → повторный запуск в `resume`.
При отмене журнал объясняет остановку, `pipeline_run.json` получает `cancelled`.
Прогресс показывается как число обработанных записей и текущий этап;
процент/оставшееся время не оцениваются.

## 5. Результаты проведённой проверки

Машиночитаемые результаты находятся в:

- `tests/differential_result.json` — сравнение функций;
- `validation/comparison.json` — полная сверка 19 файлов;
- `validation/comparison_single.json` — сверка одного гена;
- `validation/benchmark.json` — время и пиковая память полного C++-прогона;
- `validation/ctest.log` — нативные/интеграционные проверки;
- `validation/full/pipeline_run.json` — состояние последнего полного запуска;
- `validation/all_gene_profiles.jsonl` — итоговые профили всех генов.

Результаты относятся к конкретным исходникам, SHA-256 которых сохранены в
`docs/reference_sources.json`, и установленному Biopython 1.88. Полное совпадение
на этих данных и тестах не является доказательством для всех возможных входов.
Изменение алгоритма Python или версии выравнивателя требует повторной сверки.

Фактический результат проверки 12 ноября 2026 года:

| Проверка | Результат |
|---|---:|
| Полный набор, все пять этапов | 19/19 файлов совпали |
| FBgn0262810, все пять этапов | 19/19 файлов совпали |
| Дифференциальные проверки функций | 6382/6382 |
| Тесты исходного Python-стенда | 66/66 |
| Нативные/интеграционные/GUI-проверки CTest | 3/3 |
| Resume полного набора | 0 пересчитано, 5 использовано |
| Время полного C++-прогона | 107.781 с |
| Пиковая RSS полного C++-прогона | 142.85 МиБ |

Измерение времени и памяти выполнено на текущем компьютере; это не отдельный
контролируемый сравнительный benchmark Python/C++. Сводные результаты и
SHA-256 проверенного приложения: `validation/verification.json`.
