#!/usr/bin/env python3

"""
ORF Analysis — Stage 1
======================

Потоковый парсер supplementary_file_unannotated_alignments.txt.

Возможности:
    1. Чтение трёх секций: Raw, Spaln, GeneWise.
    2. Разбор заголовков последовательностей.
    3. Проверка структуры входных данных.
    4. Проверка допустимых нуклеотидов.
    5. Обработка одного гена или всего набора.
    6. Вывод прогресса в терминал.
    7. Сохранение JSONL, CSV и JSON-отчёта.

Требования:
    Python >= 3.10

Сторонние библиотеки не требуются.

ВАЖНО:
    На этом этапе поиск ORF не выполняется.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys

from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path


# ============================================================
# 1. CONSTANTS
# ============================================================

SECTION_NAMES = {
    "I": "raw",
    "II": "spaln",
    "III": "genewise",
}

SECTIONS = tuple(SECTION_NAMES.values())

SECTION_RE = re.compile(
    r"^Section\s+(I|II|III)\.\s+",
    re.IGNORECASE
)

GROUP_RE = re.compile(
    r"^(raw sequence|spaln|genewise) for dmel_(\S+)$",
    re.IGNORECASE
)

LOCUS_RE = re.compile(
    r"^(?P<contig>.+)(?P<strand>[+-]):"
    r"(?P<start>\d+):(?P<end>\d+)$"
)

# Допустимые IUPAC-символы ДНК
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")

# Однозначные нуклеотиды
UNAMBIGUOUS_DNA = frozenset("ACGT")

CSV_FIELDS = (
    "gene_id",
    "section",
    "species",
    "locus",
    "sequence_kind",
    "length_nt",
    "ambiguous_nt",
    "invalid_nt",
    "score_label",
    "shift_label",
    "exon_count",
    "contig",
    "strand",
    "start",
    "end",
    "header_line",
    "issues",
)


# ============================================================
# 2. EXCEPTIONS
# ============================================================

class InputFormatError(ValueError):
    """Ошибка структуры входного файла."""
    pass


# ============================================================
# 3. DATA MODEL
# ============================================================

@dataclass(slots=True)
class SequenceRecord:
    """
    Унифицированное представление одной последовательности.
    """

    gene_id: str
    section: str

    species: str
    locus: str
    sequence_kind: str

    sequence: str
    length_nt: int

    ambiguous_nt: int
    invalid_nt: int

    score_label: str | None
    shift_label: str | None
    exon_count: int | None

    contig: str | None
    strand: str | None
    start: int | None
    end: int | None

    original_header: str
    header_line: int

    issues: list[str]


# ============================================================
# 4. FASTA RECORD PARSER
# ============================================================

def parse_fasta_record(
    section: str,
    gene_id: str,
    header: str,
    header_line: int,
    fragments: list[str],
) -> SequenceRecord:

    """
    Преобразует FASTA-запись в SequenceRecord.

    Пример:

    >dmel|FBgn0262810|score_1E10|shift_No|1
    ATGAAATAA

    Также поддерживаются координатные заголовки:

    >dper|NW_020825322.1-:709401:710983
    ATGAAATAA
    """

    fields = header.removeprefix(">").split("|")

    issues = []

    # --------------------------------------------------------
    # Проверка заголовка
    # --------------------------------------------------------

    if len(fields) not in (2, 5):
        issues.append("header_wrong_field_count")

    if section != "raw" and len(fields) != 5:
        issues.append("prediction_missing_metadata")

    species = fields[0] if fields else ""
    locus = fields[1] if len(fields) >= 2 else ""

    if not species or not locus:
        issues.append("header_missing_species_or_locus")

    # --------------------------------------------------------
    # Геномные координаты
    # --------------------------------------------------------

    location = LOCUS_RE.fullmatch(locus)

    contig = None
    strand = None
    start = None
    end = None

    if location:

        contig = location.group("contig")
        strand = location.group("strand")

        start = int(location.group("start"))
        end = int(location.group("end"))

        if start >= end:
            issues.append(
                "genomic_locus_coordinates_invalid"
            )

    # --------------------------------------------------------
    # Тип последовательности
    # --------------------------------------------------------

    if section == "raw" and len(fields) == 2:

        sequence_kind = "genomic_region"

        if location is None:
            issues.append("genomic_locus_format_invalid")

    elif len(fields) == 5:

        if location:
            sequence_kind = "predicted_cds"
        else:
            sequence_kind = "annotated_cds"

    else:
        sequence_kind = "unknown"

    # --------------------------------------------------------
    # Метаданные предсказаний
    # --------------------------------------------------------

    score_label = None
    shift_label = None
    exon_count = None

    if len(fields) == 5:

        # Score
        if fields[2].startswith("score_"):

            score_label = fields[2][len("score_"):]

            try:
                float(score_label)
            except ValueError:
                issues.append("score_not_numeric")

        else:
            issues.append("score_field_invalid")

        # Shift
        if fields[3].startswith("shift_"):

            shift_label = fields[3][len("shift_"):]

            if shift_label not in ("Yes", "No"):
                issues.append("shift_field_invalid")

        else:
            issues.append("shift_field_invalid")

        # Exon count
        try:

            exon_count = int(fields[4])

            if exon_count < 1:
                issues.append("exon_count_invalid")

        except ValueError:
            issues.append("exon_count_invalid")

    # --------------------------------------------------------
    # Последовательность ДНК
    # --------------------------------------------------------

    sequence = "".join(fragments).upper()

    length_nt = len(sequence)

    invalid_nt = sum(
        nucleotide not in IUPAC_DNA
        for nucleotide in sequence
    )

    ambiguous_nt = sum(
        nucleotide in IUPAC_DNA - UNAMBIGUOUS_DNA
        for nucleotide in sequence
    )

    if not sequence:

        if section == "raw":
            issues.append("empty_sequence")
        else:
            issues.append("empty_prediction")

    if invalid_nt:
        issues.append("invalid_dna_characters")

    # --------------------------------------------------------
    # Результат
    # --------------------------------------------------------

    return SequenceRecord(
        gene_id=gene_id,
        section=section,
        species=species,
        locus=locus,
        sequence_kind=sequence_kind,
        sequence=sequence,
        length_nt=length_nt,
        ambiguous_nt=ambiguous_nt,
        invalid_nt=invalid_nt,
        score_label=score_label,
        shift_label=shift_label,
        exon_count=exon_count,
        contig=contig,
        strand=strand,
        start=start,
        end=end,
        original_header=header,
        header_line=header_line,
        issues=issues,
    )


# ============================================================
# 5. STREAMING FILE READER
# ============================================================

def iter_file_events(input_path: Path):

    """
    Потоковое чтение файла.

    Генерирует события:

        ("section", section_name)

        ("gene", (section_name, gene_id))

        ("record", SequenceRecord)

    В оперативной памяти находится только текущая
    последовательность, а не весь исходный файл.
    """

    section = None
    gene_id = None

    header = None
    header_line = 0

    fragments = []

    with Path(input_path).open(
        "r",
        encoding="utf-8"
    ) as file:

        for line_no, raw_line in enumerate(file, start=1):

            line = raw_line.strip()

            if not line:
                continue

            # ------------------------------------------------
            # Новый раздел
            # ------------------------------------------------

            if line.startswith("Section "):

                if header is not None:

                    yield (
                        "record",
                        parse_fasta_record(
                            section,
                            gene_id,
                            header,
                            header_line,
                            fragments
                        )
                    )

                    header = None
                    fragments = []

                match = SECTION_RE.match(line)

                if match is None:
                    raise InputFormatError(
                        f"Строка {line_no}: "
                        f"неизвестный раздел: {line}"
                    )

                section = SECTION_NAMES[
                    match.group(1).upper()
                ]

                gene_id = None

                yield ("section", section)

                continue

            # ------------------------------------------------
            # Новая группа гена
            # ------------------------------------------------

            group_match = GROUP_RE.fullmatch(line)

            if group_match:

                if section is None:
                    raise InputFormatError(
                        f"Строка {line_no}: "
                        "группа гена находится вне секции"
                    )

                if header is not None:

                    yield (
                        "record",
                        parse_fasta_record(
                            section,
                            gene_id,
                            header,
                            header_line,
                            fragments
                        )
                    )

                    header = None
                    fragments = []

                expected_group = {
                    "raw": "raw sequence",
                    "spaln": "spaln",
                    "genewise": "genewise",
                }[section]

                if group_match.group(1).lower() != expected_group:

                    raise InputFormatError(
                        f"Строка {line_no}: "
                        f"неверная группа для {section}"
                    )

                gene_id = group_match.group(2)

                yield (
                    "gene",
                    (section, gene_id)
                )

                continue

            # ------------------------------------------------
            # FASTA-заголовок
            # ------------------------------------------------

            if line.startswith(">"):

                if section is None or gene_id is None:

                    raise InputFormatError(
                        f"Строка {line_no}: "
                        "FASTA-запись вне группы гена"
                    )

                if header is not None:

                    yield (
                        "record",
                        parse_fasta_record(
                            section,
                            gene_id,
                            header,
                            header_line,
                            fragments
                        )
                    )

                header = line
                header_line = line_no

                fragments = []

                continue

            # ------------------------------------------------
            # Нуклеотидная последовательность
            # ------------------------------------------------

            if header is None:

                raise InputFormatError(
                    f"Строка {line_no}: "
                    f"неожиданный текст: {line[:100]}"
                )

            fragments.append(line)

    # --------------------------------------------------------
    # Последняя запись файла
    # --------------------------------------------------------

    if header is not None:

        yield (
            "record",
            parse_fasta_record(
                section,
                gene_id,
                header,
                header_line,
                fragments
            )
        )


# ============================================================
# 6. CSV CONVERTER
# ============================================================

def record_to_csv_row(record: SequenceRecord) -> dict:

    """
    Преобразует SequenceRecord в CSV-строку.

    Полная ДНК в CSV не записывается.
    Она сохраняется в JSONL.
    """

    values = asdict(record)

    values["issues"] = ";".join(record.issues)

    return {
        field: values[field]
        for field in CSV_FIELDS
    }


# ============================================================
# 7. MAIN PROCESSING FUNCTION
# ============================================================

def run(
    input_path: Path,
    records_path: Path,
    summary_path: Path,
    report_path: Path,
    gene_filter: str | None = None,
    progress_every: int = 100,
    verbose: bool = True,
) -> dict:

    """
    Основная обработка.

    Результаты:

        records.jsonl
        summary.csv
        validation.json

    Если указан gene_filter, сохраняются только записи
    выбранного гена.

    При этом входной файл просматривается полностью,
    поскольку один ген встречается в трёх секциях.
    """

    # --------------------------------------------------------
    # Подготовка путей
    # --------------------------------------------------------

    input_path = Path(input_path)

    records_path = Path(records_path)
    summary_path = Path(summary_path)
    report_path = Path(report_path)

    if progress_every < 0:
        raise ValueError(
            "progress_every не может быть отрицательным"
        )

    all_paths = [
        path.resolve()
        for path in (
            input_path,
            records_path,
            summary_path,
            report_path
        )
    ]

    if len(set(all_paths)) != 4:
        raise ValueError(
            "Входной и выходные файлы должны иметь разные пути"
        )

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Входной файл не найден: {input_path.resolve()}\n"
            f"Текущая директория: {Path.cwd()}"
        )

    # --------------------------------------------------------
    # Логирование
    # --------------------------------------------------------

    def log(message: str):

        if verbose:
            print(
                f"[ORF Parser] {message}",
                flush=True
            )

    # --------------------------------------------------------
    # Создание выходных директорий
    # --------------------------------------------------------

    for output_path in (
        records_path,
        summary_path,
        report_path
    ):

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True
        )

    log("=" * 55)

    log("Начало обработки:")
    log(str(input_path.resolve()))

    if gene_filter:
        log(f"Режим: один ген ({gene_filter})")
    else:
        log("Режим: все гены")

    log("=" * 55)

    # --------------------------------------------------------
    # Статистика
    # --------------------------------------------------------

    seen_sections = Counter()

    gene_groups = defaultdict(Counter)

    hits = defaultdict(
        lambda: defaultdict(Counter)
    )

    record_counts = Counter()

    warning_counts = Counter()
    error_counts = Counter()

    warning_examples = []
    error_examples = []

    groups_scanned = 0
    selected_records = 0

    ambiguous_records = 0
    ambiguous_bases = 0

    # --------------------------------------------------------
    # Обработка
    # --------------------------------------------------------

    with (
        records_path.open(
            "w",
            encoding="utf-8"
        ) as records_file,

        summary_path.open(
            "w",
            newline="",
            encoding="utf-8"
        ) as summary_file
    ):

        writer = csv.DictWriter(
            summary_file,
            fieldnames=CSV_FIELDS
        )

        writer.writeheader()

        for event, payload in iter_file_events(input_path):

            # ------------------------------------------------
            # Секция
            # ------------------------------------------------

            if event == "section":

                seen_sections[payload] += 1

                log(f"Чтение секции: {payload}")

                continue

            # ------------------------------------------------
            # Группа гена
            # ------------------------------------------------

            if event == "gene":

                section, gene_id = payload

                groups_scanned += 1

                if gene_filter is None or gene_filter == gene_id:

                    gene_groups[section][gene_id] += 1

                    # Регистрируем группу даже без записей
                    hits[gene_id][section]

                    if gene_filter is not None:

                        log(
                            f"Найден ген {gene_id}, "
                            f"секция {section}"
                        )

                if (
                    progress_every > 0
                    and groups_scanned % progress_every == 0
                ):

                    log(
                        f"Прогресс: обработано групп "
                        f"{groups_scanned}, "
                        f"выбрано записей {selected_records}"
                    )

                continue

            # ------------------------------------------------
            # Последовательность
            # ------------------------------------------------

            record = payload

            if (
                gene_filter is not None
                and record.gene_id != gene_filter
            ):
                continue

            selected_records += 1

            record_counts[record.section] += 1

            hits[
                record.gene_id
            ][
                record.section
            ][
                record.species
            ] += 1

            # ------------------------------------------------
            # Неоднозначные нуклеотиды
            # ------------------------------------------------

            if record.ambiguous_nt:

                ambiguous_records += 1
                ambiguous_bases += record.ambiguous_nt

            # ------------------------------------------------
            # Ошибки и предупреждения
            # ------------------------------------------------

            for issue in record.issues:

                is_warning = issue == "empty_prediction"

                if is_warning:

                    warning_counts[issue] += 1

                    if len(warning_examples) < 15:

                        warning_examples.append({
                            "line": record.header_line,
                            "gene": record.gene_id,
                            "section": record.section,
                            "species": record.species,
                            "issue": issue,
                        })

                else:

                    error_counts[issue] += 1

                    if len(error_examples) < 15:

                        error_examples.append({
                            "line": record.header_line,
                            "gene": record.gene_id,
                            "section": record.section,
                            "species": record.species,
                            "issue": issue,
                        })

            # ------------------------------------------------
            # JSONL
            # ------------------------------------------------

            records_file.write(
                json.dumps(
                    asdict(record),
                    ensure_ascii=False
                ) + "\n"
            )

            # ------------------------------------------------
            # CSV
            # ------------------------------------------------

            writer.writerow(
                record_to_csv_row(record)
            )

    # ========================================================
    # 8. CROSS-SECTION VALIDATION
    # ========================================================

    selected_genes = set(hits)

    for gene_id, by_section in hits.items():

        # Проверка наличия всех секций для гена
        for section in SECTIONS:

            if section not in by_section:

                error_counts["gene_missing_in_section"] += 1

                if len(error_examples) < 15:

                    error_examples.append({
                        "gene": gene_id,
                        "section": section,
                        "issue": "gene_missing_in_section",
                    })

        # Сравнение количества записей каждого вида
        if all(
            section in by_section
            for section in SECTIONS
        ):

            for method in ("spaln", "genewise"):

                if by_section["raw"] != by_section[method]:

                    error_counts[
                        "species_hit_count_mismatch"
                    ] += 1

                    if len(error_examples) < 15:

                        error_examples.append({
                            "gene": gene_id,
                            "section": method,
                            "issue": "species_hit_count_mismatch",
                        })

    # --------------------------------------------------------
    # Проверка секций файла
    # --------------------------------------------------------

    for section in SECTIONS:

        if seen_sections[section] != 1:

            error_counts[
                "missing_or_repeated_section"
            ] += 1

        for count in gene_groups[section].values():

            if count > 1:

                error_counts[
                    "repeated_gene_group"
                ] += 1

    # --------------------------------------------------------
    # Проверка выбранного гена
    # --------------------------------------------------------

    if (
        gene_filter is not None
        and gene_filter not in selected_genes
    ):

        error_counts[
            "requested_gene_not_found"
        ] += 1

    # ========================================================
    # 9. VALIDATION REPORT
    # ========================================================

    report = {

        "input_file": str(input_path.resolve()),

        "gene_filter": gene_filter,

        "sections_seen": {
            section: seen_sections[section]
            for section in SECTIONS
        },

        "genes_per_section": {
            section: len(gene_groups[section])
            for section in SECTIONS
        },

        "unique_genes": len(selected_genes),

        "records_per_section": {
            section: record_counts[section]
            for section in SECTIONS
        },

        "selected_records": selected_records,

        "groups_scanned": groups_scanned,

        "records_with_iupac_ambiguity": ambiguous_records,

        "total_ambiguous_bases": ambiguous_bases,

        "warning_count": sum(warning_counts.values()),

        "warning_counts": dict(warning_counts),

        "warning_examples": warning_examples,

        "error_count": sum(error_counts.values()),

        "error_counts": dict(error_counts),

        "error_examples": error_examples,
    }

    # --------------------------------------------------------
    # Сохранение отчёта
    # --------------------------------------------------------

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ) + "\n",
        encoding="utf-8"
    )

    # ========================================================
    # 10. FINAL LOG
    # ========================================================

    log("=" * 55)

    log(
        f"Готово: {len(selected_genes)} генов, "
        f"{selected_records} записей"
    )

    log(
        f"Ошибок: {report['error_count']}"
    )

    log(
        f"Предупреждений: {report['warning_count']}"
    )

    log("Выходные файлы:")

    log(
        f"JSONL: {records_path.resolve()}"
    )

    log(
        f"CSV:   {summary_path.resolve()}"
    )

    log(
        f"JSON:  {report_path.resolve()}"
    )

    log("=" * 55)

    return report


# ============================================================
# 11. COMMAND LINE INTERFACE
# ============================================================

def main(argv: list[str] | None = None) -> int:

    parser = argparse.ArgumentParser(
        description=(
            "ORF Analysis Stage 1: "
            "парсер и проверка исходных данных"
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Путь к входному TXT-файлу"
    )

    parser.add_argument(
        "--records-out",
        required=True,
        type=Path,
        help="Путь к выходному JSONL"
    )

    parser.add_argument(
        "--summary-out",
        required=True,
        type=Path,
        help="Путь к выходному CSV"
    )

    parser.add_argument(
        "--report-out",
        required=True,
        type=Path,
        help="Путь к JSON-отчёту"
    )

    parser.add_argument(
        "--gene",
        default=None,
        help="ID гена; без параметра обрабатываются все гены"
    )

    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help=(
            "Показывать прогресс каждые N групп генов; "
            "0 — отключить периодические сообщения"
        )
    )

    args = parser.parse_args(argv)

    try:

        report = run(
            input_path=args.input,
            records_path=args.records_out,
            summary_path=args.summary_out,
            report_path=args.report_out,
            gene_filter=args.gene,
            progress_every=args.progress_every,
            verbose=True,
        )

    except (
        OSError,
        InputFormatError,
        ValueError
    ) as error:

        print(
            f"[ORF Parser] ОШИБКА: {error}",
            file=sys.stderr,
            flush=True
        )

        return 2

    if report["error_count"] > 0:
        return 1

    return 0


# ============================================================
# 12. ENTRY POINT
# ============================================================

if __name__ == "__main__":
    raise SystemExit(main())
