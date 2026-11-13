#!/usr/bin/env python3

"""
ORF Analysis — Stage 2
======================

Поиск ORF и проверка CDS.

Вход:
    JSONL, созданный orf_parser.py.

Выход:
    orf_records.jsonl
    orf_summary.csv
    orf_report.json

Особенности:
    - Потоковая обработка.
    - Поиск ORF в шести рамках.
    - Поиск вложенных ORF.
    - Проверка предсказанных CDS.
    - Обработка неоднозначных нуклеотидов.
    - Вывод прогресса.
    - Ограничение количества сохраняемых ORF.

Python >= 3.10
Biopython
"""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import sys

from collections import Counter
from pathlib import Path

from Bio.Seq import Seq


# ============================================================
# 1. CONSTANTS
# ============================================================

START_CODON = "ATG"

STOP_CODONS = frozenset({
    "TAA",
    "TAG",
    "TGA",
})

DNA_BASES = frozenset("ACGT")

IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")

COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)

ORF_TYPES = (
    "complete",
    "open_end",
    "ambiguous_boundary",
)

CDS_KINDS = {
    "annotated_cds",
    "predicted_cds",
}

SUMMARY_FIELDS = [
    "gene_id",
    "species",
    "section",
    "locus",
    "sequence_kind",
    "length_nt",
    "ambiguous_nt",
    "invalid_nt",
    "reported_shift",
    "reported_score",
    "cds_status",
    "cds_start_atg",
    "cds_terminal_stop",
    "cds_length_mod3",
    "cds_internal_stop_count",
    "cds_ambiguous_codons",
    "complete_orfs_all",
    "complete_orfs_min_aa",
    "open_end_orfs_min_aa",
    "ambiguous_boundary_min_aa",
    "longest_complete_aa",
    "top_orf_frame",
    "top_orf_type",
    "top_orf_aa",
    "top_orf_start_nt",
    "top_orf_end_nt",
]


# ============================================================
# 2. SEQUENCE UTILITIES
# ============================================================

def reverse_complement(sequence: str) -> str:
    """
    Возвращает обратную комплементарную последовательность.

    Поддерживает неоднозначные IUPAC-нуклеотиды.
    """

    sequence = sequence.upper()

    return sequence.translate(
        COMPLEMENT_TABLE
    )[::-1]


def count_ambiguous_bases(sequence: str) -> int:
    """
    Количество допустимых, но неоднозначных
    IUPAC-нуклеотидов.
    """

    return sum(
        base in IUPAC_DNA and base not in DNA_BASES
        for base in sequence
    )


def count_invalid_bases(sequence: str) -> int:
    """
    Количество символов вне IUPAC DNA.
    """

    return sum(
        base not in IUPAC_DNA
        for base in sequence
    )


# ============================================================
# 3. ORF SCANNER
# ============================================================

def scan_orfs(
    sequence: str,
    min_aa: int = 25,
    top_k: int = 10,
    max_peptide_aa: int = 120,
) -> dict:

    """
    Ищет ORF в шести рамках считывания.

    Правила:
        - Старт: ATG.
        - Стоп: TAA, TAG, TGA.
        - Каждый ATG рассматривается как возможный старт.
        - ORF завершается первым стоп-кодоном
          в соответствующей рамке.
        - Неоднозначный кодон прерывает поиск.
        - Неполные ORF учитываются отдельно.

    Координаты:
        0-based, half-open.

    Знак рамки относится к входной последовательности.

    min_aa:
        Минимальное количество аминокислот.
        Терминальный стоп-кодон не учитывается.

    top_k:
        Максимальное количество сохраняемых ORF.

    Все ORF подсчитываются, но сохраняются только
    top_k наиболее длинных, удовлетворяющих min_aa.
    """

    if min_aa < 1:
        raise ValueError("min_aa должен быть >= 1")

    if top_k < 0:
        raise ValueError("top_k не может быть отрицательным")

    if max_peptide_aa < 0:
        raise ValueError(
            "max_peptide_aa не может быть отрицательным"
        )

    sequence = sequence.upper()

    sequence_length = len(sequence)

    # --------------------------------------------------------
    # Статистика
    # --------------------------------------------------------

    counts_all = Counter({
        kind: 0
        for kind in ORF_TYPES
    })

    counts_passing = Counter({
        kind: 0
        for kind in ORF_TYPES
    })

    frame_counts = {
        f"{strand}{offset + 1}": {
            kind: 0
            for kind in ORF_TYPES
        }
        for strand in ("+", "-")
        for offset in range(3)
    }

    longest_complete_aa = 0

    # Куча хранит только лучшие ORF.
    top_heap = []

    serial = 0

    # --------------------------------------------------------
    # Регистрация найденной ORF
    # --------------------------------------------------------

    def register_orf(
        start: int,
        end: int,
        frame: str,
        orf_type: str,
        stop_codon: str | None,
        strand: str,
    ):

        nonlocal serial
        nonlocal longest_complete_aa

        # Полная ORF содержит терминальный стоп.
        # Он не кодирует аминокислоту.

        total_codons = (end - start) // 3

        if orf_type == "complete":
            aa_length = total_codons - 1
        else:
            aa_length = total_codons

        if aa_length < 0:
            return

        counts_all[orf_type] += 1

        if orf_type == "complete":
            longest_complete_aa = max(
                longest_complete_aa,
                aa_length,
            )

        if aa_length < min_aa:
            return

        counts_passing[orf_type] += 1

        frame_counts[frame][orf_type] += 1

        # ----------------------------------------------------
        # Преобразование координат
        # ----------------------------------------------------

        if strand == "+":

            original_start = start
            original_end = end

        else:

            original_start = sequence_length - end
            original_end = sequence_length - start

        candidate = {
            "frame": frame,
            "strand_relative_to_input": strand,
            "type": orf_type,
            "start_nt": original_start,
            "end_nt": original_end,
            "length_nt": end - start,
            "aa_length": aa_length,
            "start_codon": START_CODON,
            "stop_codon": stop_codon,
        }

        # ----------------------------------------------------
        # Сохраняем только top_k ORF
        # ----------------------------------------------------

        if top_k == 0:
            return

        serial += 1

        # При равной длине предпочитаем полную ORF.
        # Это правило сортировки, а не биологическая
        # оценка вероятности кодирования.

        priority = (
            aa_length,
            1 if orf_type == "complete" else 0,
            -original_start,
            -original_end,
        )

        item = (
            priority,
            serial,
            candidate,
        )

        if len(top_heap) < top_k:

            heapq.heappush(
                top_heap,
                item,
            )

        elif priority > top_heap[0][0]:

            heapq.heapreplace(
                top_heap,
                item,
            )

    # ========================================================
    # 4. SIX READING FRAMES
    # ========================================================

    orientations = (
        ("+", sequence),
        ("-", reverse_complement(sequence)),
    )

    for strand, oriented_sequence in orientations:

        n = len(oriented_sequence)

        for offset in range(3):

            frame = f"{strand}{offset + 1}"

            # Все ATG после последнего стоп-кодона
            # или неоднозначного кодона.

            active_starts = []

            position = offset

            while position + 3 <= n:

                codon = oriented_sequence[
                    position:position + 3
                ]

                # --------------------------------------------
                # Неоднозначный кодон
                # --------------------------------------------

                if any(
                    base not in DNA_BASES
                    for base in codon
                ):

                    for start in active_starts:

                        register_orf(
                            start=start,
                            end=position,
                            frame=frame,
                            orf_type="ambiguous_boundary",
                            stop_codon=None,
                            strand=strand,
                        )

                    active_starts.clear()

                # --------------------------------------------
                # Стоп-кодон
                # --------------------------------------------

                elif codon in STOP_CODONS:

                    for start in active_starts:

                        register_orf(
                            start=start,
                            end=position + 3,
                            frame=frame,
                            orf_type="complete",
                            stop_codon=codon,
                            strand=strand,
                        )

                    active_starts.clear()

                # --------------------------------------------
                # Старт-кодон
                # --------------------------------------------

                elif codon == START_CODON:

                    active_starts.append(position)

                position += 3

            # -----------------------------------------------
            # ATG без последующего стоп-кодона
            # -----------------------------------------------

            for start in active_starts:

                register_orf(
                    start=start,
                    end=position,
                    frame=frame,
                    orf_type="open_end",
                    stop_codon=None,
                    strand=strand,
                )

    # ========================================================
    # 5. TOP ORFS
    # ========================================================

    sorted_candidates = sorted(
        top_heap,
        key=lambda item: (item[0], item[1]),
        reverse=True,
    )

    top_orfs = []

    for _, _, candidate in sorted_candidates:

        start = candidate["start_nt"]
        end = candidate["end_nt"]

        strand = candidate[
            "strand_relative_to_input"
        ]

        dna = sequence[start:end]

        if strand == "-":
            dna = reverse_complement(dna)

        # Все отобранные участки кратны трём.
        peptide = str(
            Seq(dna).translate(table=1)
        )

        if candidate["type"] == "complete":
            peptide = peptide.removesuffix("*")

        candidate["peptide_preview"] = (
            peptide[:max_peptide_aa]
        )

        candidate["peptide_truncated"] = (
            len(peptide) > max_peptide_aa
        )

        top_orfs.append(candidate)

    return {
        "min_aa": min_aa,
        "counts_all": dict(counts_all),
        "counts_passing_min_aa": dict(counts_passing),
        "counts_by_frame_passing_min_aa": frame_counts,
        "longest_complete_aa": longest_complete_aa,
        "top_orfs": top_orfs,
    }


# ============================================================
# 6. CDS VALIDATOR
# ============================================================

def validate_cds(
    sequence: str,
    sequence_kind: str,
) -> dict:

    """
    Проверяет структуру предполагаемой CDS.

    Проверка производится в рамке, начинающейся
    с первого нуклеотида переданной последовательности.

    Не пытается искать альтернативный старт.

    Не доказывает, что последовательность реально
    транслируется в клетке.
    """

    if sequence_kind not in CDS_KINDS:

        return {
            "status": "not_applicable",
            "reason_codes": [],
        }

    sequence = sequence.upper()

    n = len(sequence)

    if n == 0:

        return {
            "status": "empty",
            "reason_codes": ["empty_cds"],
            "start_atg": False,
            "terminal_stop": False,
            "length_mod3": 0,
            "internal_stop_count": 0,
            "internal_stop_positions": [],
            "ambiguous_codon_count": 0,
            "invalid_nt": 0,
        }

    # --------------------------------------------------------
    # Общие характеристики
    # --------------------------------------------------------

    remainder = n % 3

    complete_codon_end = n - remainder

    start_atg = (
        n >= 3
        and sequence[:3] == START_CODON
    )

    terminal_stop = (
        remainder == 0
        and n >= 3
        and sequence[-3:] in STOP_CODONS
    )

    internal_stop_positions = []

    ambiguous_codon_count = 0

    # --------------------------------------------------------
    # Анализ кодонов
    # --------------------------------------------------------

    for position in range(
        0,
        complete_codon_end,
        3,
    ):

        codon = sequence[
            position:position + 3
        ]

        if any(
            base not in DNA_BASES
            for base in codon
        ):

            ambiguous_codon_count += 1

        if codon in STOP_CODONS:

            is_terminal = (
                terminal_stop
                and position == n - 3
            )

            if not is_terminal:

                internal_stop_positions.append(
                    position
                )

    # --------------------------------------------------------
    # Формирование предупреждений
    # --------------------------------------------------------

    reasons = []

    if remainder != 0:
        reasons.append("length_not_multiple_of_3")

    if not start_atg:
        reasons.append("start_not_atg")

    if not terminal_stop:
        reasons.append("terminal_stop_missing")

    if internal_stop_positions:
        reasons.append("internal_stop_codons")

    if ambiguous_codon_count:
        reasons.append("ambiguous_codons")

    invalid_nt = count_invalid_bases(sequence)

    if invalid_nt:
        reasons.append("invalid_characters")

    # --------------------------------------------------------
    # Итоговый статус
    # --------------------------------------------------------

    if internal_stop_positions:

        status = "interrupted"

    elif invalid_nt or ambiguous_codon_count:

        status = "uncertain"

    elif remainder != 0 or not start_atg or not terminal_stop:

        status = "partial_or_noncanonical"

    else:

        status = "complete_candidate"

    return {
        "status": status,
        "reason_codes": reasons,
        "start_atg": start_atg,
        "terminal_stop": terminal_stop,
        "length_mod3": remainder,
        "internal_stop_count": len(
            internal_stop_positions
        ),
        "internal_stop_positions": (
            internal_stop_positions[:50]
        ),
        "internal_stop_positions_truncated": (
            len(internal_stop_positions) > 50
        ),
        "ambiguous_codon_count": (
            ambiguous_codon_count
        ),
        "invalid_nt": invalid_nt,
    }


# ============================================================
# 7. SUMMARY ROW
# ============================================================

def build_summary_row(
    record: dict,
    scan: dict,
    cds: dict,
) -> dict:

    """
    Преобразует результат анализа одной последовательности
    в плоскую CSV-строку.
    """

    top = (
        scan["top_orfs"][0]
        if scan["top_orfs"]
        else {}
    )

    passed = scan["counts_passing_min_aa"]

    all_counts = scan["counts_all"]

    return {
        "gene_id": record["gene_id"],
        "species": record["species"],
        "section": record["section"],
        "locus": record["locus"],
        "sequence_kind": record["sequence_kind"],
        "length_nt": len(record["sequence"]),
        "ambiguous_nt": count_ambiguous_bases(
            record["sequence"]
        ),
        "invalid_nt": count_invalid_bases(
            record["sequence"]
        ),
        "reported_shift": record.get("shift_label"),
        "reported_score": record.get("score_label"),
        "cds_status": cds["status"],
        "cds_start_atg": cds.get("start_atg"),
        "cds_terminal_stop": cds.get("terminal_stop"),
        "cds_length_mod3": cds.get("length_mod3"),
        "cds_internal_stop_count": cds.get(
            "internal_stop_count"
        ),
        "cds_ambiguous_codons": cds.get(
            "ambiguous_codon_count"
        ),
        "complete_orfs_all": all_counts["complete"],
        "complete_orfs_min_aa": passed["complete"],
        "open_end_orfs_min_aa": passed["open_end"],
        "ambiguous_boundary_min_aa": passed[
            "ambiguous_boundary"
        ],
        "longest_complete_aa": scan[
            "longest_complete_aa"
        ],
        "top_orf_frame": top.get("frame"),
        "top_orf_type": top.get("type"),
        "top_orf_aa": top.get("aa_length"),
        "top_orf_start_nt": top.get("start_nt"),
        "top_orf_end_nt": top.get("end_nt"),
    }


# ============================================================
# 8. MAIN PROCESSING
# ============================================================

def run(
    input_path: Path,
    records_path: Path,
    summary_path: Path,
    report_path: Path,
    gene_filter: str | None = None,
    min_aa: int = 25,
    top_k: int = 10,
    max_peptide_aa: int = 120,
    progress_every: int = 1000,
    verbose: bool = True,
) -> dict:

    """
    Потоково обрабатывает JSONL первого этапа.
    """

    input_path = Path(input_path)

    records_path = Path(records_path)
    summary_path = Path(summary_path)
    report_path = Path(report_path)

    if min_aa < 1:
        raise ValueError("--min-aa должен быть >= 1")

    if top_k < 0:
        raise ValueError("--top-k должен быть >= 0")

    if max_peptide_aa < 0:
        raise ValueError(
            "--max-peptide-aa должен быть >= 0"
        )

    if progress_every < 0:
        raise ValueError(
            "--progress-every должен быть >= 0"
        )

    paths = [
        input_path.resolve(),
        records_path.resolve(),
        summary_path.resolve(),
        report_path.resolve(),
    ]

    if len(set(paths)) != 4:
        raise ValueError(
            "Входной и выходные файлы должны иметь разные пути"
        )

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Входной файл не найден: {input_path.resolve()}"
        )

    def log(message: str):

        if verbose:
            print(
                f"[ORF Scanner] {message}",
                flush=True,
            )

    # --------------------------------------------------------
    # Выходные директории
    # --------------------------------------------------------

    for path in (
        records_path,
        summary_path,
        report_path,
    ):

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    # --------------------------------------------------------
    # Статистика
    # --------------------------------------------------------

    processed = 0
    scanned_lines = 0

    genes = set()

    section_counts = Counter()
    cds_status_counts = Counter()

    records_with_complete_orf = 0
    records_with_ambiguous_bases = 0
    records_with_invalid_bases = 0

    total_complete_orfs = 0
    total_complete_orfs_passing = 0

    log("=" * 55)

    log("Начало ORF Analysis")
    log(f"Вход: {input_path.resolve()}")
    log(f"Минимальная ORF: {min_aa} aa")
    log(f"Сохраняем top {top_k} ORF на запись")

    if gene_filter:
        log(f"Ген: {gene_filter}")
    else:
        log("Режим: все гены")

    # ========================================================
    # 9. STREAMING PROCESSING
    # ========================================================

    with (
        input_path.open(
            "r",
            encoding="utf-8",
        ) as input_file,

        records_path.open(
            "w",
            encoding="utf-8",
        ) as records_file,

        summary_path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as summary_file
    ):

        writer = csv.DictWriter(
            summary_file,
            fieldnames=SUMMARY_FIELDS,
        )

        writer.writeheader()

        for line_number, line in enumerate(
            input_file,
            start=1,
        ):

            if not line.strip():
                continue

            scanned_lines += 1

            try:
                record = json.loads(line)

            except json.JSONDecodeError as error:

                raise ValueError(
                    f"Некорректный JSON в строке "
                    f"{line_number}: {error}"
                ) from error

            if not isinstance(record, dict):
                raise ValueError(
                    f"Строка {line_number}: "
                    "ожидался JSON-объект"
                )

            required = (
                "gene_id",
                "species",
                "section",
                "locus",
                "sequence_kind",
                "sequence",
            )

            for field in required:

                if (
                    field not in record
                    or not isinstance(record[field], str)
                ):

                    raise ValueError(
                        f"Строка {line_number}: "
                        f"нет строкового поля {field}"
                    )

            if (
                gene_filter is not None
                and record["gene_id"] != gene_filter
            ):
                continue

            sequence = record["sequence"].upper()

            record["sequence"] = sequence

            processed += 1

            genes.add(record["gene_id"])

            section_counts[
                record["section"]
            ] += 1

            # -----------------------------------------------
            # Поиск ORF
            # -----------------------------------------------

            scan = scan_orfs(
                sequence=sequence,
                min_aa=min_aa,
                top_k=top_k,
                max_peptide_aa=max_peptide_aa,
            )

            # -----------------------------------------------
            # Проверка CDS
            # -----------------------------------------------

            cds = validate_cds(
                sequence=sequence,
                sequence_kind=record["sequence_kind"],
            )

            cds_status_counts[
                cds["status"]
            ] += 1

            # -----------------------------------------------
            # Статистика
            # -----------------------------------------------

            complete_all = scan[
                "counts_all"
            ]["complete"]

            complete_passing = scan[
                "counts_passing_min_aa"
            ]["complete"]

            total_complete_orfs += complete_all

            total_complete_orfs_passing += (
                complete_passing
            )

            if complete_passing:
                records_with_complete_orf += 1

            if count_ambiguous_bases(sequence):
                records_with_ambiguous_bases += 1

            if count_invalid_bases(sequence):
                records_with_invalid_bases += 1

            # -----------------------------------------------
            # Полный результат JSONL
            # -----------------------------------------------

            result = {
                "gene_id": record["gene_id"],
                "species": record["species"],
                "section": record["section"],
                "locus": record["locus"],
                "sequence_kind": record["sequence_kind"],
                "sequence": sequence,
                "source_metadata": {
                    "original_header": record.get(
                        "original_header"
                    ),
                    "header_line": record.get(
                        "header_line"
                    ),
                    "score_label": record.get(
                        "score_label"
                    ),
                    "shift_label": record.get(
                        "shift_label"
                    ),
                    "exon_count": record.get(
                        "exon_count"
                    ),
                    "strand": record.get(
                        "strand"
                    ),
                    "start": record.get(
                        "start"
                    ),
                    "end": record.get(
                        "end"
                    ),
                    "issues": record.get(
                        "issues", []
                    ),
                },
                "orf_scan": scan,
                "cds_validation": cds,
            }

            records_file.write(
                json.dumps(
                    result,
                    ensure_ascii=False,
                ) + "\n"
            )

            # -----------------------------------------------
            # CSV
            # -----------------------------------------------

            writer.writerow(
                build_summary_row(
                    record,
                    scan,
                    cds,
                )
            )

            # -----------------------------------------------
            # Прогресс
            # -----------------------------------------------

            if (
                progress_every > 0
                and processed % progress_every == 0
            ):

                log(
                    f"Обработано записей: {processed}; "
                    f"генов: {len(genes)}"
                )

    # ========================================================
    # 10. FINAL REPORT
    # ========================================================

    if processed == 0:

        raise ValueError(
            "Не найдено записей для обработки. "
            "Проверьте входной файл и --gene."
        )

    report = {
        "stage": "ORF Analysis Stage 2",
        "input_file": str(input_path.resolve()),
        "gene_filter": gene_filter,
        "settings": {
            "min_aa": min_aa,
            "top_k": top_k,
            "max_peptide_aa": max_peptide_aa,
            "start_codon": START_CODON,
            "stop_codons": sorted(STOP_CODONS),
            "genetic_code": 1,
        },
        "scanned_jsonl_records": scanned_lines,
        "processed_records": processed,
        "unique_genes": len(genes),
        "records_per_section": dict(section_counts),
        "cds_status_counts": dict(cds_status_counts),
        "records_with_complete_orf_min_aa": (
            records_with_complete_orf
        ),
        "total_complete_orfs_all_lengths": (
            total_complete_orfs
        ),
        "total_complete_orfs_min_aa": (
            total_complete_orfs_passing
        ),
        "records_with_ambiguous_bases": (
            records_with_ambiguous_bases
        ),
        "records_with_invalid_bases": (
            records_with_invalid_bases
        ),
        "coordinate_system": "0-based_half-open",
        "strand_definition": (
            "relative_to_input_sequence"
        ),
    }

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    log("=" * 55)

    log("Анализ завершён")

    log(f"Обработано последовательностей: {processed}")
    log(f"Уникальных генов: {len(genes)}")

    log(
        "Записей с полной ORF не короче порога: "
        f"{records_with_complete_orf}"
    )

    log(
        f"Всего полных ORF не короче порога: "
        f"{total_complete_orfs_passing}"
    )

    log("CDS validation:")

    for status, count in sorted(
        cds_status_counts.items()
    ):

        log(f"  {status}: {count}")

    log("Выходные файлы:")
    log(str(records_path.resolve()))
    log(str(summary_path.resolve()))
    log(str(report_path.resolve()))

    log("=" * 55)

    return report


# ============================================================
# 11. COMMAND LINE INTERFACE
# ============================================================

def main(argv: list[str] | None = None) -> int:

    parser = argparse.ArgumentParser(
        description=(
            "ORF Analysis Stage 2: "
            "поиск ORF и проверка CDS"
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="JSONL первого этапа",
    )

    parser.add_argument(
        "--records-out",
        type=Path,
        required=True,
        help="Выходной JSONL",
    )

    parser.add_argument(
        "--summary-out",
        type=Path,
        required=True,
        help="Выходной CSV",
    )

    parser.add_argument(
        "--report-out",
        type=Path,
        required=True,
        help="Выходной JSON-отчёт",
    )

    parser.add_argument(
        "--gene",
        default=None,
        help="ID одного гена (необязательно)",
    )

    parser.add_argument(
        "--min-aa",
        type=int,
        default=25,
        help="Минимальная длина ORF в аминокислотах",
    )

    parser.add_argument(
        "--top-k",
        type=int,
        default=10,
        help="Сколько наиболее длинных ORF сохранять",
    )

    parser.add_argument(
        "--max-peptide-aa",
        type=int,
        default=120,
        help="Максимальная длина сохраняемого превью белка",
    )

    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help=(
            "Периодичность сообщений о прогрессе; "
            "0 отключает периодические сообщения"
        ),
    )

    args = parser.parse_args(argv)

    try:

        run(
            input_path=args.input,
            records_path=args.records_out,
            summary_path=args.summary_out,
            report_path=args.report_out,
            gene_filter=args.gene,
            min_aa=args.min_aa,
            top_k=args.top_k,
            max_peptide_aa=args.max_peptide_aa,
            progress_every=args.progress_every,
        )

    except (
        OSError,
        ValueError,
    ) as error:

        print(
            f"[ORF Scanner] ОШИБКА: {error}",
            file=sys.stderr,
            flush=True,
        )

        return 2

    return 0


# ============================================================
# 12. ENTRY POINT
# ============================================================

if __name__ == "__main__":
    raise SystemExit(main())
