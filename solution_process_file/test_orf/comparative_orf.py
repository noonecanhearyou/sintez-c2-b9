#!/usr/bin/env python3

"""
ORF Analysis — Stage 3
======================

Comparative CDS Analysis.

Вход:
    JSONL из orf_scanner.py (этап 2).

Выход:
    comparisons.jsonl
    comparison_summary.csv
    tool_agreement.csv
    comparison_report.json

Возможности:
    - Сравнение CDS с референсом D. melanogaster.
    - Независимый анализ Spaln и GeneWise.
    - Глобальное попарное выравнивание.
    - Подсчёт совпадений и замен.
    - Поиск разрывов выравнивания.
    - Выявление внутренних разрывов длиной не кратной 3.
    - Проверка качества выравнивания.
    - Сравнение предсказаний Spaln и GeneWise.
    - Двухпроходное чтение JSONL.
    - Ограничение вычислительной сложности.
    - CLI и вывод прогресса.

Важно:
    Обнаруженные различия не являются доказанными
    эволюционными мутациями.

Требования:
    Python >= 3.10
    Biopython
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys

from collections import Counter, defaultdict
from pathlib import Path

from Bio.Align import PairwiseAligner


# ============================================================
# 1. CONSTANTS
# ============================================================

REFERENCE_SPECIES = "dmel"

METHODS = ("spaln", "genewise")

DNA_BASES = frozenset("ACGT")

CDS_KINDS = frozenset({
    "annotated_cds",
    "predicted_cds",
})

SUMMARY_FIELDS = [
    "gene_id",
    "species",
    "method",
    "locus",
    "reference_length_nt",
    "target_length_nt",
    "reference_cds_status",
    "target_cds_status",
    "reported_shift",
    "alignment_status",
    "alignment_score",
    "identity",
    "reference_coverage",
    "query_coverage",
    "matched_nt",
    "mismatch_nt",
    "gap_in_query_nt",
    "gap_in_reference_nt",
    "gap_event_count",
    "internal_nontriplet_gap_events",
    "internal_triplet_gap_events",
    "terminal_gap_events",
    "quality_ok",
    "reason",
]

AGREEMENT_FIELDS = [
    "gene_id",
    "species",
    "locus",
    "agreement_type",
    "spaln_count",
    "genewise_count",
    "spaln_length_nt",
    "genewise_length_nt",
    "spaln_cds_status",
    "genewise_cds_status",
    "spaln_shift_label",
    "genewise_shift_label",
    "same_sequence",
    "same_cds_status",
    "same_shift_label",
]


# ============================================================
# 2. INPUT READER
# ============================================================

def iter_records(path: Path):
    """
    Последовательно читает JSONL второго этапа.

    Проверяет наличие основных полей.
    Не загружает весь файл в оперативную память.
    """

    with Path(path).open(
        "r",
        encoding="utf-8",
    ) as file:

        for line_number, line in enumerate(
            file,
            start=1,
        ):

            if not line.strip():
                continue

            try:
                record = json.loads(line)

            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Строка {line_number}: "
                    f"некорректный JSON: {error}"
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
                        f"отсутствует строковое поле {field}"
                    )

            if record["section"] not in (
                "raw",
                "spaln",
                "genewise",
            ):
                raise ValueError(
                    f"Строка {line_number}: "
                    "неизвестная секция"
                )

            record["sequence"] = (
                record["sequence"].upper()
            )

            yield record


# ============================================================
# 3. ALIGNER
# ============================================================

def create_aligner() -> PairwiseAligner:
    """
    Создаёт глобальный нуклеотидный выравниватель.

    Схема оценивания:

        match       +2
        mismatch    -1
        gap open    -5
        gap extend  -0.5

    Это рабочие параметры прототипа,
    а не специальная модель эволюции кодонов.
    """

    aligner = PairwiseAligner()

    aligner.mode = "global"

    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0

    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5

    return aligner


# ============================================================
# 4. ALIGNMENT ANALYSIS
# ============================================================

def analyze_alignment(
    reference: str,
    query: str,
    alignment,
    max_events: int = 25,
    preview_columns: int = 90,
) -> dict:
    """
    Извлекает измеримые характеристики выравнивания.

    Координаты:
        0-based, half-open.

    Gap in query:
        Нуклеотиды есть в референсе,
        но соответствующие позиции являются
        пропуском в query.

    Gap in reference:
        Нуклеотиды есть в query,
        но отсутствуют в референсном выравнивании.

    Направление эволюционного события
    по этим данным не устанавливается.
    """

    if max_events < 0:
        raise ValueError("max_events должен быть >= 0")

    if preview_columns < 0:
        raise ValueError(
            "preview_columns должен быть >= 0"
        )

    coordinates = alignment.coordinates

    matched_nt = 0
    mismatch_nt = 0
    aligned_pairs = 0

    gap_in_query_nt = 0
    gap_in_reference_nt = 0

    gap_event_count = 0

    internal_nontriplet = 0
    internal_triplet = 0
    terminal_gaps = 0

    mismatch_examples = []
    gap_examples = []

    # --------------------------------------------------------
    # Последовательно рассматриваем блоки выравнивания
    # --------------------------------------------------------

    for index in range(
        coordinates.shape[1] - 1
    ):

        ref_start = int(coordinates[0, index])
        ref_end = int(coordinates[0, index + 1])

        query_start = int(coordinates[1, index])
        query_end = int(coordinates[1, index + 1])

        ref_delta = ref_end - ref_start
        query_delta = query_end - query_start

        # ----------------------------------------------------
        # Блок сопоставленных нуклеотидов
        # ----------------------------------------------------

        if ref_delta > 0 and query_delta > 0:

            if ref_delta != query_delta:
                raise ValueError(
                    "Некорректный блок выравнивания"
                )

            for offset in range(ref_delta):

                ref_position = ref_start + offset
                query_position = query_start + offset

                ref_base = reference[ref_position]
                query_base = query[query_position]

                aligned_pairs += 1

                if ref_base == query_base:

                    matched_nt += 1

                else:

                    mismatch_nt += 1

                    if len(mismatch_examples) < max_events:

                        mismatch_examples.append({
                            "reference_position": ref_position,
                            "query_position": query_position,
                            "reference_base": ref_base,
                            "query_base": query_base,
                        })

        # ----------------------------------------------------
        # Разрыв в query
        # ----------------------------------------------------

        elif ref_delta > 0 and query_delta == 0:

            length_nt = ref_delta

            gap_in_query_nt += length_nt
            gap_event_count += 1

            # Краевые разрывы могут отражать
            # неполные CDS или разницу границ.

            is_internal = (
                ref_start > 0
                and ref_end < len(reference)
                and query_start > 0
                and query_start < len(query)
            )

            if is_internal:

                if length_nt % 3 == 0:
                    internal_triplet += 1
                else:
                    internal_nontriplet += 1

            else:
                terminal_gaps += 1

            if len(gap_examples) < max_events:

                gap_examples.append({
                    "kind": "gap_in_query",
                    "reference_start": ref_start,
                    "reference_end": ref_end,
                    "query_position": query_start,
                    "length_nt": length_nt,
                    "internal": is_internal,
                    "nontriplet": (
                        length_nt % 3 != 0
                    ),
                })

        # ----------------------------------------------------
        # Разрыв в reference
        # ----------------------------------------------------

        elif ref_delta == 0 and query_delta > 0:

            length_nt = query_delta

            gap_in_reference_nt += length_nt
            gap_event_count += 1

            is_internal = (
                query_start > 0
                and query_end < len(query)
                and ref_start > 0
                and ref_start < len(reference)
            )

            if is_internal:

                if length_nt % 3 == 0:
                    internal_triplet += 1
                else:
                    internal_nontriplet += 1

            else:
                terminal_gaps += 1

            if len(gap_examples) < max_events:

                gap_examples.append({
                    "kind": "gap_in_reference",
                    "reference_position": ref_start,
                    "query_start": query_start,
                    "query_end": query_end,
                    "length_nt": length_nt,
                    "internal": is_internal,
                    "nontriplet": (
                        length_nt % 3 != 0
                    ),
                })

        elif ref_delta == 0 and query_delta == 0:
            continue

        else:
            raise ValueError(
                "Неожиданная структура выравнивания"
            )

    # ========================================================
    # 5. METRICS
    # ========================================================

    identity = (
        matched_nt / aligned_pairs
        if aligned_pairs
        else 0.0
    )

    reference_coverage = (
        aligned_pairs / len(reference)
        if reference
        else 0.0
    )

    query_coverage = (
        aligned_pairs / len(query)
        if query
        else 0.0
    )

    # --------------------------------------------------------
    # Короткий фрагмент выравнивания
    # --------------------------------------------------------

    preview = None

    if preview_columns > 0:

        aligned_reference = alignment[0]
        aligned_query = alignment[1]

        preview = {
            "reference": aligned_reference[
                :preview_columns
            ],
            "query": aligned_query[
                :preview_columns
            ],
            "truncated": (
                alignment.length > preview_columns
            ),
        }

    return {
        "alignment_score": float(alignment.score),
        "alignment_columns": int(alignment.length),
        "matched_nt": matched_nt,
        "mismatch_nt": mismatch_nt,
        "aligned_pairs": aligned_pairs,
        "identity": round(identity, 6),
        "reference_coverage": round(
            reference_coverage, 6
        ),
        "query_coverage": round(
            query_coverage, 6
        ),
        "gap_in_query_nt": gap_in_query_nt,
        "gap_in_reference_nt": gap_in_reference_nt,
        "gap_event_count": gap_event_count,
        "internal_nontriplet_gap_events": (
            internal_nontriplet
        ),
        "internal_triplet_gap_events": (
            internal_triplet
        ),
        "terminal_gap_events": terminal_gaps,
        "mismatch_examples": mismatch_examples,
        "gap_examples": gap_examples,
        "mismatch_examples_truncated": (
            mismatch_nt > len(mismatch_examples)
        ),
        "gap_examples_truncated": (
            gap_event_count > len(gap_examples)
        ),
        "alignment_preview": preview,
    }


# ============================================================
# 6. SAFE CDS COMPARISON
# ============================================================

def compare_cds(
    reference: str,
    query: str,
    aligner: PairwiseAligner,
    max_align_length: int = 4000,
    max_align_cells: int = 4_000_000,
    min_identity: float = 0.60,
    min_coverage: float = 0.60,
    max_events: int = 25,
    preview_columns: int = 90,
) -> dict:
    """
    Безопасное сравнение двух CDS.

    Проверки выполняются до запуска выравнивания.

    Возвращает status и, при успешном выравнивании,
    количественные характеристики.

    min_coverage применяется независимо к референсу
    и сравниваемой последовательности.
    """

    reference = reference.upper()
    query = query.upper()

    if not reference:
        return {
            "status": "skipped_empty_reference",
            "reason": "Пустая референсная CDS",
        }

    if not query:
        return {
            "status": "skipped_empty_query",
            "reason": "Пустая сравниваемая CDS",
        }

    # --------------------------------------------------------
    # Проверка нуклеотидного алфавита
    # --------------------------------------------------------

    if (
        any(base not in DNA_BASES for base in reference)
        or any(base not in DNA_BASES for base in query)
    ):

        return {
            "status": "skipped_non_acgt",
            "reason": (
                "Последовательность содержит неоднозначные "
                "или недопустимые нуклеотиды"
            ),
        }

    # --------------------------------------------------------
    # Ограничение вычислительной сложности
    # --------------------------------------------------------

    if (
        len(reference) > max_align_length
        or len(query) > max_align_length
    ):

        return {
            "status": "skipped_size_limit",
            "reason": (
                "Превышена максимальная длина "
                "последовательности"
            ),
        }

    if (
        len(reference) * len(query)
        > max_align_cells
    ):

        return {
            "status": "skipped_size_limit",
            "reason": (
                "Превышен лимит произведения "
                "длин последовательностей"
            ),
        }

    # --------------------------------------------------------
    # Глобальное выравнивание
    # --------------------------------------------------------

    try:

        alignments = aligner.align(
            reference,
            query,
        )

        # Используем одно оптимальное выравнивание.
        # Возможны другие выравнивания с тем же score.

        alignment = alignments[0]

        metrics = analyze_alignment(
            reference=reference,
            query=query,
            alignment=alignment,
            max_events=max_events,
            preview_columns=preview_columns,
        )

    except (ValueError, RuntimeError) as error:

        return {
            "status": "alignment_failed",
            "reason": str(error),
        }

    # --------------------------------------------------------
    # Проверка качества
    # --------------------------------------------------------

    identity = (
        metrics["matched_nt"] / metrics["aligned_pairs"]
        if metrics["aligned_pairs"]
        else 0.0
    )

    reference_coverage = (
        metrics["aligned_pairs"] / len(reference)
    )

    query_coverage = (
        metrics["aligned_pairs"] / len(query)
    )

    quality_ok = (
        identity >= min_identity
        and reference_coverage >= min_coverage
        and query_coverage >= min_coverage
    )

    metrics["quality_ok"] = quality_ok

    if quality_ok:

        metrics["status"] = "aligned"
        metrics["reason"] = None

    else:

        metrics["status"] = "low_quality_alignment"
        metrics["reason"] = (
            "Identity или покрытие ниже заданного порога"
        )

    # Эти события — кандидаты на изменения
    # кодирующей структуры, но не доказанные
    # исторические мутации.

    metrics["interpretation"] = {
        "internal_nontriplet_gaps_are": (
            "potential_frame_disruption_candidates"
        ),
        "evolutionary_direction_inferred": False,
        "codon_aware_alignment": False,
        "quality_passed": quality_ok,
    }

    return metrics


# ============================================================
# 7. TOOL AGREEMENT
# ============================================================

def prediction_signature(record: dict) -> dict:
    """
    Краткие сведения о предсказании Spaln/GeneWise.

    Полную последовательность не сохраняем,
    чтобы экономить память.
    """

    sequence = record["sequence"]

    cds = record.get("cds_validation") or {}
    metadata = record.get("source_metadata") or {}

    return {
        "length_nt": len(sequence),
        "cds_status": cds.get("status"),
        "shift_label": metadata.get("shift_label"),
        "sha256": hashlib.sha256(
            sequence.encode("ascii")
        ).hexdigest(),
    }


def build_agreement_row(
    gene_id: str,
    species: str,
    locus: str,
    profiles: dict,
) -> dict:
    """
    Сравнивает два предсказания только при
    однозначном соответствии:

        gene_id + species + locus

    Если одного предсказания нет или существует
    несколько записей одного метода, это отмечается.
    """

    spaln = profiles.get("spaln", [])
    genewise = profiles.get("genewise", [])

    row = {
        "gene_id": gene_id,
        "species": species,
        "locus": locus,
        "agreement_type": None,
        "spaln_count": len(spaln),
        "genewise_count": len(genewise),
        "spaln_length_nt": None,
        "genewise_length_nt": None,
        "spaln_cds_status": None,
        "genewise_cds_status": None,
        "spaln_shift_label": None,
        "genewise_shift_label": None,
        "same_sequence": None,
        "same_cds_status": None,
        "same_shift_label": None,
    }

    if len(spaln) > 1 or len(genewise) > 1:

        row["agreement_type"] = (
            "ambiguous_multiple_hits"
        )

        return row

    if not spaln:

        row["agreement_type"] = "missing_spaln"

        return row

    if not genewise:

        row["agreement_type"] = "missing_genewise"

        return row

    a = spaln[0]
    b = genewise[0]

    same_sequence = (
        a["sha256"] == b["sha256"]
        and a["length_nt"] == b["length_nt"]
    )

    same_status = (
        a["cds_status"] == b["cds_status"]
    )

    same_shift = (
        a["shift_label"] == b["shift_label"]
    )

    if same_sequence:

        agreement_type = "exact_sequence"

    elif same_status:

        agreement_type = "same_cds_status_only"

    else:

        agreement_type = "different_cds_status"

    row.update({
        "agreement_type": agreement_type,
        "spaln_length_nt": a["length_nt"],
        "genewise_length_nt": b["length_nt"],
        "spaln_cds_status": a["cds_status"],
        "genewise_cds_status": b["cds_status"],
        "spaln_shift_label": a["shift_label"],
        "genewise_shift_label": b["shift_label"],
        "same_sequence": same_sequence,
        "same_cds_status": same_status,
        "same_shift_label": same_shift,
    })

    return row


# ============================================================
# 8. MAIN PROCESSING
# ============================================================

def run(
    input_path: Path,
    comparisons_path: Path,
    summary_path: Path,
    agreement_path: Path,
    report_path: Path,
    gene_filter: str | None = None,
    max_align_length: int = 4000,
    max_align_cells: int = 4_000_000,
    min_identity: float = 0.60,
    min_coverage: float = 0.60,
    max_events: int = 25,
    preview_columns: int = 90,
    progress_every: int = 200,
    verbose: bool = True,
) -> dict:

    input_path = Path(input_path)

    comparisons_path = Path(comparisons_path)
    summary_path = Path(summary_path)
    agreement_path = Path(agreement_path)
    report_path = Path(report_path)

    # --------------------------------------------------------
    # Проверка параметров
    # --------------------------------------------------------

    if max_align_length < 1:
        raise ValueError(
            "max_align_length должен быть >= 1"
        )

    if max_align_cells < 1:
        raise ValueError(
            "max_align_cells должен быть >= 1"
        )

    if not 0 <= min_identity <= 1:
        raise ValueError(
            "min_identity должен быть от 0 до 1"
        )

    if not 0 <= min_coverage <= 1:
        raise ValueError(
            "min_coverage должен быть от 0 до 1"
        )

    if max_events < 0:
        raise ValueError(
            "max_events должен быть >= 0"
        )

    if preview_columns < 0:
        raise ValueError(
            "preview_columns должен быть >= 0"
        )

    if progress_every < 0:
        raise ValueError(
            "progress_every должен быть >= 0"
        )

    paths = [
        input_path.resolve(),
        comparisons_path.resolve(),
        summary_path.resolve(),
        agreement_path.resolve(),
        report_path.resolve(),
    ]

    if len(set(paths)) != 5:

        raise ValueError(
            "Входной и выходные пути должны различаться"
        )

    if not input_path.is_file():

        raise FileNotFoundError(
            f"Входной файл не найден: {input_path.resolve()}"
        )

    # --------------------------------------------------------
    # Logging
    # --------------------------------------------------------

    def log(message: str):

        if verbose:
            print(
                f"[Comparative ORF] {message}",
                flush=True,
            )

    for path in (
        comparisons_path,
        summary_path,
        agreement_path,
        report_path,
    ):

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    log("=" * 55)

    log("Stage 3: Comparative CDS Analysis")
    log(f"Вход: {input_path.resolve()}")

    if gene_filter:
        log(f"Ген: {gene_filter}")
    else:
        log("Режим: все гены")

    log(
        f"Лимит выравнивания: {max_align_length} nt"
    )

    log(
        f"Лимит сложности: {max_align_cells} клеток"
    )

    # ========================================================
    # 9. FIRST PASS: REFERENCES
    # ========================================================

    log("Проход 1/2: поиск референсных CDS dmel")

    references = defaultdict(list)

    reference_records = 0

    for record in iter_records(input_path):

        if record["section"] not in METHODS:
            continue

        if (
            gene_filter is not None
            and record["gene_id"] != gene_filter
        ):
            continue

        if record["species"] != REFERENCE_SPECIES:
            continue

        key = (
            record["gene_id"],
            record["section"],
        )

        # Для определения неоднозначности
        # достаточно хранить максимум две записи.

        if len(references[key]) < 2:

            references[key].append({
                "sequence": record["sequence"],
                "sequence_kind": record["sequence_kind"],
                "cds_status": (
                    record.get("cds_validation") or {}
                ).get("status"),
                "locus": record["locus"],
            })

        reference_records += 1

    log(
        f"Найдено референсных записей: "
        f"{reference_records}"
    )

    log(
        f"Уникальных сочетаний ген + метод: "
        f"{len(references)}"
    )

    # ========================================================
    # 10. SECOND PASS: COMPARISONS
    # ========================================================

    log("Проход 2/2: сравнительный анализ CDS")

    aligner = create_aligner()

    processed = 0

    genes = set()

    status_counts = Counter()
    method_counts = Counter()

    accepted_nontriplet_gaps = 0

    # Храним только краткие сведения о предсказаниях,
    # без полных последовательностей.

    tool_profiles = defaultdict(
        lambda: {
            "spaln": [],
            "genewise": [],
        }
    )

    with (
        comparisons_path.open(
            "w",
            encoding="utf-8",
        ) as comparisons_file,

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

        for record in iter_records(input_path):

            method = record["section"]

            if method not in METHODS:
                continue

            if (
                gene_filter is not None
                and record["gene_id"] != gene_filter
            ):
                continue

            # Референсные записи не сравниваем
            # сами с собой.

            if record["species"] == REFERENCE_SPECIES:
                continue

            gene_id = record["gene_id"]
            species = record["species"]
            locus = record["locus"]
            sequence = record["sequence"]

            processed += 1

            genes.add(gene_id)
            method_counts[method] += 1

            # -----------------------------------------------
            # Профиль Spaln/GeneWise
            # -----------------------------------------------

            profile_key = (
                gene_id,
                species,
                locus,
            )

            tool_profiles[
                profile_key
            ][method].append(
                prediction_signature(record)
            )

            # -----------------------------------------------
            # Референс
            # -----------------------------------------------

            ref_key = (
                gene_id,
                method,
            )

            ref_list = references.get(
                ref_key, []
            )

            reference = None

            if not ref_list:

                comparison = {
                    "status": "skipped_missing_reference",
                    "reason": (
                        "Референсная CDS dmel отсутствует"
                    ),
                }

            elif len(ref_list) > 1:

                comparison = {
                    "status": "skipped_multiple_references",
                    "reason": (
                        "Найдено несколько референсных "
                        "последовательностей"
                    ),
                }

            else:

                reference = ref_list[0]

                if (
                    reference["sequence_kind"]
                    != "annotated_cds"
                    or reference["cds_status"]
                    != "complete_candidate"
                ):

                    comparison = {
                        "status": (
                            "skipped_reference_noncanonical"
                        ),
                        "reason": (
                            "Референсная CDS не прошла "
                            "проверку структуры"
                        ),
                    }

                elif record["sequence_kind"] not in CDS_KINDS:

                    comparison = {
                        "status": "skipped_non_cds_type",
                        "reason": (
                            "Последовательность не является "
                            "подготовленной CDS"
                        ),
                    }

                else:

                    comparison = compare_cds(
                        reference=reference["sequence"],
                        query=sequence,
                        aligner=aligner,
                        max_align_length=max_align_length,
                        max_align_cells=max_align_cells,
                        min_identity=min_identity,
                        min_coverage=min_coverage,
                        max_events=max_events,
                        preview_columns=preview_columns,
                    )

            # -----------------------------------------------
            # Статистика
            # -----------------------------------------------

            status = comparison["status"]

            status_counts[status] += 1

            if status == "aligned":

                accepted_nontriplet_gaps += (
                    comparison.get(
                        "internal_nontriplet_gap_events",
                        0,
                    )
                )

            metadata = (
                record.get("source_metadata") or {}
            )

            target_cds_status = (
                record.get("cds_validation") or {}
            ).get("status")

            ref_sequence = (
                reference["sequence"]
                if reference is not None
                else None
            )

            ref_cds_status = (
                reference["cds_status"]
                if reference is not None
                else None
            )

            # -----------------------------------------------
            # JSONL
            # -----------------------------------------------

            result = {
                "gene_id": gene_id,
                "species": species,
                "method": method,
                "locus": locus,
                "reference": {
                    "species": REFERENCE_SPECIES,
                    "length_nt": (
                        len(ref_sequence)
                        if ref_sequence is not None
                        else None
                    ),
                    "cds_status": ref_cds_status,
                },
                "target": {
                    "length_nt": len(sequence),
                    "sequence_kind": record[
                        "sequence_kind"
                    ],
                    "cds_status": target_cds_status,
                    "reported_shift": metadata.get(
                        "shift_label"
                    ),
                    "reported_score": metadata.get(
                        "score_label"
                    ),
                },
                "comparison": comparison,
            }

            comparisons_file.write(
                json.dumps(
                    result,
                    ensure_ascii=False,
                ) + "\n"
            )

            # -----------------------------------------------
            # CSV
            # -----------------------------------------------

            row = {
                "gene_id": gene_id,
                "species": species,
                "method": method,
                "locus": locus,
                "reference_length_nt": (
                    len(ref_sequence)
                    if ref_sequence is not None
                    else None
                ),
                "target_length_nt": len(sequence),
                "reference_cds_status": ref_cds_status,
                "target_cds_status": target_cds_status,
                "reported_shift": metadata.get(
                    "shift_label"
                ),
                "alignment_status": status,
                "alignment_score": comparison.get(
                    "alignment_score"
                ),
                "identity": comparison.get(
                    "identity"
                ),
                "reference_coverage": comparison.get(
                    "reference_coverage"
                ),
                "query_coverage": comparison.get(
                    "query_coverage"
                ),
                "matched_nt": comparison.get(
                    "matched_nt"
                ),
                "mismatch_nt": comparison.get(
                    "mismatch_nt"
                ),
                "gap_in_query_nt": comparison.get(
                    "gap_in_query_nt"
                ),
                "gap_in_reference_nt": comparison.get(
                    "gap_in_reference_nt"
                ),
                "gap_event_count": comparison.get(
                    "gap_event_count"
                ),
                "internal_nontriplet_gap_events": (
                    comparison.get(
                        "internal_nontriplet_gap_events"
                    )
                ),
                "internal_triplet_gap_events": (
                    comparison.get(
                        "internal_triplet_gap_events"
                    )
                ),
                "terminal_gap_events": comparison.get(
                    "terminal_gap_events"
                ),
                "quality_ok": comparison.get(
                    "quality_ok", False
                ),
                "reason": comparison.get(
                    "reason"
                ),
            }

            writer.writerow(row)

            # -----------------------------------------------
            # Progress
            # -----------------------------------------------

            if (
                progress_every > 0
                and processed % progress_every == 0
            ):

                log(
                    f"Обработано: {processed}; "
                    f"генов: {len(genes)}; "
                    f"успешных выравниваний: "
                    f"{status_counts['aligned']}"
                )

    if processed == 0:
        raise ValueError(
            "Нет последовательностей для сравнения. "
            "Проверьте входной файл и --gene."
        )

    # ========================================================
    # 11. SPALN / GENEWISE AGREEMENT
    # ========================================================

    log("Сравнение предсказаний Spaln и GeneWise")

    agreement_counts = Counter()

    with agreement_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as agreement_file:

        writer = csv.DictWriter(
            agreement_file,
            fieldnames=AGREEMENT_FIELDS,
        )

        writer.writeheader()

        for key in sorted(tool_profiles):

            gene_id, species, locus = key

            agreement_row = build_agreement_row(
                gene_id,
                species,
                locus,
                tool_profiles[key],
            )

            writer.writerow(agreement_row)

            agreement_counts[
                agreement_row["agreement_type"]
            ] += 1

    # ========================================================
    # 12. FINAL REPORT
    # ========================================================

    report = {
        "stage": "Comparative ORF Analysis Stage 3",
        "input_file": str(input_path.resolve()),
        "gene_filter": gene_filter,
        "reference_species": REFERENCE_SPECIES,
        "methods": list(METHODS),
        "settings": {
            "max_align_length": max_align_length,
            "max_align_cells": max_align_cells,
            "min_identity": min_identity,
            "min_coverage": min_coverage,
            "max_events": max_events,
            "preview_columns": preview_columns,
            "alignment_mode": "global",
            "match_score": 2.0,
            "mismatch_score": -1.0,
            "open_gap_score": -5.0,
            "extend_gap_score": -0.5,
        },
        "reference_records": reference_records,
        "reference_gene_method_pairs": len(references),
        "processed_comparisons": processed,
        "unique_genes": len(genes),
        "comparisons_per_method": dict(method_counts),
        "alignment_status_counts": dict(status_counts),
        "accepted_internal_nontriplet_gap_events": (
            accepted_nontriplet_gaps
        ),
        "tool_agreement_counts": dict(agreement_counts),
        "coordinate_system": "0-based_half-open",
        "limitations": [
            "Pairwise nucleotide alignment is not codon-aware.",
            "Alignment gaps are not verified evolutionary indels.",
            "Evolutionary direction is not inferred.",
            "Predicted CDS differences may reflect gene prediction errors.",
            "Sequence homology and orthology are not established here.",
            "Frameshift-like gaps are preliminary candidates.",
        ],
    }

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    # ========================================================
    # 13. FINAL LOG
    # ========================================================

    log("=" * 55)

    log("Comparative CDS Analysis завершён")

    log(f"Сравнений: {processed}")
    log(f"Уникальных генов: {len(genes)}")

    log("Статусы выравниваний:")

    for status, count in sorted(
        status_counts.items()
    ):

        log(f"  {status}: {count}")

    log(
        "Кандидатных внутренних некратных 3 разрывов "
        "в выравниваниях, прошедших пороги качества: "
        f"{accepted_nontriplet_gaps}"
    )

    log("Согласованность Spaln / GeneWise:")

    for status, count in sorted(
        agreement_counts.items()
    ):

        log(f"  {status}: {count}")

    log("Результаты:")

    log(str(comparisons_path.resolve()))
    log(str(summary_path.resolve()))
    log(str(agreement_path.resolve()))
    log(str(report_path.resolve()))

    log("=" * 55)

    return report


# ============================================================
# 14. COMMAND LINE INTERFACE
# ============================================================

def main(argv: list[str] | None = None) -> int:

    parser = argparse.ArgumentParser(
        description=(
            "ORF Analysis Stage 3: "
            "межвидовое сравнение CDS"
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="JSONL второго этапа",
    )

    parser.add_argument(
        "--comparisons-out",
        required=True,
        type=Path,
        help="Выходной JSONL сравнений",
    )

    parser.add_argument(
        "--summary-out",
        required=True,
        type=Path,
        help="Выходной CSV сравнений",
    )

    parser.add_argument(
        "--agreement-out",
        required=True,
        type=Path,
        help="CSV согласованности Spaln/GeneWise",
    )

    parser.add_argument(
        "--report-out",
        required=True,
        type=Path,
        help="Итоговый JSON-отчёт",
    )

    parser.add_argument(
        "--gene",
        default=None,
        help="ID одного гена; иначе весь набор",
    )

    parser.add_argument(
        "--max-align-length",
        type=int,
        default=4000,
        help="Максимальная длина каждой CDS",
    )

    parser.add_argument(
        "--max-align-cells",
        type=int,
        default=4_000_000,
        help="Максимальное произведение длин CDS",
    )

    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.60,
        help="Минимальная доля совпадений (0..1)",
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.60,
        help="Минимальное покрытие обеих CDS (0..1)",
    )

    parser.add_argument(
        "--max-events",
        type=int,
        default=25,
        help="Сколько примеров замен и разрывов сохранять",
    )

    parser.add_argument(
        "--preview-columns",
        type=int,
        default=90,
        help="Длина превью выравнивания",
    )

    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
        help="Сообщение о прогрессе каждые N сравнений",
    )

    args = parser.parse_args(argv)

    try:

        run(
            input_path=args.input,
            comparisons_path=args.comparisons_out,
            summary_path=args.summary_out,
            agreement_path=args.agreement_out,
            report_path=args.report_out,
            gene_filter=args.gene,
            max_align_length=args.max_align_length,
            max_align_cells=args.max_align_cells,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            max_events=args.max_events,
            preview_columns=args.preview_columns,
            progress_every=args.progress_every,
        )

    except (
        OSError,
        ValueError,
    ) as error:

        print(
            f"[Comparative ORF] ОШИБКА: {error}",
            file=sys.stderr,
            flush=True,
        )

        return 2

    return 0


# ============================================================
# 15. ENTRY POINT
# ============================================================

if __name__ == "__main__":
    raise SystemExit(main())
