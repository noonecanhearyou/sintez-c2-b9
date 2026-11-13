"""
ORF Analysis Stage 3 — Tests.

Запуск:
    python -m pytest -v test_comparative_orf.py

Тестовые данные создаются во временных каталогах.
Реальные исходные файлы не изменяются.
"""

import csv
import json
import subprocess
import sys

from pathlib import Path

import comparative_orf

from comparative_orf import (
    create_aligner,
    compare_cds,
    prediction_signature,
    build_agreement_row,
)


# ============================================================
# 1. TEST HELPERS
# ============================================================

def compare(
    reference: str,
    query: str,
    **kwargs,
):

    aligner = create_aligner()

    return compare_cds(
        reference=reference,
        query=query,
        aligner=aligner,
        **kwargs,
    )


def make_record(
    gene_id: str,
    species: str,
    section: str,
    sequence: str,
    locus: str,
    cds_status: str = "complete_candidate",
    shift: str = "No",
) -> dict:

    return {
        "gene_id": gene_id,
        "species": species,
        "section": section,
        "locus": locus,
        "sequence_kind": (
            "annotated_cds"
            if species == "dmel"
            else "predicted_cds"
        ),
        "sequence": sequence,
        "cds_validation": {
            "status": cds_status,
        },
        "source_metadata": {
            "shift_label": shift,
            "score_label": "10.0",
        },
    }


# ============================================================
# 2. IDENTICAL CDS
# ============================================================

def test_identical_cds():

    sequence = "ATGAAACCCGGGTAA"

    result = compare(
        sequence,
        sequence,
    )

    assert result["status"] == "aligned"

    assert result["identity"] == 1.0

    assert result["reference_coverage"] == 1.0

    assert result["query_coverage"] == 1.0

    assert result["matched_nt"] == 15

    assert result["mismatch_nt"] == 0

    assert result["gap_event_count"] == 0

    assert result["quality_ok"] is True


# ============================================================
# 3. SINGLE NUCLEOTIDE SUBSTITUTION
# ============================================================

def test_single_substitution():

    reference = "ATGAAACCCGGGTAA"

    query = "ATGAAGCCCGGGTAA"

    result = compare(
        reference,
        query,
    )

    assert result["status"] == "aligned"

    assert result["matched_nt"] == 14

    assert result["mismatch_nt"] == 1

    assert result["gap_event_count"] == 0

    assert len(result["mismatch_examples"]) == 1

    event = result["mismatch_examples"][0]

    assert event["reference_position"] == 5

    assert event["reference_base"] == "A"

    assert event["query_base"] == "G"


# ============================================================
# 4. IN-FRAME GAP
# ============================================================

def test_inframe_gap():

    reference = "ATGAAACCCGGGTAA"

    query = "ATGAAAGGGTAA"

    result = compare(
        reference,
        query,
    )

    assert result["status"] == "aligned"

    assert result["gap_in_query_nt"] == 3

    assert (
        result["internal_triplet_gap_events"]
        == 1
    )

    assert (
        result["internal_nontriplet_gap_events"]
        == 0
    )


# ============================================================
# 5. NONTRIPLET INTERNAL GAP
# ============================================================

def test_nontriplet_gap():

    reference = "ATGAAACCCGGGTAA"

    # Отсутствует один C внутри последовательности.
    query = "ATGAAACCGGGTAA"

    result = compare(
        reference,
        query,
    )

    assert result["status"] == "aligned"

    assert result["gap_in_query_nt"] == 1

    assert (
        result["internal_nontriplet_gap_events"]
        == 1
    )

    event = result["gap_examples"][0]

    assert event["length_nt"] == 1

    assert event["internal"] is True

    assert event["nontriplet"] is True


# ============================================================
# 6. TERMINAL GAP
# ============================================================

def test_terminal_gap():

    reference = "ATGAAACCCGGGTAA"

    query = "ATGAAACCCGGGT"

    result = compare(
        reference,
        query,
    )

    assert result["status"] == "aligned"

    assert result["gap_event_count"] >= 1

    assert result["terminal_gap_events"] >= 1

    assert (
        result["internal_nontriplet_gap_events"]
        == 0
    )


# ============================================================
# 7. AMBIGUOUS NUCLEOTIDE
# ============================================================

def test_ambiguous_sequence_skipped():

    result = compare(
        "ATGAAACCCGGGTAA",
        "ATGAAANCCGGGTAA",
    )

    assert result["status"] == "skipped_non_acgt"


# ============================================================
# 8. COMPUTATIONAL LIMIT
# ============================================================

def test_alignment_size_limit():

    result = compare(
        "A" * 100,
        "A" * 100,
        max_align_cells=5000,
    )

    assert result["status"] == "skipped_size_limit"


# ============================================================
# 9. LOW QUALITY ALIGNMENT
# ============================================================

def test_low_quality_alignment():

    reference = "ATGAAACCCGGGTAA"

    query = "ATGAAGCCCGGGTAA"

    result = compare(
        reference,
        query,
        min_identity=1.0,
    )

    assert (
        result["status"]
        == "low_quality_alignment"
    )

    assert result["quality_ok"] is False


# ============================================================
# 10. TOOL AGREEMENT
# ============================================================

def test_tool_agreement():

    spaln_record = make_record(
        gene_id="FBgn0000001",
        species="dsim",
        section="spaln",
        sequence="ATGAAACCCGGGTAA",
        locus="chr2L+:100:115",
    )

    genewise_record = make_record(
        gene_id="FBgn0000001",
        species="dsim",
        section="genewise",
        sequence="ATGAAGCCCGGGTAA",
        locus="chr2L+:100:115",
    )

    profiles = {
        "spaln": [
            prediction_signature(spaln_record)
        ],
        "genewise": [
            prediction_signature(genewise_record)
        ],
    }

    row = build_agreement_row(
        gene_id="FBgn0000001",
        species="dsim",
        locus="chr2L+:100:115",
        profiles=profiles,
    )

    assert (
        row["agreement_type"]
        == "same_cds_status_only"
    )

    assert row["same_sequence"] is False

    assert row["same_cds_status"] is True

    assert row["same_shift_label"] is True


# ============================================================
# 11. FULL CLI TEST
# ============================================================

def test_cli_pipeline(tmp_path):

    input_dir = (
        tmp_path
        / "3_results_second"
    )

    input_dir.mkdir()

    input_path = (
        input_dir
        / "orf_records.jsonl"
    )

    gene = "FBgn0000001"

    reference = "ATGAAACCCGGGTAA"

    records = [
        # Сначала target — проверяем двухпроходное чтение.
        make_record(
            gene,
            "dsim",
            "spaln",
            "ATGAAGCCCGGGTAA",
            "chr2L+:100:115",
        ),

        make_record(
            gene,
            "dmel",
            "spaln",
            reference,
            gene,
        ),

        make_record(
            gene,
            "dsim",
            "genewise",
            reference,
            "chr2L+:100:115",
        ),

        make_record(
            gene,
            "dmel",
            "genewise",
            reference,
            gene,
        ),

        # Другая секция должна игнорироваться.
        {
            "gene_id": gene,
            "species": "dsim",
            "section": "raw",
            "locus": "chr2L+:100:115",
            "sequence_kind": "genomic_region",
            "sequence": reference,
        },
    ]

    with input_path.open(
        "w",
        encoding="utf-8",
    ) as file:

        for record in records:

            file.write(
                json.dumps(record) + "\n"
            )

    command = [
        sys.executable,
        str(
            Path(comparative_orf.__file__).resolve()
        ),

        "--input",
        "3_results_second/orf_records.jsonl",

        "--comparisons-out",
        "4_results_third/comparisons.jsonl",

        "--summary-out",
        "4_results_third/comparison_summary.csv",

        "--agreement-out",
        "4_results_third/tool_agreement.csv",

        "--report-out",
        "4_results_third/comparison_report.json",

        "--gene",
        gene,

        "--progress-every",
        "1",
    ]

    completed = subprocess.run(
        command,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr

    assert (
        "Comparative CDS Analysis завершён"
        in completed.stdout
    )

    output_dir = (
        tmp_path
        / "4_results_third"
    )

    expected_files = [
        "comparisons.jsonl",
        "comparison_summary.csv",
        "tool_agreement.csv",
        "comparison_report.json",
    ]

    for filename in expected_files:

        assert (
            output_dir / filename
        ).is_file()

    # --------------------------------------------------------
    # JSONL
    # --------------------------------------------------------

    with (
        output_dir / "comparisons.jsonl"
    ).open(
        encoding="utf-8"
    ) as file:

        comparisons = [
            json.loads(line)
            for line in file
            if line.strip()
        ]

    assert len(comparisons) == 2

    assert all(
        record["comparison"]["status"] == "aligned"
        for record in comparisons
    )

    # --------------------------------------------------------
    # Summary CSV
    # --------------------------------------------------------

    with (
        output_dir / "comparison_summary.csv"
    ).open(
        newline="",
        encoding="utf-8",
    ) as file:

        rows = list(
            csv.DictReader(file)
        )

    assert len(rows) == 2

    # --------------------------------------------------------
    # Tool agreement
    # --------------------------------------------------------

    with (
        output_dir / "tool_agreement.csv"
    ).open(
        newline="",
        encoding="utf-8",
    ) as file:

        agreement = list(
            csv.DictReader(file)
        )

    assert len(agreement) == 1

    assert (
        agreement[0]["agreement_type"]
        == "same_cds_status_only"
    )

    # --------------------------------------------------------
    # Report
    # --------------------------------------------------------

    report = json.loads(
        (
            output_dir / "comparison_report.json"
        ).read_text(
            encoding="utf-8"
        )
    )

    assert report["processed_comparisons"] == 2

    assert report["unique_genes"] == 1

    assert (
        report["alignment_status_counts"]["aligned"]
        == 2
    )


# ============================================================
# 12. MISSING REFERENCE
# ============================================================

def test_missing_reference(tmp_path):

    input_dir = (
        tmp_path
        / "3_results_second"
    )

    input_dir.mkdir()

    input_path = (
        input_dir / "orf_records.jsonl"
    )

    record = make_record(
        gene_id="FBgn0000001",
        species="dsim",
        section="spaln",
        sequence="ATGAAACCCGGGTAA",
        locus="chr2L+:100:115",
    )

    input_path.write_text(
        json.dumps(record) + "\n",
        encoding="utf-8",
    )

    report = comparative_orf.run(
        input_path=input_path,
        comparisons_path=(
            tmp_path
            / "4_results_third"
            / "comparisons.jsonl"
        ),
        summary_path=(
            tmp_path
            / "4_results_third"
            / "summary.csv"
        ),
        agreement_path=(
            tmp_path
            / "4_results_third"
            / "agreement.csv"
        ),
        report_path=(
            tmp_path
            / "4_results_third"
            / "report.json"
        ),
        verbose=False,
    )

    assert (
        report["alignment_status_counts"][
            "skipped_missing_reference"
        ] == 1
    )


# ============================================================
# 13. MISSING INPUT
# ============================================================

def test_cli_missing_input(tmp_path):

    command = [
        sys.executable,
        str(
            Path(comparative_orf.__file__).resolve()
        ),

        "--input",
        "3_results_second/missing.jsonl",

        "--comparisons-out",
        "4_results_third/comparisons.jsonl",

        "--summary-out",
        "4_results_third/summary.csv",

        "--agreement-out",
        "4_results_third/agreement.csv",

        "--report-out",
        "4_results_third/report.json",
    ]

    completed = subprocess.run(
        command,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 2

    assert (
        "Входной файл не найден"
        in completed.stderr
    )
