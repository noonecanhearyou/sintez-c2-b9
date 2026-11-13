"""
ORF Scanner — Stage 2 Tests

Запуск:

    python -m pytest -v test_orf_scanner.py
"""

import csv
import json
import subprocess
import sys

from pathlib import Path

from orf_scanner import (
    reverse_complement,
    scan_orfs,
    validate_cds,
)

import orf_scanner


# ============================================================
# 1. REVERSE COMPLEMENT
# ============================================================

def test_reverse_complement():

    assert reverse_complement(
        "ATGCCG"
    ) == "CGGCAT"

    assert reverse_complement(
        "ATGN"
    ) == "NCAT"


# ============================================================
# 2. COMPLETE ORF
# ============================================================

def test_complete_orf():

    sequence = "ATGAAATAA"

    result = scan_orfs(
        sequence,
        min_aa=2,
    )

    assert result["counts_all"]["complete"] >= 1

    candidates = [
        orf
        for orf in result["top_orfs"]
        if (
            orf["frame"] == "+1"
            and orf["type"] == "complete"
        )
    ]

    assert len(candidates) == 1

    orf = candidates[0]

    assert orf["start_nt"] == 0
    assert orf["end_nt"] == 9
    assert orf["aa_length"] == 2

    assert orf["stop_codon"] == "TAA"

    assert orf["peptide_preview"] == "MK"


# ============================================================
# 3. REVERSE STRAND ORF
# ============================================================

def test_reverse_strand_orf():

    # Обратная комплементарная:
    # ATGAAATAA

    sequence = "TTATTTCAT"

    result = scan_orfs(
        sequence,
        min_aa=2,
    )

    candidates = [
        orf
        for orf in result["top_orfs"]
        if (
            orf["frame"] == "-1"
            and orf["type"] == "complete"
        )
    ]

    assert len(candidates) == 1

    orf = candidates[0]

    assert orf["start_nt"] == 0
    assert orf["end_nt"] == 9
    assert orf["aa_length"] == 2

    assert orf["peptide_preview"] == "MK"


# ============================================================
# 4. NESTED ORFS
# ============================================================

def test_nested_orfs():

    # ATG | ATG | AAA | TAA

    sequence = "ATGATGAAATAA"

    result = scan_orfs(
        sequence,
        min_aa=2,
    )

    candidates = [
        orf
        for orf in result["top_orfs"]
        if (
            orf["frame"] == "+1"
            and orf["type"] == "complete"
        )
    ]

    lengths = sorted(
        orf["aa_length"]
        for orf in candidates
    )

    assert lengths == [2, 3]


# ============================================================
# 5. ORF WITHOUT STOP
# ============================================================

def test_open_end_orf():

    sequence = "ATGAAAAAA"

    result = scan_orfs(
        sequence,
        min_aa=2,
    )

    candidates = [
        orf
        for orf in result["top_orfs"]
        if (
            orf["frame"] == "+1"
            and orf["type"] == "open_end"
        )
    ]

    assert len(candidates) == 1

    assert candidates[0]["aa_length"] == 3

    assert candidates[0]["stop_codon"] is None


# ============================================================
# 6. AMBIGUOUS NUCLEOTIDES
# ============================================================

def test_ambiguous_boundary():

    sequence = "ATGNNNTAA"

    result = scan_orfs(
        sequence,
        min_aa=1,
    )

    assert result["counts_all"]["complete"] == 0

    assert (
        result["counts_all"]["ambiguous_boundary"]
        >= 1
    )

    assert not any(
        orf["type"] == "complete"
        for orf in result["top_orfs"]
    )


# ============================================================
# 7. COMPLETE CDS VALIDATION
# ============================================================

def test_complete_cds():

    result = validate_cds(
        "ATGAAATAA",
        "predicted_cds",
    )

    assert result["status"] == "complete_candidate"

    assert result["start_atg"] is True

    assert result["terminal_stop"] is True

    assert result["internal_stop_count"] == 0

    assert result["length_mod3"] == 0


# ============================================================
# 8. INTERNAL STOP
# ============================================================

def test_internal_stop():

    # ATG | TAG | TAA

    result = validate_cds(
        "ATGTAGTAA",
        "annotated_cds",
    )

    assert result["status"] == "interrupted"

    assert result["internal_stop_count"] == 1

    assert result["internal_stop_positions"] == [3]


# ============================================================
# 9. PARTIAL CDS
# ============================================================

def test_partial_cds():

    result = validate_cds(
        "ATGAAAAAA",
        "predicted_cds",
    )

    assert (
        result["status"]
        == "partial_or_noncanonical"
    )

    assert result["terminal_stop"] is False


# ============================================================
# 10. UNCERTAIN CDS
# ============================================================

def test_uncertain_cds():

    result = validate_cds(
        "ATGNNNTAA",
        "predicted_cds",
    )

    assert result["status"] == "uncertain"

    assert result["ambiguous_codon_count"] == 1


# ============================================================
# 11. GENOMIC REGION
# ============================================================

def test_genomic_region_not_cds():

    result = validate_cds(
        "ATGAAATAA",
        "genomic_region",
    )

    assert result["status"] == "not_applicable"


# ============================================================
# 12. EMPTY CDS
# ============================================================

def test_empty_cds():

    result = validate_cds(
        "",
        "predicted_cds",
    )

    assert result["status"] == "empty"


# ============================================================
# 13. FULL CLI TEST
# ============================================================

def test_cli_pipeline(tmp_path):

    input_dir = tmp_path / "2_results_first"

    input_dir.mkdir()

    input_path = input_dir / "records.jsonl"

    records = [
        {
            "gene_id": "FBgn0000001",
            "species": "dmel",
            "section": "raw",
            "locus": "FBgn0000001",
            "sequence_kind": "annotated_cds",
            "sequence": "ATGAAATAA",
            "score_label": "1E10",
            "shift_label": "No",
        },
        {
            "gene_id": "FBgn0000001",
            "species": "dsim",
            "section": "spaln",
            "locus": "chr2L+:100:109",
            "sequence_kind": "predicted_cds",
            "sequence": "ATGTAGTAA",
            "score_label": "10.5",
            "shift_label": "Yes",
        },
        {
            "gene_id": "FBgn0000001",
            "species": "dsim",
            "section": "genewise",
            "locus": "chr2L+:100:109",
            "sequence_kind": "predicted_cds",
            "sequence": "",
            "score_label": "9.2",
            "shift_label": "No",
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
        str(Path(orf_scanner.__file__).resolve()),

        "--input",
        "2_results_first/records.jsonl",

        "--records-out",
        "3_results_second/orf_records.jsonl",

        "--summary-out",
        "3_results_second/orf_summary.csv",

        "--report-out",
        "3_results_second/orf_report.json",

        "--gene",
        "FBgn0000001",

        "--min-aa",
        "2",

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

    assert "Начало ORF Analysis" in completed.stdout

    assert "Обработано записей:" in completed.stdout

    output_dir = tmp_path / "3_results_second"

    assert (
        output_dir / "orf_records.jsonl"
    ).is_file()

    assert (
        output_dir / "orf_summary.csv"
    ).is_file()

    assert (
        output_dir / "orf_report.json"
    ).is_file()

    # --------------------------------------------------------
    # Проверка JSON-отчёта
    # --------------------------------------------------------

    report = json.loads(
        (output_dir / "orf_report.json").read_text(
            encoding="utf-8"
        )
    )

    assert report["processed_records"] == 3

    assert report["unique_genes"] == 1

    assert (
        report["cds_status_counts"]["complete_candidate"]
        == 1
    )

    assert (
        report["cds_status_counts"]["interrupted"]
        == 1
    )

    assert (
        report["cds_status_counts"]["empty"]
        == 1
    )

    # --------------------------------------------------------
    # Проверка CSV
    # --------------------------------------------------------

    with (
        output_dir / "orf_summary.csv"
    ).open(
        newline="",
        encoding="utf-8",
    ) as file:

        rows = list(
            csv.DictReader(file)
        )

    assert len(rows) == 3

    # --------------------------------------------------------
    # Проверка JSONL
    # --------------------------------------------------------

    with (
        output_dir / "orf_records.jsonl"
    ).open(
        encoding="utf-8",
    ) as file:

        results = [
            json.loads(line)
            for line in file
            if line.strip()
        ]

    assert len(results) == 3

    assert all(
        "orf_scan" in record
        for record in results
    )

    assert all(
        "cds_validation" in record
        for record in results
    )
