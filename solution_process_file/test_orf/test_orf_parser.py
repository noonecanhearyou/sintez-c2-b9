"""
Тесты ORF Parser — Stage 1.

Запуск:

    python -m pytest -v test_orf_parser.py

Для работы необходим pytest.
"""

import csv
import json
import subprocess
import sys

from pathlib import Path

import orf_parser

from orf_parser import (
    iter_file_events,
    run,
)


# ============================================================
# 1. TEST DATA
# ============================================================

DATA = """Section I. Raw hits extracted from progressive cactus alignments
raw sequence for dmel_FBgn0000001
>dmel|FBgn0000001|score_1E10|shift_No|1
ATGAAATAA
>dsim|chr2L-:100:112
ATGRYN

Section II. Spliced align by Spaln
spaln for dmel_FBgn0000001
>dmel|FBgn0000001|score_1E10|shift_No|1
ATGAAATAA
>dsim|chr2L-:100:112|score_10.5|shift_No|1

Section III. Spliced align by Genewise
genewise for dmel_FBgn0000001
>dmel|FBgn0000001|score_1E10|shift_No|1
ATGAAATAA
>dsim|chr2L-:100:112|score_9.2|shift_Yes|2
ATGTAA
"""


# ============================================================
# 2. TEST HELPER
# ============================================================

def execute(
    tmp_path: Path,
    data: str = DATA,
    gene: str | None = None,
):

    """
    Создаёт временный входной файл,
    запускает парсер и считывает результаты.
    """

    input_path = (
        tmp_path
        / "1_raw_data"
        / "input.txt"
    )

    input_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    input_path.write_text(
        data,
        encoding="utf-8"
    )

    records_path = (
        tmp_path
        / "2_results_first"
        / "records.jsonl"
    )

    summary_path = (
        tmp_path
        / "2_results_first"
        / "summary.csv"
    )

    report_path = (
        tmp_path
        / "2_results_first"
        / "validation.json"
    )

    report = run(
        input_path=input_path,
        records_path=records_path,
        summary_path=summary_path,
        report_path=report_path,
        gene_filter=gene,
        progress_every=0,
        verbose=False,
    )

    # Чтение JSONL
    with records_path.open(
        encoding="utf-8"
    ) as file:

        records = [
            json.loads(line)
            for line in file
            if line.strip()
        ]

    # Чтение CSV
    with summary_path.open(
        newline="",
        encoding="utf-8"
    ) as file:

        rows = list(
            csv.DictReader(file)
        )

    # Чтение отчёта
    disk_report = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    return (
        report,
        records,
        rows,
        disk_report,
    )


# ============================================================
# 3. BASIC PARSER TEST
# ============================================================

def test_parse_and_validate_with_empty_prediction(tmp_path):

    report, records, rows, disk_report = execute(
        tmp_path
    )

    assert report == disk_report

    assert len(records) == 6
    assert len(rows) == 6

    assert report["unique_genes"] == 1

    assert report["error_count"] == 0

    assert report["warning_counts"] == {
        "empty_prediction": 1
    }

    assert report["records_per_section"] == {
        "raw": 2,
        "spaln": 2,
        "genewise": 2,
    }

    raw = next(
        record
        for record in records
        if (
            record["section"] == "raw"
            and record["species"] == "dsim"
        )
    )

    assert raw["sequence_kind"] == "genomic_region"

    assert (
        raw["contig"],
        raw["strand"],
        raw["start"],
        raw["end"]
    ) == (
        "chr2L",
        "-",
        100,
        112
    )

    assert raw["sequence"] == "ATGRYN"

    assert raw["ambiguous_nt"] == 3

    spaln = next(
        record
        for record in records
        if (
            record["section"] == "spaln"
            and record["species"] == "dsim"
        )
    )

    assert spaln["sequence_kind"] == "predicted_cds"

    assert spaln["sequence"] == ""

    assert spaln["issues"] == [
        "empty_prediction"
    ]


# ============================================================
# 4. GENE FILTER TEST
# ============================================================

def test_gene_filter(tmp_path):

    data = DATA.replace(
        "Section II. Spliced align by Spaln",
        (
            "raw sequence for dmel_FBgn0000002\n"
            ">dmel|FBgn0000002|score_1E10|shift_No|1\n"
            "ATGTAA\n"
            "Section II. Spliced align by Spaln"
        )
    ).replace(
        "Section III. Spliced align by Genewise",
        (
            "spaln for dmel_FBgn0000002\n"
            ">dmel|FBgn0000002|score_1E10|shift_No|1\n"
            "ATGTAA\n"
            "Section III. Spliced align by Genewise"
        )
    ) + (
        "genewise for dmel_FBgn0000002\n"
        ">dmel|FBgn0000002|score_1E10|shift_No|1\n"
        "ATGTAA\n"
    )

    filtered, records, _, _ = execute(
        tmp_path,
        data=data,
        gene="FBgn0000001"
    )

    assert filtered["unique_genes"] == 1

    assert len(records) == 6

    assert all(
        record["gene_id"] == "FBgn0000001"
        for record in records
    )

    full, full_records, _, _ = execute(
        tmp_path,
        data=data
    )

    assert full["unique_genes"] == 2

    assert len(full_records) == 9


# ============================================================
# 5. MISSING GENE TEST
# ============================================================

def test_missing_gene(tmp_path):

    report, records, _, _ = execute(
        tmp_path,
        gene="FBgn9999999"
    )

    assert (
        report["error_counts"]["requested_gene_not_found"]
        == 1
    )

    assert records == []


# ============================================================
# 6. INVALID DNA TEST
# ============================================================

def test_invalid_dna_character(tmp_path):

    invalid_data = DATA.replace(
        "ATGRYN",
        "ATGRZYN"
    )

    report, records, _, _ = execute(
        tmp_path,
        data=invalid_data
    )

    assert (
        report["error_counts"]["invalid_dna_characters"]
        == 1
    )

    assert any(
        record["invalid_nt"] == 1
        for record in records
    )


# ============================================================
# 7. MULTILINE SEQUENCE TEST
# ============================================================

def test_multiline_sequence(tmp_path):

    data = DATA.replace(
        "ATGAAATAA\n>dsim",
        "ATGA\nAATAA\n>dsim"
    )

    report, records, _, _ = execute(
        tmp_path,
        data=data
    )

    assert report["error_count"] == 0

    assert records[0]["sequence"] == "ATGAAATAA"


# ============================================================
# 8. EVENT STREAM TEST
# ============================================================

def test_event_stream(tmp_path):

    input_path = tmp_path / "input.txt"

    input_path.write_text(
        DATA,
        encoding="utf-8"
    )

    events = list(
        iter_file_events(input_path)
    )

    sections = [
        value
        for name, value in events
        if name == "section"
    ]

    assert sections == [
        "raw",
        "spaln",
        "genewise",
    ]

    record_count = sum(
        name == "record"
        for name, _ in events
    )

    assert record_count == 6


# ============================================================
# 9. CLI TEST
# ============================================================

def test_cli_creates_outputs_and_prints_progress(tmp_path):

    source = (
        tmp_path
        / "1_raw_data"
        / "input.txt"
    )

    source.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    source.write_text(
        DATA,
        encoding="utf-8"
    )

    command = [
        sys.executable,
        str(Path(orf_parser.__file__).resolve()),

        "--input",
        "1_raw_data/input.txt",

        "--records-out",
        "2_results_first/records.jsonl",

        "--summary-out",
        "2_results_first/summary.csv",

        "--report-out",
        "2_results_first/validation.json",

        "--gene",
        "FBgn0000001",

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
        "[ORF Parser] Начало обработки:"
        in completed.stdout
    )

    assert "Прогресс:" in completed.stdout

    assert "Готово:" in completed.stdout

    assert (
        tmp_path
        / "2_results_first"
        / "records.jsonl"
    ).is_file()

    assert (
        tmp_path
        / "2_results_first"
        / "summary.csv"
    ).is_file()

    assert (
        tmp_path
        / "2_results_first"
        / "validation.json"
    ).is_file()


# ============================================================
# 10. CLI ERROR TEST
# ============================================================

def test_cli_missing_input_reports_error(tmp_path):

    command = [
        sys.executable,
        str(Path(orf_parser.__file__).resolve()),

        "--input",
        "1_raw_data/not_found.txt",

        "--records-out",
        "2_results_first/records.jsonl",

        "--summary-out",
        "2_results_first/summary.csv",

        "--report-out",
        "2_results_first/validation.json",
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
