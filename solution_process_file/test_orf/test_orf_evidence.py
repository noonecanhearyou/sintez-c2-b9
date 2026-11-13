"""
Tests for ORF Evidence Analysis — Stage 4.

Run:
    python -m pytest -v test_orf_evidence.py
"""

import csv
import json
import subprocess
import sys

from pathlib import Path

import orf_evidence

from orf_evidence import (
    analyze_gap_path,
    assess_evidence,
    load_orf_profiles,
    match_profile,
)


GENE = "FBgn0000001"
LOCUS = "chr2L+:100:115"


# ------------------------------------------------------------
# 1. Test fixtures
# ------------------------------------------------------------

def make_comparison(events=None, **overrides):

    if events is None:
        events = []

    nontriplet = sum(
        1 for event in events
        if (
            event["internal"]
            and event["length_nt"] % 3 != 0
        )
    )

    triplet = sum(
        1 for event in events
        if (
            event["internal"]
            and event["length_nt"] % 3 == 0
        )
    )

    terminal = sum(
        1 for event in events
        if not event["internal"]
    )

    comparison = {
        "status": "aligned",
        "quality_ok": True,
        "identity": 1.0,
        "reference_coverage": 1.0,
        "query_coverage": 1.0,
        "gap_event_count": len(events),
        "gap_examples": events,
        "gap_examples_truncated": False,
        "internal_nontriplet_gap_events": nontriplet,
        "internal_triplet_gap_events": triplet,
        "terminal_gap_events": terminal,
    }

    comparison.update(overrides)

    return comparison


def gap_in_query(length=1):

    return {
        "kind": "gap_in_query",
        "reference_start": 3,
        "reference_end": 3 + length,
        "query_position": 3,
        "length_nt": length,
        "internal": True,
        "nontriplet": length % 3 != 0,
    }


def gap_in_reference(length=1):

    return {
        "kind": "gap_in_reference",
        "reference_position": 9,
        "query_start": 8,
        "query_end": 8 + length,
        "length_nt": length,
        "internal": True,
        "nontriplet": length % 3 != 0,
    }


def make_stage2(method, gene=GENE):

    return {
        "gene_id": gene,
        "species": "dsim",
        "section": method,
        "locus": LOCUS,
        "sequence_kind": "predicted_cds",
        "sequence": "ATGAAACCCGGGTAA",
        "orf_scan": {
            "longest_complete_aa": 4,
        },
        "cds_validation": {
            "status": "complete_candidate",
            "internal_stop_count": 0,
        },
        "source_metadata": {
            "shift_label": "No",
        },
    }


def make_stage3(method, gene=GENE):

    return {
        "gene_id": gene,
        "species": "dsim",
        "method": method,
        "locus": LOCUS,
        "reference": {
            "length_nt": 15,
            "cds_status": "complete_candidate",
        },
        "target": {
            "length_nt": 15,
            "cds_status": "complete_candidate",
            "reported_shift": "No",
        },
        "comparison": make_comparison(),
    }


def write_jsonl(path, records):

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
    ) as handle:

        for record in records:
            handle.write(json.dumps(record) + "\n")


# ------------------------------------------------------------
# 2. No gap events
# ------------------------------------------------------------

def test_no_gaps():

    result = analyze_gap_path(
        make_comparison()
    )

    assert result["status"] == "no_gaps"
    assert result["events_complete"] is True
    assert result["phase_final_mod3"] == 0
    assert result["phase_restorations"] == 0


# ------------------------------------------------------------
# 3. In-frame triplet gap
# ------------------------------------------------------------

def test_triplet_gap():

    comparison = make_comparison([
        gap_in_query(length=3)
    ])

    result = analyze_gap_path(comparison)

    assert result["status"] == "in_frame_gaps_only"
    assert result["phase_final_mod3"] == 0

    assert result["internal_nontriplet_observed"] == 0
    assert result["codon_boundary_gap_events"] == 1


# ------------------------------------------------------------
# 4. Nontriplet gap without restoration
# ------------------------------------------------------------

def test_unresolved_phase():

    comparison = make_comparison([
        gap_in_query(length=1)
    ])

    result = analyze_gap_path(comparison)

    assert result["status"] == "phase_unresolved"
    assert result["phase_final_mod3"] == 2
    assert result["phase_restorations"] == 0

    assert result["details"][0]["phase_before"] == 0
    assert result["details"][0]["phase_after"] == 2


# ------------------------------------------------------------
# 5. Phase restoration by subsequent gap
# ------------------------------------------------------------

def test_phase_restoration():

    comparison = make_comparison([
        gap_in_query(length=1),
        gap_in_reference(length=1),
    ])

    result = analyze_gap_path(comparison)

    assert result["status"] == "phase_restored"

    assert result["phase_final_mod3"] == 0
    assert result["phase_restorations"] == 1
    assert result["phase_excursions"] == 1

    assert result["internal_nontriplet_observed"] == 2


# ------------------------------------------------------------
# 6. Truncated gap examples
# ------------------------------------------------------------

def test_truncated_gap_examples():

    comparison = make_comparison(
        [gap_in_query(length=1)],
        gap_event_count=3,
        gap_examples_truncated=True,
    )

    result = analyze_gap_path(comparison)

    assert result["status"] == "incomplete_gap_examples"

    assert result["events_complete"] is False
    assert result["phase_final_mod3"] is None

    assert result["phase_restorations"] is None


# ------------------------------------------------------------
# 7. Inconsistent gap counts
# ------------------------------------------------------------

def test_inconsistent_gap_counts():

    comparison = make_comparison(
        [gap_in_query(length=1)],
        internal_nontriplet_gap_events=2,
    )

    result = analyze_gap_path(comparison)

    assert result["status"] == "inconsistent_gap_counts"
    assert result["events_complete"] is False


# ------------------------------------------------------------
# 8. Conservative evidence assessment
# ------------------------------------------------------------

def test_evidence_assessment():

    no_gaps = analyze_gap_path(
        make_comparison()
    )

    result = assess_evidence(
        "matched",
        "complete_candidate",
        make_comparison(),
        no_gaps,
    )

    assert result[0] == "preserved_cds_candidate"

    disrupted = assess_evidence(
        "matched",
        "interrupted",
        make_comparison(),
        no_gaps,
    )

    assert disrupted[0] == "potential_cds_disruption"

    one_gap = make_comparison([
        gap_in_query(length=1)
    ])

    conflicting = assess_evidence(
        "matched",
        "complete_candidate",
        one_gap,
        analyze_gap_path(one_gap),
    )

    assert conflicting[0] == "conflicting_signals"

    low_quality = make_comparison(
        status="low_quality_alignment",
        quality_ok=False,
    )

    insufficient = assess_evidence(
        "matched",
        "complete_candidate",
        low_quality,
        no_gaps,
    )

    assert insufficient[0] == "insufficient_evidence"


# ------------------------------------------------------------
# 9. Detect mismatched Stage-2 profiles
# ------------------------------------------------------------

def test_profile_matching_and_duplicates(tmp_path):

    path = tmp_path / "orf_records.jsonl"

    write_jsonl(
        path,
        [make_stage2("spaln")],
    )

    index, count = load_orf_profiles(path)

    assert count == 1

    comparison = make_stage3("spaln")

    status, profile = match_profile(
        index,
        comparison,
    )

    assert status == "matched"
    assert profile["length_nt"] == 15

    comparison["target"]["length_nt"] = 14

    status, _ = match_profile(index, comparison)

    assert status == "length_mismatch"

    # Two records with the same key must be
    # reported as ambiguous rather than selected.

    write_jsonl(
        path,
        [
            make_stage2("spaln"),
            make_stage2("spaln"),
        ],
    )

    duplicate_index, count = load_orf_profiles(path)

    assert count == 2

    comparison["target"]["length_nt"] = 15

    status, profile = match_profile(
        duplicate_index,
        comparison,
    )

    assert status == "ambiguous_multiple_profiles"
    assert profile is None


# ------------------------------------------------------------
# 10. Full CLI pipeline
# ------------------------------------------------------------

def test_cli_pipeline(tmp_path):

    input2 = (
        tmp_path
        / "3_results_second"
        / "orf_records.jsonl"
    )

    input3 = (
        tmp_path
        / "4_results_third"
        / "comparisons.jsonl"
    )

    agreement = (
        tmp_path
        / "4_results_third"
        / "tool_agreement.csv"
    )

    write_jsonl(
        input2,
        [
            make_stage2("spaln"),
            make_stage2("genewise"),
            make_stage2("spaln", gene="FBgn9999999"),
        ],
    )

    write_jsonl(
        input3,
        [
            make_stage3("spaln"),
            make_stage3("genewise"),
            make_stage3("spaln", gene="FBgn9999999"),
        ],
    )

    agreement.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with agreement.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene_id",
                "species",
                "locus",
                "agreement_type",
            ],
        )

        writer.writeheader()

        writer.writerow({
            "gene_id": GENE,
            "species": "dsim",
            "locus": LOCUS,
            "agreement_type": "exact_sequence",
        })

    command = [
        sys.executable,
        str(Path(orf_evidence.__file__).resolve()),

        "--orf-input",
        "3_results_second/orf_records.jsonl",

        "--comparisons-input",
        "4_results_third/comparisons.jsonl",

        "--agreement-input",
        "4_results_third/tool_agreement.csv",

        "--records-out",
        "5_results_fourth/evidence_records.jsonl",

        "--summary-out",
        "5_results_fourth/evidence_summary.csv",

        "--genes-out",
        "5_results_fourth/gene_summary.csv",

        "--report-out",
        "5_results_fourth/evidence_report.json",

        "--gene",
        GENE,

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

    assert "ORF Evidence Analysis завершён" in (
        completed.stdout
    )

    output = tmp_path / "5_results_fourth"

    for filename in [
        "evidence_records.jsonl",
        "evidence_summary.csv",
        "gene_summary.csv",
        "evidence_report.json",
    ]:
        assert (output / filename).is_file()

    with (
        output / "evidence_records.jsonl"
    ).open(encoding="utf-8") as handle:

        records = [
            json.loads(line)
            for line in handle
            if line.strip()
        ]

    assert len(records) == 2

    assert all(
        record["gene_id"] == GENE
        for record in records
    )

    assert all(
        record["evidence"]["class"]
        == "preserved_cds_candidate"
        for record in records
    )

    with (
        output / "gene_summary.csv"
    ).open(
        encoding="utf-8",
        newline="",
    ) as handle:

        genes = list(csv.DictReader(handle))

    assert len(genes) == 1

    assert genes[0]["comparisons"] == "2"
    assert genes[0]["unique_species"] == "1"
    assert genes[0]["exact_agreement_loci"] == "1"

    report = json.loads(
        (output / "evidence_report.json").read_text(
            encoding="utf-8"
        )
    )

    assert report["processed_comparisons"] == 2
    assert report["unique_genes"] == 1

    assert (
        report["evidence_class_counts"][
            "preserved_cds_candidate"
        ] == 2
    )


# ------------------------------------------------------------
# 11. Missing input
# ------------------------------------------------------------

def test_cli_missing_input(tmp_path):

    command = [
        sys.executable,
        str(Path(orf_evidence.__file__).resolve()),

        "--orf-input",
        "missing_orf.jsonl",

        "--comparisons-input",
        "missing_comparisons.jsonl",

        "--agreement-input",
        "missing_agreement.csv",

        "--records-out",
        "5_results_fourth/records.jsonl",

        "--summary-out",
        "5_results_fourth/summary.csv",

        "--genes-out",
        "5_results_fourth/genes.csv",

        "--report-out",
        "5_results_fourth/report.json",
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
