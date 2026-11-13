"""
Tests for ORF Profile Aggregation — Stage 5.

Run:
    python -m pytest -v test_orf_profile.py
"""

import csv
import json
import subprocess
import sys

from pathlib import Path

import pytest
import orf_profile

from orf_profile import (
    compact_record,
    classify_locus,
    aggregate_species,
    build_gene_profile,
)


GENE = "FBgn0000001"
LOCUS = "chr2L+:100:115"


# ============================================================
# 1. TEST DATA
# ============================================================

def make_record(
    method,
    evidence="preserved_cds_candidate",
    cds_status="complete_candidate",
    agreement="exact_sequence",
    alignment_status="aligned",
    locus=LOCUS,
    gene=GENE,
    species="dsim",
    nontriplet=0,
    gap_status="no_gaps",
):
    return {
        "gene_id": gene,
        "species": species,
        "method": method,
        "locus": locus,

        "stage2_profile_match": "matched",

        "cds": {
            "status": cds_status,
            "internal_stop_count": (
                1 if cds_status == "interrupted" else 0
            ),
            "longest_complete_aa": 4,
            "reported_shift": "No",
        },

        "tool_agreement": agreement,

        "alignment": {
            "status": alignment_status,
            "identity": 0.95,
            "reference_coverage": 1.0,
            "query_coverage": 1.0,
            "internal_nontriplet_gap_events": nontriplet,
        },

        "gap_path": {
            "status": gap_status,
        },

        "evidence": {
            "class": evidence,
            "reason": "test",
            "evolutionary_direction_inferred": False,
            "confirmed_frameshift": False,
        },
    }


def classify(*records):

    compact = [
        compact_record(record)
        for record in records
    ]

    return classify_locus(
        GENE,
        "dsim",
        LOCUS,
        compact,
    )


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
            handle.write(
                json.dumps(record) + "\n"
            )


# ============================================================
# 2. IDENTICAL PRESERVED CDS
# ============================================================

def test_preserved_identical():

    row = classify(
        make_record("spaln"),
        make_record("genewise"),
    )

    assert row["status"] == "preserved_both_identical"
    assert row["spaln_records"] == 1
    assert row["genewise_records"] == 1
    assert row["comparison_records"] == 2


# ============================================================
# 3. DIFFERENT SEQUENCES, SAME CDS STATUS
# ============================================================

def test_preserved_different():

    row = classify(
        make_record(
            "spaln",
            agreement="same_cds_status_only",
        ),
        make_record(
            "genewise",
            agreement="same_cds_status_only",
        ),
    )

    assert row["status"] == "preserved_both_different"


# ============================================================
# 4. INTERRUPTED CDS FROM BOTH METHODS
# ============================================================

def test_both_interrupted():

    row = classify(
        make_record(
            "spaln",
            evidence="potential_cds_disruption",
            cds_status="interrupted",
        ),
        make_record(
            "genewise",
            evidence="potential_cds_disruption",
            cds_status="interrupted",
        ),
    )

    assert row["status"] == "predicted_interruption_both"


# ============================================================
# 5. DISCORDANT PREDICTIONS
# ============================================================

def test_discordant_methods():

    row = classify(
        make_record(
            "spaln",
            agreement="different_cds_status",
        ),
        make_record(
            "genewise",
            evidence="potential_cds_disruption",
            cds_status="interrupted",
            agreement="different_cds_status",
        ),
    )

    assert row["status"] == "discordant_methods"


# ============================================================
# 6. ALIGNMENT CONFLICT
# ============================================================

def test_alignment_conflict():

    row = classify(
        make_record(
            "spaln",
            evidence="conflicting_signals",
            nontriplet=2,
            gap_status="phase_restored",
        ),
        make_record("genewise"),
    )

    assert row["status"] == "conflicting_alignment"

    assert row["nontriplet_gap_observations"] == 2
    assert row["phase_restored_observations"] == 1


# ============================================================
# 7. PARTIAL EVIDENCE
# ============================================================

def test_partial_evidence():

    row = classify(
        make_record("spaln"),
        make_record(
            "genewise",
            evidence="insufficient_evidence",
            alignment_status="low_quality_alignment",
        ),
    )

    assert row["status"] == "partial_evidence"


# ============================================================
# 8. DUPLICATE PREDICTIONS
# ============================================================

def test_duplicate_predictions():

    row = classify(
        make_record("spaln"),
        make_record("spaln"),
        make_record("genewise"),
    )

    assert row["status"] == "multiple_predictions"


# ============================================================
# 9. INCONSISTENT METADATA
# ============================================================

def test_inconsistent_metadata():

    row = classify(
        make_record(
            "spaln",
            agreement="exact_sequence",
        ),
        make_record(
            "genewise",
            agreement="same_cds_status_only",
        ),
    )

    assert row["status"] == "inconsistent_metadata"


# ============================================================
# 10. MULTIPLE LOCI
# ============================================================

def test_multiple_loci():

    rows = [
        {
            "gene_id": GENE,
            "species": "dsim",
            "locus": "chr2L:100",
            "status": "preserved_both_identical",
            "comparison_records": 2,
            "nontriplet_gap_observations": 0,
            "phase_restored_observations": 0,
        },
        {
            "gene_id": GENE,
            "species": "dsim",
            "locus": "chr3R:500",
            "status": "preserved_both_identical",
            "comparison_records": 2,
            "nontriplet_gap_observations": 0,
            "phase_restored_observations": 0,
        },
    ]

    result = aggregate_species(rows)

    assert len(result) == 1

    assert result[0]["status"] == "multiple_loci_uncertain"
    assert result[0]["locus_count"] == 2

    profile, summary = build_gene_profile(
        GENE,
        result,
    )

    assert profile["species_assessed"] == 1
    assert profile["loci_assessed"] == 2

    assert summary[
        "species_multiple_loci_uncertain"
    ] == 1


# ============================================================
# 11. FULL CLI PIPELINE
# ============================================================

def test_cli_pipeline(tmp_path):

    input_path = (
        tmp_path
        / "5_results_fourth"
        / "evidence_records.jsonl"
    )

    records = [
        make_record("spaln"),
        make_record("genewise"),

        # Another gene to test --gene filtering.
        make_record(
            "spaln",
            gene="FBgn9999999",
        ),
    ]

    write_jsonl(input_path, records)

    command = [
        sys.executable,
        str(Path(orf_profile.__file__).resolve()),

        "--input",
        "5_results_fourth/evidence_records.jsonl",

        "--locus-out",
        "6_results_fifth/locus_profiles.csv",

        "--species-out",
        "6_results_fifth/species_profiles.csv",

        "--genes-out",
        "6_results_fifth/gene_profiles.jsonl",

        "--summary-out",
        "6_results_fifth/gene_summary.csv",

        "--report-out",
        "6_results_fifth/profile_report.json",

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

    assert (
        "ORF Profile Aggregation завершён"
        in completed.stdout
    )

    output = tmp_path / "6_results_fifth"

    for filename in (
        "locus_profiles.csv",
        "species_profiles.csv",
        "gene_profiles.jsonl",
        "gene_summary.csv",
        "profile_report.json",
    ):
        assert (output / filename).is_file()

    with (
        output / "locus_profiles.csv"
    ).open(
        encoding="utf-8",
        newline="",
    ) as handle:

        loci = list(csv.DictReader(handle))

    assert len(loci) == 1

    assert loci[0]["status"] == (
        "preserved_both_identical"
    )

    with (
        output / "species_profiles.csv"
    ).open(
        encoding="utf-8",
        newline="",
    ) as handle:

        species = list(csv.DictReader(handle))

    assert len(species) == 1
    assert species[0]["locus_count"] == "1"

    with (
        output / "gene_profiles.jsonl"
    ).open(encoding="utf-8") as handle:

        profiles = [
            json.loads(line)
            for line in handle
            if line.strip()
        ]

    assert len(profiles) == 1
    assert profiles[0]["gene_id"] == GENE
    assert profiles[0]["species_assessed"] == 1
    assert profiles[0]["origin_classification"] is None

    report = json.loads(
        (output / "profile_report.json").read_text(
            encoding="utf-8"
        )
    )

    assert report["comparison_records"] == 2
    assert report["unique_loci"] == 1
    assert report["gene_species_pairs"] == 1
    assert report["unique_genes"] == 1


# ============================================================
# 12. MISSING GENE
# ============================================================

def test_missing_gene(tmp_path):

    input_path = (
        tmp_path
        / "5_results_fourth"
        / "evidence_records.jsonl"
    )

    write_jsonl(
        input_path,
        [make_record("spaln")],
    )

    with pytest.raises(
        ValueError,
        match="Не найдено сравнений",
    ):

        orf_profile.run(
            input_path=input_path,
            locus_out=tmp_path / "out/loci.csv",
            species_out=tmp_path / "out/species.csv",
            genes_out=tmp_path / "out/genes.jsonl",
            summary_out=tmp_path / "out/summary.csv",
            report_out=tmp_path / "out/report.json",
            gene_filter="FBgn9999999",
            verbose=False,
        )


# ============================================================
# 13. MISSING INPUT FILE
# ============================================================

def test_missing_input(tmp_path):

    command = [
        sys.executable,
        str(Path(orf_profile.__file__).resolve()),

        "--input",
        "missing.jsonl",

        "--locus-out",
        "6_results_fifth/locus.csv",

        "--species-out",
        "6_results_fifth/species.csv",

        "--genes-out",
        "6_results_fifth/genes.jsonl",

        "--summary-out",
        "6_results_fifth/summary.csv",

        "--report-out",
        "6_results_fifth/report.json",
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
