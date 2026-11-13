#!/usr/bin/env python3

"""
ORF Analysis — Stage 5
======================

Aggregation of ORF evidence.

Input:
    evidence_records.jsonl from Stage 4

Outputs:
    locus_profiles.csv
    species_profiles.csv
    gene_profiles.jsonl
    gene_summary.csv
    profile_report.json

Principles:
    - Do not classify gene origin.
    - Do not infer ancestral states.
    - Do not count Spaln/GeneWise as independent evidence.
    - Do not select one of multiple loci automatically.
    - Preserve ambiguous and conflicting results.

Python >= 3.10
No third-party dependencies.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys

from collections import Counter, defaultdict
from pathlib import Path


METHODS = ("spaln", "genewise")

EVIDENCE_CLASSES = {
    "preserved_cds_candidate",
    "potential_cds_disruption",
    "conflicting_signals",
    "insufficient_evidence",
}

LOCUS_STATUSES = (
    "preserved_both_identical",
    "preserved_both_different",
    "predicted_interruption_both",
    "discordant_methods",
    "conflicting_alignment",
    "partial_evidence",
    "insufficient_evidence",
    "missing_method",
    "multiple_predictions",
    "inconsistent_metadata",
)

SPECIES_STATUSES = (
    *LOCUS_STATUSES,
    "multiple_loci_uncertain",
)

LOCUS_FIELDS = [
    "gene_id",
    "species",
    "locus",
    "status",
    "reason",
    "tool_agreement",
    "spaln_records",
    "genewise_records",
    "spaln_evidence",
    "genewise_evidence",
    "spaln_cds_status",
    "genewise_cds_status",
    "comparison_records",
    "nontriplet_gap_observations",
    "phase_restored_observations",
]

SPECIES_FIELDS = [
    "gene_id",
    "species",
    "status",
    "locus_count",
    "comparison_records",
    "nontriplet_gap_observations",
    "phase_restored_observations",
    "loci",
]

GENE_FIELDS = [
    "gene_id",
    "species_assessed",
    "loci_assessed",
    "comparison_records",
    "nontriplet_gap_observations",
    "phase_restored_observations",
] + [
    f"species_{status}"
    for status in SPECIES_STATUSES
]


# ============================================================
# 1. INPUT
# ============================================================

def iter_jsonl(path: Path):
    """Read and validate JSONL one line at a time."""

    with Path(path).open(
        "r",
        encoding="utf-8",
    ) as handle:

        for line_number, line in enumerate(
            handle,
            start=1,
        ):

            if not line.strip():
                continue

            try:
                record = json.loads(line)

            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Line {line_number}: invalid JSON"
                ) from error

            if not isinstance(record, dict):
                raise ValueError(
                    f"Line {line_number}: expected object"
                )

            yield record


def compact_record(record: dict) -> dict:
    """
    Extract only the features required for aggregation.

    No DNA sequences or full alignments are retained.
    """

    required = (
        "gene_id",
        "species",
        "method",
        "locus",
    )

    for field in required:

        if not isinstance(record.get(field), str):
            raise ValueError(
                f"Missing string field: {field}"
            )

    if record["method"] not in METHODS:
        raise ValueError(
            f"Unknown method: {record['method']}"
        )

    cds = record.get("cds")
    alignment = record.get("alignment")
    evidence = record.get("evidence")
    gap_path = record.get("gap_path")

    if not all(
        isinstance(item, dict)
        for item in (
            cds,
            alignment,
            evidence,
            gap_path,
        )
    ):
        raise ValueError(
            "Missing Stage-4 analysis objects"
        )

    evidence_class = evidence.get("class")

    if evidence_class not in EVIDENCE_CLASSES:
        raise ValueError(
            f"Unknown evidence class: {evidence_class}"
        )

    gap_count = alignment.get(
        "internal_nontriplet_gap_events"
    )

    if gap_count is None:
        gap_count = 0

    if type(gap_count) is not int or gap_count < 0:
        raise ValueError(
            "Invalid nontriplet gap count"
        )

    alignment_status = alignment.get("status")

    # This is a count of observations from individual
    # alignments, not independent mutations.

    accepted_gap_count = (
        gap_count
        if alignment_status == "aligned"
        else 0
    )

    phase_restored = int(
        alignment_status == "aligned"
        and gap_path.get("status") == "phase_restored"
    )

    return {
        "gene_id": record["gene_id"],
        "species": record["species"],
        "method": record["method"],
        "locus": record["locus"],

        "evidence_class": evidence_class,
        "cds_status": cds.get("status"),

        "profile_match": record.get(
            "stage2_profile_match"
        ),

        "tool_agreement": record.get(
            "tool_agreement"
        ),

        "alignment_status": alignment_status,

        "reported_shift": cds.get(
            "reported_shift"
        ),

        "nontriplet_gap_observations": (
            accepted_gap_count
        ),

        "phase_restored_observations": (
            phase_restored
        ),
    }


# ============================================================
# 2. CONSISTENCY CHECKS
# ============================================================

def valid_method_record(record: dict) -> bool:
    """
    Check whether Stage-4 evidence is consistent
    with the attached CDS and alignment information.
    """

    if record["profile_match"] != "matched":
        return False

    evidence = record["evidence_class"]
    cds = record["cds_status"]
    alignment = record["alignment_status"]

    if evidence == "preserved_cds_candidate":

        return (
            cds == "complete_candidate"
            and alignment == "aligned"
        )

    if evidence == "potential_cds_disruption":

        return (
            cds == "interrupted"
            and alignment == "aligned"
        )

    if evidence == "conflicting_signals":

        return (
            cds == "complete_candidate"
            and alignment == "aligned"
        )

    # Insufficient evidence can have several causes:
    # poor alignment, incomplete CDS, missing data etc.
    return True


# ============================================================
# 3. LOCUS-LEVEL AGGREGATION
# ============================================================

def classify_locus(
    gene_id: str,
    species: str,
    locus: str,
    records: list[dict],
) -> dict:
    """
    Aggregate predictions for one genomic locus.

    Identity:
        gene_id + species + locus

    This is a conservative rule-based aggregation,
    not a statistical model of gene evolution.
    """

    by_method = {
        method: []
        for method in METHODS
    }

    for record in records:
        by_method[record["method"]].append(record)

    spaln = by_method["spaln"]
    genewise = by_method["genewise"]

    agreement_values = {
        record["tool_agreement"]
        for record in records
    }

    agreement = (
        next(iter(agreement_values))
        if len(agreement_values) == 1
        else "inconsistent"
    )

    gap_observations = sum(
        record["nontriplet_gap_observations"]
        for record in records
    )

    phase_observations = sum(
        record["phase_restored_observations"]
        for record in records
    )

    row = {
        "gene_id": gene_id,
        "species": species,
        "locus": locus,
        "status": None,
        "reason": None,
        "tool_agreement": agreement,

        "spaln_records": len(spaln),
        "genewise_records": len(genewise),

        "spaln_evidence": (
            spaln[0]["evidence_class"]
            if len(spaln) == 1
            else None
        ),

        "genewise_evidence": (
            genewise[0]["evidence_class"]
            if len(genewise) == 1
            else None
        ),

        "spaln_cds_status": (
            spaln[0]["cds_status"]
            if len(spaln) == 1
            else None
        ),

        "genewise_cds_status": (
            genewise[0]["cds_status"]
            if len(genewise) == 1
            else None
        ),

        "comparison_records": len(records),

        "nontriplet_gap_observations": (
            gap_observations
        ),

        "phase_restored_observations": (
            phase_observations
        ),
    }

    def finish(status, reason):
        row["status"] = status
        row["reason"] = reason
        return row

    # --------------------------------------------------------
    # Multiple predictions for the same method and locus
    # --------------------------------------------------------

    if len(spaln) > 1 or len(genewise) > 1:

        return finish(
            "multiple_predictions",
            "More than one prediction for the same method",
        )

    # --------------------------------------------------------
    # Missing one of the methods
    # --------------------------------------------------------

    if not spaln or not genewise:

        return finish(
            "missing_method",
            "Only one prediction method is represented",
        )

    a = spaln[0]
    b = genewise[0]

    # --------------------------------------------------------
    # Validate Stage-4 and agreement metadata
    # --------------------------------------------------------

    if not valid_method_record(a):
        return finish(
            "inconsistent_metadata",
            "Inconsistent Spaln evidence",
        )

    if not valid_method_record(b):
        return finish(
            "inconsistent_metadata",
            "Inconsistent GeneWise evidence",
        )

    valid_agreements = {
        "exact_sequence",
        "same_cds_status_only",
        "different_cds_status",
    }

    if agreement not in valid_agreements:

        return finish(
            "inconsistent_metadata",
            "Missing or inconsistent method agreement",
        )

    cds_equal = (
        a["cds_status"] == b["cds_status"]
    )

    if (
        agreement in (
            "exact_sequence",
            "same_cds_status_only",
        )
        and not cds_equal
    ):

        return finish(
            "inconsistent_metadata",
            "Agreement contradicts CDS statuses",
        )

    if (
        agreement == "different_cds_status"
        and cds_equal
    ):

        return finish(
            "inconsistent_metadata",
            "Different-CDS agreement contradicts statuses",
        )

    x = a["evidence_class"]
    y = b["evidence_class"]

    # --------------------------------------------------------
    # Both methods support structurally intact CDS
    # --------------------------------------------------------

    if (
        x == "preserved_cds_candidate"
        and y == "preserved_cds_candidate"
    ):

        if agreement == "exact_sequence":

            return finish(
                "preserved_both_identical",
                "Both methods predict identical intact CDS",
            )

        if agreement == "same_cds_status_only":

            return finish(
                "preserved_both_different",
                "Both predict intact CDS but sequences differ",
            )

        return finish(
            "inconsistent_metadata",
            "Unexpected agreement for preserved CDS",
        )

    # --------------------------------------------------------
    # Both methods predict an interrupted CDS
    # --------------------------------------------------------

    if (
        x == "potential_cds_disruption"
        and y == "potential_cds_disruption"
    ):

        return finish(
            "predicted_interruption_both",
            "Internal stops predicted by both methods",
        )

    # --------------------------------------------------------
    # One method suggests preservation,
    # the other suggests a disrupted CDS.
    # --------------------------------------------------------

    if {
        x, y
    } == {
        "preserved_cds_candidate",
        "potential_cds_disruption",
    }:

        return finish(
            "discordant_methods",
            "Predictions contradict each other",
        )

    # --------------------------------------------------------
    # Frame-like alignment gaps contradict
    # structurally intact CDS predictions.
    # --------------------------------------------------------

    if "conflicting_signals" in (x, y):

        return finish(
            "conflicting_alignment",
            "Nontriplet gap signal conflicts with CDS evidence",
        )

    # --------------------------------------------------------
    # Only one method provides an informative prediction
    # --------------------------------------------------------

    if "insufficient_evidence" in (x, y):

        if (
            "preserved_cds_candidate" in (x, y)
            or "potential_cds_disruption" in (x, y)
        ):

            return finish(
                "partial_evidence",
                "Only one method provides usable evidence",
            )

        return finish(
            "insufficient_evidence",
            "Both methods lack sufficient evidence",
        )

    return finish(
        "insufficient_evidence",
        "No supported aggregation rule",
    )


# ============================================================
# 4. SPECIES-LEVEL AGGREGATION
# ============================================================

def aggregate_species(
    locus_rows: list[dict],
) -> list[dict]:
    """
    Aggregate all observed loci for each gene and species.

    Multiple loci are marked as ambiguous.
    No ortholog is selected automatically.
    """

    groups = defaultdict(list)

    for row in locus_rows:

        key = (
            row["gene_id"],
            row["species"],
        )

        groups[key].append(row)

    species_rows = []

    for key in sorted(groups):

        gene_id, species = key

        loci = groups[key]

        if len(loci) == 1:

            status = loci[0]["status"]

        else:

            status = "multiple_loci_uncertain"

        species_rows.append({
            "gene_id": gene_id,
            "species": species,
            "status": status,

            "locus_count": len(loci),

            "comparison_records": sum(
                row["comparison_records"]
                for row in loci
            ),

            "nontriplet_gap_observations": sum(
                row["nontriplet_gap_observations"]
                for row in loci
            ),

            "phase_restored_observations": sum(
                row["phase_restored_observations"]
                for row in loci
            ),

            "loci": ";".join(
                sorted(
                    row["locus"]
                    for row in loci
                )
            ),
        })

    return species_rows


# ============================================================
# 5. GENE-LEVEL AGGREGATION
# ============================================================

def build_gene_profile(
    gene_id: str,
    species_rows: list[dict],
) -> tuple[dict, dict]:
    """
    Produce a structured gene profile and a CSV summary.

    Important:
        The result is not a gene-origin classification.
    """

    statuses = Counter(
        row["status"]
        for row in species_rows
    )

    loci_total = sum(
        row["locus_count"]
        for row in species_rows
    )

    comparisons_total = sum(
        row["comparison_records"]
        for row in species_rows
    )

    gaps_total = sum(
        row["nontriplet_gap_observations"]
        for row in species_rows
    )

    phase_total = sum(
        row["phase_restored_observations"]
        for row in species_rows
    )

    profile = {
        "gene_id": gene_id,

        "species_assessed": len(species_rows),
        "loci_assessed": loci_total,
        "comparison_records": comparisons_total,

        "species_status_counts": dict(statuses),

        "nontriplet_gap_observations": gaps_total,
        "phase_restored_observations": phase_total,

        "species_profiles": [
            {
                "species": row["species"],
                "status": row["status"],
                "locus_count": row["locus_count"],
                "loci": row["loci"].split(";"),
            }
            for row in species_rows
        ],

        "origin_classification": None,

        "limitations": [
            "Only observed loci are evaluated.",
            "Missing species do not imply gene absence.",
            "Multiple loci are not resolved.",
            "Spaln and GeneWise are correlated predictions.",
            "Gap observations are not verified mutations.",
            "Ancestral coding status is not inferred.",
        ],
    }

    csv_row = {
        "gene_id": gene_id,
        "species_assessed": len(species_rows),
        "loci_assessed": loci_total,
        "comparison_records": comparisons_total,
        "nontriplet_gap_observations": gaps_total,
        "phase_restored_observations": phase_total,
    }

    for status in SPECIES_STATUSES:

        csv_row[f"species_{status}"] = (
            statuses[status]
        )

    return profile, csv_row


# ============================================================
# 6. PIPELINE
# ============================================================

def run(
    input_path: Path,
    locus_out: Path,
    species_out: Path,
    genes_out: Path,
    summary_out: Path,
    report_out: Path,
    gene_filter: str | None = None,
    progress_every: int = 1000,
    verbose: bool = True,
) -> dict:
    """
    Aggregate Stage-4 evidence into locus,
    species and gene profiles.

    Keeps compact metadata in memory.
    Does not retain sequence strings.
    """

    input_path = Path(input_path)

    output_paths = [
        Path(p)
        for p in (
            locus_out,
            species_out,
            genes_out,
            summary_out,
            report_out,
        )
    ]

    paths = [input_path] + output_paths

    if len({
        path.resolve()
        for path in paths
    }) != len(paths):

        raise ValueError(
            "All input/output paths must be different"
        )

    if progress_every < 0:
        raise ValueError(
            "progress_every must be >= 0"
        )

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Входной файл не найден: {input_path.resolve()}"
        )

    def log(message):

        if verbose:
            print(
                f"[ORF Profile] {message}",
                flush=True,
            )

    # --------------------------------------------------------
    # Load compact evidence records
    # --------------------------------------------------------

    log("Stage 5: ORF Profile Aggregation")
    log(f"Input: {input_path.resolve()}")

    groups = defaultdict(list)
    processed = 0

    for record in iter_jsonl(input_path):

        if (
            gene_filter is not None
            and record.get("gene_id") != gene_filter
        ):
            continue

        compact = compact_record(record)

        key = (
            compact["gene_id"],
            compact["species"],
            compact["locus"],
        )

        groups[key].append(compact)

        processed += 1

        if (
            progress_every > 0
            and processed % progress_every == 0
        ):

            log(
                f"Прочитано сравнений: {processed}; "
                f"локусов: {len(groups)}"
            )

    if processed == 0:

        raise ValueError(
            "Не найдено сравнений для обработки. "
            "Проверьте --input и --gene."
        )

    log(f"Всего сравнений: {processed}")
    log("Агрегация на уровне локусов...")

    # --------------------------------------------------------
    # Locus-level profiles
    # --------------------------------------------------------

    locus_rows = []

    for key in sorted(groups):

        gene_id, species, locus = key

        row = classify_locus(
            gene_id,
            species,
            locus,
            groups[key],
        )

        locus_rows.append(row)

    # --------------------------------------------------------
    # Species-level profiles
    # --------------------------------------------------------

    log("Агрегация на уровне видов...")

    species_rows = aggregate_species(
        locus_rows
    )

    # --------------------------------------------------------
    # Gene-level profiles
    # --------------------------------------------------------

    by_gene = defaultdict(list)

    for row in species_rows:
        by_gene[row["gene_id"]].append(row)

    genes = {}

    for gene_id in sorted(by_gene):

        genes[gene_id] = build_gene_profile(
            gene_id,
            by_gene[gene_id],
        )

    # --------------------------------------------------------
    # Output directories
    # --------------------------------------------------------

    for path in output_paths:

        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    # --------------------------------------------------------
    # Locus CSV
    # --------------------------------------------------------

    with output_paths[0].open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=LOCUS_FIELDS,
        )

        writer.writeheader()
        writer.writerows(locus_rows)

    # --------------------------------------------------------
    # Species CSV
    # --------------------------------------------------------

    with output_paths[1].open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=SPECIES_FIELDS,
        )

        writer.writeheader()
        writer.writerows(species_rows)

    # --------------------------------------------------------
    # Gene JSONL and CSV
    # --------------------------------------------------------

    with (
        output_paths[2].open(
            "w",
            encoding="utf-8",
        ) as genes_file,

        output_paths[3].open(
            "w",
            encoding="utf-8",
            newline="",
        ) as summary_file
    ):

        writer = csv.DictWriter(
            summary_file,
            fieldnames=GENE_FIELDS,
        )

        writer.writeheader()

        for gene_id in sorted(genes):

            profile, csv_row = genes[gene_id]

            genes_file.write(
                json.dumps(
                    profile,
                    ensure_ascii=False,
                ) + "\n"
            )

            writer.writerow(csv_row)

    # --------------------------------------------------------
    # Final report
    # --------------------------------------------------------

    locus_counts = Counter(
        row["status"]
        for row in locus_rows
    )

    species_counts = Counter(
        row["status"]
        for row in species_rows
    )

    report = {
        "stage": "ORF Profile Aggregation Stage 5",

        "input_file": str(input_path.resolve()),
        "gene_filter": gene_filter,

        "comparison_records": processed,
        "unique_loci": len(locus_rows),
        "gene_species_pairs": len(species_rows),
        "unique_genes": len(genes),

        "locus_status_counts": dict(locus_counts),
        "species_status_counts": dict(species_counts),

        "nontriplet_gap_observations": sum(
            row["nontriplet_gap_observations"]
            for row in locus_rows
        ),

        "phase_restored_observations": sum(
            row["phase_restored_observations"]
            for row in locus_rows
        ),

        "origin_classification_performed": False,

        "limitations": [
            "Prediction agreement is not independent evidence.",
            "Observed loci are not verified orthologs.",
            "Multiple loci are left unresolved.",
            "No species absence calls are made.",
            "Gap observations are not evolutionary events.",
            "No ancestral reconstruction is performed.",
        ],
    }

    output_paths[4].write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    log("=" * 55)
    log("ORF Profile Aggregation завершён")

    log(f"Генов: {len(genes)}")
    log(f"Видовых профилей: {len(species_rows)}")
    log(f"Локусов: {len(locus_rows)}")

    log("Статусы локусов:")

    for status, count in sorted(locus_counts.items()):
        log(f"  {status}: {count}")

    log("Выходные файлы:")

    for path in output_paths:
        log(str(path.resolve()))

    return report


# ============================================================
# 7. CLI
# ============================================================

def main(argv=None):

    parser = argparse.ArgumentParser(
        description=(
            "ORF Analysis Stage 5: gene-level aggregation"
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="JSONL четвёртого этапа",
    )

    parser.add_argument(
        "--locus-out",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--species-out",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--genes-out",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--summary-out",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--report-out",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--gene",
        default=None,
        help="ID одного гена; иначе все гены",
    )

    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
    )

    args = parser.parse_args(argv)

    try:

        run(
            input_path=args.input,
            locus_out=args.locus_out,
            species_out=args.species_out,
            genes_out=args.genes_out,
            summary_out=args.summary_out,
            report_out=args.report_out,
            gene_filter=args.gene,
            progress_every=args.progress_every,
        )

    except (OSError, ValueError, TypeError) as error:

        print(
            f"[ORF Profile] ОШИБКА: {error}",
            file=sys.stderr,
            flush=True,
        )

        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
