#!/usr/bin/env python3

"""
ORF Analysis — Stage 4
Evidence assessment from stages 2 and 3.

Python 3.12; standard library only.

Inputs:
    Stage 2: orf_records.jsonl
    Stage 3: comparisons.jsonl
    Stage 3: tool_agreement.csv

Outputs:
    evidence_records.jsonl
    evidence_summary.csv
    gene_summary.csv
    evidence_report.json

No evolutionary mutation direction is inferred.
No gene-origin classification is performed.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys

from collections import Counter, defaultdict
from pathlib import Path


METHODS = {"spaln", "genewise"}

SUMMARY_FIELDS = [
    "gene_id",
    "species",
    "method",
    "locus",
    "profile_match",
    "cds_status",
    "internal_stop_count",
    "longest_complete_aa",
    "reported_shift",
    "tool_agreement",
    "alignment_status",
    "identity",
    "reference_coverage",
    "query_coverage",
    "internal_nontriplet_gap_events",
    "gap_path_status",
    "phase_final_mod3",
    "phase_restorations",
    "codon_boundary_gap_events",
    "off_boundary_gap_events",
    "evidence_class",
    "reason",
]

GENE_FIELDS = [
    "gene_id",
    "comparisons",
    "unique_species",
    "unique_loci",
    "preserved_cds_candidate",
    "potential_cds_disruption",
    "conflicting_signals",
    "insufficient_evidence",
    "exact_agreement_loci",
    "nontriplet_event_observations",
]


# ------------------------------------------------------------
# 1. Streaming JSONL input
# ------------------------------------------------------------

def iter_jsonl(path: Path):
    """Read JSONL one record at a time."""

    with Path(path).open("r", encoding="utf-8") as handle:

        for line_no, line in enumerate(handle, 1):

            if not line.strip():
                continue

            try:
                obj = json.loads(line)

            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"{path}, line {line_no}: invalid JSON"
                ) from exc

            if not isinstance(obj, dict):
                raise ValueError(
                    f"{path}, line {line_no}: expected object"
                )

            yield obj


def require_strings(obj: dict, fields, source: str):
    """Validate record identity fields."""

    for field in fields:

        if not isinstance(obj.get(field), str):
            raise ValueError(
                f"{source}: missing string field '{field}'"
            )


# ------------------------------------------------------------
# 2. Compact metadata from Stage 2
# ------------------------------------------------------------

def load_orf_profiles(
    path: Path,
    gene_filter: str | None = None,
):
    """
    Retain only compact CDS/ORF features.

    Full DNA sequences are discarded after each record.
    Duplicate keys are marked ambiguous.
    """

    index = {}
    total_records = 0

    for record in iter_jsonl(path):

        require_strings(
            record,
            ("gene_id", "species", "section",
             "locus", "sequence"),
            "Stage 2",
        )

        if record["section"] not in METHODS:
            continue

        if gene_filter and record["gene_id"] != gene_filter:
            continue

        total_records += 1

        key = (
            record["gene_id"],
            record["species"],
            record["section"],
            record["locus"],
        )

        cds = record.get("cds_validation") or {}
        scan = record.get("orf_scan") or {}
        metadata = record.get("source_metadata") or {}

        if not all(
            isinstance(value, dict)
            for value in (cds, scan, metadata)
        ):
            raise ValueError(
                "Stage 2: invalid CDS/ORF metadata"
            )

        compact = {
            "length_nt": len(record["sequence"]),
            "cds_status": cds.get("status"),
            "internal_stop_count": cds.get(
                "internal_stop_count"
            ),
            "longest_complete_aa": scan.get(
                "longest_complete_aa"
            ),
            "reported_shift": metadata.get(
                "shift_label"
            ),
        }

        if key not in index:

            index[key] = {
                "count": 1,
                "profile": compact,
            }

        else:

            index[key]["count"] += 1

            # We must not arbitrarily select one
            # of multiple matching predictions.
            index[key]["profile"] = None

    return index, total_records


# ------------------------------------------------------------
# 3. Spaln/GeneWise agreement
# ------------------------------------------------------------

def load_agreements(
    path: Path,
    gene_filter: str | None = None,
):
    """
    Read compact agreement metadata.

    Keys:
        gene_id + species + locus

    Duplicate rows are not silently resolved.
    """

    index = {}

    with Path(path).open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:

        reader = csv.DictReader(handle)

        required = {
            "gene_id",
            "species",
            "locus",
            "agreement_type",
        }

        if not required.issubset(
            set(reader.fieldnames or [])
        ):
            raise ValueError(
                "Agreement CSV: missing required columns"
            )

        for row in reader:

            if gene_filter and row["gene_id"] != gene_filter:
                continue

            key = (
                row["gene_id"],
                row["species"],
                row["locus"],
            )

            if key in index:
                index[key] = "duplicate_agreement_rows"
            else:
                index[key] = row["agreement_type"]

    return index


# ------------------------------------------------------------
# 4. Join Stage 2 and Stage 3
# ------------------------------------------------------------

def match_profile(
    profiles: dict,
    comparison_record: dict,
):
    """
    Match exactly one Stage-2 prediction.

    Check length, CDS status and shift label.
    Do not infer that records are identical merely
    because their gene and species match.
    """

    key = (
        comparison_record["gene_id"],
        comparison_record["species"],
        comparison_record["method"],
        comparison_record["locus"],
    )

    entry = profiles.get(key)

    if entry is None:
        return "missing_stage2_profile", None

    if entry["count"] != 1:
        return "ambiguous_multiple_profiles", None

    profile = entry["profile"]
    target = comparison_record["target"]

    if profile["length_nt"] != target.get("length_nt"):
        return "length_mismatch", profile

    if profile["cds_status"] != target.get("cds_status"):
        return "cds_status_mismatch", profile

    if profile["reported_shift"] != target.get(
        "reported_shift"
    ):
        return "shift_label_mismatch", profile

    return "matched", profile


# ------------------------------------------------------------
# 5. Gap path and relative reading-frame offset
# ------------------------------------------------------------

def analyze_gap_path(
    comparison: dict,
    max_details: int = 10,
) -> dict:
    """
    Analyze the gap path from Stage 3.

    Each internal gap changes the relative nucleotide
    offset between aligned CDS sequences.

        gap_in_reference:  +gap_length
        gap_in_query:      -gap_length

    The offset is tracked modulo 3.

    This is NOT a codon-aware alignment and NOT
    a reconstruction of ancestral frameshifts.

    Complete gap examples are required to infer
    whether offset returns to zero.
    """

    result = {
        "status": "not_available",
        "events_complete": False,
        "phase_final_mod3": None,
        "phase_restorations": None,
        "phase_excursions": None,
        "codon_boundary_gap_events": None,
        "off_boundary_gap_events": None,
        "internal_nontriplet_observed": None,
        "details": [],
    }

    if comparison.get("status") not in (
        "aligned",
        "low_quality_alignment",
    ):
        return result

    events = comparison.get("gap_examples")
    total_events = comparison.get("gap_event_count")

    if (
        not isinstance(events, list)
        or not isinstance(total_events, int)
        or total_events < 0
    ):
        result["status"] = "missing_gap_data"
        return result

    # Stage 3 normally stores at most 25 examples.
    # Never extrapolate phase from a truncated list.

    if (
        comparison.get("gap_examples_truncated") is True
        or len(events) != total_events
    ):
        result["status"] = "incomplete_gap_examples"
        return result

    phase = 0
    restorations = 0
    excursions = 0

    boundary_count = 0
    off_boundary_count = 0

    internal_nontriplet = 0
    internal_triplet = 0
    terminal_count = 0

    details = []

    for event in events:

        if not isinstance(event, dict):
            result["status"] = "invalid_gap_data"
            return result

        kind = event.get("kind")
        length = event.get("length_nt")
        internal = event.get("internal")

        if (
            kind not in ("gap_in_query", "gap_in_reference")
            or not isinstance(length, int)
            or length <= 0
            or not isinstance(internal, bool)
        ):
            result["status"] = "invalid_gap_data"
            return result

        if not internal:
            terminal_count += 1
            continue

        # Codon-boundary checks use the coordinates
        # of the two CDS sequences, not genomic
        # coordinates and not inferred ancestral CDS.

        if kind == "gap_in_query":

            ref_position = event.get("reference_start")
            query_position = event.get("query_position")
            delta = -length

        else:

            ref_position = event.get("reference_position")
            query_position = event.get("query_start")
            delta = length

        if (
            not isinstance(ref_position, int)
            or not isinstance(query_position, int)
            or ref_position < 0
            or query_position < 0
        ):
            result["status"] = "invalid_gap_data"
            return result

        at_boundary = (
            ref_position % 3 == 0
            and query_position % 3 == 0
        )

        if at_boundary:
            boundary_count += 1
        else:
            off_boundary_count += 1

        if length % 3:
            internal_nontriplet += 1
        else:
            internal_triplet += 1

        before = phase
        phase = (phase + delta) % 3

        if before == 0 and phase != 0:
            excursions += 1

        if before != 0 and phase == 0:
            restorations += 1

        if len(details) < max_details:

            details.append({
                "kind": kind,
                "length_nt": length,
                "reference_position": ref_position,
                "query_position": query_position,
                "at_codon_boundary": at_boundary,
                "phase_before": before,
                "phase_after": phase,
            })

    # Check consistency against complete event counts
    # emitted by Stage 3.

    expected = (
        comparison.get("internal_nontriplet_gap_events"),
        comparison.get("internal_triplet_gap_events"),
        comparison.get("terminal_gap_events"),
    )

    observed = (
        internal_nontriplet,
        internal_triplet,
        terminal_count,
    )

    if expected != observed:
        result["status"] = "inconsistent_gap_counts"
        return result

    if total_events == 0:

        status = "no_gaps"

    elif internal_nontriplet == 0:

        status = "in_frame_gaps_only"

    elif phase == 0:

        status = "phase_restored"

    else:

        status = "phase_unresolved"

    result.update({
        "status": status,
        "events_complete": True,
        "phase_final_mod3": phase,
        "phase_restorations": restorations,
        "phase_excursions": excursions,
        "codon_boundary_gap_events": boundary_count,
        "off_boundary_gap_events": off_boundary_count,
        "internal_nontriplet_observed": (
            internal_nontriplet
        ),
        "details": details,
    })

    return result


# ------------------------------------------------------------
# 6. Evidence assessment
# ------------------------------------------------------------

def assess_evidence(
    profile_match: str,
    target_status: str | None,
    comparison: dict,
    gap_path: dict,
):
    """
    Conservative, rule-based evidence assessment.

    Returns:
        (evidence_class, reason)
    """

    if profile_match != "matched":

        return (
            "insufficient_evidence",
            "stage2_stage3_data_mismatch",
        )

    if (
        comparison.get("status") != "aligned"
        or comparison.get("quality_ok") is not True
    ):

        return (
            "insufficient_evidence",
            "alignment_did_not_pass_quality_thresholds",
        )

    if target_status == "interrupted":

        return (
            "potential_cds_disruption",
            "internal_stop_in_predicted_cds",
        )

    if target_status == "complete_candidate":

        if gap_path["status"] in (
            "no_gaps",
            "in_frame_gaps_only",
        ):

            return (
                "preserved_cds_candidate",
                "intact_cds_without_nontriplet_internal_gaps",
            )

        if gap_path["status"] in (
            "phase_restored",
            "phase_unresolved",
        ):

            return (
                "conflicting_signals",
                "intact_cds_but_nontriplet_alignment_gaps",
            )

        return (
            "insufficient_evidence",
            "gap_path_not_fully_evaluable",
        )

    return (
        "insufficient_evidence",
        f"target_cds_status_{target_status}",
    )


# ------------------------------------------------------------
# 7. Main pipeline
# ------------------------------------------------------------

def run(
    orf_input: Path,
    comparisons_input: Path,
    agreement_input: Path,
    records_out: Path,
    summary_out: Path,
    genes_out: Path,
    report_out: Path,
    gene_filter: str | None = None,
    max_gap_details: int = 10,
    progress_every: int = 200,
    verbose: bool = True,
) -> dict:

    paths = [
        Path(p) for p in (
            orf_input,
            comparisons_input,
            agreement_input,
            records_out,
            summary_out,
            genes_out,
            report_out,
        )
    ]

    if len({p.resolve() for p in paths}) != len(paths):
        raise ValueError(
            "Input and output paths must be distinct"
        )

    if max_gap_details < 0 or progress_every < 0:
        raise ValueError(
            "max_gap_details and progress_every must be >= 0"
        )

    for path in paths[:3]:

        if not path.is_file():
            raise FileNotFoundError(
                f"Входной файл не найден: {path.resolve()}"
            )

    for path in paths[3:]:
        path.parent.mkdir(parents=True, exist_ok=True)

    def log(message):

        if verbose:
            print(
                f"[ORF Evidence] {message}",
                flush=True,
            )

    log("Stage 4: ORF Evidence Analysis")
    log("Загрузка коротких характеристик ORF...")

    profiles, profile_records = load_orf_profiles(
        paths[0], gene_filter
    )

    agreements = load_agreements(
        paths[2], gene_filter
    )

    log(f"Профилей этапа 2: {profile_records}")
    log(f"Сопоставлений методов: {len(agreements)}")
    log("Обработка сравнений этапа 3...")

    processed = 0

    evidence_counts = Counter()
    alignment_counts = Counter()
    gap_status_counts = Counter()
    profile_match_counts = Counter()

    # These counts are observations across methods,
    # not independent evolutionary events.
    nontriplet_observations = 0

    genes = defaultdict(lambda: {
        "comparisons": 0,
        "species": set(),
        "loci": set(),
        "classes": Counter(),
        "exact_loci": set(),
        "nontriplet_observations": 0,
    })

    with (
        paths[4].open(
            "w", encoding="utf-8", newline=""
        ) as summary_file,

        paths[3].open(
            "w", encoding="utf-8"
        ) as records_file
    ):

        writer = csv.DictWriter(
            summary_file,
            fieldnames=SUMMARY_FIELDS,
        )

        writer.writeheader()

        for record in iter_jsonl(paths[1]):

            require_strings(
                record,
                ("gene_id", "species", "method", "locus"),
                "Stage 3",
            )

            if gene_filter and record["gene_id"] != gene_filter:
                continue

            if record["method"] not in METHODS:
                raise ValueError(
                    "Stage 3: unknown prediction method"
                )

            comparison = record.get("comparison")
            target = record.get("target")

            if not isinstance(comparison, dict):
                raise ValueError(
                    "Stage 3: comparison must be an object"
                )

            if not isinstance(target, dict):
                raise ValueError(
                    "Stage 3: target must be an object"
                )

            gene_id = record["gene_id"]
            species = record["species"]
            locus = record["locus"]
            method = record["method"]

            profile_match, profile = match_profile(
                profiles,
                record,
            )

            agreement_type = agreements.get(
                (gene_id, species, locus),
                "missing_agreement",
            )

            gap_path = analyze_gap_path(
                comparison,
                max_details=max_gap_details,
            )

            evidence_class, reason = assess_evidence(
                profile_match,
                target.get("cds_status"),
                comparison,
                gap_path,
            )

            # ------------------------------------------------
            # Compact detailed JSONL
            # ------------------------------------------------

            result = {
                "gene_id": gene_id,
                "species": species,
                "method": method,
                "locus": locus,

                "stage2_profile_match": profile_match,
                "cds": {
                    "status": target.get("cds_status"),
                    "internal_stop_count": (
                        profile.get("internal_stop_count")
                        if profile else None
                    ),
                    "longest_complete_aa": (
                        profile.get("longest_complete_aa")
                        if profile else None
                    ),
                    "reported_shift": target.get(
                        "reported_shift"
                    ),
                },

                "tool_agreement": agreement_type,

                "alignment": {
                    "status": comparison.get("status"),
                    "identity": comparison.get("identity"),
                    "reference_coverage": comparison.get(
                        "reference_coverage"
                    ),
                    "query_coverage": comparison.get(
                        "query_coverage"
                    ),
                    "internal_nontriplet_gap_events": (
                        comparison.get(
                            "internal_nontriplet_gap_events"
                        )
                    ),
                },

                "gap_path": gap_path,

                "evidence": {
                    "class": evidence_class,
                    "reason": reason,
                    "evolutionary_direction_inferred": False,
                    "confirmed_frameshift": False,
                },
            }

            records_file.write(
                json.dumps(
                    result,
                    ensure_ascii=False,
                ) + "\n"
            )

            # ------------------------------------------------
            # CSV
            # ------------------------------------------------

            row = {
                "gene_id": gene_id,
                "species": species,
                "method": method,
                "locus": locus,

                "profile_match": profile_match,

                "cds_status": target.get("cds_status"),

                "internal_stop_count": (
                    profile.get("internal_stop_count")
                    if profile else None
                ),

                "longest_complete_aa": (
                    profile.get("longest_complete_aa")
                    if profile else None
                ),

                "reported_shift": target.get(
                    "reported_shift"
                ),

                "tool_agreement": agreement_type,

                "alignment_status": comparison.get("status"),
                "identity": comparison.get("identity"),

                "reference_coverage": comparison.get(
                    "reference_coverage"
                ),

                "query_coverage": comparison.get(
                    "query_coverage"
                ),

                "internal_nontriplet_gap_events": (
                    comparison.get(
                        "internal_nontriplet_gap_events"
                    )
                ),

                "gap_path_status": gap_path["status"],

                "phase_final_mod3": gap_path[
                    "phase_final_mod3"
                ],

                "phase_restorations": gap_path[
                    "phase_restorations"
                ],

                "codon_boundary_gap_events": gap_path[
                    "codon_boundary_gap_events"
                ],

                "off_boundary_gap_events": gap_path[
                    "off_boundary_gap_events"
                ],

                "evidence_class": evidence_class,
                "reason": reason,
            }

            writer.writerow(row)

            # ------------------------------------------------
            # Global and gene-level statistics
            # ------------------------------------------------

            processed += 1

            evidence_counts[evidence_class] += 1
            alignment_counts[
                comparison.get("status", "unknown")
            ] += 1

            gap_status_counts[gap_path["status"]] += 1
            profile_match_counts[profile_match] += 1

            gap_count = comparison.get(
                "internal_nontriplet_gap_events"
            )

            # Avoid interpreting missing values as zero.
            if not isinstance(gap_count, int):
                gap_count = 0

            if comparison.get("status") == "aligned":
                nontriplet_observations += gap_count
            else:
                gap_count = 0

            gene = genes[gene_id]

            gene["comparisons"] += 1
            gene["species"].add(species)
            gene["loci"].add((species, locus))
            gene["classes"][evidence_class] += 1

            gene["nontriplet_observations"] += gap_count

            if agreement_type == "exact_sequence":
                gene["exact_loci"].add((species, locus))

            if (
                progress_every > 0
                and processed % progress_every == 0
            ):
                log(
                    f"Сравнений: {processed}; "
                    f"генов: {len(genes)}"
                )

    if processed == 0:
        raise ValueError(
            "Нет сравнений для анализа. "
            "Проверьте --gene и входной JSONL."
        )

    # --------------------------------------------------------
    # Gene summary
    # --------------------------------------------------------

    with paths[5].open(
        "w", encoding="utf-8", newline=""
    ) as gene_file:

        writer = csv.DictWriter(
            gene_file,
            fieldnames=GENE_FIELDS,
        )

        writer.writeheader()

        for gene_id in sorted(genes):

            gene = genes[gene_id]
            counts = gene["classes"]

            writer.writerow({
                "gene_id": gene_id,
                "comparisons": gene["comparisons"],
                "unique_species": len(gene["species"]),
                "unique_loci": len(gene["loci"]),

                "preserved_cds_candidate": counts[
                    "preserved_cds_candidate"
                ],

                "potential_cds_disruption": counts[
                    "potential_cds_disruption"
                ],

                "conflicting_signals": counts[
                    "conflicting_signals"
                ],

                "insufficient_evidence": counts[
                    "insufficient_evidence"
                ],

                "exact_agreement_loci": len(
                    gene["exact_loci"]
                ),

                "nontriplet_event_observations": gene[
                    "nontriplet_observations"
                ],
            })

    # --------------------------------------------------------
    # Validation/report
    # --------------------------------------------------------

    report = {
        "stage": "ORF Evidence Analysis Stage 4",
        "orf_input": str(paths[0].resolve()),
        "comparisons_input": str(paths[1].resolve()),
        "agreement_input": str(paths[2].resolve()),
        "gene_filter": gene_filter,

        "processed_comparisons": processed,
        "unique_genes": len(genes),

        "evidence_class_counts": dict(evidence_counts),
        "alignment_status_counts": dict(alignment_counts),
        "gap_path_status_counts": dict(gap_status_counts),
        "profile_match_counts": dict(profile_match_counts),

        "nontriplet_event_observations": (
            nontriplet_observations
        ),

        "settings": {
            "max_gap_details": max_gap_details,
        },

        "limitations": [
            "Stage 3 uses ordinary nucleotide alignment.",
            "Gap coordinates are not genomic mutation calls.",
            "Truncated gap lists do not support phase reconstruction.",
            "The phase model is relative to the selected alignment.",
            "Spaln and GeneWise are not independent experiments.",
            "Event observations may be duplicated across methods.",
            "No ancestral state or gene origin is inferred.",
        ],
    }

    paths[6].write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    log("=" * 55)
    log("ORF Evidence Analysis завершён")
    log(f"Сравнений: {processed}")
    log(f"Генов: {len(genes)}")

    for name, count in sorted(evidence_counts.items()):
        log(f"  {name}: {count}")

    log("Результаты:")
    for path in paths[3:]:
        log(str(path.resolve()))

    return report


# ------------------------------------------------------------
# 8. CLI
# ------------------------------------------------------------

def main(argv=None):

    parser = argparse.ArgumentParser(
        description="ORF Evidence Analysis — Stage 4"
    )

    parser.add_argument(
        "--orf-input",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--comparisons-input",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--agreement-input",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--records-out",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--summary-out",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--genes-out",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--report-out",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--gene",
        default=None,
    )

    parser.add_argument(
        "--max-gap-details",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
    )

    args = parser.parse_args(argv)

    try:

        run(
            orf_input=args.orf_input,
            comparisons_input=args.comparisons_input,
            agreement_input=args.agreement_input,
            records_out=args.records_out,
            summary_out=args.summary_out,
            genes_out=args.genes_out,
            report_out=args.report_out,
            gene_filter=args.gene,
            max_gap_details=args.max_gap_details,
            progress_every=args.progress_every,
        )

    except (OSError, ValueError, TypeError) as exc:

        print(
            f"[ORF Evidence] ОШИБКА: {exc}",
            file=sys.stderr,
            flush=True,
        )

        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
