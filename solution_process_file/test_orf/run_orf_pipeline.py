#!/usr/bin/env python3

"""
Unified ORF Pipeline — Stage 6

Runs:
    1. orf_parser.py
    2. orf_scanner.py
    3. comparative_orf.py
    4. orf_evidence.py
    5. orf_profile.py

Requirements:
    Python >= 3.10
    Linux (fcntl file locking)

Existing stage scripts require their original dependencies.

Features:
    - One gene or all genes
    - Fresh and resume modes
    - Independent input, output root and final output paths
    - SHA256 verification of input and output files
    - Stage parameter and source-code change detection
    - Progress displayed directly in terminal
    - Failure propagation
    - Atomic writing of manifest and final result
    - One pipeline execution per output directory
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


MANIFEST_VERSION = 1

SCRIPT_NAMES = (
    "orf_parser.py",
    "orf_scanner.py",
    "comparative_orf.py",
    "orf_evidence.py",
    "orf_profile.py",
)

STAGE_DIRS = (
    "1_raw_data",
    "2_results_first",
    "3_results_second",
    "4_results_third",
    "5_results_fourth",
    "6_results_fifth",
)


class PipelineError(RuntimeError):
    """Error in pipeline orchestration."""


@dataclass(frozen=True)
class Stage:
    name: str
    script: Path
    inputs: tuple[Path, ...]
    outputs: tuple[Path, ...]
    command: tuple[str, ...]


# ============================================================
# 1. Logging and timestamps
# ============================================================

def now_utc() -> str:
    return datetime.now(
        timezone.utc
    ).isoformat(timespec="seconds")


def log(message: str):
    print(
        f"[ORF Pipeline] {message}",
        flush=True,
    )


# ============================================================
# 2. File signatures
# ============================================================

def file_signature(path: Path) -> dict:
    """
    Calculate SHA256 without loading the complete file.

    Detects changes even if file size is unchanged.
    """

    path = Path(path).resolve()

    if not path.is_file():
        raise PipelineError(
            f"Отсутствует необходимый файл: {path}"
        )

    digest = hashlib.sha256()

    with path.open("rb") as handle:

        while True:

            chunk = handle.read(1024 * 1024)

            if not chunk:
                break

            digest.update(chunk)

    return {
        "path": str(path),
        "size": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def object_hash(value: dict) -> str:
    """Hash a JSON-serializable object."""

    data = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")

    return hashlib.sha256(data).hexdigest()


# ============================================================
# 3. Atomic JSON writing and final copy
# ============================================================

def atomic_write_json(path: Path, value: dict):
    """
    Save JSON through a temporary file.

    Prevents a partially written manifest.
    """

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    temporary = path.with_name(
        path.name + ".tmp"
    )

    with temporary.open(
        "w",
        encoding="utf-8",
    ) as handle:

        json.dump(
            value,
            handle,
            ensure_ascii=False,
            indent=2,
        )

        handle.write("\n")

    os.replace(temporary, path)


def atomic_copy(source: Path, destination: Path):
    """
    Copy the final result in chunks.

    Replace destination only after copying succeeds.
    """

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary_path = None

    try:

        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=".orf_final_",
            dir=destination.parent,
            delete=False,
        ) as temporary:

            temporary_path = Path(temporary.name)

            with source.open("rb") as input_file:

                shutil.copyfileobj(
                    input_file,
                    temporary,
                    length=1024 * 1024,
                )

            temporary.flush()
            os.fsync(temporary.fileno())

        os.replace(
            temporary_path,
            destination,
        )

    finally:

        if (
            temporary_path is not None
            and temporary_path.exists()
        ):
            temporary_path.unlink()


# ============================================================
# 4. Output validation
# ============================================================

def validate_output(path: Path):
    """
    Basic validation of a produced file.

    This validates the output format, not biological
    correctness of the content.
    """

    if not path.is_file():
        raise PipelineError(
            f"Этап не создал файл: {path}"
        )

    if path.stat().st_size == 0:
        raise PipelineError(
            f"Этап создал пустой файл: {path}"
        )

    if path.suffix == ".json":

        try:

            with path.open(
                "r",
                encoding="utf-8",
            ) as handle:

                data = json.load(handle)

            if not isinstance(data, dict):
                raise ValueError(
                    "Expected JSON object"
                )

        except (ValueError, UnicodeError) as exc:

            raise PipelineError(
                f"Некорректный JSON: {path}: {exc}"
            ) from exc

    elif path.suffix == ".jsonl":

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            first_line = next(
                (
                    line
                    for line in handle
                    if line.strip()
                ),
                None,
            )

        if first_line is None:
            raise PipelineError(
                f"JSONL не содержит записей: {path}"
            )

        try:

            first_record = json.loads(first_line)

            if not isinstance(first_record, dict):
                raise ValueError(
                    "Expected JSON object"
                )

        except ValueError as exc:

            raise PipelineError(
                f"Некорректный JSONL: {path}: {exc}"
            ) from exc

    elif path.suffix == ".csv":

        with path.open(
            "r",
            encoding="utf-8",
            newline="",
        ) as handle:

            reader = csv.reader(handle)
            header = next(reader, None)

        if not header:
            raise PipelineError(
                f"CSV не содержит заголовка: {path}"
            )


# ============================================================
# 5. Exclusive pipeline lock
# ============================================================

@contextmanager
def pipeline_lock(output_root: Path):
    """
    Prevent simultaneous runs into the same directory.

    Linux automatically releases the lock when the
    process terminates.
    """

    lock_path = output_root / ".pipeline.lock"

    with lock_path.open(
        "a+",
        encoding="utf-8",
    ) as handle:

        try:

            fcntl.flock(
                handle.fileno(),
                fcntl.LOCK_EX | fcntl.LOCK_NB,
            )

        except BlockingIOError as exc:

            raise PipelineError(
                "В этом выходном каталоге уже "
                "выполняется другой конвейер."
            ) from exc

        try:
            yield

        finally:

            fcntl.flock(
                handle.fileno(),
                fcntl.LOCK_UN,
            )


# ============================================================
# 6. Stage construction
# ============================================================

def build_stages(args) -> list[Stage]:
    """
    Construct commands and paths for all five stages.

    All paths are absolute.
    """

    root = args.output_root.resolve()
    scripts = args.scripts_dir.resolve()
    source = args.input.resolve()

    d2 = root / "2_results_first"
    d3 = root / "3_results_second"
    d4 = root / "4_results_third"
    d5 = root / "5_results_fourth"
    d6 = root / "6_results_fifth"

    # Stage 1
    p2_records = d2 / "records.jsonl"
    p2_summary = d2 / "summary.csv"
    p2_report = d2 / "validation.json"

    # Stage 2
    p3_records = d3 / "orf_records.jsonl"
    p3_summary = d3 / "orf_summary.csv"
    p3_report = d3 / "orf_report.json"

    # Stage 3
    p4_records = d4 / "comparisons.jsonl"
    p4_summary = d4 / "comparison_summary.csv"
    p4_agreement = d4 / "tool_agreement.csv"
    p4_report = d4 / "comparison_report.json"

    # Stage 4
    p5_records = d5 / "evidence_records.jsonl"
    p5_summary = d5 / "evidence_summary.csv"
    p5_genes = d5 / "gene_summary.csv"
    p5_report = d5 / "evidence_report.json"

    # Stage 5
    p6_locus = d6 / "locus_profiles.csv"
    p6_species = d6 / "species_profiles.csv"
    p6_genes = d6 / "gene_profiles.jsonl"
    p6_summary = d6 / "gene_summary.csv"
    p6_report = d6 / "profile_report.json"

    gene_options = (
        ["--gene", args.gene]
        if args.gene
        else []
    )

    progress = [
        "--progress-every",
        str(args.progress_every),
    ]

    def command(script_name: str, options: list):
        return tuple([
            sys.executable,
            "-u",
            str(scripts / script_name),
            *map(str, options),
            *gene_options,
            *progress,
        ])

    return [
        Stage(
            name="parser",
            script=scripts / "orf_parser.py",
            inputs=(source,),
            outputs=(
                p2_records,
                p2_summary,
                p2_report,
            ),
            command=command(
                "orf_parser.py",
                [
                    "--input", source,
                    "--records-out", p2_records,
                    "--summary-out", p2_summary,
                    "--report-out", p2_report,
                ],
            ),
        ),

        Stage(
            name="scanner",
            script=scripts / "orf_scanner.py",
            inputs=(p2_records,),
            outputs=(
                p3_records,
                p3_summary,
                p3_report,
            ),
            command=command(
                "orf_scanner.py",
                [
                    "--input", p2_records,
                    "--records-out", p3_records,
                    "--summary-out", p3_summary,
                    "--report-out", p3_report,
                    "--min-aa", args.min_aa,
                    "--top-k", args.top_k,
                    "--max-peptide-aa", args.max_peptide_aa,
                ],
            ),
        ),

        Stage(
            name="comparative",
            script=scripts / "comparative_orf.py",
            inputs=(p3_records,),
            outputs=(
                p4_records,
                p4_summary,
                p4_agreement,
                p4_report,
            ),
            command=command(
                "comparative_orf.py",
                [
                    "--input", p3_records,
                    "--comparisons-out", p4_records,
                    "--summary-out", p4_summary,
                    "--agreement-out", p4_agreement,
                    "--report-out", p4_report,
                    "--max-align-length",
                    args.max_align_length,
                    "--max-align-cells",
                    args.max_align_cells,
                    "--min-identity", args.min_identity,
                    "--min-coverage", args.min_coverage,
                    "--max-events", args.max_events,
                    "--preview-columns",
                    args.preview_columns,
                ],
            ),
        ),

        Stage(
            name="evidence",
            script=scripts / "orf_evidence.py",
            inputs=(
                p3_records,
                p4_records,
                p4_agreement,
            ),
            outputs=(
                p5_records,
                p5_summary,
                p5_genes,
                p5_report,
            ),
            command=command(
                "orf_evidence.py",
                [
                    "--orf-input", p3_records,
                    "--comparisons-input", p4_records,
                    "--agreement-input", p4_agreement,
                    "--records-out", p5_records,
                    "--summary-out", p5_summary,
                    "--genes-out", p5_genes,
                    "--report-out", p5_report,
                    "--max-gap-details",
                    args.max_gap_details,
                ],
            ),
        ),

        Stage(
            name="profile",
            script=scripts / "orf_profile.py",
            inputs=(p5_records,),
            outputs=(
                p6_locus,
                p6_species,
                p6_genes,
                p6_summary,
                p6_report,
            ),
            command=command(
                "orf_profile.py",
                [
                    "--input", p5_records,
                    "--locus-out", p6_locus,
                    "--species-out", p6_species,
                    "--genes-out", p6_genes,
                    "--summary-out", p6_summary,
                    "--report-out", p6_report,
                ],
            ),
        ),
    ]


# ============================================================
# 7. Manifest
# ============================================================

def new_manifest() -> dict:
    return {
        "version": MANIFEST_VERSION,
        "stages": {},
    }


def load_manifest(path: Path) -> dict:

    if not path.is_file():
        return new_manifest()

    try:

        with path.open(
            "r",
            encoding="utf-8",
        ) as handle:

            manifest = json.load(handle)

    except (OSError, ValueError) as exc:

        raise PipelineError(
            "Невозможно прочитать pipeline_manifest.json. "
            "Для нового расчёта используйте --mode fresh."
        ) from exc

    if (
        not isinstance(manifest, dict)
        or manifest.get("version") != MANIFEST_VERSION
        or not isinstance(manifest.get("stages"), dict)
    ):

        raise PipelineError(
            "Неподдерживаемый формат манифеста. "
            "Используйте --mode fresh."
        )

    return manifest


def stage_signature(stage: Stage) -> str:
    """
    Stage signature includes:
        - command-line parameters
        - script contents
        - input file contents

    SHA256 is computed in chunks.
    """

    state = {
        "stage": stage.name,
        "command": list(stage.command),
        "script": file_signature(stage.script),
        "inputs": [
            file_signature(path)
            for path in stage.inputs
        ],
    }

    return object_hash(state)


def may_reuse(
    stage: Stage,
    previous: dict | None,
    signature: str,
) -> bool:

    if not isinstance(previous, dict):
        return False

    if previous.get("status") != "completed":
        return False

    if previous.get("signature") != signature:
        return False

    saved_outputs = previous.get("outputs")

    if not isinstance(saved_outputs, dict):
        return False

    for path in stage.outputs:

        key = str(path.resolve())

        if key not in saved_outputs:
            return False

        if not path.is_file():
            return False

        if file_signature(path) != saved_outputs[key]:
            return False

        try:
            validate_output(path)

        except (PipelineError, OSError, ValueError):
            return False

    return True


# ============================================================
# 8. Configuration validation
# ============================================================

def validate_configuration(args, stages):

    if not args.input.is_file():
        raise PipelineError(
            f"Исходный файл не найден: "
            f"{args.input.resolve()}"
        )

    for script_name in SCRIPT_NAMES:

        script = args.scripts_dir / script_name

        if not script.is_file():
            raise PipelineError(
                f"Не найден скрипт: {script.resolve()}"
            )

    if args.final_out.suffix.lower() not in (
        ".jsonl",
        ".csv",
    ):
        raise PipelineError(
            "--final-out должен оканчиваться "
            "на .jsonl или .csv"
        )

    numeric_checks = (
        ("min-aa", args.min_aa, 1),
        ("top-k", args.top_k, 0),
        ("max-peptide-aa", args.max_peptide_aa, 0),
        ("max-align-length", args.max_align_length, 1),
        ("max-align-cells", args.max_align_cells, 1),
        ("max-events", args.max_events, 0),
        ("preview-columns", args.preview_columns, 0),
        ("max-gap-details", args.max_gap_details, 0),
        ("progress-every", args.progress_every, 0),
    )

    for name, value, minimum in numeric_checks:

        if value < minimum:
            raise PipelineError(
                f"--{name} должен быть >= {minimum}"
            )

    for name, value in (
        ("min-identity", args.min_identity),
        ("min-coverage", args.min_coverage),
    ):

        if not 0 <= value <= 1:
            raise PipelineError(
                f"--{name} должен быть в диапазоне [0, 1]"
            )

    root = args.output_root.resolve()

    reserved_paths = {
        args.input.resolve(),
        root / "pipeline_manifest.json",
        root / "pipeline_run.json",
        root / "1_raw_data" / "source_reference.json",
    }

    for stage in stages:

        reserved_paths.add(stage.script.resolve())

        for path in stage.outputs:
            reserved_paths.add(path.resolve())

    if args.final_out.resolve() in reserved_paths:

        raise PipelineError(
            "--final-out совпадает с входным, служебным "
            "или промежуточным файлом"
        )


# ============================================================
# 9. Pipeline execution
# ============================================================

def execute(args) -> dict:

    args.input = args.input.expanduser().resolve()

    args.output_root = (
        args.output_root.expanduser().resolve()
    )

    args.final_out = (
        args.final_out.expanduser().resolve()
    )

    args.scripts_dir = (
        args.scripts_dir.expanduser().resolve()
    )

    stages = build_stages(args)

    validate_configuration(args, stages)

    root = args.output_root

    for directory in STAGE_DIRS:
        (root / directory).mkdir(
            parents=True,
            exist_ok=True,
        )

    manifest_path = root / "pipeline_manifest.json"
    run_report_path = root / "pipeline_run.json"

    source_reference_path = (
        root / "1_raw_data" / "source_reference.json"
    )

    with pipeline_lock(root):

        if args.mode == "resume":

            manifest = load_manifest(manifest_path)

        else:

            manifest = new_manifest()

            # In fresh mode, never reuse previous
            # stage completion metadata.
            atomic_write_json(
                manifest_path,
                manifest,
            )

        source_signature = file_signature(args.input)

        atomic_write_json(
            source_reference_path,
            {
                "input": source_signature,
                "note": (
                    "The original input file is not copied."
                ),
            },
        )

        started_at = now_utc()

        run_report = {
            "status": "running",
            "started_at": started_at,
            "finished_at": None,
            "mode": args.mode,
            "gene": args.gene,
            "input": str(args.input),
            "output_root": str(root),
            "final_output": str(args.final_out),
            "stages": [],
        }

        atomic_write_json(
            run_report_path,
            run_report,
        )

        current_stage = None

        try:

            log("=" * 60)
            log("Запуск единого ORF-конвейера")
            log(f"Режим: {args.mode}")
            log(f"Исходный файл: {args.input}")
            log(f"Рабочий каталог: {root}")
            log(f"Итоговый файл: {args.final_out}")
            log(f"Ген: {args.gene or 'все гены'}")
            log("=" * 60)

            for stage_number, stage in enumerate(
                stages,
                start=1,
            ):

                current_stage = stage.name

                log(
                    f"Этап {stage_number}/{len(stages)}: "
                    f"{stage.name}"
                )

                signature = stage_signature(stage)

                previous = manifest["stages"].get(
                    stage.name
                )

                if (
                    args.mode == "resume"
                    and may_reuse(
                        stage,
                        previous,
                        signature,
                    )
                ):

                    log(
                        "Результаты проверены — "
                        "этап пропущен."
                    )

                    run_report["stages"].append({
                        "stage": stage.name,
                        "action": "reused",
                        "return_code": 0,
                    })

                    atomic_write_json(
                        run_report_path,
                        run_report,
                    )

                    continue

                # A recomputed stage invalidates its
                # own manifest entry and later entries.

                for later_stage in stages[
                    stage_number - 1:
                ]:

                    manifest["stages"].pop(
                        later_stage.name,
                        None,
                    )

                atomic_write_json(
                    manifest_path,
                    manifest,
                )

                log("Выполнение...")

                log(
                    "Команда: "
                    + " ".join(stage.command)
                )

                completed = subprocess.run(
                    list(stage.command),
                    cwd=args.scripts_dir,
                    check=False,
                )

                if completed.returncode != 0:

                    raise PipelineError(
                        f"Этап '{stage.name}' завершился "
                        f"с кодом {completed.returncode}"
                    )

                for path in stage.outputs:
                    validate_output(path)

                output_signatures = {
                    str(path.resolve()): file_signature(path)
                    for path in stage.outputs
                }

                manifest["stages"][stage.name] = {
                    "status": "completed",
                    "signature": signature,
                    "completed_at": now_utc(),
                    "outputs": output_signatures,
                }

                atomic_write_json(
                    manifest_path,
                    manifest,
                )

                run_report["stages"].append({
                    "stage": stage.name,
                    "action": "computed",
                    "return_code": 0,
                })

                atomic_write_json(
                    run_report_path,
                    run_report,
                )

                log(f"Этап '{stage.name}' завершён.")

            # ------------------------------------------------
            # Final output
            # ------------------------------------------------

            if args.final_out.suffix.lower() == ".jsonl":

                final_source = (
                    root
                    / "6_results_fifth"
                    / "gene_profiles.jsonl"
                )

            else:

                final_source = (
                    root
                    / "6_results_fifth"
                    / "gene_summary.csv"
                )

            validate_output(final_source)

            log("Копирование итогового результата...")

            atomic_copy(
                final_source,
                args.final_out,
            )

            run_report["status"] = "completed"
            run_report["finished_at"] = now_utc()

            run_report["final_source"] = str(
                final_source
            )

            run_report["final_signature"] = (
                file_signature(args.final_out)
            )

            run_report["computed_stages"] = sum(
                item["action"] == "computed"
                for item in run_report["stages"]
            )

            run_report["reused_stages"] = sum(
                item["action"] == "reused"
                for item in run_report["stages"]
            )

            atomic_write_json(
                run_report_path,
                run_report,
            )

            log("=" * 60)
            log("Конвейер успешно завершён.")
            log(
                f"Пересчитано этапов: "
                f"{run_report['computed_stages']}"
            )
            log(
                f"Использовано готовых этапов: "
                f"{run_report['reused_stages']}"
            )
            log(
                f"Итоговый результат: {args.final_out}"
            )
            log(
                f"Отчёт запуска: {run_report_path}"
            )
            log("=" * 60)

            return run_report

        except BaseException as exc:

            run_report["status"] = "failed"
            run_report["finished_at"] = now_utc()
            run_report["failed_stage"] = current_stage
            run_report["error"] = str(exc)

            atomic_write_json(
                run_report_path,
                run_report,
            )

            raise


# ============================================================
# 10. CLI
# ============================================================

def main(argv=None) -> int:

    parser = argparse.ArgumentParser(
        description=(
            "Единый запуск конвейера анализа ORF"
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Исходный TXT с тремя секциями",
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Корневая папка результатов всех этапов",
    )

    parser.add_argument(
        "--final-out",
        type=Path,
        required=True,
        help="Отдельный итоговый файл .jsonl или .csv",
    )

    parser.add_argument(
        "--scripts-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help=(
            "Каталог пяти Python-модулей; "
            "по умолчанию каталог управляющего скрипта"
        ),
    )

    parser.add_argument(
        "--mode",
        choices=("fresh", "resume"),
        default="fresh",
        help="Полный пересчёт или продолжение",
    )

    parser.add_argument(
        "--gene",
        default=None,
        help="ID одного гена; иначе весь набор",
    )

    # Stage 2 parameters

    parser.add_argument(
        "--min-aa",
        type=int,
        default=25,
    )

    parser.add_argument(
        "--top-k",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--max-peptide-aa",
        type=int,
        default=120,
    )

    # Stage 3 parameters

    parser.add_argument(
        "--max-align-length",
        type=int,
        default=4000,
    )

    parser.add_argument(
        "--max-align-cells",
        type=int,
        default=4_000_000,
    )

    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.60,
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.60,
    )

    parser.add_argument(
        "--max-events",
        type=int,
        default=200,
        help=(
            "Количество сохраняемых примеров "
            "замен и разрывов на выравнивание"
        ),
    )

    parser.add_argument(
        "--preview-columns",
        type=int,
        default=90,
    )

    # Stage 4 parameters

    parser.add_argument(
        "--max-gap-details",
        type=int,
        default=10,
    )

    # All stages

    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
    )

    args = parser.parse_args(argv)

    try:

        execute(args)

    except KeyboardInterrupt:

        print(
            "\n[ORF Pipeline] Прервано пользователем.",
            file=sys.stderr,
            flush=True,
        )

        return 130

    except (
        OSError,
        ValueError,
        PipelineError,
    ) as exc:

        print(
            f"[ORF Pipeline] ОШИБКА: {exc}",
            file=sys.stderr,
            flush=True,
        )

        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
