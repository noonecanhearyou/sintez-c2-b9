"""
Tests for the unified ORF pipeline.

Run:
    python -m pytest -v test_run_orf_pipeline.py

No real biological processing is performed.
Temporary scripts simulate the five stages.
"""

import json
import sys

from pathlib import Path

import run_orf_pipeline as pipeline


# ============================================================
# 1. Fake stage script
# ============================================================

FAKE_STAGE = r'''
import json
import sys

from pathlib import Path

script_name = Path(__file__).name
script_dir = Path(__file__).parent

args = sys.argv[1:]

if len(args) % 2 != 0:
    sys.exit(3)

options = {}

for index in range(0, len(args), 2):
    options[args[index]] = args[index + 1]

with (script_dir / "calls.log").open(
    "a",
    encoding="utf-8",
) as log:
    log.write(script_name + "\n")

# Deliberate failure for failure-handling tests.
if (script_dir / "fail_stage.txt").is_file():

    failing = (
        script_dir / "fail_stage.txt"
    ).read_text(encoding="utf-8").strip()

    if failing == script_name:
        sys.exit(17)

for flag, filename in options.items():

    if not flag.endswith("-out"):
        continue

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.suffix == ".jsonl":

        path.write_text(
            json.dumps({"gene_id": "FBgn0000001"})
            + "\n",
            encoding="utf-8",
        )

    elif path.suffix == ".json":

        path.write_text(
            json.dumps({"ok": True}) + "\n",
            encoding="utf-8",
        )

    elif path.suffix == ".csv":

        path.write_text(
            "gene_id\nFBgn0000001\n",
            encoding="utf-8",
        )

print("Fake stage completed:", script_name)
'''


def prepare_environment(tmp_path):

    scripts = tmp_path / "scripts"
    scripts.mkdir()

    for name in pipeline.SCRIPT_NAMES:

        (scripts / name).write_text(
            FAKE_STAGE,
            encoding="utf-8",
        )

    input_file = (
        tmp_path
        / "raw"
        / "input.txt"
    )

    input_file.parent.mkdir()

    input_file.write_text(
        "FAKE INPUT\n",
        encoding="utf-8",
    )

    root = tmp_path / "results"

    final = (
        tmp_path
        / "final"
        / "profiles.jsonl"
    )

    args = [
        "--input", str(input_file),
        "--output-root", str(root),
        "--final-out", str(final),
        "--scripts-dir", str(scripts),
        "--gene", "FBgn0000001",
        "--progress-every", "1",
    ]

    return scripts, input_file, root, final, args


def call_log(scripts):

    path = scripts / "calls.log"

    if not path.exists():
        return []

    return path.read_text(
        encoding="utf-8"
    ).splitlines()


# ============================================================
# 2. First full run
# ============================================================

def test_fresh_pipeline(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    assert call_log(scripts) == list(
        pipeline.SCRIPT_NAMES
    )

    for folder in pipeline.STAGE_DIRS:
        assert (root / folder).is_dir()

    assert final.is_file()

    assert (
        root / "pipeline_manifest.json"
    ).is_file()

    report = json.loads(
        (root / "pipeline_run.json").read_text(
            encoding="utf-8"
        )
    )

    assert report["status"] == "completed"
    assert report["computed_stages"] == 5
    assert report["reused_stages"] == 0


# ============================================================
# 3. Resume without changes
# ============================================================

def test_resume_reuses_all_stages(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    assert pipeline.main(
        args + ["--mode", "resume"]
    ) == 0

    assert len(call_log(scripts)) == 5

    report = json.loads(
        (root / "pipeline_run.json").read_text(
            encoding="utf-8"
        )
    )

    assert report["computed_stages"] == 0
    assert report["reused_stages"] == 5


# ============================================================
# 4. Changed source input
# ============================================================

def test_input_change_forces_recompute(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    source.write_text(
        "CHANGED INPUT\n",
        encoding="utf-8",
    )

    assert pipeline.main(
        args + ["--mode", "resume"]
    ) == 0

    assert len(call_log(scripts)) == 10

    report = json.loads(
        (root / "pipeline_run.json").read_text(
            encoding="utf-8"
        )
    )

    assert report["computed_stages"] == 5


# ============================================================
# 5. Changed comparative script
# ============================================================

def test_changed_script_recomputes_downstream(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    comparative = scripts / "comparative_orf.py"

    comparative.write_text(
        FAKE_STAGE + "\n# Modified version\n",
        encoding="utf-8",
    )

    assert pipeline.main(
        args + ["--mode", "resume"]
    ) == 0

    calls = call_log(scripts)

    assert len(calls) == 8

    assert calls[-3:] == [
        "comparative_orf.py",
        "orf_evidence.py",
        "orf_profile.py",
    ]


# ============================================================
# 6. Changed scanner settings
# ============================================================

def test_parameter_change_recomputes_downstream(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    assert pipeline.main(
        args + [
            "--mode", "resume",
            "--min-aa", "40",
        ]
    ) == 0

    calls = call_log(scripts)

    assert len(calls) == 9

    assert calls[-4:] == [
        "orf_scanner.py",
        "comparative_orf.py",
        "orf_evidence.py",
        "orf_profile.py",
    ]


# ============================================================
# 7. Deleted intermediate output
# ============================================================

def test_missing_output_recomputed(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    missing = (
        root
        / "4_results_third"
        / "comparison_summary.csv"
    )

    missing.unlink()

    assert pipeline.main(
        args + ["--mode", "resume"]
    ) == 0

    calls = call_log(scripts)

    assert len(calls) == 8

    assert calls[-3:] == [
        "comparative_orf.py",
        "orf_evidence.py",
        "orf_profile.py",
    ]


# ============================================================
# 8. Failure propagation
# ============================================================

def test_failed_stage_stops_pipeline(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    (scripts / "fail_stage.txt").write_text(
        "orf_scanner.py",
        encoding="utf-8",
    )

    result = pipeline.main(
        args + ["--mode", "fresh"]
    )

    assert result == 2

    assert call_log(scripts) == [
        "orf_parser.py",
        "orf_scanner.py",
    ]

    assert not final.exists()

    report = json.loads(
        (root / "pipeline_run.json").read_text(
            encoding="utf-8"
        )
    )

    assert report["status"] == "failed"
    assert report["failed_stage"] == "scanner"


# ============================================================
# 9. Final CSV
# ============================================================

def test_final_csv_output(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    final_csv = (
        tmp_path
        / "external"
        / "summary.csv"
    )

    index = args.index("--final-out")

    args[index + 1] = str(final_csv)

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    internal_csv = (
        root
        / "6_results_fifth"
        / "gene_summary.csv"
    )

    assert final_csv.is_file()

    assert final_csv.read_bytes() == (
        internal_csv.read_bytes()
    )


# ============================================================
# 10. Invalid final path
# ============================================================

def test_invalid_final_extension(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    index = args.index("--final-out")

    args[index + 1] = str(
        tmp_path / "final" / "results.txt"
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 2

    assert call_log(scripts) == []


# ============================================================
# 11. Missing source
# ============================================================

def test_missing_input(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    source.unlink()

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 2

    assert call_log(scripts) == []


# ============================================================
# 12. Gene filter changes stage signature
# ============================================================

def test_gene_change_recomputes_all(tmp_path):

    scripts, source, root, final, args = (
        prepare_environment(tmp_path)
    )

    assert pipeline.main(
        args + ["--mode", "fresh"]
    ) == 0

    gene_index = args.index("--gene")

    args[gene_index + 1] = "FBgn0000002"

    assert pipeline.main(
        args + ["--mode", "resume"]
    ) == 0

    assert len(call_log(scripts)) == 10
