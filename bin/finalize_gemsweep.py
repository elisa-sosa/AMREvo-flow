#!/usr/bin/env python3
"""Validate and register pinned GEMSWEEP results as provisional sweep outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

from build_reports import TERMINAL_ORDER, collect_cohort_versions, render_plots, write_tsv
from completion_state import atomic_json, build_record, snapshot_paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--resolved", required=True, type=Path)
    parser.add_argument("--gemsweep-output", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--index-prefix", required=True, type=Path)
    parser.add_argument("--ref-groups", required=True, type=Path)
    parser.add_argument("--commit", default="d06885f")
    parser.add_argument("--run-id", required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_index_files(prefix: Path) -> list[Path]:
    parent = prefix.parent.resolve()
    matches = sorted(path.resolve() for path in parent.glob(prefix.name + "*") if path.is_file())
    if not matches:
        raise FileNotFoundError(f"No Themisto index files match prefix: {prefix}")
    return matches


def aggregate_hash(paths: list[Path]) -> str:
    rows = snapshot_paths(paths)
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def abundance_result(root: Path, sample_id: str) -> Path:
    patterns = (
        f"**/{sample_id}_mSWEEP_abundances.txt",
        f"**/{sample_id}.abundances.tsv",
        f"**/{sample_id}*abundances*.txt",
        f"**/{sample_id}*abundances*.tsv",
    )
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(path.resolve() for path in root.glob(pattern) if path.is_file() and path.stat().st_size > 0)
    matches = sorted(set(matches))
    if len(matches) != 1:
        raise FileNotFoundError(f"{sample_id}: expected exactly one GEMSWEEP abundance result, found {len(matches)}")
    return matches[0]


def main() -> int:
    args = parse_args()
    if args.commit != "d06885f":
        raise SystemExit(f"Refusing unreviewed GEMSWEEP commit: {args.commit}")
    if not args.ref_groups.is_file() or args.ref_groups.stat().st_size == 0:
        raise SystemExit(f"Reference group file missing or empty: {args.ref_groups}")
    index_files = find_index_files(args.index_prefix)
    index_hash = aggregate_hash(index_files)
    groups_hash = aggregate_hash([args.ref_groups])

    resolved = {row["sample_id"]: row for row in read_tsv(args.resolved)}
    with args.manifest.open(newline="") as handle:
        manifest_rows = list(csv.DictReader(handle))
    if not manifest_rows:
        raise SystemExit("GEMSWEEP manifest contains no samples")
    if list(manifest_rows[0]) != ["ID", "R1", "R2"]:
        raise SystemExit("GEMSWEEP manifest header must be ID,R1,R2")

    status_path = args.results_root / "cohort" / "current" / "sample_status.tsv"
    cohort_rows = read_tsv(status_path)
    cohort_by_id = {row["sample_id"]: row for row in cohort_rows}
    completed = []
    for manifest_row in manifest_rows:
        sample_id = manifest_row["ID"]
        row = resolved.get(sample_id)
        if row is None or row.get("material") != "sweep":
            raise SystemExit(f"{sample_id}: absent from resolved sweep manifest")
        abundance = abundance_result(args.gemsweep_output, sample_id)
        sample_dir = args.results_root / "samples" / "sweeps" / sample_id
        result_dir = sample_dir / "02_gemsweep"
        result_dir.mkdir(parents=True, exist_ok=True)
        destination = result_dir / f"{sample_id}.PROVISIONAL_GEMSWEEP.abundances.tsv"
        shutil.copy2(abundance, destination)
        provenance = {
            "classification": "PROVISIONAL_GEMSWEEP",
            "gemsweep_commit": args.commit,
            "index_prefix": str(args.index_prefix.resolve()),
            "index_sha256": index_hash,
            "ref_groups": str(args.ref_groups.resolve()),
            "ref_groups_sha256": groups_hash,
            "run_id": args.run_id,
            "upstream_output": str(abundance),
        }
        provenance_path = result_dir / "provenance.json"
        atomic_json(provenance_path, provenance)
        versions = sample_dir / "software_versions.txt"
        with versions.open("a") as handle:
            handle.write(f"gemsweep_commit\t{args.commit}\n")
            handle.write(f"themisto_index_sha256\t{index_hash}\n")
            handle.write(f"reference_groups_sha256\t{groups_hash}\n")
        sample_status = cohort_by_id.get(sample_id)
        if sample_status is None:
            raise SystemExit(f"{sample_id}: absent from canonical cohort status")
        sample_status["status"] = "PROVISIONAL_SWEEP_COMPLETE"
        sample_status["reason"] = "pinned frozen-index GEMSWEEP outputs validated"
        sample_status["run_id"] = args.run_id
        write_tsv(sample_dir / "sample_status.tsv", [sample_status], list(sample_status))

        input_paths = [
            Path(value)
            for key in ("illumina_r1_json", "illumina_r2_json")
            for value in json.loads(row[key])
        ]
        input_paths.extend(Path(value) for value in json.loads(row["trim_metrics_json"]) if value)
        input_paths.extend(index_files)
        input_paths.append(args.ref_groups.resolve())
        record = build_record(
            sample_id=sample_id,
            sample_class="sweep",
            stage="FINAL",
            status="PROVISIONAL_SWEEP_COMPLETE",
            inputs=input_paths,
            outputs=[destination, provenance_path, versions],
            parameter_hash=row["parameter_hash"],
            software_hash=row["software_hash"],
            run_id=args.run_id,
            metadata=provenance,
        )
        atomic_json(sample_dir / "FINAL_COMPLETE.json", record)
        completed.append(sample_id)

    fields = list(cohort_rows[0])
    write_tsv(status_path, cohort_rows, fields)
    counts = Counter(row["status"] for row in cohort_rows)
    retried_count = sum(row.get("retried_previous_failure") == "true" for row in cohort_rows)
    summary = args.results_root / "cohort" / "current" / "run_summary.txt"
    with summary.open("w") as handle:
        handle.write(f"AMREvoFlow run: {args.run_id}\n")
        handle.write(f"Samples: {len(cohort_rows)}\n")
        for status in TERMINAL_ORDER:
            value = retried_count if status == "RETRIED_PREVIOUS_FAILURE" else counts.get(status, 0)
            handle.write(f"{status}\t{value}\n")
        for status in sorted(set(counts) - set(TERMINAL_ORDER)):
            handle.write(f"{status}\t{counts[status]}\n")
        handle.write(f"STATUS_TABLE\t{status_path.resolve()}\n")
        handle.write(f"PLOTS\t{(args.results_root / 'cohort/current/plots').resolve()}\n")
        handle.write(f"DASHBOARD\t{(args.results_root / 'cohort/current/plots/pipeline_status.png').resolve()}\n")
    render_plots(cohort_rows, args.results_root / "cohort" / "current" / "plots")
    cohort_versions = args.results_root / "cohort" / "current" / "software_versions.tsv"
    collect_cohort_versions(args.results_root, cohort_versions)

    snapshot = args.results_root / "cohort" / "snapshots" / args.run_id
    snapshot.mkdir(parents=True, exist_ok=True)
    shutil.copy2(status_path, snapshot / "sample_status.tsv")
    shutil.copy2(summary, snapshot / "run_summary.txt")
    shutil.copy2(cohort_versions, snapshot / "software_versions.tsv")
    plots = snapshot / "plots"
    if plots.exists():
        shutil.rmtree(plots)
    shutil.copytree(args.results_root / "cohort" / "current" / "plots", plots)
    print(f"Validated {len(completed)} PROVISIONAL_GEMSWEEP samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
