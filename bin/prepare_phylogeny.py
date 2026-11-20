#!/usr/bin/env python3
"""Prepare overall and patient/Achtman-ST phylogeny jobs from current results."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--minimum-samples", type=int, default=3)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def timepoint_key(value: str) -> tuple[int, int | str]:
    try:
        return (0, int(value))
    except ValueError:
        match = re.search(r"-?\d+", value)
        return (0, int(match.group())) if match else (1, value)


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return cleaned.strip("_") or "group"


def assembly_path(results_root: Path, sample_id: str) -> Path:
    root = results_root / "samples" / "isolates" / sample_id
    matches = sorted(root.glob("02_hybrid_assembly/04_pypolca/*.polished.fasta"))
    if len(matches) != 1 or not matches[0].is_file() or matches[0].stat().st_size == 0:
        raise FileNotFoundError(f"{sample_id}: expected exactly one polished assembly, found {len(matches)}")
    return matches[0].resolve()


def reference_key(row: dict[str, str]) -> tuple[tuple[int, int | str], float, float, str]:
    completeness = float(row["completeness"] or 0)
    contamination = float(row["contamination"] or 999)
    return (timepoint_key(row["timepoint"]), -completeness, contamination, row["sample_id"])


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(args.status)
    eligible = []
    for row in rows:
        if row.get("material") != "isolate":
            continue
        if row.get("species") != "Escherichia coli":
            continue
        if row.get("status") not in {"COMPLETED", "SKIPPED_COMPLETE"}:
            continue
        if not row.get("sequence_type"):
            continue
        row = dict(row)
        row["assembly"] = str(assembly_path(args.results_root, row["sample_id"]))
        eligible.append(row)

    groups: dict[str, tuple[str, list[dict[str, str]]]] = {}
    if eligible:
        groups["overall"] = ("overall", eligible)
    events: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in eligible:
        events[(row["patient_id"], row["sequence_type"])].append(row)
    for (patient, sequence_type), members in sorted(events.items()):
        group_id = safe_id(f"patient_{patient}__ST{sequence_type}")
        groups[group_id] = ("patient_st_event", members)

    jobs = []
    references = []
    statuses = []
    for group_id, (scope, members) in groups.items():
        reference = min(members, key=reference_key)
        references.append(
            {
                "group_id": group_id,
                "scope": scope,
                "reference_sample_id": reference["sample_id"],
                "reference_assembly": reference["assembly"],
                "timepoint": reference["timepoint"],
                "completeness": reference["completeness"],
                "contamination": reference["contamination"],
                "selection_rule": "earliest_timepoint_then_highest_completeness_then_lowest_contamination_then_sample_id",
            }
        )
        if len(members) < args.minimum_samples:
            statuses.append(
                {
                    "group_id": group_id,
                    "scope": scope,
                    "sample_count": str(len(members)),
                    "status": "INSUFFICIENT_SAMPLES",
                    "reason": f"requires at least {args.minimum_samples} samples",
                }
            )
            continue
        statuses.append(
            {
                "group_id": group_id,
                "scope": scope,
                "sample_count": str(len(members)),
                "status": "READY",
                "reason": "",
            }
        )
        for row in sorted(members, key=lambda item: item["sample_id"]):
            jobs.append(
                {
                    "group_id": group_id,
                    "scope": scope,
                    "sample_id": row["sample_id"],
                    "assembly": row["assembly"],
                    "reference": reference["assembly"],
                    "reference_sample_id": reference["sample_id"],
                }
            )

    outputs = (
        ("phylogeny_jobs.tsv", jobs, ["group_id", "scope", "sample_id", "assembly", "reference", "reference_sample_id"]),
        ("reference_manifest.tsv", references, ["group_id", "scope", "reference_sample_id", "reference_assembly", "timepoint", "completeness", "contamination", "selection_rule"]),
        ("phylogeny_group_status.tsv", statuses, ["group_id", "scope", "sample_count", "status", "reason"]),
    )
    for name, values, fields in outputs:
        with (args.output_dir / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(values)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
