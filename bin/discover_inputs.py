#!/usr/bin/env python3
"""Discover review candidates without selecting or changing sequencing data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path


FIELDS = [
    "sample_id",
    "read_unit_id",
    "platform",
    "source_run",
    "r1",
    "r2",
    "ont_fastq",
    "trim_metrics",
    "decision",
    "facility_qc_status",
    "merge_order",
    "comment",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--illumina-root", action="append", default=[], type=Path)
    parser.add_argument("--ont-root", action="append", default=[], type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def stable_id(platform: str, path: Path) -> str:
    digest = hashlib.sha256(str(path).encode()).hexdigest()[:10]
    return f"candidate-{platform}-{digest}"


def illumina_candidates(root: Path) -> list[dict[str, str]]:
    rows = []
    patterns = ("*_r1.fastq.gz", "*_R1.fastq.gz", "*_1.fastq.gz")
    seen: set[Path] = set()
    for pattern in patterns:
        for r1 in root.rglob(pattern):
            r1 = r1.resolve()
            if r1 in seen or r1.parent.name != "0_trimming":
                continue
            seen.add(r1)
            replacements = (
                ("_r1.fastq.gz", "_r2.fastq.gz"),
                ("_R1.fastq.gz", "_R2.fastq.gz"),
                ("_1.fastq.gz", "_2.fastq.gz"),
            )
            r2 = None
            for left, right in replacements:
                if r1.name.endswith(left):
                    candidate = r1.with_name(r1.name[: -len(left)] + right)
                    if candidate.is_file():
                        r2 = candidate.resolve()
                    break
            sample_id = r1.parent.parent.name
            trim_options = sorted(r1.parent.glob("*.quality_read_trimm_info"))
            trim = trim_options[0].resolve() if len(trim_options) == 1 else None
            source_run = r1.parent.parent.parent.name
            comment = []
            if r2 is None:
                comment.append("MISSING_MATE")
            if len(trim_options) > 1:
                comment.append("AMBIGUOUS_TRIM_METRICS")
            rows.append(
                {
                    "sample_id": sample_id,
                    "read_unit_id": stable_id("illumina", r1),
                    "platform": "illumina",
                    "source_run": source_run,
                    "r1": str(r1),
                    "r2": "" if r2 is None else str(r2),
                    "ont_fastq": "",
                    "trim_metrics": "" if trim is None else str(trim),
                    "decision": "candidate",
                    "facility_qc_status": "accepted" if r2 else "",
                    "merge_order": "",
                    "comment": ";".join(comment) or "REVIEW_REQUIRED",
                }
            )
    return rows


def infer_ont_sample(path: Path) -> str:
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    name = re.sub(r"(?:_pass|_ont|_nanopore|_merged)$", "", name, flags=re.IGNORECASE)
    return name


def ont_candidates(root: Path) -> list[dict[str, str]]:
    rows = []
    seen: set[Path] = set()
    for pattern in ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"):
        for fastq in root.rglob(pattern):
            fastq = fastq.resolve()
            if fastq in seen:
                continue
            seen.add(fastq)
            sample_id = infer_ont_sample(fastq)
            rows.append(
                {
                    "sample_id": sample_id,
                    "read_unit_id": stable_id("ont", fastq),
                    "platform": "ont",
                    "source_run": fastq.parent.name,
                    "r1": "",
                    "r2": "",
                    "ont_fastq": str(fastq),
                    "trim_metrics": "",
                    "decision": "candidate",
                    "facility_qc_status": "",
                    "merge_order": "",
                    "comment": "REVIEW_REQUIRED",
                }
            )
    return rows


def main() -> int:
    args = parse_args()
    rows: list[dict[str, str]] = []
    for root in args.illumina_root:
        if not root.is_dir():
            print(f"ERROR: Illumina root is not a directory: {root}", file=sys.stderr)
            return 2
        rows.extend(illumina_candidates(root))
    for root in args.ont_root:
        if not root.is_dir():
            print(f"ERROR: ONT root is not a directory: {root}", file=sys.stderr)
            return 2
        rows.extend(ont_candidates(root))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["sample_id"], row["platform"], row["source_run"], row["r1"] or row["ont_fastq"])))
    print(f"Wrote {len(rows)} unapproved candidates to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
