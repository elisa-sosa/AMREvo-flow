#!/usr/bin/env python3
"""Validate explicit Illumina and ONT file mappings without modifying inputs."""

import argparse
import csv
import sys
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_id",
    "illumina_r1_relpath",
    "illumina_r2_relpath",
    "ont_relpath",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--illumina-root", required=True, type=Path)
    parser.add_argument("--ont-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def safe_path(root: Path, relative_path: str, label: str, sample_id: str) -> Path:
    candidate = (root / relative_path).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError(f"{sample_id}: {label} escapes its configured root: {relative_path}") from error
    if not candidate.is_file() or candidate.stat().st_size == 0:
        raise ValueError(f"{sample_id}: missing or empty {label}: {candidate}")
    return candidate


def main():
    args = parse_args()
    if not args.illumina_root.is_dir():
        sys.exit(f"Illumina root is not a directory: {args.illumina_root}")
    if not args.ont_root.is_dir():
        sys.exit(f"ONT root is not a directory: {args.ont_root}")

    with args.samples.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames or []))
            sys.exit(f"Sample sheet is missing required columns: {', '.join(missing)}")

        rows = list(reader)

    if not rows:
        sys.exit("Sample sheet has no data rows")

    seen = set()
    validated_rows = []
    try:
        for row in rows:
            sample_id = row["sample_id"].strip()
            if not sample_id:
                raise ValueError("Sample sheet contains an empty sample_id")
            if sample_id in seen:
                raise ValueError(f"Duplicate sample_id in sample sheet: {sample_id}")
            seen.add(sample_id)

            r1 = safe_path(args.illumina_root, row["illumina_r1_relpath"], "Illumina R1", sample_id)
            r2 = safe_path(args.illumina_root, row["illumina_r2_relpath"], "Illumina R2", sample_id)
            ont = safe_path(args.ont_root, row["ont_relpath"], "ONT FASTQ", sample_id)
            validated_rows.append({"sample_id": sample_id, "illumina_r1": r1, "illumina_r2": r2, "ont_fastq": ont})
    except ValueError as error:
        sys.exit(str(error))

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=validated_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(validated_rows)

    print(f"Validated {len(validated_rows)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
