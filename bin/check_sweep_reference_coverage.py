#!/usr/bin/env python3
"""Warn when reviewed carriage-event references are absent from a frozen sweep index."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--carriage-references", required=True, type=Path)
    parser.add_argument("--ref-groups", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with args.carriage_references.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        id_column = next(
            (
                name
                for name in (
                    "bin_id", "sweep_group_id", "group_id", "reference_sample_id",
                    "sample_id", "Unique_ID", "reference_id",
                )
                if name in fields
            ),
            None,
        )
        if id_column is None:
            raise SystemExit(f"Cannot identify reference ID column in {args.carriage_references}: {fields}")
        expected = sorted({(row.get(id_column) or "").strip() for row in reader if (row.get(id_column) or "").strip()})
    represented = {
        line.strip().split("\t", 1)[0]
        for line in args.ref_groups.read_text(errors="replace").splitlines()
        if line.strip()
    }
    rows = [
        {
            "reference_id": reference_id,
            "status": "PRESENT" if reference_id in represented else "ABSENT_FROM_FROZEN_INDEX",
            "action": "" if reference_id in represented else "PROVISIONAL_RESULT_WARNING",
        }
        for reference_id in expected
    ]
    output = args.output.open("w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(output, fieldnames=["reference_id", "status", "action"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if args.output:
            output.close()
    missing = sum(row["status"] == "ABSENT_FROM_FROZEN_INDEX" for row in rows)
    print(
        f"WARNING: {missing} of {len(rows)} reviewed carriage-event references are absent from the frozen sweep index."
        if missing else f"All {len(rows)} reviewed carriage-event references are represented in the frozen index.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
