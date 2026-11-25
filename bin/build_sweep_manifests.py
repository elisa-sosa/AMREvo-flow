#!/usr/bin/env python3
"""Build sweep manifests exclusively from the reconstructed all/ manifest."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from build_isolate_manifests import READ_FIELDS, SAMPLE_FIELDS, write_tsv
from discover_inputs import stable_id


SWEEP_ID = re.compile(r"^(P\d{3})-(MCA)-(-?\d{2})$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--production-manifest", type=Path)
    source.add_argument("--all-root", type=Path)
    parser.add_argument("--samples-output", required=True, type=Path)
    parser.add_argument("--read-units-output", required=True, type=Path)
    return parser.parse_args()


def load_source_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    if args.production_manifest:
        if not args.production_manifest.is_file():
            raise ValueError(f"Production manifest does not exist: {args.production_manifest}")
        with args.production_manifest.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ["ID", "R1", "R2"]:
                raise ValueError("Sweep production manifest header must be ID,R1,R2")
            return list(reader)
    if not args.all_root.is_dir():
        raise ValueError(f"Reconstructed sweep directory does not exist: {args.all_root}")
    rows: list[dict[str, str]] = []
    for r1 in sorted(args.all_root.glob("*_R1.fastq.gz")):
        sample_id = r1.name[: -len("_R1.fastq.gz")]
        r2 = args.all_root / f"{sample_id}_R2.fastq.gz"
        if not r2.is_file():
            raise ValueError(f"{sample_id}: reconstructed R2 mate is missing")
        rows.append({"ID": sample_id, "R1": str(r1.resolve()), "R2": str(r2.resolve())})
    orphan_r2 = [
        path
        for path in args.all_root.glob("*_R2.fastq.gz")
        if not (args.all_root / path.name.replace("_R2.fastq.gz", "_R1.fastq.gz")).is_file()
    ]
    if orphan_r2:
        raise ValueError(f"Reconstructed sweep R2 has no R1 mate: {orphan_r2[0]}")
    return rows


def main() -> int:
    args = parse_args()
    try:
        source_rows = load_source_rows(args)
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    if not source_rows:
        print("ERROR: Sweep production manifest is empty", file=sys.stderr)
        return 2

    samples: list[dict[str, str]] = []
    reads: list[dict[str, str]] = []
    seen: set[str] = set()
    for source in source_rows:
        sample_id = source["ID"].strip()
        match = SWEEP_ID.fullmatch(sample_id)
        if not match:
            print(f"ERROR: Invalid canonical sweep ID: {sample_id}", file=sys.stderr)
            return 2
        if sample_id in seen:
            print(f"ERROR: Duplicate sweep ID: {sample_id}", file=sys.stderr)
            return 2
        seen.add(sample_id)
        r1 = Path(source["R1"].strip())
        r2 = Path(source["R2"].strip())
        if r1.parent.name != "all" or r2.parent.name != "all":
            print(f"ERROR: {sample_id}: sweep inputs must come from sweeps_fastq/all/", file=sys.stderr)
            return 2
        if r1.name != f"{sample_id}_R1.fastq.gz" or r2.name != f"{sample_id}_R2.fastq.gz":
            print(f"ERROR: {sample_id}: reconstructed pair filenames do not match the ID", file=sys.stderr)
            return 2
        timepoint = str(int(match.group(3)))
        samples.append(
            {
                "sample_id": sample_id,
                "patient_id": match.group(1),
                "timepoint": timepoint,
                "material": "sweep",
                "agar_type": match.group(2),
                "expect_illumina": "true",
                "expect_ont": "false",
                "active": "true",
            }
        )
        reads.append(
            {
                "sample_id": sample_id,
                "read_unit_id": stable_id("illumina", r1),
                "platform": "illumina",
                "source_run": "reconstructed_all",
                "r1": str(r1),
                "r2": str(r2),
                "ont_fastq": "",
                "trim_metrics": "",
                "decision": "approved",
                "facility_qc_status": "accepted",
                "merge_order": "1",
                "comment": "APPROVED_RECONSTRUCTED_SWEEP_ALL",
            }
        )

    write_tsv(args.samples_output, SAMPLE_FIELDS, sorted(samples, key=lambda row: row["sample_id"]))
    write_tsv(args.read_units_output, READ_FIELDS, sorted(reads, key=lambda row: row["sample_id"]))
    print(f"Wrote {len(samples)} reconstructed sweep samples", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
