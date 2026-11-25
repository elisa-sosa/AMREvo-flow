#!/usr/bin/env python3
"""Build AMREvo isolate manifests from completed LIMS, pending v3 and reads."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from discover_inputs import FIELDS as READ_FIELDS
from discover_inputs import illumina_candidates, ont_candidates


SAMPLE_FIELDS = [
    "sample_id",
    "patient_id",
    "timepoint",
    "material",
    "agar_type",
    "expect_illumina",
    "expect_ont",
    "active",
]
RUN_PATTERN = re.compile(r"^esccol_(\d{6})$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--completed-lims", required=True, type=Path)
    parser.add_argument("--pending-v3", required=True, type=Path)
    parser.add_argument("--illumina-root", required=True, type=Path)
    parser.add_argument("--ont-root", required=True, type=Path)
    parser.add_argument("--samples-output", required=True, type=Path)
    parser.add_argument("--read-units-output", required=True, type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"Input table does not exist: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def first(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = (row.get(name) or "").strip()
        if value:
            return value
    return ""


def patient_id(value: str, sample_id: str) -> str:
    value = value.strip()
    if value.upper().startswith("P"):
        return value.upper()
    if value.isdigit():
        return f"P{int(value):03d}"
    match = re.match(r"^(P\d{3})-", sample_id)
    if match:
        return match.group(1)
    raise ValueError(f"Cannot determine patient_id for {sample_id!r}")


def agar_type(value: str, sample_id: str) -> str:
    if value.strip():
        return value.strip().upper()
    match = re.match(r"^P\d{3}-([A-Za-z]+)-", sample_id)
    if match:
        return match.group(1).upper()
    raise ValueError(f"Cannot determine agar_type for {sample_id!r}")


def sample_row(row: dict[str, str], source: str) -> dict[str, str]:
    sample_id = first(row, "sample_id", "Unique_ID")
    if not sample_id:
        raise ValueError(f"{source}: row has no sample_id or Unique_ID")
    return {
        "sample_id": sample_id,
        "patient_id": patient_id(first(row, "Patient_Nr", "PatientNr"), sample_id),
        "timepoint": first(row, "Timepoint"),
        "material": "isolate",
        "agar_type": agar_type(first(row, "Agar_Type", "AgarType"), sample_id),
        "expect_illumina": "true",
        "expect_ont": "true",
        "active": "true",
    }


def build_registry(completed: Path, pending: Path) -> dict[str, dict[str, str]]:
    registry: dict[str, dict[str, str]] = {}
    for source, path in (("completed LIMS", completed), ("pending v3", pending)):
        for raw in read_csv(path):
            row = sample_row(raw, source)
            sample_id = row["sample_id"]
            previous = registry.get(sample_id)
            if previous and previous != row:
                raise ValueError(
                    f"Conflicting metadata for {sample_id} between completed LIMS and pending v3"
                )
            registry[sample_id] = row
    if not registry:
        raise ValueError("Combined isolate registry is empty")
    return registry


def append_comment(row: dict[str, str], value: str) -> None:
    existing = row.get("comment", "")
    row["comment"] = ";".join(part for part in (existing, value) if part)


def run_date(source_run: str) -> datetime | None:
    match = RUN_PATTERN.fullmatch(source_run)
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%y%m%d")
    except ValueError:
        return None


def select_illumina(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample_id"]].append(dict(row))

    selected: list[dict[str, str]] = []
    for sample_id, candidates in grouped.items():
        dated = [(run_date(row["source_run"]), row) for row in candidates]
        invalid = [row for date, row in dated if date is None]
        if invalid:
            for row in candidates:
                row["decision"] = "candidate"
                row["merge_order"] = ""
                append_comment(row, "INVALID_OR_AMBIGUOUS_RUN_DATE")
            selected.extend(candidates)
            continue

        latest_date = max(date for date, _ in dated if date is not None)
        latest = [row for date, row in dated if date == latest_date]
        for row in candidates:
            if row not in latest:
                row["decision"] = "rejected"
                row["merge_order"] = ""
                append_comment(row, "SUPERSEDED_BY_LATEST_DATED_RUN")

        if len(latest) != 1:
            for row in latest:
                row["decision"] = "candidate"
                row["merge_order"] = ""
                append_comment(row, "AMBIGUOUS_LATEST_DATED_RUN")
        else:
            row = latest[0]
            if row["r1"] and row["r2"] and "AMBIGUOUS_TRIM_METRICS" not in row["comment"]:
                row["decision"] = "approved"
                row["facility_qc_status"] = "accepted"
                row["merge_order"] = "1"
                append_comment(row, "AUTO_APPROVED_LATEST_DATED_TRIMMED_PAIR")
            else:
                row["decision"] = "candidate"
                row["merge_order"] = ""
                append_comment(row, "LATEST_DATED_RUN_INCOMPLETE_NO_FALLBACK")
        selected.extend(candidates)
    return selected


def select_ont(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["sample_id"]].append(dict(row))
    selected: list[dict[str, str]] = []
    for candidates in grouped.values():
        if len(candidates) == 1:
            candidates[0]["decision"] = "approved"
            candidates[0]["merge_order"] = "1"
            append_comment(candidates[0], "AUTO_APPROVED_UNIQUE_PASS_FASTQ")
        else:
            for row in candidates:
                row["decision"] = "candidate"
                row["merge_order"] = ""
                append_comment(row, "AMBIGUOUS_MULTIPLE_ONT_FASTQ")
        selected.extend(candidates)
    return selected


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        registry = build_registry(args.completed_lims, args.pending_v3)
        if not args.illumina_root.is_dir() or not args.ont_root.is_dir():
            raise ValueError("Illumina and ONT roots must both exist")
        illumina = select_illumina(illumina_candidates(args.illumina_root))
        ont = select_ont(ont_candidates(args.ont_root))
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    known_reads = [row for row in illumina + ont if row["sample_id"] in registry]
    unknown = sorted({row["sample_id"] for row in illumina + ont if row["sample_id"] not in registry})
    samples = [registry[key] for key in sorted(registry)]
    reads = sorted(
        known_reads,
        key=lambda row: (row["sample_id"], row["platform"], row["source_run"], row["read_unit_id"]),
    )
    write_tsv(args.samples_output, SAMPLE_FIELDS, samples)
    write_tsv(args.read_units_output, READ_FIELDS, reads)

    decisions: dict[str, int] = defaultdict(int)
    for row in reads:
        decisions[row["decision"]] += 1
    print(f"Wrote {len(samples)} isolate samples to {args.samples_output}", file=sys.stderr)
    print(
        f"Wrote {len(reads)} read units to {args.read_units_output}: "
        + ", ".join(f"{key}={decisions[key]}" for key in sorted(decisions)),
        file=sys.stderr,
    )
    if unknown:
        print(
            f"WARNING: ignored {len(unknown)} read sample IDs absent from the LIMS+v3 registry; "
            f"first IDs: {', '.join(unknown[:10])}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
