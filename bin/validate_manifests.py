#!/usr/bin/env python3
"""Validate AMREvoFlow sample/read-unit manifests and resolve approved inputs.

The program is intentionally strict about approved inputs and intentionally
passive about discovery: it never selects a newest run and never modifies raw
data.  Full FASTQ validation is I/O intensive and belongs in a scheduled
Nextflow validation process, not an HPC login shell.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator, TextIO

from completion_state import fingerprint, sha256_file, validate_record


SAMPLE_COLUMNS = [
    "sample_id",
    "patient_id",
    "timepoint",
    "material",
    "agar_type",
    "expect_illumina",
    "expect_ont",
    "active",
]
READ_UNIT_COLUMNS = [
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
OVERRIDE_COLUMNS = ["sample_id", "stage", "reviewer", "date", "reason", "approval"]
VALID_DECISIONS = {"approved", "candidate", "rejected"}
VALID_FACILITY = {"accepted", "rerun_pending", "rejected"}
VALID_OVERRIDE_STAGES = {"ILLUMINA_QC", "ONT_QC", "ASSEMBLY_QC"}
FAILED_SAMPLE_STATUSES = {
    "ILLUMINA_FACILITY_RERUN_PENDING",
    "BLOCKED_ONT_QC",
    "BLOCKED_ASSEMBLY_QC",
    "PROCESS_FAILURE",
}
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class ValidationError(ValueError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-class", required=True, choices=("isolate", "sweep"))
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--read-units", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--allowed-root", action="append", default=[], type=Path)
    parser.add_argument(
        "--allowed-roots-json",
        help="JSON array alternative to repeated --allowed-root (convenient from Nextflow)",
    )
    parser.add_argument("--qc-overrides", type=Path)
    parser.add_argument("--results-root", type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--parameter-hash", required=True)
    parser.add_argument("--software-hash", required=True)
    parser.add_argument("--dependency", action="append", default=[], type=Path)
    parser.add_argument("--dependency-prefix", type=Path)
    parser.add_argument(
        "--fastq-validation",
        choices=("none", "headers", "full"),
        default="full",
        help="headers checks the first --header-records pairs; full scans every record",
    )
    parser.add_argument("--header-records", type=int, default=10_000)
    return parser.parse_args()


def read_tsv(path: Path, required_columns: list[str]) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValidationError(f"Manifest is missing or empty: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        actual = reader.fieldnames or []
        if actual != required_columns:
            raise ValidationError(
                f"{path}: header must be exactly {required_columns!r}; observed {actual!r}"
            )
        rows = []
        for number, row in enumerate(reader, start=2):
            clean = {key: (value or "").strip() for key, value in row.items()}
            clean["_line"] = str(number)
            rows.append(clean)
    return rows


def parse_bool(value: str, *, field: str, path: Path, line: str) -> bool:
    lowered = value.lower()
    if lowered not in {"true", "false"}:
        raise ValidationError(f"{path}:{line}: {field} must be true or false")
    return lowered == "true"


def within(path: Path, roots: list[Path]) -> bool:
    return any(path == root or root in path.parents for root in roots)


def checked_path(raw: str, *, label: str, roots: list[Path], context: str) -> Path:
    if not raw:
        raise ValidationError(f"{context}: missing {label} path")
    path = Path(raw).expanduser()
    if not path.is_absolute():
        raise ValidationError(f"{context}: {label} must be an absolute path: {raw}")
    resolved = path.resolve()
    if not within(resolved, roots):
        raise ValidationError(f"{context}: {label} is outside all --allowed-root paths: {resolved}")
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise ValidationError(f"{context}: missing or empty {label}: {resolved}")
    return resolved


def open_fastq(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def read_id(header: str) -> str:
    token = header[1:].strip().split()[0]
    return re.sub(r"/[12]$", "", token)


def fastq_records(path: Path, limit: int | None) -> Iterator[tuple[str, int]]:
    with open_fastq(path) as handle:
        record = 0
        while limit is None or record < limit:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record += 1
            if not (sequence and plus and quality):
                raise ValidationError(f"{path}: truncated FASTQ record {record}")
            if not header.startswith("@"):
                raise ValidationError(f"{path}: record {record} header does not start with @")
            if not plus.startswith("+"):
                raise ValidationError(f"{path}: record {record} separator does not start with +")
            if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                raise ValidationError(f"{path}: sequence/quality length mismatch at record {record}")
            yield read_id(header), record


def validate_pair(
    r1: Path,
    r2: Path,
    *,
    mode: str,
    limit: int,
    cumulative_ids: set[str],
) -> int | None:
    if mode == "none":
        return None
    maximum = None if mode == "full" else limit
    left = fastq_records(r1, maximum)
    right = fastq_records(r2, maximum)
    count = 0
    while True:
        one = next(left, None)
        two = next(right, None)
        if one is None and two is None:
            break
        if one is None or two is None:
            raise ValidationError(f"Paired FASTQ record count differs: {r1} vs {r2}")
        count += 1
        if one[0] != two[0]:
            raise ValidationError(
                f"Paired FASTQ identifiers differ at record {count}: {one[0]!r} vs {two[0]!r}"
            )
        if one[0] in cumulative_ids:
            raise ValidationError(f"Duplicate approved read identifier for merged units: {one[0]}")
        cumulative_ids.add(one[0])
    if count == 0:
        raise ValidationError(f"FASTQ contains no records: {r1}")
    return count


def validate_single_fastq(
    path: Path, *, mode: str, limit: int, cumulative_ids: set[str]
) -> int | None:
    if mode == "none":
        return None
    maximum = None if mode == "full" else limit
    count = 0
    for identifier, _ in fastq_records(path, maximum):
        count += 1
        if identifier in cumulative_ids:
            raise ValidationError(f"Duplicate approved read identifier for merged units: {identifier}")
        cumulative_ids.add(identifier)
    if count == 0:
        raise ValidationError(f"FASTQ contains no records: {path}")
    return count


def file_record(path: Path) -> dict[str, Any]:
    return {"path": str(path), "size": path.stat().st_size, "sha256": sha256_file(path)}


def load_overrides(path: Path | None, sample_ids: set[str]) -> dict[str, set[str]]:
    approved: dict[str, set[str]] = defaultdict(set)
    if path is None or not path.exists():
        return approved
    for row in read_tsv(path, OVERRIDE_COLUMNS):
        context = f"{path}:{row['_line']}"
        if row["sample_id"] not in sample_ids:
            raise ValidationError(f"{context}: override sample_id is absent from sample manifest")
        if row["stage"] not in VALID_OVERRIDE_STAGES:
            raise ValidationError(f"{context}: unsupported override stage {row['stage']!r}")
        if not row["reviewer"] or not row["date"] or not row["reason"]:
            raise ValidationError(f"{context}: reviewer, date and reason are required")
        if row["approval"] not in {"approved", "withdrawn"}:
            raise ValidationError(f"{context}: approval must be approved or withdrawn")
        if row["approval"] == "approved":
            approved[row["sample_id"]].add(row["stage"])
        else:
            approved[row["sample_id"]].discard(row["stage"])
    return approved


def previous_status(sample_dir: Path) -> str:
    status_file = sample_dir / "sample_status.tsv"
    if not status_file.is_file():
        return ""
    try:
        with status_file.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        return (rows[-1].get("status") or "").strip() if rows else ""
    except (OSError, csv.Error):
        return ""


def json_cell(value: Any) -> str:
    return json.dumps(value, separators=(",", ":"), sort_keys=True)


def main() -> int:
    args = parse_args()
    try:
        configured_roots = list(args.allowed_root)
        if args.allowed_roots_json:
            try:
                configured_roots.extend(Path(value) for value in json.loads(args.allowed_roots_json))
            except (json.JSONDecodeError, TypeError) as error:
                raise ValidationError("--allowed-roots-json must be a JSON array of paths") from error
        if not configured_roots:
            raise ValidationError("At least one explicit --allowed-root is required")
        roots = [root.expanduser().resolve() for root in configured_roots]
        missing_roots = [root for root in roots if not root.is_dir()]
        if missing_roots:
            raise ValidationError(f"Allowed root does not exist: {missing_roots[0]}")

        dependency_paths = [path.expanduser().resolve() for path in args.dependency]
        if args.dependency_prefix:
            prefix = args.dependency_prefix.expanduser()
            prefix_matches = sorted(path.resolve() for path in prefix.parent.glob(prefix.name + "*") if path.is_file())
            if not prefix_matches:
                raise ValidationError(f"No dependency files match prefix: {args.dependency_prefix}")
            dependency_paths.extend(prefix_matches)
        dependency_paths = sorted(set(dependency_paths))
        for dependency in dependency_paths:
            if not dependency.is_file() or dependency.stat().st_size == 0:
                raise ValidationError(f"Dependency is missing or empty: {dependency}")
        dependency_records = [file_record(path) for path in dependency_paths]

        sample_rows = read_tsv(args.samples, SAMPLE_COLUMNS)
        if not sample_rows:
            raise ValidationError(f"Sample manifest contains no rows: {args.samples}")
        samples: dict[str, dict[str, Any]] = {}
        for row in sample_rows:
            context = f"{args.samples}:{row['_line']}"
            sample_id = row["sample_id"]
            if not sample_id:
                raise ValidationError(f"{context}: empty sample_id")
            if not SAFE_SAMPLE_ID.fullmatch(sample_id):
                raise ValidationError(
                    f"{context}: sample_id contains unsafe path characters: {sample_id!r}"
                )
            if sample_id in samples:
                raise ValidationError(f"{context}: duplicate sample_id {sample_id!r}")
            if row["material"] != args.sample_class:
                raise ValidationError(
                    f"{context}: material must be {args.sample_class!r}, observed {row['material']!r}"
                )
            samples[sample_id] = {
                **{column: row[column] for column in SAMPLE_COLUMNS},
                "expect_illumina": parse_bool(
                    row["expect_illumina"], field="expect_illumina", path=args.samples, line=row["_line"]
                ),
                "expect_ont": parse_bool(
                    row["expect_ont"], field="expect_ont", path=args.samples, line=row["_line"]
                ),
                "active": parse_bool(row["active"], field="active", path=args.samples, line=row["_line"]),
            }

        overrides = load_overrides(args.qc_overrides, set(samples))
        units_by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
        seen_unit_ids: set[str] = set()
        for row in read_tsv(args.read_units, READ_UNIT_COLUMNS):
            context = f"{args.read_units}:{row['_line']}"
            sample_id = row["sample_id"]
            if sample_id not in samples:
                raise ValidationError(f"{context}: sample_id is absent from sample manifest: {sample_id}")
            if not row["read_unit_id"] or row["read_unit_id"] in seen_unit_ids:
                raise ValidationError(f"{context}: empty or duplicate read_unit_id {row['read_unit_id']!r}")
            seen_unit_ids.add(row["read_unit_id"])
            if row["platform"].lower() not in {"illumina", "ont"}:
                raise ValidationError(f"{context}: platform must be illumina or ont")
            if row["decision"].lower() not in VALID_DECISIONS:
                raise ValidationError(f"{context}: decision must be approved, candidate or rejected")
            if not row["source_run"]:
                raise ValidationError(f"{context}: source_run is required")
            row["platform"] = row["platform"].lower()
            row["decision"] = row["decision"].lower()
            row["_context"] = context
            units_by_sample[sample_id].append(row)

        args.output_dir.mkdir(parents=True, exist_ok=True)
        resolved_rows: list[dict[str, str]] = []
        checksum_rows: list[dict[str, str]] = []
        seen_approved_paths: dict[str, str] = {}
        summary_counts: dict[str, int] = defaultdict(int)

        for sample_id in sorted(samples):
            sample = samples[sample_id]
            approved = [row for row in units_by_sample[sample_id] if row["decision"] == "approved"]
            platform_units: dict[str, list[dict[str, Any]]] = {}
            input_records: list[dict[str, Any]] = list(dependency_records if args.sample_class == "sweep" else [])
            illumina_input_records: list[dict[str, Any]] = []
            ont_input_records: list[dict[str, Any]] = []
            record_counts: dict[str, list[int | None]] = {"illumina": [], "ont": []}
            for platform in ("illumina", "ont"):
                units = [row for row in approved if row["platform"] == platform]
                if len(units) > 1 and any(not row["merge_order"] for row in units):
                    raise ValidationError(
                        f"{sample_id}: multiple approved {platform} units require explicit merge_order values"
                    )
                orders: set[int] = set()
                for row in units:
                    raw_order = row["merge_order"] or "1"
                    try:
                        order = int(raw_order)
                    except ValueError as error:
                        raise ValidationError(f"{row['_context']}: merge_order must be a positive integer") from error
                    if order < 1 or order in orders:
                        raise ValidationError(f"{row['_context']}: duplicate or invalid merge_order {order}")
                    orders.add(order)
                    row["_merge_order"] = order
                platform_units[platform] = sorted(units, key=lambda row: row.get("_merge_order", 1))

            illumina_paths: list[tuple[Path, Path, Path | None]] = []
            illumina_ids: set[str] = set()
            facility_states: list[str] = []
            for unit in platform_units["illumina"]:
                context = unit["_context"]
                if unit["ont_fastq"]:
                    raise ValidationError(f"{context}: Illumina unit must not declare ont_fastq")
                r1 = checked_path(unit["r1"], label="r1", roots=roots, context=context)
                r2 = checked_path(unit["r2"], label="r2", roots=roots, context=context)
                if r1 == r2:
                    raise ValidationError(f"{context}: r1 and r2 resolve to the same path")
                trim = (
                    checked_path(unit["trim_metrics"], label="trim_metrics", roots=roots, context=context)
                    if unit["trim_metrics"]
                    else None
                )
                facility = unit["facility_qc_status"].lower() or "accepted"
                if facility not in VALID_FACILITY:
                    raise ValidationError(f"{context}: invalid facility_qc_status {facility!r}")
                facility_states.append(facility)
                count = validate_pair(
                    r1,
                    r2,
                    mode=args.fastq_validation,
                    limit=args.header_records,
                    cumulative_ids=illumina_ids,
                )
                record_counts["illumina"].append(count)
                illumina_paths.append((r1, r2, trim))
                for kind, path in (("r1", r1), ("r2", r2), ("trim_metrics", trim)):
                    if path is None:
                        continue
                    owner = seen_approved_paths.get(str(path))
                    if owner:
                        raise ValidationError(f"Approved path is reused by {owner} and {sample_id}: {path}")
                    seen_approved_paths[str(path)] = sample_id
                    file_row = file_record(path)
                    input_records.append(file_row)
                    illumina_input_records.append(file_row)
                    checksum_rows.append(
                        {"sample_id": sample_id, "read_unit_id": unit["read_unit_id"], "kind": kind, **file_row}
                    )

            ont_paths: list[Path] = []
            ont_ids: set[str] = set()
            for unit in platform_units["ont"]:
                context = unit["_context"]
                if unit["r1"] or unit["r2"] or unit["trim_metrics"]:
                    raise ValidationError(f"{context}: ONT unit must only declare ont_fastq")
                if unit["facility_qc_status"]:
                    raise ValidationError(f"{context}: ONT unit must leave facility_qc_status blank")
                ont = checked_path(unit["ont_fastq"], label="ont_fastq", roots=roots, context=context)
                count = validate_single_fastq(
                    ont,
                    mode=args.fastq_validation,
                    limit=args.header_records,
                    cumulative_ids=ont_ids,
                )
                record_counts["ont"].append(count)
                owner = seen_approved_paths.get(str(ont))
                if owner:
                    raise ValidationError(f"Approved path is reused by {owner} and {sample_id}: {ont}")
                seen_approved_paths[str(ont)] = sample_id
                file_row = file_record(ont)
                input_records.append(file_row)
                ont_input_records.append(file_row)
                checksum_rows.append(
                    {"sample_id": sample_id, "read_unit_id": unit["read_unit_id"], "kind": "ont_fastq", **file_row}
                )
                ont_paths.append(ont)

            facility_status = ""
            if facility_states:
                if "rejected" in facility_states:
                    facility_status = "rejected"
                elif "rerun_pending" in facility_states:
                    facility_status = "rerun_pending"
                else:
                    facility_status = "accepted"
            has_illumina = bool(illumina_paths)
            has_ont = bool(ont_paths)
            if not sample["active"]:
                status = "INACTIVE"
            elif facility_status in {"rerun_pending", "rejected"} and "ILLUMINA_QC" not in overrides[sample_id]:
                status = "ILLUMINA_FACILITY_RERUN_PENDING"
            elif has_illumina and has_ont:
                status = "READY_FOR_QC"
            elif has_illumina:
                status = "WAITING_FOR_ONT" if sample["expect_ont"] else "ILLUMINA_QC_ONLY"
            elif has_ont:
                status = "WAITING_FOR_ILLUMINA" if sample["expect_illumina"] else "ONT_QC_ONLY"
            else:
                status = "WAITING_FOR_READS"

            input_fingerprint = fingerprint(input_records) if input_records else ""
            illumina_fingerprint = fingerprint(illumina_input_records) if illumina_input_records else ""
            ont_fingerprint = fingerprint(ont_input_records) if ont_input_records else ""
            sample_root = (
                args.results_root / "samples" / ("isolates" if args.sample_class == "isolate" else "sweeps") / sample_id
                if args.results_root
                else None
            )
            prior = previous_status(sample_root) if sample_root else ""
            final_record = sample_root / "FINAL_COMPLETE.json" if sample_root else None
            skip_complete = False
            invalid_completion_reasons: list[str] = []
            if final_record and final_record.exists():
                skip_complete, invalid_completion_reasons, _ = validate_record(
                    final_record,
                    expected_sample_id=sample_id,
                    expected_parameter_hash=args.parameter_hash,
                    expected_software_hash=args.software_hash,
                    expected_input_fingerprint=input_fingerprint,
                    require_final=True,
                )
                if skip_complete:
                    status = "SKIPPED_COMPLETE"
            previously_failed = prior in FAILED_SAMPLE_STATUSES
            retry_token = args.run_id if previously_failed else input_fingerprint or "no-input"
            if previously_failed and status != "SKIPPED_COMPLETE":
                status_before_retry = status
                status = "RETRIED_PREVIOUS_FAILURE"
            else:
                status_before_retry = ""

            row = {
                "sample_id": sample_id,
                "patient_id": str(sample["patient_id"]),
                "timepoint": str(sample["timepoint"]),
                "material": str(sample["material"]),
                "agar_type": str(sample["agar_type"]),
                "expect_illumina": str(sample["expect_illumina"]).lower(),
                "expect_ont": str(sample["expect_ont"]).lower(),
                "active": str(sample["active"]).lower(),
                "illumina_r1_json": json_cell([str(item[0]) for item in illumina_paths]),
                "illumina_r2_json": json_cell([str(item[1]) for item in illumina_paths]),
                "trim_metrics_json": json_cell([str(item[2]) if item[2] else "" for item in illumina_paths]),
                "ont_fastq_json": json_cell([str(item) for item in ont_paths]),
                "source_runs_json": json_cell([row["source_run"] for row in approved]),
                "read_unit_ids_json": json_cell([row["read_unit_id"] for row in approved]),
                "facility_qc_status": facility_status,
                "input_fingerprint": input_fingerprint,
                "illumina_fingerprint": illumina_fingerprint,
                "ont_fingerprint": ont_fingerprint,
                "parameter_hash": args.parameter_hash,
                "software_hash": args.software_hash,
                "run_id": args.run_id,
                "retry_token": retry_token,
                "prior_status": prior,
                "status_before_retry": status_before_retry,
                "status": status,
                "skip_complete": str(skip_complete).lower(),
                "illumina_qc_override": str("ILLUMINA_QC" in overrides[sample_id]).lower(),
                "ont_qc_override": str("ONT_QC" in overrides[sample_id]).lower(),
                "assembly_qc_override": str("ASSEMBLY_QC" in overrides[sample_id]).lower(),
                "invalid_completion_reasons_json": json_cell(invalid_completion_reasons),
                "fastq_validation": args.fastq_validation,
                "illumina_records_checked_json": json_cell(record_counts["illumina"]),
                "ont_records_checked_json": json_cell(record_counts["ont"]),
            }
            resolved_rows.append(row)
            summary_counts[status] += 1

        resolved_path = args.output_dir / "resolved_samples.tsv"
        with resolved_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(resolved_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(resolved_rows)
        checksum_path = args.output_dir / "input_checksums.tsv"
        checksum_fields = ["sample_id", "read_unit_id", "kind", "path", "size", "sha256"]
        with checksum_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=checksum_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(checksum_rows)
        summary = {
            "sample_class": args.sample_class,
            "samples": len(resolved_rows),
            "approved_files": len(checksum_rows),
            "fastq_validation": args.fastq_validation,
            "run_id": args.run_id,
            "status_counts": dict(sorted(summary_counts.items())),
        }
        with (args.output_dir / "validation_summary.json").open("w") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(json.dumps(summary, sort_keys=True), file=sys.stderr)
        return 0
    except (ValidationError, OSError, gzip.BadGzipFile) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
