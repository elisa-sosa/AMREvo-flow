#!/usr/bin/env python3
"""Create and validate checksummed AMREvoFlow completion records.

File existence and Nextflow's cache are useful implementation details, but
neither is proof that a biological result is still current.  This helper makes
the skip contract explicit and is intentionally usable both from Nextflow
processes and from tests/recovery commands.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1
CACHE_VALID_FINAL_STATUSES = {
    "COMPLETE",
    "NON_ECOLI_COMPLETE",
    "PROVISIONAL_SWEEP_COMPLETE",
}


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def iter_files(paths: Iterable[Path]) -> list[Path]:
    files: list[Path] = []
    for original in paths:
        path = original.expanduser().resolve()
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files.extend(sorted(item.resolve() for item in path.rglob("*") if item.is_file()))
        else:
            raise FileNotFoundError(f"Required path is missing: {path}")
    return sorted(set(files), key=lambda item: str(item))


def snapshot_paths(paths: Iterable[Path]) -> list[dict[str, Any]]:
    records = []
    for path in iter_files(paths):
        stat = path.stat()
        records.append(
            {
                "path": str(path),
                "size": stat.st_size,
                "sha256": sha256_file(path),
            }
        )
    return records


def fingerprint(records: Iterable[dict[str, Any]]) -> str:
    stable = sorted([
        {"path": row["path"], "size": int(row["size"]), "sha256": row["sha256"]}
        for row in records
    ], key=lambda row: row["path"])
    return sha256_text(json.dumps(stable, sort_keys=True, separators=(",", ":")))


def build_record(
    *,
    sample_id: str,
    sample_class: str,
    stage: str,
    status: str,
    inputs: Iterable[Path],
    outputs: Iterable[Path],
    parameter_hash: str,
    software_hash: str,
    run_id: str,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    input_rows = snapshot_paths(inputs)
    output_rows = snapshot_paths(outputs)
    if not output_rows:
        raise ValueError("A completion record must contain at least one output file")
    return {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "sample_class": sample_class,
        "stage": stage,
        "status": status,
        "run_id": run_id,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "parameter_hash": parameter_hash,
        "software_hash": software_hash,
        "input_fingerprint": fingerprint(input_rows),
        "output_fingerprint": fingerprint(output_rows),
        "inputs": input_rows,
        "outputs": output_rows,
        "metadata": metadata or {},
    }


def _check_rows(rows: Any, label: str) -> list[str]:
    problems: list[str] = []
    if not isinstance(rows, list) or not rows:
        return [f"{label} list is empty or malformed"]
    seen: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            problems.append(f"{label} entry is not an object")
            continue
        try:
            path = Path(row["path"]).expanduser().resolve()
            expected_size = int(row["size"])
            expected_sha = str(row["sha256"])
        except (KeyError, TypeError, ValueError):
            problems.append(f"Malformed {label} entry: {row!r}")
            continue
        if str(path) in seen:
            problems.append(f"Duplicate {label} path in record: {path}")
            continue
        seen.add(str(path))
        if not path.is_file():
            problems.append(f"Missing {label}: {path}")
            continue
        size = path.stat().st_size
        if size != expected_size:
            problems.append(f"Changed {label} size: {path} ({size} != {expected_size})")
            continue
        observed_sha = sha256_file(path)
        if observed_sha != expected_sha:
            problems.append(f"Changed {label} checksum: {path}")
    return problems


def validate_record(
    record_path: Path,
    *,
    expected_sample_id: str | None = None,
    expected_parameter_hash: str | None = None,
    expected_software_hash: str | None = None,
    expected_input_fingerprint: str | None = None,
    require_final: bool = False,
) -> tuple[bool, list[str], dict[str, Any] | None]:
    problems: list[str] = []
    try:
        with record_path.open() as handle:
            record = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        return False, [f"Cannot read completion record {record_path}: {error}"], None

    if record.get("schema_version") != SCHEMA_VERSION:
        problems.append(f"Unsupported schema_version: {record.get('schema_version')!r}")
    if expected_sample_id and record.get("sample_id") != expected_sample_id:
        problems.append("sample_id does not match")
    if expected_parameter_hash and record.get("parameter_hash") != expected_parameter_hash:
        problems.append("parameter/configuration hash changed")
    if expected_software_hash and record.get("software_hash") != expected_software_hash:
        problems.append("software/database fingerprint changed")
    if expected_input_fingerprint and record.get("input_fingerprint") != expected_input_fingerprint:
        problems.append("approved input fingerprint changed")
    if require_final:
        if record.get("stage") != "FINAL":
            problems.append("record is not a final completion record")
        if record.get("status") not in CACHE_VALID_FINAL_STATUSES:
            problems.append(f"final status is not cache-valid: {record.get('status')!r}")
    elif record.get("status") == "FAIL":
        problems.append("failed stages are deliberately not cache-valid")

    problems.extend(_check_rows(record.get("inputs"), "input"))
    problems.extend(_check_rows(record.get("outputs"), "output"))
    if not problems:
        if fingerprint(record["inputs"]) != record.get("input_fingerprint"):
            problems.append("stored input fingerprint is inconsistent")
        if fingerprint(record["outputs"]) != record.get("output_fingerprint"):
            problems.append("stored output fingerprint is inconsistent")
    return not problems, problems, record


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    create = subparsers.add_parser("create", help="write a checksummed completion record")
    create.add_argument("--record", required=True, type=Path)
    create.add_argument("--sample-id", required=True)
    create.add_argument("--sample-class", required=True, choices=("isolate", "sweep"))
    create.add_argument("--stage", required=True)
    create.add_argument("--status", required=True)
    create.add_argument("--input", action="append", default=[], type=Path)
    create.add_argument("--output", action="append", default=[], type=Path)
    create.add_argument("--parameter-hash", required=True)
    create.add_argument("--software-hash", required=True)
    create.add_argument("--run-id", required=True)
    create.add_argument("--metadata-json", default="{}")

    validate = subparsers.add_parser("validate", help="validate a completion record and its files")
    validate.add_argument("--record", required=True, type=Path)
    validate.add_argument("--sample-id")
    validate.add_argument("--parameter-hash")
    validate.add_argument("--software-hash")
    validate.add_argument("--input-fingerprint")
    validate.add_argument("--final", action="store_true")
    validate.add_argument("--json-output", type=Path)

    hash_parser = subparsers.add_parser("hash", help="print a deterministic fingerprint for files")
    hash_parser.add_argument("paths", nargs="+", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "hash":
        print(fingerprint(snapshot_paths(args.paths)))
        return 0
    if args.command == "create":
        try:
            metadata = json.loads(args.metadata_json)
            record = build_record(
                sample_id=args.sample_id,
                sample_class=args.sample_class,
                stage=args.stage,
                status=args.status,
                inputs=args.input,
                outputs=args.output,
                parameter_hash=args.parameter_hash,
                software_hash=args.software_hash,
                run_id=args.run_id,
                metadata=metadata,
            )
            atomic_json(args.record, record)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return 2
        return 0

    valid, problems, record = validate_record(
        args.record,
        expected_sample_id=args.sample_id,
        expected_parameter_hash=args.parameter_hash,
        expected_software_hash=args.software_hash,
        expected_input_fingerprint=args.input_fingerprint,
        require_final=args.final,
    )
    payload = {
        "record": str(args.record),
        "sample_id": None if record is None else record.get("sample_id"),
        "valid": valid,
        "problems": problems,
    }
    if args.json_output:
        atomic_json(args.json_output, payload)
    print(json.dumps(payload, sort_keys=True))
    return 0 if valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
