#!/usr/bin/env python3
"""Archive a sample's previous final results when its approved state changed."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-class", required=True, choices=("isolate", "sweep"))
    parser.add_argument("--resolved", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def safe_component(value: str) -> str:
    return "".join(character if character.isalnum() or character in "._-" else "_" for character in value)


def main() -> int:
    args = parse_args()
    sample_group = "isolates" if args.sample_class == "isolate" else "sweeps"
    with args.resolved.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
        fields = list(handle and (rows[0].keys() if rows else []))
    if not rows:
        raise SystemExit("Resolved manifest contains no rows")

    archived: list[dict[str, str]] = []
    for row in rows:
        if row.get("skip_complete") == "true":
            continue
        reasons = json.loads(row.get("invalid_completion_reasons_json") or "[]")
        sample_dir = args.results_root / "samples" / sample_group / row["sample_id"]
        final = sample_dir / "FINAL_COMPLETE.json"
        if not final.exists() or not reasons:
            continue
        try:
            previous = json.loads(final.read_text())
            old_fingerprint = str(previous.get("input_fingerprint") or "unknown")[:16]
        except (OSError, json.JSONDecodeError):
            old_fingerprint = "unreadable-marker"
        history = sample_dir / "history" / f"{safe_component(old_fingerprint)}-{safe_component(args.run_id)}"
        if history.exists():
            raise SystemExit(f"Refusing to overwrite existing history location: {history}")
        history.mkdir(parents=True)
        for child in list(sample_dir.iterdir()):
            if child.name == "history":
                continue
            shutil.move(str(child), str(history / child.name))
        archived.append(
            {
                "sample_id": row["sample_id"],
                "history_path": str(history),
                "reasons": "; ".join(reasons),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.resolved, args.output)
    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "history_path", "reasons"], delimiter="\t")
        writer.writeheader()
        writer.writerows(archived)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
