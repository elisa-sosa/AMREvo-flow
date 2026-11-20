#!/usr/bin/env python3
"""Apply the AMREvoFlow Illumina QC policy to one accepted read pair."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def paired_survival_percentage(path: Path) -> float | None:
    """Return the paired-read survival percentage when a supplied metric has it."""
    if not path.is_file() or path.stat().st_size == 0:
        return None
    text = path.read_text(errors="replace")
    patterns = (
        r"Both Surviving:\s*[0-9,]+\s*\(([0-9.]+)%\)",
        r"both[_ ]surviving(?:[_ ]percentage)?[\t,: ]+([0-9.]+)",
        r"percentage[\t,: ]+([0-9.]+)",
    )
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return float(match.group(1))
    return None


def metric_flag(percentage: float | None) -> str:
    if percentage is None:
        return "NOT_AVAILABLE"
    if percentage >= 95:
        return "GOOD"
    if percentage >= 90:
        return "REVIEW"
    return "STRONG_REVIEW"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--trim-metrics", required=True, type=Path)
    parser.add_argument("--facility-status", default="accepted")
    parser.add_argument("--override", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ready-marker", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    facility = args.facility_status.strip().lower() or "accepted"
    percentage = paired_survival_percentage(args.trim_metrics)
    accepted = facility == "accepted" or args.override
    override_applied = args.override and facility != "accepted"

    if override_applied:
        reason = f"override:facility_qc_status={facility}"
    elif accepted:
        reason = ""
    else:
        reason = f"facility_qc_status={facility}"

    args.output.write_text(
        "sample_id\tstage\tstatus\tfacility_qc_status\tpaired_survival_pct\t"
        "metric_flag\toverride_applied\treason\n"
        f"{args.sample_id}\t{args.stage}\t{'PASS' if accepted else 'FAIL'}\t"
        f"{facility}\t{'' if percentage is None else f'{percentage:.6g}'}\t"
        f"{metric_flag(percentage)}\t{str(override_applied).lower()}\t{reason}\n"
    )
    if accepted:
        args.ready_marker.write_text("PASS\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
