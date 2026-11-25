#!/usr/bin/env python3
"""Convert a CheckM2 report into the AMREvoFlow assembly-QC decision."""

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--override", required=True)
    parser.add_argument("--pass-completeness", required=True, type=float)
    parser.add_argument("--min-completeness", required=True, type=float)
    parser.add_argument("--pass-contamination", required=True, type=float)
    parser.add_argument("--max-contamination", required=True, type=float)
    args = parser.parse_args()
    with args.report.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    completeness = float(row["Completeness"])
    contamination = float(row["Contamination"])
    if completeness < args.min_completeness or contamination >= args.max_contamination:
        qc = "FAIL"
    elif completeness < args.pass_completeness or contamination > args.pass_contamination:
        qc = "WARN"
    else:
        qc = "PASS"
    approved_override = args.override.lower() == "true"
    override_applied = qc == "FAIL" and approved_override
    reasons = []
    if completeness < args.min_completeness:
        reasons.append(f"completeness<{args.min_completeness:g}")
    if contamination >= args.max_contamination:
        reasons.append(f"contamination>={args.max_contamination:g}")
    args.output.write_text(
        "sample_id\tstage\tstatus\tcompleteness\tcontamination\toverride_applied\treason\n"
        f"{args.sample_id}\tASSEMBLY_QC\t{qc}\t{completeness:.6g}\t"
        f"{contamination:.6g}\t{str(override_applied).lower()}\t{';'.join(reasons)}\n"
    )
    if qc != "FAIL" or approved_override:
        Path(f"{args.sample_id}.ASSEMBLY_READY").write_text(f"{qc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
