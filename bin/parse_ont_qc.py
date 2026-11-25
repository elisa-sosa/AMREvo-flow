#!/usr/bin/env python3
"""Convert NanoPlot NanoStats into the AMREvoFlow ONT-QC decision."""

import argparse
import re
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--override", required=True)
    parser.add_argument("--min-bases", required=True, type=float)
    parser.add_argument("--min-n50", required=True, type=float)
    parser.add_argument("--min-mean-qual", required=True, type=float)
    args = parser.parse_args()
    text = args.stats.read_text(errors="replace")

    def value(*labels: str) -> float:
        for label in labels:
            match = re.search(rf"(?im)^\s*{label}\s*:\s*([0-9,.]+)", text)
            if match:
                return float(match.group(1).replace(",", ""))
        raise SystemExit(f"Could not parse {labels[0]} from {args.stats}")

    bases = value("Number of bases", "Total bases")
    n50 = value("Read length N50", "N50")
    mean_qual = value("Mean read quality", "Mean quality")
    failures = []
    if bases < args.min_bases:
        failures.append(f"number_of_bases<{args.min_bases:g}")
    if n50 < args.min_n50:
        failures.append(f"n50<{args.min_n50:g}")
    if mean_qual < args.min_mean_qual:
        failures.append(f"mean_qual<{args.min_mean_qual:g}")

    approved_override = args.override.lower() == "true"
    override_applied = bool(failures) and approved_override
    passed = not failures or approved_override
    args.output.write_text(
        "sample_id\tstage\tstatus\tnumber_of_bases\tn50\tmean_qual\toverride_applied\treason\n"
        f"{args.sample_id}\tONT_QC\t{'PASS' if passed else 'FAIL'}\t"
        f"{bases:.0f}\t{n50:.0f}\t{mean_qual:.4g}\t"
        f"{str(override_applied).lower()}\t{';'.join(failures)}\n"
    )
    if passed:
        Path(f"{args.sample_id}.ONT_READY").write_text("PASS\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
