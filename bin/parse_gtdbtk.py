#!/usr/bin/env python3
"""Extract the exact GTDB species label and AMREvoFlow E. coli gate."""

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    args = parser.parse_args()
    with args.summary.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    classification = row.get("classification", "")
    species = next(
        (rank[3:] for rank in classification.split(";") if rank.startswith("s__")),
        "",
    )
    status = "ECOLI" if species == "Escherichia coli" else "NON_ECOLI"
    args.output.write_text(
        "sample_id\tstage\tstatus\tspecies\tclassification\n"
        f"{args.sample_id}\tGTDBTK\t{status}\t{species}\t{classification}\n"
    )
    if status == "ECOLI":
        Path(f"{args.sample_id}.ECOLI_READY").write_text("Escherichia coli\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
