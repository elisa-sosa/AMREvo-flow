#!/usr/bin/env python3
"""Count variable sites in a FASTA alignment, excluding N, gaps, and ?."""

import argparse
from pathlib import Path


def read_sequences(path: Path) -> list[str]:
    sequences = []
    current = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
            current = []
        elif line:
            current.append(line.upper())
    if current:
        sequences.append("".join(current))
    return sequences


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("alignment", type=Path)
    args = parser.parse_args()
    sequences = read_sequences(args.alignment)
    if sequences and len({len(sequence) for sequence in sequences}) != 1:
        raise SystemExit("Alignment sequences have different lengths")
    count = sum(
        len({base for base in column if base not in {"N", "-", "?"}}) > 1
        for column in zip(*sequences)
    )
    print(count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
