#!/usr/bin/env python3
"""Regression checks for TSV parsing at the Nextflow workflow boundary."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class NextflowTsvContractTests(unittest.TestCase):
    def test_prepared_sample_tsv_readers_enable_csv_quoting(self) -> None:
        """JSON arrays written by csv.DictWriter must remain complete cells."""
        workflow = (ROOT / "main.nf").read_text()
        parser = ".splitCsv(header: true, sep: '\\t', quote: '\"')"
        self.assertEqual(workflow.count(parser), 2)


if __name__ == "__main__":
    unittest.main()
