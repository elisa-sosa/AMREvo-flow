#!/usr/bin/env python3
"""Tests for helper programs executed inside analytical task containers."""

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"


class TaskHelperTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def run_helper(self, name: str, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(BIN / name), *map(str, args)],
            cwd=self.root,
            text=True,
            capture_output=True,
            check=True,
        )

    def test_ont_qc_pass_and_fail(self) -> None:
        stats = self.root / "NanoStats.txt"
        stats.write_text(
            "Number of bases: 6,000,000\n"
            "Read length N50: 1,500\n"
            "Mean read quality: 10.2\n"
        )
        output = self.root / "ont.tsv"
        self.run_helper(
            "parse_ont_qc.py", "--stats", stats, "--output", output,
            "--sample-id", "PASS", "--override", "false",
            "--min-bases", "5000000", "--min-n50", "1000",
            "--min-mean-qual", "9",
        )
        self.assertIn("\tPASS\t", output.read_text())
        self.assertTrue((self.root / "PASS.ONT_READY").exists())

        stats.write_text(
            "Number of bases: 4,000,000\nRead length N50: 900\nMean read quality: 8\n"
        )
        self.run_helper(
            "parse_ont_qc.py", "--stats", stats, "--output", output,
            "--sample-id", "FAIL", "--override", "false",
            "--min-bases", "5000000", "--min-n50", "1000",
            "--min-mean-qual", "9",
        )
        self.assertIn("\tFAIL\t", output.read_text())
        self.assertFalse((self.root / "FAIL.ONT_READY").exists())

    def test_checkm2_policy(self) -> None:
        report = self.root / "quality_report.tsv"
        report.write_text("Name\tCompleteness\tContamination\nS\t97\t2\n")
        output = self.root / "checkm2.tsv"
        self.run_helper(
            "parse_checkm2_qc.py", "--report", report, "--output", output,
            "--sample-id", "S", "--override", "false",
            "--pass-completeness", "99", "--min-completeness", "95",
            "--pass-contamination", "1", "--max-contamination", "5",
        )
        self.assertIn("\tWARN\t", output.read_text())
        self.assertEqual((self.root / "S.ASSEMBLY_READY").read_text(), "WARN\n")

    def test_exact_gtdb_ecoli_gate(self) -> None:
        summary = self.root / "summary.tsv"
        summary.write_text(
            "user_genome\tclassification\n"
            "S\td__Bacteria;p__Pseudomonadota;s__Escherichia coli\n"
        )
        output = self.root / "gtdb.tsv"
        self.run_helper(
            "parse_gtdbtk.py", "--summary", summary, "--output", output,
            "--sample-id", "S",
        )
        with output.open(newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["status"], "ECOLI")
        self.assertTrue((self.root / "S.ECOLI_READY").exists())

    def test_variable_site_count(self) -> None:
        alignment = self.root / "core.fa"
        alignment.write_text(">a\nACNT\n>b\nAG-T\n>c\nAC?T\n")
        result = self.run_helper("count_variable_sites.py", alignment)
        self.assertEqual(result.stdout.strip(), "1")


if __name__ == "__main__":
    unittest.main()
