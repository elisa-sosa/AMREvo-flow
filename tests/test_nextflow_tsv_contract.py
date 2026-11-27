#!/usr/bin/env python3
"""Regression checks for TSV parsing at the Nextflow workflow boundary."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class NextflowTsvContractTests(unittest.TestCase):
    def test_prepared_sample_tsv_readers_use_unquoted_tsv(self) -> None:
        """Resolved JSON cells are deliberately parsed without CSV quoting."""
        workflow = (ROOT / "main.nf").read_text()
        parser = ".splitCsv(header: true, sep: '\\t')"
        self.assertEqual(workflow.count(parser), 2)

    def test_retry_overwrites_run_scoped_technical_reports(self) -> None:
        config = (ROOT / "nextflow.config").read_text()
        self.assertEqual(config.count("overwrite = true"), 4)

    def test_reporting_container_supports_nextflow_metrics(self) -> None:
        definition = (ROOT / "containers" / "reporting-multiqc-1.33.def").read_text()
        cluster = (ROOT / "conf" / "slurm.config").read_text()
        submitter = (ROOT / "submit_pilot.sh").read_text()
        image = "reporting_multiqc-1.33_matplotlib-3.11.2_procps.sif"
        self.assertIn("procps", definition)
        self.assertIn("command -v ps", definition)
        self.assertEqual(cluster.count(image), 2)
        self.assertIn(image, submitter)


if __name__ == "__main__":
    unittest.main()
