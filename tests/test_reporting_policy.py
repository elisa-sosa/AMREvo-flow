import base64
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from build_reports import classify  # noqa: E402


def row(**changes):
    base = {
        "skip_complete": "false", "active": "true", "material": "isolate",
        "illumina_r1_json": '["r1"]', "ont_fastq_json": '["ont"]',
        "facility_qc_status": "accepted", "illumina_qc_override": "false",
        "assembly_qc_override": "false",
    }
    base.update(changes)
    return base


class ReportingPolicyTests(unittest.TestCase):
    def test_ont_failure_blocks(self):
        status, _ = classify(row(), {"ILLUMINA_QC": {"status": "PASS"}, "ONT_QC": {"status": "FAIL", "reason": "n50<1000"}})
        self.assertEqual(status, "BLOCKED_ONT_QC")

    def test_non_ecoli_is_complete(self):
        stages = {
            "ILLUMINA_QC": {"status": "PASS"}, "ONT_QC": {"status": "PASS"},
            "ASSEMBLY_QC": {"status": "WARN"}, "GTDBTK": {"status": "NON_ECOLI"},
        }
        status, reason = classify(row(), stages)
        self.assertEqual(status, "NON_ECOLI_COMPLETE")
        self.assertIn("NOT_APPLICABLE", reason)

    def test_ecoli_requires_all_specific_stages(self):
        stages = {
            "ILLUMINA_QC": {"status": "PASS"}, "ONT_QC": {"status": "PASS"},
            "ASSEMBLY_QC": {"status": "PASS"}, "GTDBTK": {"status": "ECOLI"},
            "MLST_ACHTMAN": {"status": "PASS"}, "AMRFINDERPLUS": {"status": "PASS"},
            "KTYPR": {"status": "PASS"},
        }
        self.assertEqual(classify(row(), stages)[0], "COMPLETED")
        del stages["KTYPR"]
        self.assertEqual(classify(row(), stages)[0], "PROCESS_FAILURE")


class IlluminaPolicyHelperTests(unittest.TestCase):
    def run_helper(self, facility_status, *, override=False):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        metrics = root / "trim.txt"
        output = root / "qc.tsv"
        marker = root / "ready"
        metrics.write_text("Both Surviving: 100 (82.0%)\n")
        command = [
            sys.executable,
            str(BIN / "assess_illumina_qc.py"),
            "--sample-id", "S1",
            "--stage", "ILLUMINA_QC",
            "--trim-metrics", str(metrics),
            "--facility-status", facility_status,
            "--output", str(output),
            "--ready-marker", str(marker),
        ]
        if override:
            command.append("--override")
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open(newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        return temporary, row, marker.exists()

    def test_facility_accepted_strong_review_is_usable(self):
        temporary, result, marker_exists = self.run_helper("accepted")
        self.addCleanup(temporary.cleanup)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["metric_flag"], "STRONG_REVIEW")
        self.assertTrue(marker_exists)

    def test_facility_rerun_pending_blocks_without_override(self):
        temporary, result, marker_exists = self.run_helper("rerun_pending")
        self.addCleanup(temporary.cleanup)
        self.assertEqual(result["status"], "FAIL")
        self.assertFalse(marker_exists)


class ReportArtifactTests(unittest.TestCase):
    def test_compact_snapshot_and_retry_count(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            results = root / "results"
            output = results / "cohort" / "current"
            snapshot = results / "cohort" / "snapshots" / "test-run"
            artifacts = root / "artifacts"
            artifacts.mkdir()
            (artifacts / "isolate.input_checksums.tsv").write_text(
                "sample_id\tpath\tsize\tsha256\n"
            )
            (artifacts / "isolate.validation_summary.json").write_text("{}\n")
            resolved = root / "isolate.prepared_samples.tsv"
            fields = [
                "sample_id", "patient_id", "timepoint", "material", "agar_type",
                "active", "skip_complete", "illumina_r1_json", "illumina_r2_json",
                "ont_fastq_json", "trim_metrics_json", "facility_qc_status",
                "illumina_qc_override", "assembly_qc_override", "prior_status",
                "parameter_hash", "software_hash",
            ]
            with resolved.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow({
                    "sample_id": "PENDING", "patient_id": "P1", "timepoint": "1",
                    "material": "isolate", "agar_type": "MCA", "active": "true",
                    "skip_complete": "false", "illumina_r1_json": "[]",
                    "illumina_r2_json": "[]", "ont_fastq_json": "[]",
                    "trim_metrics_json": "[]", "facility_qc_status": "",
                    "illumina_qc_override": "false", "assembly_qc_override": "false",
                    "prior_status": "BLOCKED_ONT_QC", "parameter_hash": "p",
                    "software_hash": "s",
                })
            configuration = base64.b64encode(
                json.dumps({"parameters": {"ont_min_bases": 5_000_000}}).encode()
            ).decode()
            result = subprocess.run([
                sys.executable, str(BIN / "build_reports.py"),
                "--resolved", str(resolved),
                "--artifact-dir", str(artifacts),
                "--results-root", str(results),
                "--output-dir", str(output),
                "--snapshot-dir", str(snapshot),
                "--run-id", "test-run",
                "--configuration-base64", configuration,
                "--no-plots",
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("RETRIED_PREVIOUS_FAILURE\t1", (output / "run_summary.txt").read_text())
            self.assertTrue((output / "software_versions.tsv").is_file())
            self.assertTrue((snapshot / "manifests" / resolved.name).is_file())
            self.assertTrue((snapshot / "provenance" / "run_configuration.json").is_file())
            self.assertTrue((snapshot / "provenance" / "isolate.input_checksums.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
