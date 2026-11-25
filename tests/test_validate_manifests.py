import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "bin" / "validate_manifests.py"
sys.path.insert(0, str(ROOT / "bin"))

from completion_state import atomic_json, build_record  # noqa: E402


SAMPLE_FIELDS = [
    "sample_id", "patient_id", "timepoint", "material", "agar_type",
    "expect_illumina", "expect_ont", "active",
]
UNIT_FIELDS = [
    "sample_id", "read_unit_id", "platform", "source_run", "r1", "r2",
    "ont_fastq", "trim_metrics", "decision", "facility_qc_status",
    "merge_order", "comment",
]


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_fastq(path, identifiers, mate=None):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as handle:
        for identifier in identifiers:
            suffix = f"/{mate}" if mate else ""
            handle.write(f"@{identifier}{suffix}\nACGT\n+\nIIII\n")


class ManifestValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.samples = self.root / "samples.tsv"
        self.units = self.root / "units.tsv"
        self.results = self.root / "results"
        self.samples_rows = [
            {"sample_id": "BOTH", "patient_id": "P1", "timepoint": "1", "material": "isolate", "agar_type": "MCA", "expect_illumina": "true", "expect_ont": "true", "active": "true"},
            {"sample_id": "ILL", "patient_id": "P1", "timepoint": "2", "material": "isolate", "agar_type": "MCA", "expect_illumina": "true", "expect_ont": "true", "active": "true"},
            {"sample_id": "ONT", "patient_id": "P2", "timepoint": "1", "material": "isolate", "agar_type": "URI", "expect_illumina": "true", "expect_ont": "true", "active": "true"},
            {"sample_id": "NONE", "patient_id": "P2", "timepoint": "2", "material": "isolate", "agar_type": "MCA", "expect_illumina": "true", "expect_ont": "true", "active": "true"},
        ]
        write_tsv(self.samples, SAMPLE_FIELDS, self.samples_rows)
        self.unit_rows = []
        self._add_illumina("BOTH", "B-I1", ["b1", "b2"])
        self._add_ont("BOTH", "B-O1", ["bo1", "bo2"])
        self._add_illumina("ILL", "I-I1", ["i1", "i2"], survival=82)
        self._add_ont("ONT", "O-O1", ["o1", "o2"])
        write_tsv(self.units, UNIT_FIELDS, self.unit_rows)

    def tearDown(self):
        self.temp.cleanup()

    def _add_illumina(self, sample, unit, identifiers, survival=98, order="1"):
        r1 = self.root / f"{unit}_R1.fastq.gz"
        r2 = self.root / f"{unit}_R2.fastq.gz"
        metric = self.root / f"{unit}.quality_read_trimm_info"
        write_fastq(r1, identifiers, mate=1)
        write_fastq(r2, identifiers, mate=2)
        metric.write_text(f"Both Surviving: 100 ({survival}.0%)\n")
        self.unit_rows.append({
            "sample_id": sample, "read_unit_id": unit, "platform": "illumina",
            "source_run": "run1", "r1": r1, "r2": r2, "ont_fastq": "",
            "trim_metrics": metric, "decision": "approved",
            "facility_qc_status": "accepted", "merge_order": order, "comment": "",
        })

    def _add_ont(self, sample, unit, identifiers, order="1"):
        ont = self.root / f"{unit}.fastq.gz"
        write_fastq(ont, identifiers)
        self.unit_rows.append({
            "sample_id": sample, "read_unit_id": unit, "platform": "ont",
            "source_run": "run2", "r1": "", "r2": "", "ont_fastq": ont,
            "trim_metrics": "", "decision": "approved", "facility_qc_status": "",
            "merge_order": order, "comment": "",
        })

    def run_validator(self, run_id="run1", expect=0):
        output = self.root / f"validation-{run_id}"
        command = [
            sys.executable, str(VALIDATOR), "--sample-class", "isolate",
            "--samples", str(self.samples), "--read-units", str(self.units),
            "--output-dir", str(output), "--allowed-root", str(self.root),
            "--results-root", str(self.results), "--run-id", run_id,
            "--parameter-hash", "parameters", "--software-hash", "software",
            "--fastq-validation", "full",
        ]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, expect, result.stderr)
        return output, result

    def test_both_one_and_neither_routes(self):
        output, _ = self.run_validator()
        physical_lines = (output / "resolved_samples.tsv").read_text().splitlines()
        self.assertIn('["', physical_lines[1])
        self.assertNotIn('"[""', physical_lines[1])
        rows = {row["sample_id"]: row for row in read_tsv(output / "resolved_samples.tsv")}
        self.assertEqual(rows["BOTH"]["status"], "READY_FOR_QC")
        self.assertEqual(rows["ILL"]["status"], "WAITING_FOR_ONT")
        self.assertEqual(rows["ONT"]["status"], "WAITING_FOR_ILLUMINA")
        self.assertEqual(rows["NONE"]["status"], "WAITING_FOR_READS")
        self.assertEqual(rows["ILL"]["facility_qc_status"], "accepted")

    def test_ordered_units_and_duplicate_read_ids(self):
        self._add_illumina("ILL", "I-I2", ["i3", "i4"], order="2")
        write_tsv(self.units, UNIT_FIELDS, self.unit_rows)
        output, _ = self.run_validator()
        rows = {row["sample_id"]: row for row in read_tsv(output / "resolved_samples.tsv")}
        self.assertEqual(len(json.loads(rows["ILL"]["illumina_r1_json"])), 2)

        write_fastq(self.root / "I-I2_R1.fastq.gz", ["i1", "i4"], mate=1)
        write_fastq(self.root / "I-I2_R2.fastq.gz", ["i1", "i4"], mate=2)
        _, result = self.run_validator(run_id="duplicate", expect=2)
        self.assertIn("Duplicate approved read identifier", result.stderr)

    def test_mismatched_mates_fail(self):
        write_fastq(self.root / "B-I1_R2.fastq.gz", ["wrong1", "wrong2"], mate=2)
        _, result = self.run_validator(expect=2)
        self.assertIn("Paired FASTQ identifiers differ", result.stderr)

    def test_facility_rerun_pending_blocks_hybrid_eligibility(self):
        self.unit_rows[0]["facility_qc_status"] = "rerun_pending"
        write_tsv(self.units, UNIT_FIELDS, self.unit_rows)
        output, _ = self.run_validator()
        rows = {row["sample_id"]: row for row in read_tsv(output / "resolved_samples.tsv")}
        self.assertEqual(rows["BOTH"]["status"], "ILLUMINA_FACILITY_RERUN_PENDING")

    def test_new_ont_keeps_illumina_branch_fingerprint(self):
        first, _ = self.run_validator(run_id="illumina-only")
        before = {row["sample_id"]: row for row in read_tsv(first / "resolved_samples.tsv")}["ILL"]
        self._add_ont("ILL", "I-O1", ["io1", "io2"])
        write_tsv(self.units, UNIT_FIELDS, self.unit_rows)
        second, _ = self.run_validator(run_id="ont-arrived")
        after = {row["sample_id"]: row for row in read_tsv(second / "resolved_samples.tsv")}["ILL"]
        self.assertEqual(before["illumina_fingerprint"], after["illumina_fingerprint"])
        self.assertNotEqual(before["input_fingerprint"], after["input_fingerprint"])
        self.assertTrue(after["ont_fingerprint"])

    def test_complete_skip_corruption_and_failed_retry(self):
        first, _ = self.run_validator()
        rows = {row["sample_id"]: row for row in read_tsv(first / "resolved_samples.tsv")}
        sample_dir = self.results / "samples" / "isolates" / "BOTH"
        sample_dir.mkdir(parents=True)
        final_output = sample_dir / "result.tsv"
        final_output.write_text("result\n")
        raw_inputs = []
        for key in ("illumina_r1_json", "illumina_r2_json", "ont_fastq_json", "trim_metrics_json"):
            raw_inputs.extend(Path(value) for value in json.loads(rows["BOTH"][key]) if value)
        record = build_record(
            sample_id="BOTH", sample_class="isolate", stage="FINAL", status="COMPLETE",
            inputs=raw_inputs, outputs=[final_output], parameter_hash="parameters",
            software_hash="software", run_id="first",
        )
        atomic_json(sample_dir / "FINAL_COMPLETE.json", record)
        second, _ = self.run_validator(run_id="second")
        current = {row["sample_id"]: row for row in read_tsv(second / "resolved_samples.tsv")}
        self.assertEqual(current["BOTH"]["status"], "SKIPPED_COMPLETE")

        final_output.write_text("corrupt\n")
        third, _ = self.run_validator(run_id="third")
        current = {row["sample_id"]: row for row in read_tsv(third / "resolved_samples.tsv")}
        self.assertNotEqual(current["BOTH"]["status"], "SKIPPED_COMPLETE")
        self.assertTrue(json.loads(current["BOTH"]["invalid_completion_reasons_json"]))

        ill_dir = self.results / "samples" / "isolates" / "ILL"
        ill_dir.mkdir(parents=True)
        write_tsv(ill_dir / "sample_status.tsv", ["sample_id", "status"], [{"sample_id": "ILL", "status": "BLOCKED_ONT_QC"}])
        fourth, _ = self.run_validator(run_id="retry-token")
        current = {row["sample_id"]: row for row in read_tsv(fourth / "resolved_samples.tsv")}
        self.assertEqual(current["ILL"]["status"], "RETRIED_PREVIOUS_FAILURE")
        self.assertEqual(current["ILL"]["retry_token"], "retry-token")


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


if __name__ == "__main__":
    unittest.main()
