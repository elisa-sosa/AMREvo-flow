import json
import sys
import tempfile
import unittest
from pathlib import Path


BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

from completion_state import atomic_json, build_record, validate_record  # noqa: E402


class CompletionStateTests(unittest.TestCase):
    def test_valid_final_and_corrupt_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "reads.fastq"
            output = root / "result.tsv"
            source.write_text("@r1\nAC\n+\nII\n")
            output.write_text("sample\tvalue\nS1\tPASS\n")
            record_path = root / "FINAL_COMPLETE.json"
            record = build_record(
                sample_id="S1",
                sample_class="isolate",
                stage="FINAL",
                status="COMPLETE",
                inputs=[source],
                outputs=[output],
                parameter_hash="parameters",
                software_hash="software",
                run_id="test",
            )
            atomic_json(record_path, record)

            valid, problems, _ = validate_record(
                record_path,
                expected_sample_id="S1",
                expected_parameter_hash="parameters",
                expected_software_hash="software",
                expected_input_fingerprint=record["input_fingerprint"],
                require_final=True,
            )
            self.assertTrue(valid, problems)

            output.write_text("corrupted\n")
            valid, problems, _ = validate_record(record_path, require_final=True)
            self.assertFalse(valid)
            self.assertTrue(any("output" in problem for problem in problems))

    def test_failed_stage_is_not_cache_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input"
            output = root / "output"
            source.write_text("input")
            output.write_text("output")
            record = build_record(
                sample_id="S1",
                sample_class="isolate",
                stage="ONT_QC",
                status="FAIL",
                inputs=[source],
                outputs=[output],
                parameter_hash="p",
                software_hash="s",
                run_id="run",
            )
            path = root / "stage.json"
            atomic_json(path, record)
            valid, problems, _ = validate_record(path)
            self.assertFalse(valid)
            self.assertIn("failed stages are deliberately not cache-valid", problems)


if __name__ == "__main__":
    unittest.main()
