import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PhylogenyAndCoverageTests(unittest.TestCase):
    def test_reference_selection_and_insufficient_event(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            results = root / "results"
            status = root / "status.tsv"
            rows = []
            specifications = [
                ("S1", "P1", "2", "10", "99", "0.2"),
                ("S2", "P1", "1", "10", "98", "0.1"),
                ("S3", "P1", "1", "10", "99", "0.5"),
                ("S4", "P2", "1", "69", "99", "0.1"),
                ("S5", "P2", "2", "69", "99", "0.1"),
            ]
            for sample_id, patient, timepoint, sequence_type, completeness, contamination in specifications:
                assembly_dir = results / "samples" / "isolates" / sample_id / "02_hybrid_assembly" / "04_pypolca"
                assembly_dir.mkdir(parents=True)
                (assembly_dir / f"{sample_id}.polished.fasta").write_text(">contig\nACGT\n")
                rows.append({
                    "sample_id": sample_id, "patient_id": patient, "timepoint": timepoint,
                    "material": "isolate", "species": "Escherichia coli",
                    "sequence_type": sequence_type, "status": "COMPLETED",
                    "completeness": completeness, "contamination": contamination,
                })
            write_tsv(status, list(rows[0]), rows)
            output = root / "phylogeny"
            result = subprocess.run([
                sys.executable, str(ROOT / "bin" / "prepare_phylogeny.py"),
                "--status", str(status), "--results-root", str(results),
                "--output-dir", str(output), "--minimum-samples", "3",
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            references = {row["group_id"]: row for row in read_tsv(output / "reference_manifest.tsv")}
            self.assertEqual(references["patient_P1__ST10"]["reference_sample_id"], "S3")
            group_status = {row["group_id"]: row for row in read_tsv(output / "phylogeny_group_status.tsv")}
            self.assertEqual(group_status["patient_P2__ST69"]["status"], "INSUFFICIENT_SAMPLES")
            jobs = read_tsv(output / "phylogeny_jobs.tsv")
            self.assertFalse(any(row["group_id"] == "patient_P2__ST69" for row in jobs))

    def test_frozen_reference_warning(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            references = root / "references.tsv"
            groups = root / "groups.txt"
            output = root / "coverage.tsv"
            write_tsv(references, ["reference_sample_id"], [
                {"reference_sample_id": "REF_A"}, {"reference_sample_id": "REF_AB"},
            ])
            groups.write_text("REF_A\tcluster_1\n")
            result = subprocess.run([
                sys.executable, str(ROOT / "bin" / "check_sweep_reference_coverage.py"),
                "--carriage-references", str(references), "--ref-groups", str(groups),
                "--output", str(output),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = {row["reference_id"]: row for row in read_tsv(output)}
            self.assertEqual(rows["REF_A"]["status"], "PRESENT")
            self.assertEqual(rows["REF_AB"]["status"], "ABSENT_FROM_FROZEN_INDEX")


if __name__ == "__main__":
    unittest.main()
