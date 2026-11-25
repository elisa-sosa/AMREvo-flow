import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "build_sweep_manifests.py"


class SweepManifestBuilderTests(unittest.TestCase):
    def test_only_reconstructed_all_paths_are_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            all_dir = root / "sweeps_fastq" / "all"
            all_dir.mkdir(parents=True)
            (all_dir / "P004-MCA--01_R1.fastq.gz").write_bytes(b"r1")
            (all_dir / "P004-MCA--01_R2.fastq.gz").write_bytes(b"r2")
            samples = root / "sweeps.samples.local.tsv"
            reads = root / "sweeps.read_units.local.tsv"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--all-root",
                    str(all_dir),
                    "--samples-output",
                    str(samples),
                    "--read-units-output",
                    str(reads),
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with samples.open(newline="") as handle:
                sample_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(sample_row["timepoint"], "-1")
            with reads.open(newline="") as handle:
                read_row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(read_row["decision"], "approved")
            self.assertEqual(read_row["source_run"], "reconstructed_all")

    def test_legacy_cat_path_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            legacy = root / "all_fastq_files_withCATunclearWhatThisIs"
            legacy.mkdir()
            manifest = root / "manifest.csv"
            with manifest.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["ID", "R1", "R2"])
                writer.writerow(
                    [
                        "P004-MCA--01",
                        legacy / "P004-MCA--01_R1.fastq.gz",
                        legacy / "P004-MCA--01_R2.fastq.gz",
                    ]
                )
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--production-manifest",
                    str(manifest),
                    "--samples-output",
                    str(root / "samples.tsv"),
                    "--read-units-output",
                    str(root / "reads.tsv"),
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("sweeps_fastq/all", result.stderr)
