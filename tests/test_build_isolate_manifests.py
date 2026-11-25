import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "build_isolate_manifests.py"


def write_csv(path: Path, fields: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def prepare_metadata(tmp_path: Path) -> tuple[Path, Path]:
    completed = tmp_path / "lims.csv"
    pending = tmp_path / "v3.csv"
    write_csv(
        completed,
        ["sample_id", "PatientNr", "Timepoint", "AgarType"],
        [["P001-MCA--01-01", "001", "-1", "MCA"]],
    )
    write_csv(
        pending,
        ["Unique_ID", "Patient_Nr", "Timepoint", "Agar_Type"],
        [["P014-MCA-02-02", "P014", "2", "MCA"]],
    )
    return completed, pending


def add_pair(root: Path, run: str, sample: str, with_r2: bool = True) -> None:
    trim = root / run / sample / "0_trimming"
    trim.mkdir(parents=True)
    (trim / f"{sample}_R1.fastq.gz").write_bytes(b"r1")
    if with_r2:
        (trim / f"{sample}_R2.fastq.gz").write_bytes(b"r2")
    (trim / f"{sample}.quality_read_trimm_info").write_text("ok\n")


class IsolateManifestBuilderTests(unittest.TestCase):
    def run_builder(self, tmp_path: Path, newest_complete: bool) -> tuple[Path, Path]:
        completed, pending = prepare_metadata(tmp_path)
        illumina = tmp_path / "illumina"
        ont = tmp_path / "ont"
        add_pair(illumina, "esccol_250101", "P001-MCA--01-01")
        add_pair(
            illumina,
            "esccol_250201",
            "P001-MCA--01-01",
            with_r2=newest_complete,
        )
        ont.mkdir()
        if newest_complete:
            (ont / "P001-MCA--01-01.fastq").write_text("@r\nA\n+\nI\n")
        samples = tmp_path / "samples.tsv"
        reads = tmp_path / "reads.tsv"
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--completed-lims",
                str(completed),
                "--pending-v3",
                str(pending),
                "--illumina-root",
                str(illumina),
                "--ont-root",
                str(ont),
                "--samples-output",
                str(samples),
                "--read-units-output",
                str(reads),
            ],
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return samples, reads

    def test_latest_complete_run_is_approved(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            samples, reads = self.run_builder(Path(directory), newest_complete=True)
            sample_rows = read_tsv(samples)
            self.assertEqual(
                {row["sample_id"] for row in sample_rows},
                {"P001-MCA--01-01", "P014-MCA-02-02"},
            )
            illumina_rows = [
                row for row in read_tsv(reads) if row["platform"] == "illumina"
            ]
            by_run = {row["source_run"]: row for row in illumina_rows}
            self.assertEqual(by_run["esccol_250101"]["decision"], "rejected")
            self.assertEqual(by_run["esccol_250201"]["decision"], "approved")
            ont_row = next(row for row in read_tsv(reads) if row["platform"] == "ont")
            self.assertEqual(ont_row["decision"], "approved")

    def test_latest_incomplete_run_does_not_fall_back(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            _, reads = self.run_builder(Path(directory), newest_complete=False)
            by_run = {
                row["source_run"]: row
                for row in read_tsv(reads)
                if row["platform"] == "illumina"
            }
            self.assertEqual(by_run["esccol_250101"]["decision"], "rejected")
            self.assertEqual(by_run["esccol_250201"]["decision"], "candidate")
            self.assertIn(
                "LATEST_DATED_RUN_INCOMPLETE_NO_FALLBACK",
                by_run["esccol_250201"]["comment"],
            )
