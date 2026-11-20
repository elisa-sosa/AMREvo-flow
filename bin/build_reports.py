#!/usr/bin/env python3
"""Build the canonical AMREvoFlow status table, summaries and PNG dashboard."""

from __future__ import annotations

import argparse
import base64
import csv
import json
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from completion_state import atomic_json, build_record


TERMINAL_ORDER = [
    "COMPLETED",
    "SKIPPED_COMPLETE",
    "RETRIED_PREVIOUS_FAILURE",
    "WAITING_FOR_ILLUMINA",
    "WAITING_FOR_ONT",
    "WAITING_FOR_READS",
    "ILLUMINA_FACILITY_RERUN_PENDING",
    "BLOCKED_ONT_QC",
    "BLOCKED_ASSEMBLY_QC",
    "NON_ECOLI_COMPLETE",
    "PROCESS_FAILURE",
    "PROVISIONAL_SWEEP_COMPLETE",
    "READY_FOR_GEMSWEEP",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolved", action="append", required=True, type=Path)
    parser.add_argument("--artifact-dir", type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--snapshot-dir", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--configuration-base64", required=True)
    parser.add_argument("--no-plots", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def iter_artifacts(artifact_dir: Path | None, results_root: Path) -> Iterable[Path]:
    if artifact_dir and artifact_dir.exists():
        yield from (path for path in artifact_dir.rglob("*") if path.is_file())
    for group in ("isolates", "sweeps"):
        root = results_root / "samples" / group
        if root.exists():
            yield from (path for path in root.rglob("*.tsv") if "history" not in path.parts)


def parse_stage_artifacts(
    paths: Iterable[Path],
) -> tuple[dict[str, dict[str, dict[str, str]]], dict[str, str], dict[str, dict[str, str]]]:
    stages: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    sequence_types: dict[str, str] = {}
    quast_metrics: dict[str, dict[str, str]] = {}
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen or path.stat().st_size == 0:
            continue
        seen.add(resolved)
        if path.name.endswith(".mlst.tsv") and not path.name.endswith(".stage.tsv"):
            try:
                fields = path.read_text(errors="replace").splitlines()[0].split("\t")
                if len(fields) >= 3:
                    sample_id = path.name[: -len(".mlst.tsv")]
                    sequence_types[sample_id] = fields[2]
            except (OSError, IndexError):
                pass
            continue
        if path.name.endswith(".quast.tsv") and not path.name.endswith(".stage.tsv"):
            sample_id = path.name[: -len(".quast.tsv")]
            observed: dict[str, str] = {}
            try:
                for line in path.read_text(errors="replace").splitlines():
                    cells = line.split("\t")
                    if len(cells) >= 2:
                        observed[cells[0].strip()] = cells[1].strip()
            except OSError:
                observed = {}
            quast_metrics[sample_id] = {
                "assembly_contigs": observed.get("# contigs (>= 0 bp)", observed.get("# contigs", "")),
                "assembly_n50": observed.get("N50", ""),
                "assembly_total_length": observed.get("Total length (>= 0 bp)", observed.get("Total length", "")),
            }
            continue
        try:
            rows = read_tsv(path)
        except (OSError, csv.Error, UnicodeDecodeError):
            continue
        for row in rows:
            sample_id = (row.get("sample_id") or "").strip()
            stage = (row.get("stage") or "").strip()
            if sample_id and stage and row.get("status"):
                stages[sample_id][stage] = row
    return stages, sequence_types, quast_metrics


def availability(row: dict[str, str]) -> tuple[bool, bool]:
    return bool(json.loads(row["illumina_r1_json"])), bool(json.loads(row["ont_fastq_json"]))


def classify(row: dict[str, str], stages: dict[str, dict[str, str]]) -> tuple[str, str]:
    if row["skip_complete"] == "true":
        return "SKIPPED_COMPLETE", "validated FINAL_COMPLETE.json reused"
    if row["active"] != "true":
        return "INACTIVE", "sample marked inactive"
    has_illumina, has_ont = availability(row)
    facility = row["facility_qc_status"]
    if facility in {"rerun_pending", "rejected"} and row["illumina_qc_override"] != "true":
        return "ILLUMINA_FACILITY_RERUN_PENDING", f"facility_qc_status={facility}"

    if row["material"] == "sweep":
        if not has_illumina:
            return "WAITING_FOR_READS", "no approved sweep read pair"
        qc = stages.get("SWEEP_ILLUMINA_QC", {}).get("status")
        if qc == "FAIL":
            return "ILLUMINA_FACILITY_RERUN_PENDING", stages["SWEEP_ILLUMINA_QC"].get("reason", "")
        if qc == "PASS":
            marker = Path(row.get("_sample_dir", "")) / "FINAL_COMPLETE.json"
            return ("PROVISIONAL_SWEEP_COMPLETE", "validated GEMSWEEP outputs") if marker.is_file() else (
                "READY_FOR_GEMSWEEP",
                "read QC passed; separate pinned GEMSWEEP invocation is pending",
            )
        return "PROCESS_FAILURE", "sweep QC did not produce a terminal stage record"

    if not has_illumina and not has_ont:
        return "WAITING_FOR_READS", "no approved reads"
    illumina_qc = stages.get("ILLUMINA_QC", {}).get("status")
    ont_qc = stages.get("ONT_QC", {}).get("status")
    if has_illumina and illumina_qc != "PASS":
        return "PROCESS_FAILURE", "Illumina QC process did not finish"
    if has_ont and ont_qc == "FAIL":
        return "BLOCKED_ONT_QC", stages["ONT_QC"].get("reason", "ONT thresholds failed")
    if has_ont and ont_qc != "PASS":
        return "PROCESS_FAILURE", "ONT QC process did not finish"
    if has_illumina and not has_ont:
        return "WAITING_FOR_ONT", "valid Illumina QC retained"
    if has_ont and not has_illumina:
        return "WAITING_FOR_ILLUMINA", "valid ONT QC retained"

    assembly_qc = stages.get("ASSEMBLY_QC", {}).get("status")
    if assembly_qc == "FAIL" and row["assembly_qc_override"] != "true":
        return "BLOCKED_ASSEMBLY_QC", stages["ASSEMBLY_QC"].get("reason", "assembly QC failed")
    if assembly_qc not in {"PASS", "WARN", "FAIL"}:
        return "PROCESS_FAILURE", "hybrid assembly or CheckM2 did not finish"
    gtdb = stages.get("GTDBTK", {}).get("status")
    if gtdb == "NON_ECOLI":
        return "NON_ECOLI_COMPLETE", "E. coli-specific stages are NOT_APPLICABLE"
    if gtdb != "ECOLI":
        return "PROCESS_FAILURE", "GTDB-Tk did not produce an exact species outcome"
    required = ("MLST_ACHTMAN", "AMRFINDERPLUS", "KTYPR")
    missing = [stage for stage in required if stages.get(stage, {}).get("status") != "PASS"]
    if missing:
        return "PROCESS_FAILURE", "missing required E. coli stages: " + ",".join(missing)
    return "COMPLETED", "all required isolate stages validated"


def suggestion(status: str, row: dict[str, str], stages: dict[str, dict[str, str]]) -> str:
    if status == "ILLUMINA_FACILITY_RERUN_PENDING":
        return "ILLUMINA_FACILITY_RERUN_PENDING"
    if status == "BLOCKED_ONT_QC":
        return "RESEQUENCE_ONT_REVIEW"
    if status == "BLOCKED_ASSEMBLY_QC":
        illumina = stages.get("ILLUMINA_QC", {}).get("status")
        ont = stages.get("ONT_QC", {}).get("status")
        return "ASSEMBLY_OR_SOURCE_READS_REVIEW" if illumina == ont == "PASS" else "REEXTRACTION_OR_BOTH_LIBRARIES_REVIEW"
    if status == "PROCESS_FAILURE":
        return "TECHNICAL_RECOVERY_REVIEW"
    return ""


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def collect_versions(sample_dir: Path) -> Path:
    output = sample_dir / "software_versions.txt"
    sources = sorted(
        path for path in sample_dir.rglob("*.software_versions.txt")
        if "history" not in path.parts and path.resolve() != output.resolve()
    )
    with output.open("w") as handle:
        for source in sources:
            handle.write(f"## {source.relative_to(sample_dir)}\n")
            handle.write(source.read_text(errors="replace"))
            if not source.read_text(errors="replace").endswith("\n"):
                handle.write("\n")
    return output


def collect_cohort_versions(results_root: Path, output: Path) -> None:
    rows: list[dict[str, str]] = []
    for group, material in (("isolates", "isolate"), ("sweeps", "sweep")):
        root = results_root / "samples" / group
        if not root.exists():
            continue
        for sample_dir in sorted(path for path in root.iterdir() if path.is_dir()):
            for source in sorted(sample_dir.rglob("*.software_versions.txt")):
                if "history" in source.parts:
                    continue
                for line in source.read_text(errors="replace").splitlines():
                    if not line or line.startswith("##") or "\t" not in line:
                        continue
                    tool, value = line.split("\t", 1)
                    rows.append(
                        {
                            "sample_id": sample_dir.name,
                            "material": material,
                            "source": str(source.relative_to(sample_dir)),
                            "tool_or_key": tool,
                            "version_or_value": value,
                        }
                    )
    write_tsv(
        output,
        rows,
        ["sample_id", "material", "source", "tool_or_key", "version_or_value"],
    )


def write_provenance(
    artifact_dir: Path | None,
    output_dir: Path,
    results_root: Path,
    configuration_base64: str,
    run_id: str,
) -> Path:
    provenance = output_dir / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    try:
        configuration = json.loads(base64.b64decode(configuration_base64, validate=True))
    except (ValueError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid encoded run configuration: {error}") from error
    configuration["run_id"] = run_id
    atomic_json(provenance / "run_configuration.json", configuration)

    wanted = (".input_checksums.tsv", ".validation_summary.json", ".archive_report.tsv")
    copied: dict[str, Path] = {}
    if artifact_dir and artifact_dir.exists():
        for source in sorted(path for path in artifact_dir.rglob("*") if path.is_file()):
            if not source.name.endswith(wanted):
                continue
            destination = provenance / source.name
            previous = copied.get(source.name)
            if previous and previous.read_bytes() != source.read_bytes():
                raise ValueError(f"Conflicting provenance artifacts named {source.name}")
            shutil.copy2(source, destination)
            copied[source.name] = source
    coverage = results_root / "cohort" / "current" / "sweeps" / "frozen_reference_coverage.tsv"
    if coverage.is_file():
        shutil.copy2(coverage, provenance / coverage.name)
    return provenance


def required_outputs(sample_dir: Path, status: str, sample_id: str) -> list[Path]:
    # sample_status.tsv intentionally remains outside the final checksum set:
    # its operational status changes from COMPLETED to SKIPPED_COMPLETE on a
    # valid reuse, while the scientific products remain identical.
    candidates = [sample_dir / "software_versions.txt"]
    patterns = [
        "01_read_qc/**/*.tsv",
        "02_hybrid_assembly/04_pypolca/*.polished.fasta",
        "03_assembly_qc/quast/*.quast.tsv",
        "03_assembly_qc/checkm2/*.checkm2.tsv",
        "03_assembly_qc/gtdbtk/*.gtdbtk.tsv",
    ]
    if status == "COMPLETED":
        patterns.extend(["04_typing/mlst/*.mlst.tsv", "05_amr/amrfinderplus/*.tsv", "06_capsule/ktypr/*.ktypr.tsv"])
    for pattern in patterns:
        matches = sorted(sample_dir.glob(pattern))
        if not matches:
            raise FileNotFoundError(f"{sample_id}: required final output pattern is absent: {pattern}")
        candidates.extend(matches)
    return candidates


def create_final_record(row: dict[str, str], sample_dir: Path, status: str, run_id: str) -> None:
    final_status = "COMPLETE" if status == "COMPLETED" else status
    inputs = [Path(value) for key in ("illumina_r1_json", "illumina_r2_json", "ont_fastq_json") for value in json.loads(row[key])]
    inputs.extend(Path(value) for value in json.loads(row["trim_metrics_json"]) if value)
    outputs = required_outputs(sample_dir, status, row["sample_id"])
    record = build_record(
        sample_id=row["sample_id"],
        sample_class=row["material"],
        stage="FINAL",
        status=final_status,
        inputs=inputs,
        outputs=outputs,
        parameter_hash=row["parameter_hash"],
        software_hash=row["software_hash"],
        run_id=run_id,
        metadata={"workflow_status": status},
    )
    atomic_json(sample_dir / "FINAL_COMPLETE.json", record)


def no_data(ax: Any, title: str) -> None:
    ax.text(0.5, 0.5, "No eligible data", ha="center", va="center", transform=ax.transAxes)
    ax.set_title(title)
    ax.set_axis_off()


def render_plots(rows: list[dict[str, str]], plot_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir.mkdir(parents=True, exist_ok=True)
    palette = plt.get_cmap("tab20")

    def bar_counts(path: str, title: str, labels: Iterable[str], color: str = "#3b82f6") -> None:
        counts = Counter(labels)
        fig, ax = plt.subplots(figsize=(11, 6))
        if counts:
            names, values = zip(*sorted(counts.items(), key=lambda item: (-item[1], item[0])))
            ax.bar(range(len(names)), values, color=color)
            ax.set_xticks(range(len(names)), names, rotation=45, ha="right")
            ax.set_ylabel("Samples")
            ax.set_title(title)
            for index, value in enumerate(values):
                ax.text(index, value, str(value), ha="center", va="bottom", fontsize=8)
        else:
            no_data(ax, title)
        fig.tight_layout()
        fig.savefig(plot_dir / path, dpi=180)
        plt.close(fig)

    bar_counts("pipeline_status.png", "Samples by current pipeline status", (row["status"] for row in rows))
    bar_counts("availability.png", "Approved read availability", (row["availability"] for row in rows), "#10b981")
    operational = [row["status"] for row in rows if row["status"] not in {"COMPLETED", "NON_ECOLI_COMPLETE"}]
    operational.extend(
        "RETRIED_PREVIOUS_FAILURE"
        for row in rows
        if row.get("retried_previous_failure") == "true"
    )
    bar_counts(
        "pending_blocked_failed_retried_skipped.png",
        "Pending, blocked, failed, retried and reused samples",
        operational,
        "#f59e0b",
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    metric_specs = (("ont_bases", "ONT total bases"), ("ont_n50", "ONT read N50"), ("ont_mean_q", "ONT mean Q-score"))
    for ax, (field, title) in zip(axes, metric_specs):
        values = [float(row[field]) for row in rows if row.get(field)]
        if values:
            ax.hist(values, bins=min(30, max(5, len(values) // 5)), color="#6366f1")
            ax.set_title(title)
        else:
            no_data(ax, title)
    fig.tight_layout()
    fig.savefig(plot_dir / "qc_distributions_failures.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6))
    points = [(float(row["completeness"]), float(row["contamination"]), row["status"]) for row in rows if row.get("completeness") and row.get("contamination")]
    if points:
        for status in sorted({point[2] for point in points}):
            selected = [point for point in points if point[2] == status]
            ax.scatter([point[0] for point in selected], [point[1] for point in selected], label=status, alpha=0.75)
        ax.axvline(95, color="grey", linestyle="--", linewidth=1)
        ax.axvline(99, color="grey", linestyle=":", linewidth=1)
        ax.axhline(1, color="grey", linestyle=":", linewidth=1)
        ax.axhline(5, color="grey", linestyle="--", linewidth=1)
        ax.set_xlabel("CheckM2 completeness (%)")
        ax.set_ylabel("CheckM2 contamination (%)")
        ax.legend(fontsize=7)
        ax.set_title("Hybrid assembly quality")
    else:
        no_data(ax, "Hybrid assembly quality")
    fig.tight_layout()
    fig.savefig(plot_dir / "assembly_quality.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    contiguity_specs = (("assembly_n50", "Assembly N50 (bp)"), ("assembly_contigs", "Assembly contig count"))
    for ax, (field, title) in zip(axes, contiguity_specs):
        values = [float(row[field].replace(",", "")) for row in rows if row.get(field)]
        if values:
            ax.hist(values, bins=min(30, max(5, len(values) // 5)), color="#2563eb")
            ax.set_title(title)
            ax.set_ylabel("Samples")
        else:
            no_data(ax, title)
    fig.tight_layout()
    fig.savefig(plot_dir / "assembly_contiguity.png", dpi=180)
    plt.close(fig)

    species_st = [f"{row.get('species') or 'unresolved'} | ST{row.get('sequence_type') or 'NA'}" for row in rows if row.get("species") or row.get("sequence_type")]
    bar_counts("species_st_summary.png", "GTDB species and Achtman sequence types", species_st, "#8b5cf6")

    st_rows = [row for row in rows if row.get("sequence_type")]
    bar_counts(
        "figure6_hybrid_st_frequencies.png",
        "Figure 6. Achtman sequence-type frequencies among hybrid assemblies",
        (f"ST{row['sequence_type']}" for row in st_rows),
        "#0f766e",
    )

    fig, ax = plt.subplots(figsize=(max(12, len({row['timepoint'] for row in st_rows}) * 1.2), max(5, len({row['patient_id'] for row in st_rows}) * 0.45)))
    if st_rows:
        patients = sorted({row["patient_id"] for row in st_rows})
        timepoints = sorted({row["timepoint"] for row in st_rows})
        sts = sorted({row["sequence_type"] for row in st_rows})
        colors = {st: palette(index % 20) for index, st in enumerate(sts)}
        grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in sorted(st_rows, key=lambda item: item["sample_id"]):
            grouped[(row["patient_id"], row["timepoint"])].append(row)
        for y, patient in enumerate(patients):
            for x, timepoint in enumerate(timepoints):
                for slot, row in enumerate(grouped[(patient, timepoint)][:10]):
                    x_offset = (slot % 5 - 2) * 0.10
                    y_offset = (slot // 5 - 0.5) * 0.20
                    ax.scatter(x + x_offset, y + y_offset, s=45, color=colors[row["sequence_type"]], edgecolor="black", linewidth=0.2)
        ax.set_xticks(range(len(timepoints)), timepoints, rotation=45, ha="right")
        ax.set_yticks(range(len(patients)), patients)
        ax.set_xlabel("Timepoint")
        ax.set_ylabel("Patient")
        ax.set_title("Figure 7. Condensed longitudinal Achtman-ST composition (up to 10 isolates/timepoint)")
        handles = [plt.Line2D([0], [0], marker="o", linestyle="", color=colors[st], label=f"ST{st}") for st in sts[:20]]
        ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        ax.grid(axis="x", color="#e5e7eb", linewidth=0.5)
    else:
        no_data(ax, "Figure 7. Condensed longitudinal Achtman-ST composition")
    fig.tight_layout()
    fig.savefig(plot_dir / "figure7_longitudinal_st_composition.png", dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.snapshot_dir.mkdir(parents=True, exist_ok=True)
    resolved_rows = [row for path in args.resolved for row in read_tsv(path)]
    if not resolved_rows:
        raise SystemExit("No resolved samples were supplied")
    stages_by_sample, sequence_types, quast_by_sample = parse_stage_artifacts(
        iter_artifacts(args.artifact_dir, args.results_root)
    )
    output_rows: list[dict[str, str]] = []
    finalize_errors: list[str] = []

    for row in resolved_rows:
        sample_id = row["sample_id"]
        group = "isolates" if row["material"] == "isolate" else "sweeps"
        sample_dir = args.results_root / "samples" / group / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        row["_sample_dir"] = str(sample_dir)
        stages = stages_by_sample.get(sample_id, {})
        status, reason = classify(row, stages)
        has_illumina, has_ont = availability(row)
        availability_label = "BOTH" if has_illumina and has_ont else "ILLUMINA_ONLY" if has_illumina else "ONT_ONLY" if has_ont else "NEITHER"
        gtdb = stages.get("GTDBTK", {})
        checkm2 = stages.get("ASSEMBLY_QC", {})
        ont = stages.get("ONT_QC", {})
        quast = quast_by_sample.get(sample_id, {})
        retry = row["prior_status"] in {
            "ILLUMINA_FACILITY_RERUN_PENDING", "BLOCKED_ONT_QC", "BLOCKED_ASSEMBLY_QC", "PROCESS_FAILURE"
        }
        output = {
            "sample_id": sample_id,
            "patient_id": row["patient_id"],
            "timepoint": row["timepoint"],
            "material": row["material"],
            "agar_type": row["agar_type"],
            "availability": availability_label,
            "facility_qc_status": row["facility_qc_status"],
            "illumina_qc": stages.get("ILLUMINA_QC", stages.get("SWEEP_ILLUMINA_QC", {})).get("status", ""),
            "ont_qc": ont.get("status", ""),
            "assembly_qc": checkm2.get("status", ""),
            "completeness": checkm2.get("completeness", ""),
            "contamination": checkm2.get("contamination", ""),
            "assembly_contigs": quast.get("assembly_contigs", ""),
            "assembly_n50": quast.get("assembly_n50", ""),
            "assembly_total_length": quast.get("assembly_total_length", ""),
            "ont_bases": ont.get("number_of_bases", ""),
            "ont_n50": ont.get("n50", ""),
            "ont_mean_q": ont.get("mean_qual", ""),
            "species": gtdb.get("species", ""),
            "sequence_type": sequence_types.get(sample_id, ""),
            "status": status,
            "retried_previous_failure": str(retry).lower(),
            "review_suggestion": suggestion(status, row, stages),
            "reason": reason,
            "run_id": args.run_id,
        }
        output_rows.append(output)
        write_tsv(sample_dir / "sample_status.tsv", [output], list(output))
        collect_versions(sample_dir)
        if status in {"COMPLETED", "NON_ECOLI_COMPLETE"}:
            try:
                create_final_record(row, sample_dir, status, args.run_id)
            except (OSError, ValueError) as error:
                output["status"] = "PROCESS_FAILURE"
                output["reason"] = f"final completion validation failed: {error}"
                finalize_errors.append(f"{sample_id}: {error}")
                write_tsv(sample_dir / "sample_status.tsv", [output], list(output))

    fields = list(output_rows[0])
    status_path = args.output_dir / "sample_status.tsv"
    write_tsv(status_path, output_rows, fields)
    cohort_versions = args.output_dir / "software_versions.tsv"
    collect_cohort_versions(args.results_root, cohort_versions)
    provenance = write_provenance(
        args.artifact_dir,
        args.output_dir,
        args.results_root,
        args.configuration_base64,
        args.run_id,
    )
    counts = Counter(row["status"] for row in output_rows)
    retried_count = sum(row["retried_previous_failure"] == "true" for row in output_rows)
    summary_path = args.output_dir / "run_summary.txt"
    with summary_path.open("w") as handle:
        handle.write(f"AMREvoFlow run: {args.run_id}\n")
        handle.write(f"Samples: {len(output_rows)}\n")
        for status in TERMINAL_ORDER:
            value = retried_count if status == "RETRIED_PREVIOUS_FAILURE" else counts.get(status, 0)
            handle.write(f"{status}\t{value}\n")
        for status in sorted(set(counts) - set(TERMINAL_ORDER)):
            handle.write(f"{status}\t{counts[status]}\n")
        if finalize_errors:
            handle.write("WARNING\tSome final markers could not be validated\n")
            for error in finalize_errors:
                handle.write(f"FINALIZE_ERROR\t{error}\n")
        handle.write(f"STATUS_TABLE\t{status_path.resolve()}\n")
        handle.write(f"PLOTS\t{(args.output_dir / 'plots').resolve()}\n")
        handle.write(f"DASHBOARD\t{(args.output_dir / 'plots/pipeline_status.png').resolve()}\n")

    write_tsv(args.snapshot_dir / "sample_status.tsv", output_rows, fields)
    shutil.copy2(summary_path, args.snapshot_dir / "run_summary.txt")
    shutil.copy2(cohort_versions, args.snapshot_dir / "software_versions.tsv")
    snapshot_manifests = args.snapshot_dir / "manifests"
    snapshot_manifests.mkdir(parents=True, exist_ok=True)
    for resolved in args.resolved:
        shutil.copy2(resolved, snapshot_manifests / resolved.name)
    snapshot_provenance = args.snapshot_dir / "provenance"
    snapshot_provenance.mkdir(parents=True, exist_ok=True)
    for source in provenance.iterdir():
        if source.is_file():
            shutil.copy2(source, snapshot_provenance / source.name)
    if not args.no_plots:
        render_plots(output_rows, args.output_dir / "plots")
        snapshot_plots = args.snapshot_dir / "plots"
        if snapshot_plots.exists():
            shutil.rmtree(snapshot_plots)
        shutil.copytree(args.output_dir / "plots", snapshot_plots)

    print(summary_path.read_text(), end="")
    # Per-sample finalization problems are represented as PROCESS_FAILURE rows;
    # preserving the canonical report is more useful than discarding it.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
