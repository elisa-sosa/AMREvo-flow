# AMREvoFlow

AMREvoFlow is the production Nextflow workflow for the AMREvo isolate-WGS and
plate-sweep sequencing programme. It is deliberately separate from the
exploratory AMREvo Data repository: the latter remains a research notebook and
decision source, while this repository contains the reproducible analysis that
can ultimately accompany a publication.

The workflow is designed for staggered delivery. An isolate may have Illumina
reads, ONT reads, both, or neither. Available reads receive QC immediately;
hybrid assembly and downstream E. coli analysis start only when the required
inputs and QC gates are satisfied.

## What it does

```text
explicit manifests
  ├─ Illumina present ── facility record + FastQC/MultiQC ─┐
  ├─ ONT present ─────── NanoPlot QC ──────────────────────┤
  └─ both usable ─────── Unicycler → polishing → assembly QC
                                             ├─ GTDB-Tk species gate
                                             ├─ Achtman MLST
                                             ├─ AMRFinderPlus
                                             ├─ kTYPr
                                             └─ Snippy/Gubbins cohort analyses

known sweep pairs → input/QC audit → provisional frozen-index GEMSWEEP
```

Partial samples remain visible as `WAITING_FOR_ILLUMINA`, `WAITING_FOR_ONT`, or
`WAITING_FOR_READS`. A validated complete sample is skipped on future runs.
Changed inputs invalidate only the affected branch and its downstream work.

## Manifests

Real manifests are ignored and live under `inputs/`. Four files separate the
biological sample registry from individual sequencing deliveries:

```text
inputs/isolates.samples.local.tsv
inputs/isolates.read_units.local.tsv
inputs/sweeps.samples.local.tsv
inputs/sweeps.read_units.local.tsv
```

`sample_id` is the AMREvo `Unique_ID`. Sample rows may exist before reads
arrive. The AMREvo isolate registry is the deduplicated union of the completed
LIMS sequencing interface and the pending v3 WGS handover. Pending v3 isolates
expect both Illumina and ONT.

For repeated isolate Illumina deliveries, the project-approved rule is to use
the greatest parsed `esccol_YYMMDD` run containing the sample. Only its
provider-trimmed `0_trimming` pair can be approved; an incomplete latest run is
reported and never falls back silently to an older run. A unique flat ONT pass
FASTQ is approved directly. Generate these project manifests with:

```bash
python3 bin/build_isolate_manifests.py --help
```

For completed isolates, manifest generation additionally requires the selected
latest Illumina and unique ONT paths to match the corresponding explicit paths
in the LIMS sequencing interface. Any disagreement stops generation for review.

Sweeps use only the read-ID-deduplicated reconstruction under
`sweeps_fastq/all/`; historical `cat`, `cat-new`, delivery and test folders are
not production inputs. Convert the saved `ID,R1,R2` production manifest with:

```bash
python3 bin/build_sweep_manifests.py --help
```

See `samples/` and `schemas/` for synthetic examples and machine-readable
definitions.

Use the discovery helper to generate candidates for review:

```bash
python3 bin/discover_inputs.py --help
```

Validate without launching biological analysis:

```bash
python3 bin/validate_manifests.py --help
```

## Source layout

The workflow source is intentionally shallow:

```text
main.nf                     # top-level input and branch routing
modules/                    # one file per atomic process family
  validate.nf
  isolate.nf
  sweep.nf
  reporting.nf
  phylogeny.nf
subworkflows/               # one file per biological/cohort branch
  isolates.nf
  sweeps.nf
  reporting.nf
  phylogeny.nf
bin/                        # small validation, state and plotting helpers
conf/                       # cluster and GEMSWEEP executor settings
```

There is no `modules/local/reporting/...` hierarchy and no extra `workflows/`
layer. The distinction that remains is functional: a module is one runnable
tool/process, while a subworkflow joins those processes into an isolate,
sweep, reporting, or phylogeny route. Keeping those two levels makes routing
testable without burying the code in folders.

## Running

On the UZH Science Cluster, the intended checkout is:

```text
/shares/amr.imm.uzh/data/projects/AMREvo_Flow
```

The complete proposed command must be reviewed before submission. A
validation-only isolate run is:

```bash
nextflow run main.nf -profile slurm \
  --isolate_samples inputs/isolates.samples.local.tsv \
  --isolate_read_units inputs/isolates.read_units.local.tsv \
  --allowed_roots '["/shares/amr.imm.uzh/data/projects/NGS00401-Elisa-AMREvo/05_sequencing"]' \
  --enable_sweeps false \
  --validate_only true \
  --outdir results \
  --run_id REVIEWED_VALIDATION_ID
```

For routine cluster use, provide the common raw-read root. The launcher derives
the LIMS, v3, frozen-index and reference locations, regenerates all four ignored
manifests from currently available inputs, and prints the proposed commands:

```bash
bin/run_amrevoflow.sh \
  --raw-root /shares/amr.imm.uzh/data/projects/NGS00401-Elisa-AMREvo/05_sequencing/01_data/raw_reads
```

Illumina-only, ONT-only, both-platform and sweep-only availability are normal.
Existing exploratory `hybrid_assemblies/` files are not imported: samples with
both usable read types receive a new AMREvoFlow hybrid assembly.

The reviewed two-isolate plus one-sweep cluster pilot has a tracked controller
script. Its wrapper checks the pinned GEMSWEEP checkout, pilot IDs and container
images, and prints the exact submission by default:

```bash
./submit_pilot.sh
./submit_pilot.sh --run-id pilot_YYYYMMDD --execute
```

After pilot acceptance, the same dry-run/`--execute` pattern will be used by a
routine `submit.sh` that regenerates manifests from `--raw-root` on every run.

`bin/run_amrevoflow.sh` is the single user-facing launcher for the main
workflow followed by the separately pinned GEMSWEEP pipeline. It does not call
`sbatch`; both Nextflow invocations use their configured Slurm executors. With
no `--execute` flag it only prints the three complete commands for review:

```bash
bin/run_amrevoflow.sh \
  --isolate-samples inputs/isolates.samples.local.tsv \
  --isolate-read-units inputs/isolates.read_units.local.tsv \
  --sweep-samples inputs/sweeps.samples.local.tsv \
  --sweep-read-units inputs/sweeps.read_units.local.tsv \
  --allowed-roots '["/shares/amr.imm.uzh/data/projects/NGS00401-Elisa-AMREvo/05_sequencing"]' \
  --index-prefix /reviewed/path/to/frozen_642_reference_index \
  --ref-groups /reviewed/path/to/frozen_642_reference_groups.tsv \
  --carriage-references /reviewed/path/to/current_carriage_references.tsv
```

Only add `--execute` after reviewing those printed commands and obtaining the
required approval.

For the frozen-index coverage warning, the carriage-reference table should
provide `bin_id` (preferred), `sweep_group_id`, or `group_id` values matching
the one-column GEMSWEEP reference-group file. Reference/sample IDs are accepted
only when that same identifier is literally present in the group file.

## Completion and retry semantics

`-resume` alone is not considered proof of completion. Every stage records its
inputs, outputs, checksums, parameters, and software. A sample is
`SKIPPED_COMPLETE` only when `FINAL_COMPLETE.json` and all required results
still validate.

- Valid QC is reused while a partial sample waits for its other read type.
- Replacement reads rerun their branch and every dependent stage.
- QC-failed samples are retried on every new invocation, even with `-resume`.
- Technical failures retry once and remain recoverable.
- Cohort failure does not invalidate completed per-sample work.

## QC policy

Facility acceptance is authoritative for delivered Illumina data. Trimmomatic
survival bands and FastQC findings are visible warnings rather than hard gates.
ONT requires at least 5 Mb total bases, 1 kb N50, and mean Q-score 9. Hybrid
assemblies use CheckM2 PASS/WARN/FAIL bands documented in `ai-prompt.md`.

Only an exact GTDB-Tk species call of `Escherichia coli` enters Achtman MLST,
AMRFinderPlus, kTYPr, and E. coli phylogeny.

## Results

Current per-sample outputs use an IMMense-inspired numbered tree under
`results/samples/`. Cohort status, software versions, PNG plots, and compact
snapshots live under `results/cohort/`. Results are never committed.

The durable user entry points are:

```text
results/cohort/current/sample_status.tsv
results/cohort/current/run_summary.txt
results/cohort/current/plots/pipeline_status.png
```

The complete tree and the reason for its limited nesting are documented in
`docs/result-layout.md`.

AMREvoFlow does not publish a Quarto website. MultiQC HTML is retained only as
a technical QC artifact.

## Development and safety

- Never edit, rename, move, or delete raw sequencing data.
- Never run computational analysis on a login node.
- Never submit Slurm jobs without explicit approval of the full command.
- Never commit real sample sheets, patient-linked metadata, credentials, reads,
  results, work directories, or logs.
- Do not commit or push unless explicitly requested.

The first source implementation has passed local unit tests, Nextflow 25.10.9
linting, a no-reads validation-only execution, and a no-reads control-plane run
through status/snapshot generation. It has not yet been validated with the
configured tools, databases, frozen index, real reads, PNG runtime, or Slurm.
The full implementation contract is in `ai-prompt.md`; staged work and
acceptance criteria are in `docs/development-roadmap.md`.
