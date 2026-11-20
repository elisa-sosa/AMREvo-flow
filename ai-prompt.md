# Implementation prompt for AMREvoFlow

## Objective

Implement AMREvoFlow as the clean, publishable sequencing workflow derived
from the exploratory `00_AMREvo_Data` repository.

Preserve the existing roadmap additions and executable mode changes. Inspect
`git status` before editing. Do not commit, push, copy raw data, or submit
Slurm jobs without explicit approval.

## Repository and cluster layout

Clone AMREvoFlow on the cluster into:

```text
/shares/amr.imm.uzh/data/projects/AMREvo_Flow
```

Use:

```text
AMREvo_Flow/
├── main.nf
├── conf/
├── modules/
├── subworkflows/
├── bin/
├── assets/
├── schemas/
├── samples/                 # synthetic examples only
├── inputs/                  # ignored real manifests/config
├── work/                    # ignored Nextflow cache
├── logs/                    # ignored
└── results/                 # ignored
    ├── samples/
    │   ├── isolates/<sample_id>/
    │   └── sweeps/<sample_id>/
    └── cohort/
        ├── current/
        └── snapshots/<run_id>/
```

Reference legacy raw data read-only using explicit paths. Future deliveries
may be placed in ignored `inputs/raw/` or referenced externally. Never place
raw reads, real manifests, patient-linked metadata, results, or logs in Git.

## Input interfaces

Use four project manifests:

```text
inputs/isolates.samples.local.tsv
inputs/isolates.read_units.local.tsv
inputs/sweeps.samples.local.tsv
inputs/sweeps.read_units.local.tsv
```

`sample_id` is exactly the existing AMREvo `Unique_ID`.

Sample fields:

```text
sample_id
patient_id
timepoint
material
agar_type
expect_illumina
expect_ont
active
```

Read-unit fields:

```text
sample_id
read_unit_id
platform
source_run
r1
r2
ont_fastq
trim_metrics
decision
facility_qc_status
merge_order
comment
```

Rules:

- The isolate table contains the complete expected cohort of approximately
  1,614 isolates, including pending samples.
- Sweep manifests contain only known samples; do not infer identities for the
  missing 63.
- `facility_qc_status` defaults to `accepted` for a complete delivered
  Illumina pair. Other values are `rerun_pending` and `rejected`.
- Never choose the newest sequencing run automatically.
- Combine multiple approved read units only in declared `merge_order`.
- Validate existence, size, checksums, paired-read integrity, duplicate paths,
  and duplicate read identifiers without modifying raw data.
- Missing expected inputs produce pending statuses. Missing declared inputs or
  ambiguous approved units are validation errors.
- Provide an input-discovery utility that produces candidates for review but
  never silently selects one.
- Provide an ignored `qc_overrides.local.tsv` containing sample, stage,
  reviewer, date, reason, and approval.

## Completion, skipping, and rerunning

Do not rely only on file existence or `-resume`.

Every successfully completed stage writes a structured completion record
containing sample/stage, input checksums, configuration and software/database
hashes, required output paths/sizes/checksums, timestamp, and status. Every
sample receives `FINAL_COMPLETE.json` only when all results required for that
sample class are present and validated.

Before scheduling a sample:

1. Validate `FINAL_COMPLETE.json`.
2. Confirm current approved inputs, parameters, software/databases, and final
   outputs still match.
3. If they match, emit `SKIPPED_COMPLETE` and reuse stored outputs in cohort
   reporting.
4. Otherwise invalidate the affected state and schedule required work.

Additional behavior:

- Retain valid Illumina-only or ONT-only QC while the other modality is pending.
- When the missing modality arrives, reuse valid QC and run the new branch plus
  downstream hybrid work.
- Changed replacement reads rerun the affected branch and downstream stages;
  valid QC from the unchanged modality is reused.
- A QC `FAIL` is never cache-valid and is retried on every invocation.
- Include run ID in the retry token for previous failures so `-resume` cannot
  silently reuse a failed analytical result.
- Technical process failures retry once, then remain recoverable with
  `-resume`.
- Cohort failure does not invalidate completed per-sample results.
- Non-*E. coli* isolates may finish after assembly QC and GTDB-Tk, with E.
  coli-specific outputs marked `NOT_APPLICABLE`.
- Sweeps become final-complete only after provisional GEMSWEEP outputs and
  provenance validate.

## Isolate workflow

Import provider-trimmed reads from `0_trimming/`, parse
`.quality_read_trimm_info`, and rerun FastQC/MultiQC. Provider short-read
assemblies, typing, annotation, and resistance outputs are not production
inputs.

```text
Illumina only -> Illumina QC -> WAITING_FOR_ONT
ONT only      -> NanoPlot QC -> WAITING_FOR_ILLUMINA
Both usable   -> hybrid assembly -> complete isolate workflow
Neither       -> WAITING_FOR_READS
```

```text
Unicycler
-> Medaka
-> BWA-MEM2/Polypolish
-> PyPOLCA
-> QUAST
-> CheckM2
-> GTDB-Tk
-> E. coli-only MLST, AMRFinderPlus, kTYPr and phylogeny
```

Use checksummed Apptainer images and record versions, beginning with the
cluster-validated exploratory versions. Pin kTYPr to v1.0.1, validate known
samples, and use documented 30 kb flanking mode after the E. coli gate. See
the [official kTYPr repository](https://github.com/SushiLab/kTYPr).

GTDB-Tk eligibility requires the exact species label `Escherichia coli`.
Other calls do not enter Achtman MLST, kTYPr, E. coli phylogeny, or the AMREvo
reference library. Use `ecoli_achtman_4` MLST and AMRFinderPlus with the
established `Escherichia` configuration.

## QC policy

### Illumina

Facility acceptance is authoritative. A complete delivered trimmed pair is
usable unless marked `rerun_pending` or `rejected`. Paired-read survival >=95%
is `GOOD`, 90-<95% is `REVIEW`, and <90% is `STRONG_REVIEW`; these and FastQC
findings are warnings and do not block a facility-accepted pair.
`rerun_pending` and `rejected` block hybrid work until replacement data or a
reviewed override appears.

### ONT

PASS requires total bases >=5 Mb, N50 >=1 kb, and mean Q-score >=9. Any failure
blocks hybrid work unless reviewed and overridden.

### Hybrid assembly

- PASS: CheckM2 completeness >=99% and contamination <=1%.
- WARN: completeness 95-<99% or contamination >1-<5%.
- FAIL: completeness <95% or contamination >=5%.
- WARN continues; FAIL blocks without an override.

Review suggestions remain non-definitive:

```text
ILLUMINA_FACILITY_RERUN_PENDING
RESEQUENCE_ONT_REVIEW
DNA_INTEGRITY_OR_LIBRARY_PREP_REVIEW
REEXTRACTION_OR_BOTH_LIBRARIES_REVIEW
ASSEMBLY_OR_SOURCE_READS_REVIEW
```

## Cohort phylogeny and plots

Rebuild cohort analyses from all current eligible E. coli hybrids, including
reused `SKIPPED_COMPLETE` samples. Implement overall and patient-Achtman-ST
Snippy/Gubbins phylogenies; choose references by earliest timepoint, highest
completeness, lowest contamination, then sample ID; report
`INSUFFICIENT_SAMPLES` and `NO_VARIABLE_SITES`; and save PNG tree plots.

Save plots under `results/cohort/current/plots/` and do not publish a Quarto
website. MultiQC HTML may remain as a technical artifact. Required PNGs cover
pipeline stage, availability, QC distributions/failures, assembly quality and
contiguity, species/ST, pending/blocked/failed/retried/skipped states, Figure 6
hybrid ST frequencies, and Figure 7 longitudinal ST composition generalized
to all patients with up to ten isolates per timepoint.

## Sweep workflow

Keep sweeps separate. Version 1 validates pairs, provenance, duplicate paths
and identifiers; runs FastQC and parses supplied trimming metrics; uses
facility acceptance; has no species gate; uses the frozen checksummed
preliminary 642-reference PopPUNK/Themisto index; never rebuilds it; warns when
new carriage-event references are absent; labels outputs
`PROVISIONAL_GEMSWEEP`; and tracks only known identities.

Pin separate upstream GEMSWEEP commit `d06885f`, retain the Nextflow-25 patch,
and provide one launcher that runs AMREvoFlow, GEMSWEEP, and the finalizer
sequentially. Never launch Nextflow inside a Nextflow process. Retain upstream
licensing/provenance from the
[GEMSWEEP repository](https://github.com/sanger-pathogens/gemsweep).

## Results and status reporting

Use an IMMense-inspired numbered organization based on the
[official IMMense pipeline](https://gitlab.uzh.ch/appliedmicrobiologyresearch/immense).

```text
isolates: 01_read_qc/ 02_hybrid_assembly/ 03_assembly_qc/ 04_typing/
          05_amr/ 06_capsule/ sample_status.tsv software_versions.txt
          FINAL_COMPLETE.json
sweeps:   01_read_qc/ 02_gemsweep/ sample_status.tsv software_versions.txt
          FINAL_COMPLETE.json
```

Maintain one canonical cohort `sample_status.tsv`, one `run_summary.txt`, one
pipeline-status PNG, per-sample completion/version files, one cohort
`software_versions.tsv`, and compact snapshots containing manifests,
checksums, configuration, status, provenance, and PNGs. Archive changed sample
results under a checksum-derived history location; never overwrite raw input.

Print counts for:

```text
COMPLETED
SKIPPED_COMPLETE
RETRIED_PREVIOUS_FAILURE
WAITING_FOR_ILLUMINA
WAITING_FOR_ONT
ILLUMINA_FACILITY_RERUN_PENDING
BLOCKED_ONT_QC
BLOCKED_ASSEMBLY_QC
NON_ECOLI_COMPLETE
PROCESS_FAILURE
PROVISIONAL_SWEEP_COMPLETE
```

Print exact status, summary, and dashboard paths. If cohort analysis fails,
preserve successful samples, mark the run partially complete, warn
prominently, and allow `-resume`. No module may call `sbatch`.

## Testing and acceptance

Test Illumina-only, ONT-only, both, and neither; validated skip; corrupt or
missing final output; reuse of partial QC; forced retry after QC failure;
branch-specific replacement invalidation; ordered units; duplicate IDs, paths,
reads, or mismatched mates; accepted <90% Illumina survival; facility
`rerun_pending`; ONT/assembly failure with and without overrides; non-E. coli;
and event groups with too few samples or no variants.

Cluster validation:

1. Validate manifests without computation.
2. Test a known E. coli and non-E. coli `P006-MCA-07-03`.
3. Confirm unchanged skip.
4. Confirm failed QC reruns with `-resume`.
5. Add ONT to an Illumina-only sample and reuse Illumina QC.
6. Validate kTYPr v1.0.1.
7. Validate GEMSWEEP against `P001-MCA-19`.
8. Reprocess the original approximately 550 isolates.
9. Scale to the full expected isolate manifest.

Show the complete proposed command and obtain approval before any Slurm
submission.

## Documentation and grouped commits

Update `README.md`, both roadmaps, schemas, result-tree documentation,
recovery instructions, completion semantics, and the distinction between
exploratory AMREvo Data and production AMREvoFlow.

Do not execute these commits without explicit approval:

```bash
git status --short
git add ai-prompt.md .gitignore README.md roadmap.md docs schemas samples conf nextflow.config
git commit -m "docs: define the AMREvoFlow production workflow"
```

```bash
git status --short
git add main.nf gemsweep_finalize.nf bin/*.py modules subworkflows assets containers
git commit -m "feat: add resumable isolate processing and cohort reports"
```

```bash
git status --short
git add bin/bootstrap_gemsweep.sh bin/run_amrevoflow.sh conf/gemsweep.config vendor/gemsweep patches
git commit -m "feat: integrate provisional PopPUNK GEMSWEEP analysis"
```

```bash
git status --short
git add tests
git commit -m "test: cover completion validation and partial inputs"
```

## Fixed assumptions

- The isolate target is approximately 1,614 and may grow slightly.
- `sample_id` equals `Unique_ID`.
- Incoming Illumina data are facility-processed IMMense-style outputs.
- Facility acceptance supersedes internal numeric thresholds.
- Partial read availability is normal.
- Complete validated samples are skipped; QC failures are retried every run.
- Replacement data rerun the changed branch and downstream dependencies.
- Existing approximately 550 isolates are reprocessed from reads.
- Only isolates with both usable modalities receive hybrid assembly.
- Sweeps remain provisional with the frozen current index.
- Real manifests, raw data, results, logs, and patient-linked metadata remain
  outside Git.
