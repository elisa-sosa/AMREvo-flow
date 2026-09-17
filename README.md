# AMREvoFlow

Reusable Nextflow workflow for AMREvo sequencing data. It will cover hybrid
isolate assembly first and sweep analysis later. It deliberately keeps code
separate from project data: raw reads, patient-linked sample sheets, results
and cluster-specific local settings must not be committed.

## Status

The initial scaffold validates an explicit sample manifest and its input files.
The first implementation stages are read QC, hybrid assembly, polishing,
assembly QC, MLST and metadata merging. Sweep-analysis modules will be added
after the corresponding AMREvo analysis is established.

## Inputs

The workflow accepts two arbitrary raw-data roots. Their internal run-folder
structure can vary freely. The sample sheet supplies paths relative to those
roots, so it explicitly selects the intended files when a sample has multiple
sequencing runs.

Required TSV columns:

```text
sample_id
illumina_r1_relpath
illumina_r2_relpath
ont_relpath
```

See `samples/samples.example.tsv`. Do not commit project sample sheets; use a
file ending in `.local.tsv` instead.

## First Run

```bash
nextflow run main.nf -profile slurm \\
  --samples /path/to/project-samples.local.tsv \\
  --illumina_root /path/to/illumina_reads \\
  --ont_root /path/to/nanopore_fastq \\
  --outdir /path/to/nextflow-results
```

Use `-resume` after an interrupted or extended run. Nextflow reuses completed
tasks whose inputs, command and software environment are unchanged.

## QC Plan

Each sample will produce an auditable QC row. The workflow will include:

- Illumina FastQC before and after trimming, plus Trimmomatic retention.
- ONT NanoPlot summaries.
- MultiQC reports for cohort-level read QC.
- QUAST or MetaQUAST assembly metrics.
- BUSCO and CheckM2 completeness/contamination checks.
- MLST scheme, ST and failure-reason audit fields.

QC should warn or fail explicitly; it must never silently select another run,
exclude a sample or overwrite an earlier result.

## Repository Rules

- Commit workflow code, generic configuration, documentation and synthetic examples.
- Do not commit raw reads, real sample sheets, results, logs or local cluster paths.
- Keep environment-specific overrides in ignored `conf/local.config`.
- Before public release, choose a licence, add a citation file and tag a version.
