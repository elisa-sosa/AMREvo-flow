# AMREvoFlow production roadmap

AMREvoFlow is the clean, publication-oriented sequencing workflow derived from
the exploratory `00_AMREvo_Data` repository. The two repositories have
different jobs:

- AMREvo Data remains the exploratory record for sampling, culturomics, disc
  diffusion, method development and one-off bioinformatics investigations.
- AMREvoFlow owns stable sequencing inputs, reproducible production analysis,
  status reporting and publication-ready cohort outputs.

The intended cluster checkout is
`/shares/amr.imm.uzh/data/projects/AMREvo_Flow`. Raw data remain in their
existing read-only locations or ignored delivery locations; they are never
copied into Git.

## Implemented first version

### Inputs and state

- Separate sample and read-unit manifests for isolates and sweeps.
- `sample_id` is exactly the existing AMREvo `Unique_ID`.
- Expected samples may exist before either read type arrives.
- Explicit approved read units replace newest-file discovery.
- Multiple approved units require unique `merge_order` values.
- Root containment, non-empty files, compressed FASTQ parsing, paired IDs,
  duplicate paths, duplicate read IDs and checksums are validated.
- Discovery writes unapproved candidates only.
- Reviewed QC overrides use an ignored local TSV.
- Valid `FINAL_COMPLETE.json` records are reused only after current input,
  parameter, software and output validation.
- Changed final state is archived under a checksum-derived history path before
  replacement results are written.
- Waiting samples keep a stable retry token; previously failed samples receive
  the current run ID so `-resume` cannot reuse the failed analytical result.

### Isolates

The current DSL2 route is:

```text
facility-trimmed Illumina → FastQC + trimming summary ─┐
ONT → NanoPlot + threshold gate ──────────────────────┤
both usable → Unicycler → Medaka → BWA-MEM2/Polypolish → PyPOLCA
            → QUAST → CheckM2 → GTDB-Tk
            → exact E. coli gate → Achtman MLST + AMRFinderPlus + kTYPr
```

Illumina-only and ONT-only samples finish their available QC branch and wait.
Neither state is treated as an error. Non-*E. coli* isolates finish after
assembly QC and GTDB-Tk with species-specific work marked not applicable.

QC policy is encoded as agreed: facility acceptance controls Illumina
eligibility; ONT requires 5 Mb, 1 kb N50 and Q9; assembly PASS/WARN/FAIL uses
the CheckM2 99/95% completeness and 1/5% contamination boundaries. WARN
continues and FAIL requires a reviewed override.

### Cohort products

- One canonical `sample_status.tsv`, terminal `run_summary.txt`, technical
  MultiQC report, cohort software table and PNG dashboard.
- Compact snapshots retain normalized manifests, input checksums, run
  configuration, validation/archival provenance, software versions and PNGs.
- Figure 6 reports Achtman-ST frequencies among completed hybrids.
- Figure 7 shows up to ten isolates per patient/timepoint and generalizes to
  the full cohort.
- Current eligible and reused *E. coli* hybrids feed overall and
  patient–Achtman-ST Snippy/Gubbins groups.
- Reference selection is earliest timepoint, highest completeness, lowest
  contamination, then sample ID.
- Groups with fewer than three samples or no variable sites receive explicit
  statuses rather than opaque failures.
- Cohort failure does not invalidate successful per-sample results.

### Sweeps

- Sweep identities and Illumina read units are validated separately from
  isolates; unknown missing identities are never invented.
- Facility-accepted sweep pairs receive FastQC and a pinned GEMSWEEP manifest.
- One dry-run-first launcher executes AMREvoFlow, upstream GEMSWEEP and the
  unified finalizer sequentially. No Nextflow process starts another Nextflow
  execution and no module calls `sbatch`.
- GEMSWEEP is pinned to commit `d06885f`; the tracked Nextflow-25 compatibility
  patch is applied to an ignored upstream checkout with its original license.
- The existing PopPUNK/Themisto reference index is supplied read-only and is
  never rebuilt automatically.
- Finalized outputs are explicitly labelled `PROVISIONAL_GEMSWEEP`.

## Validation gates still required

Code presence is not cluster validation. Before production use:

1. Verify all configured Apptainer paths, digests, commands and database
   releases on the Science Cluster.
2. Build and checksum the kTYPr v1.0.1 image, then validate 30 kb flanking mode
   on known isolates.
3. Run validation-only on synthetic and reviewed local manifests.
4. Run one known *E. coli* and non-*E. coli* `P006-MCA-07-03`.
5. Demonstrate unchanged skip, corrupt-output invalidation, failed-QC retry and
   addition of ONT to a previously Illumina-only isolate.
6. Validate overall and event phylogenies, including too-small and no-variant
   groups.
7. Validate the pinned frozen-index sweep route on `P001-MCA-19` and compare it
   with the reviewed exploratory output.
8. Reprocess the original approximately 550 isolates.
9. Scale to the complete approximately 1,614-isolate registry, then ingest new
   deliveries only through reviewed read units.

Local implementation validation completed with Nextflow 25.10.9: all workflow
files lint without errors, unit tests pass, a no-reads validation-only run
stops before preparation, and a no-reads control-plane run produces
`WAITING_FOR_READS`, provenance, snapshots, technical-QC placeholders and
empty phylogeny manifests. No biological tool or Slurm job was invoked. The
remaining gates above require the Science Cluster.

Every proposed cluster command must be displayed in full and explicitly
approved before execution.

## Historical decisions retained

- Negative timepoints such as `P004-MCA--01-02` must be preserved literally.
- The corrected P004 negative-timepoint hybrids are Achtman ST69.
- `P006-MCA-26-09` is Achtman ST10.
- Historical rMLST fallback groups and obsolete event trees must not drive the
  production phylogeny.
- The provisional frozen sweep index remains a deliberately versioned interim
  method; future reference-library work will be a separate reviewed change.

See `docs/legacy-script-mapping.md` for the exploratory scripts carried forward
and `docs/development-roadmap.md` for test-level progress.
