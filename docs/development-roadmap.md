# Development and validation checklist

The first implementation exists locally. Checkmarks below distinguish code
that is present and synthetically exercised from cluster validation, which is
deliberately still pending.

## 1. Input and completion contracts

- [x] Four-manifest interface for isolate/sweep samples and read units.
- [x] Complete expected isolate registry can include pending samples.
- [x] Explicit run decisions and ordered multi-unit merging.
- [x] Path, size, checksum, FASTQ structure, mate, duplicate path and duplicate
      identifier checks.
- [x] Discovery-only candidate generator.
- [x] Facility acceptance and reviewed override inputs.
- [x] Checksummed stage/final record helper.
- [x] Valid final skip, corrupt/missing output invalidation and changed-input
      archival.
- [x] Run-ID retry token for previous QC/process failures.
- [x] Synthetic tests for partial inputs, ordered units, duplicates, mates,
      skip, corruption and retry.
- [ ] Full-registry validation performance and memory benchmark on Slurm.

## 2. Isolate branch

- [x] Independent Illumina and ONT QC branches.
- [x] Facility acceptance with informational 95/90% survival bands.
- [x] ONT 5 Mb / 1 kb N50 / Q9 threshold gate.
- [x] Unicycler, Medaka, BWA-MEM2/Polypolish and PyPOLCA modules.
- [x] QUAST, CheckM2 and GTDB-Tk modules.
- [x] Exact GTDB *Escherichia coli* gate.
- [x] `ecoli_achtman_4` MLST, AMRFinderPlus `Escherichia` and kTYPr 30 kb mode.
- [x] Numbered per-sample publish locations and version artifacts.
- [ ] Validate each tool CLI and output parser against the configured cluster
      image/database pair.
- [ ] Validate kTYPr v1.0.1 image digest and known-sample result.

## 3. Status, reporting and phylogeny

- [x] Canonical cohort and per-sample status tables.
- [x] Required operational/QC/Figure 6/Figure 7 PNGs.
- [x] Current plus dated snapshot layout; no Quarto website publication.
- [x] Snapshot manifests, checksums, configuration, provenance, cohort software
      table and PNG retention.
- [x] Overall and patient-ST Snippy/Gubbins orchestration.
- [x] Deterministic reference selection.
- [x] `INSUFFICIENT_SAMPLES` and `NO_VARIABLE_SITES` handling.
- [x] Tree PNG renderer.
- [x] Nextflow 25.10.9 lint, local no-reads validation-only execution, and
      full no-reads control-plane/status/snapshot execution.
- [ ] End-to-end process/runtime validation with configured images on the
      cluster.
- [ ] Replace the still-supported but deprecated Nextflow `shell:` blocks with
      an equivalently safe interpolation pattern before a future Nextflow
      release removes them.
- [ ] Cohort comparison against the current reviewed exploratory trees/plots.

## 4. Provisional sweep branch

- [x] Separate sweep validation and QC.
- [x] Pinned upstream checkout bootstrap at `d06885f`.
- [x] Tracked Nextflow-25 compatibility patch and upstream license check.
- [x] Frozen-index GEMSWEEP hand-off and dry-run-first sequential launcher.
- [x] Provisional output registration and final completion records.
- [ ] Record and verify the reviewed 642-reference index/group-file checksums in
      the cluster-local configuration.
- [ ] Validate `P001-MCA-19` before any broader sweep execution.
- [x] Implement explicit frozen-index coverage warning table for current
      carriage-event references.
- [ ] Run that warning comparison with the reviewed cluster reference paths.

## 5. Rollout evidence

- [ ] Manifest-only cluster run.
- [ ] Known *E. coli* plus `P006-MCA-07-03` non-*E. coli* run.
- [ ] Unchanged second run emits `SKIPPED_COMPLETE`.
- [ ] Failed QC reruns despite `-resume`.
- [ ] New ONT delivery reuses valid Illumina QC.
- [ ] Replacement reads invalidate only the changed branch and downstream work.
- [ ] Historical approximately 550-isolate rerun.
- [ ] Full approximately 1,614-isolate rollout.

No Slurm submission, commit, push, raw-data copy or result deletion is part of
these checkmarks unless separately shown and explicitly approved.
