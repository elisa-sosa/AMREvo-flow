# AMREvoFlow Phase 2 roadmap

This roadmap describes the future sample-list-driven Nextflow workflow for AMREvo. It should be kept as the implementation guide while the current shell/R/Python workflow remains the validated source of truth.

## Current status

`AMREvoFlow` already contains a placeholder Nextflow scaffold. For now, do not expand the scaffold into full biological processing. The immediate requirement is documentation and a clear implementation path. The current production workflow still lives in the AMREvo sequencing repo under:

```text
/shares/amr.imm.uzh/data/projects/NGS00401-Elisa-AMREvo/05_sequencing
```

The current authoritative recovery endpoint is:

```text
05_results/phylogeny/sample_tables/26-09-22_phylogeny_metadata.tsv
```

That snapshot uses Achtman-only `st_analysis` values. Historical rMLST fallback groups are not valid drivers for the Phase 2 phylogeny or Nextflow design.

## Target user interface

The desired interface is:

```bash
nextflow run main.nf --samples selected_samples.tsv -resume
```

The user supplies a list of new, rebuilt, or reanalysis samples. The workflow processes only those samples while reusing valid existing outputs for unchanged samples.

Example:

```text
sample_id
P004-MCA--01-02
P004-MCA--01-06
P006-MCA-26-09
```

Negative timepoint IDs such as `P004-MCA--01-02` are valid and must be preserved exactly.

## Core design

Separate per-sample processing from cohort integration.

Per-sample processing should include:

- manifest validation and input discovery
- Illumina assembly where required
- hybrid Unicycler assembly
- Medaka
- Polypolish
- pypolca
- CheckM2
- Achtman MLST on Illumina assemblies
- Achtman MLST on hybrid/polished assemblies
- per-isolate QC and typing

Cohort integration should include:

- publishing an approved cohort metadata snapshot
- reference selection
- Snippy/core alignment updates
- Gubbins
- per-ST trees
- patient-ST carriage-event trees
- downstream Quarto reports and plots

Adding or replacing one isolate can change cohort membership, reference choice, alignments, Gubbins, and trees. These steps must therefore remain deliberate cohort-level actions, not automatic side effects of per-sample processing.

## Canonical manifest

Create one canonical manifest layer. A user-facing sample list can contain only `sample_id`, but a preflight step should resolve it to explicit inputs:

```text
sample_id
illumina_r1
illumina_r2
nanopore_fastq
illumina_assembly
hybrid_assembly
polished_assembly
```

Preflight must fail on:

- duplicate sample IDs
- ambiguous Illumina inputs
- ambiguous Nanopore inputs
- missing required inputs
- multiple candidate assemblies without an explicit precedence rule
- malformed sample IDs

Do not reimplement sample discovery independently in each process.

## State and provenance

The workflow must be state-aware, not merely file-existence-aware. A stage is valid only when expected outputs, sample identity, input provenance, completion markers, and validation checks agree.

The recent CheckM2 recovery showed why this matters: `quality_report.tsv` can exist even when a wrapper exits nonzero after CheckM2 removes its checksum manifest with `--force`.

Each selected run should produce an immutable run directory:

```text
runs/
  2026-09-22_p004_p006/
    samples.tsv
    resolved_input_manifest.tsv
    checksums.tsv
    qc/
    typing/
    commands/
    tool_versions.tsv
    provenance.json
    COMPLETE
```

Cohort publication should create separate immutable snapshots:

```text
cohort_snapshots/
  2026-09-22/
    metadata.tsv
    checkm2.tsv
    reference_manifest.tsv
    provenance.json
```

Never silently overwrite an earlier run or cohort snapshot.

## Implementation order

1. Finish and validate the current Phase 1 recovery endpoint.
2. Inventory existing assembly, QC, typing, and phylogeny scripts with exact inputs, outputs, resources, containers, and dependencies.
3. Create the canonical sample manifest and input-discovery/preflight layer.
4. Refactor tested scripts into atomic modules with explicit inputs.
5. Remove duplicated sample-discovery logic and special-purpose rerun scripts where possible.
6. Add immutable run directories, checksums, provenance, validation, and `COMPLETE` semantics.
7. Add a generic cohort publish/snapshot mechanism.
8. Port stable orchestration to Nextflow.
9. Add explicit cohort workflows for reference selection, Snippy/core alignment, Gubbins, per-ST trees, event trees, and reports.
10. Test with a small known sample set, then test adding one sample and confirm `-resume` does not recompute unaffected per-sample stages.

## Mapping from current scripts to future modules

Current mature scripts should be wrapped or refactored before being ported:

```text
03_scripts/03_typing_ST/assembly_mlst_compare.py
03_scripts/03_typing_ST/03_compare_assembly_mlst.sh
03_scripts/03_typing_ST/prepare_downstream_snapshot.py
03_scripts/04_phylo_tree/00_phylogeny_config_ES.sh
03_scripts/04_phylo_tree/01_prep_data_ES.sh
03_scripts/04_phylo_tree/02_snippy_on_reads.sh
03_scripts/04_phylo_tree/04_gubbins_per_ST.sh
03_scripts/04_phylo_tree/06_select_reference_per_event.R
03_scripts/04_phylo_tree/07_tree_per_carriage_event.sh
03_scripts/04_phylo_tree/08_plot_per_carriage_event.R
```

The first stable Nextflow implementation should not duplicate all biology. It should call explicit, tested one-sample or one-batch modules.

## Publish gate

Processing a selected sample does not automatically make it authoritative. Use explicit modes:

```bash
# Inspect only
nextflow run main.nf --samples new_samples.tsv --mode inspect

# Process selected samples
nextflow run main.nf --samples new_samples.tsv -resume

# Publish approved results into a cohort snapshot
nextflow run main.nf --samples new_samples.tsv --publish
```

Publication must include provenance and should make the selected cohort snapshot easy to cite from downstream QMD reports.

## Immediate next steps

- Keep this file as the Nextflow implementation roadmap.
- Keep `26-09-23_repo_open_tasks_for_ai.md` in the AMREvo data repo as the cross-AI handoff and open-task tracker.
- Do not implement more scaffold code until the current Phase 1 phylogeny outputs are confirmed.
- Once Phase 1 is confirmed, start with manifest validation and explicit path resolution, not with biological process rewrites.

## Current Phase 1 notes relevant to Phase 2

- `26-09-22_phylogeny_metadata.tsv` is the authoritative current metadata snapshot.
- The corrected P004 `TP-01` polished hybrid assemblies are Achtman `ST69`.
- `P006-MCA-26-09` is Achtman `ST10`.
- Current eligible per-ST trees should be driven by `st_analysis`.
- Current event trees should be driven by the `26-09-22` carriage-event reference audit.
- Historical rMLST fallback labels and old event trees such as obsolete `P004_ST2998` outputs should not influence the future workflow logic.
