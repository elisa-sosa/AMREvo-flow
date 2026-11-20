# Mapping from exploratory scripts to AMREvoFlow

This document records which parts of `00_AMREvo_Data/05_sequencing` inform the
production workflow. The exploratory scripts remain unchanged and are not
executed by AMREvoFlow.

## Logic carried forward

| Exploratory source | Production responsibility | Important change |
|---|---|---|
| `03_scripts/01_qc_reads/02a_ES_nanoplot_qc.sh` | ONT NanoPlot QC | Explicit manifest inputs, machine-readable thresholds, no file-existence-only skip |
| `03_scripts/01_qc_reads/04_Illumina_FastQC.sh` | Illumina FastQC | Runs per approved read pair; facility acceptance remains authoritative |
| `03_scripts/01_qc_reads/05_Illumina_MultiQC.sh` | Cohort technical QC report | Built from current and reused sample QC outputs |
| `03_scripts/02_assembly/01_unicycler.sh` | Hybrid Unicycler assembly | No automatic newest-run selection; only explicitly approved read units are used |
| `03_scripts/02_assembly/03_polish_medaka_polypolish_pypolca.sh` | Medaka, BWA-MEM2/Polypolish and PyPOLCA | Split into resumable Nextflow processes with checksummed state |
| `03_scripts/02_assembly/05_checkM2.sh` | Assembly completeness and contamination | PASS/WARN/FAIL policy is centralized and override-aware |
| `03_scripts/02_assembly/06_gtdbtk_species_qc.sh` | GTDB-Tk classification | Exact `Escherichia coli` species label gates species-specific work |
| `03_scripts/03_typing_ST/02_ST_run_mlst_on_hybrids.sh` | Achtman MLST | Runs only after the GTDB-Tk gate using `ecoli_achtman_4` |
| `03_scripts/04_resistance_virulence/01_amrfinderplus_polished_hybrids.sh` | AMRFinderPlus | Uses the established `Escherichia` organism configuration |
| `03_scripts/04_resistance_virulence/02_ktypr_polished_hybrids.sh` | kTYPr | Pinned to v1.0.1 and 30 kb flanking mode |
| `03_scripts/05_phylo_tree/*` | Snippy/Gubbins cohort and carriage-event analyses | Cohorts are rebuilt from all current eligible samples; small/no-variant groups receive explicit statuses |
| `03_scripts/07_sweeps/gemsweep/*` and `02_configs/gemsweep/*` | Provisional sweep analysis | GEMSWEEP is pinned and launched separately with the frozen index; its upstream provenance is retained |
| `03_scripts/08_results_summary/01_build_sequencing_allpatients.R` and `03_EDAandPlots/03_Sequencing_EDA_Plots.qmd` | Status summaries and PNG figures | Production scripts consume canonical workflow tables and do not publish a website |

## Deliberately not carried forward

- Provider short-read assemblies, annotations, typing, resistance and
  virulence outputs are not treated as production results. Only provider
  trimmed reads, their trimming record, and facility QC decision are imported.
- Scripts that search the filesystem and silently choose the newest run are
  replaced by reviewed read-unit manifests and a discovery-only helper.
- File-existence checks are not completion proof. Structured records validate
  input, parameter, software/database and output checksums.
- Manual Slurm array logic and embedded `sbatch`/`srun` orchestration are not
  copied. Nextflow's configured executor owns scheduling.
- Data-transfer, renaming, backup, cleanup, and one-off repair scripts remain
  outside the analytical workflow because AMREvoFlow must never mutate raw
  data.
- Culturomics, disc diffusion, sampling-design, MALDI and other work packages
  are intentionally deferred; they may later consume stable AMREvoFlow
  products without being coupled to sequencing execution.

## Provenance rule

The exploratory repository is a decision and validation source, not a runtime
dependency. Production code, schemas, patches and tests live in AMREvoFlow;
real reads and manifests are referenced read-only from ignored locations.
