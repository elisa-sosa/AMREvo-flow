#!/usr/bin/env bash
#SBATCH --job-name=amrevo_pilot
#SBATCH --account=amr.imm.uzh
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=72:00:00

# Nextflow controller for the reviewed two-isolate plus one-sweep pilot.
set -euo pipefail

PROJECT_DIR=${AMREVOFLOW_ROOT:-/shares/amr.imm.uzh/data/projects/AMREvo_Flow}
RAW_ROOT=${AMREVO_RAW_ROOT:-/shares/amr.imm.uzh/data/projects/NGS00401-Elisa-AMREvo/05_sequencing/01_data/raw_reads}
DATA_PROJECT=$(cd "$RAW_ROOT/../../.." && pwd)
RUN_ID=${AMREVO_RUN_ID:?AMREVO_RUN_ID must be supplied by submit_pilot.sh}

module load miniforge3
source /home/elsosa/data/miniconda3/bin/activate \
  /home/elsosa/data/miniconda3/envs/nextflow_25_s3it
module load apptainer/1.5.0

export NXF_OPTS="-Xms128m -Xmx512m -XX:MaxMetaspaceSize=256m -XX:ActiveProcessorCount=2"

cd "$PROJECT_DIR"
mkdir -p logs results/logs

bin/run_amrevoflow.sh \
  --isolate-samples inputs/pilot.isolates.samples.local.tsv \
  --isolate-read-units inputs/pilot.isolates.read_units.local.tsv \
  --sweep-samples inputs/pilot.sweeps.samples.local.tsv \
  --sweep-read-units inputs/pilot.sweeps.read_units.local.tsv \
  --allowed-roots "[\"$RAW_ROOT\"]" \
  --index-prefix "$DATA_PROJECT/05_sequencing/04_jobs/ecoli_reference_curation/06_themisto_index/poppunk_628_plus_14_preliminary/index" \
  --ref-groups "$DATA_PROJECT/05_sequencing/04_jobs/ecoli_reference_curation/05_final_database/poppunk_628_plus_14_preliminary/reference_groups.txt" \
  --carriage-references "$DATA_PROJECT/05_sequencing/05_results/ecoli_reference_curation/preliminary_poppunk_628_plus_14_manifest.tsv" \
  --run-id "$RUN_ID" \
  --execute
