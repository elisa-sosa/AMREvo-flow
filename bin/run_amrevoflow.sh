#!/usr/bin/env bash
# Run AMREvoFlow, pinned GEMSWEEP, then the unified sweep finalizer.
set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
RUN_ID=$(date +%Y%m%d_%H%M%S)
EXECUTE=false
ISOLATE_SAMPLES=""
ISOLATE_UNITS=""
SWEEP_SAMPLES=""
SWEEP_UNITS=""
ALLOWED_ROOTS=""
INDEX_PREFIX=""
REF_GROUPS=""
CARRIAGE_REFERENCES=""
RAW_ROOT=""
COMPLETED_LIMS=""
PENDING_V3=""

usage() {
  cat <<'EOF'
Usage: bin/run_amrevoflow.sh [options]
  --raw-root PATH
      Root containing illumina_reads/, nanopore_fastq/, hybrid_assemblies/
      and sweeps_fastq/. Generates all four ignored manifests automatically.
  [--completed-lims PATH]
  [--pending-v3 PATH]
  --isolate-samples PATH
  --isolate-read-units PATH
  --sweep-samples PATH
  --sweep-read-units PATH
  --allowed-roots JSON_ARRAY
  --index-prefix PATH
  --ref-groups PATH
  --carriage-references PATH
  [--run-id ID]
  [--execute]

With --raw-root, manifest paths and established AMREvo reference paths are
derived automatically. Without --execute the launcher generates the ignored
manifests and prints every Nextflow command. It never calls sbatch directly.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-root) RAW_ROOT=$2; shift 2 ;;
    --completed-lims) COMPLETED_LIMS=$2; shift 2 ;;
    --pending-v3) PENDING_V3=$2; shift 2 ;;
    --isolate-samples) ISOLATE_SAMPLES=$2; shift 2 ;;
    --isolate-read-units) ISOLATE_UNITS=$2; shift 2 ;;
    --sweep-samples) SWEEP_SAMPLES=$2; shift 2 ;;
    --sweep-read-units) SWEEP_UNITS=$2; shift 2 ;;
    --allowed-roots) ALLOWED_ROOTS=$2; shift 2 ;;
    --index-prefix) INDEX_PREFIX=$2; shift 2 ;;
    --ref-groups) REF_GROUPS=$2; shift 2 ;;
    --carriage-references) CARRIAGE_REFERENCES=$2; shift 2 ;;
    --run-id) RUN_ID=$2; shift 2 ;;
    --execute) EXECUTE=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -n "$RAW_ROOT" ]]; then
  [[ -d "$RAW_ROOT" ]] || { echo "ERROR: Raw root is not a directory: $RAW_ROOT" >&2; exit 2; }
  RAW_ROOT=$(cd "$RAW_ROOT" && pwd)
  DATA_PROJECT=$(cd "$RAW_ROOT/../../.." && pwd)
  COMPLETED_LIMS=${COMPLETED_LIMS:-$DATA_PROJECT/06_Overview_Data_Freezer/03_Overall_Data_sheets_LIMS/02_results/04_AMREvo_sequencing.csv}
  PENDING_V3=${PENDING_V3:-$DATA_PROJECT/05_sequencing/03_scripts/06_sampling_design/04_pending_wgs_v3_1177_with_locations.csv}
  ISOLATE_SAMPLES=${ISOLATE_SAMPLES:-$PROJECT_DIR/inputs/isolates.samples.local.tsv}
  ISOLATE_UNITS=${ISOLATE_UNITS:-$PROJECT_DIR/inputs/isolates.read_units.local.tsv}
  SWEEP_SAMPLES=${SWEEP_SAMPLES:-$PROJECT_DIR/inputs/sweeps.samples.local.tsv}
  SWEEP_UNITS=${SWEEP_UNITS:-$PROJECT_DIR/inputs/sweeps.read_units.local.tsv}
  ALLOWED_ROOTS=${ALLOWED_ROOTS:-[\"$RAW_ROOT\"]}
  INDEX_PREFIX=${INDEX_PREFIX:-$DATA_PROJECT/05_sequencing/04_jobs/ecoli_reference_curation/06_themisto_index/poppunk_628_plus_14_preliminary/index}
  REF_GROUPS=${REF_GROUPS:-$DATA_PROJECT/05_sequencing/04_jobs/ecoli_reference_curation/05_final_database/poppunk_628_plus_14_preliminary/reference_groups.txt}
  CARRIAGE_REFERENCES=${CARRIAGE_REFERENCES:-$DATA_PROJECT/05_sequencing/05_results/ecoli_reference_curation/preliminary_poppunk_628_plus_14_manifest.tsv}

  for required_dir in "$RAW_ROOT/illumina_reads" "$RAW_ROOT/nanopore_fastq" "$RAW_ROOT/sweeps_fastq/all"; do
    [[ -d "$required_dir" ]] || { echo "ERROR: Expected input directory missing: $required_dir" >&2; exit 2; }
  done
  if [[ -d "$RAW_ROOT/hybrid_assemblies" ]]; then
    echo "INFO: Existing hybrid_assemblies/ is legacy input and will not be imported; hybrids are rebuilt from reads." >&2
  fi

  python3 "$PROJECT_DIR/bin/build_isolate_manifests.py" \
    --completed-lims "$COMPLETED_LIMS" \
    --pending-v3 "$PENDING_V3" \
    --illumina-root "$RAW_ROOT/illumina_reads" \
    --ont-root "$RAW_ROOT/nanopore_fastq" \
    --samples-output "$ISOLATE_SAMPLES" \
    --read-units-output "$ISOLATE_UNITS"

  python3 "$PROJECT_DIR/bin/build_sweep_manifests.py" \
    --all-root "$RAW_ROOT/sweeps_fastq/all" \
    --samples-output "$SWEEP_SAMPLES" \
    --read-units-output "$SWEEP_UNITS"
fi

for value_name in ISOLATE_SAMPLES ISOLATE_UNITS SWEEP_SAMPLES SWEEP_UNITS ALLOWED_ROOTS INDEX_PREFIX REF_GROUPS CARRIAGE_REFERENCES; do
  [[ -n "${!value_name}" ]] || { echo "ERROR: Missing $value_name" >&2; exit 2; }
done

GEMSWEEP_DIR="$PROJECT_DIR/vendor/gemsweep/checkout"
RESULTS_ROOT="$PROJECT_DIR/results"
GEMSWEEP_MANIFEST="$RESULTS_ROOT/cohort/current/sweeps/gemsweep_manifest.csv"
RESOLVED_SWEEPS="$RESULTS_ROOT/cohort/current/provenance/validation/sweep/sweep.resolved_samples.tsv"
GEMSWEEP_OUTPUT="$RESULTS_ROOT/cohort/current/sweeps/gemsweep/$RUN_ID"

MAIN_CMD=(nextflow run "$PROJECT_DIR/main.nf" -profile slurm
  --isolate_samples "$ISOLATE_SAMPLES"
  --isolate_read_units "$ISOLATE_UNITS"
  --sweep_samples "$SWEEP_SAMPLES"
  --sweep_read_units "$SWEEP_UNITS"
  --allowed_roots "$ALLOWED_ROOTS"
  --gemsweep_index "$INDEX_PREFIX"
  --gemsweep_ref_groups "$REF_GROUPS"
  --carriage_references "$CARRIAGE_REFERENCES"
  --outdir "$RESULTS_ROOT"
  --run_id "$RUN_ID"
  -resume)

GEMSWEEP_CMD=(nextflow run "$GEMSWEEP_DIR/main.nf" -profile singularity
  -c "$PROJECT_DIR/conf/gemsweep.config"
  --manifest "$GEMSWEEP_MANIFEST"
  --ref_mode index
  --themisto_index "$INDEX_PREFIX"
  --ref_groups "$REF_GROUPS"
  --outdir "$GEMSWEEP_OUTPUT"
  -resume)

FINALIZE_CMD=(nextflow run "$PROJECT_DIR/gemsweep_finalize.nf" -profile slurm
  --manifest "$GEMSWEEP_MANIFEST"
  --resolved "$RESOLVED_SWEEPS"
  --gemsweep_output "$GEMSWEEP_OUTPUT"
  --index_prefix "$INDEX_PREFIX"
  --ref_groups "$REF_GROUPS"
  --outdir "$RESULTS_ROOT"
  --run_id "$RUN_ID"
  -resume)

printf 'Proposed AMREvoFlow command:'; printf ' %q' "${MAIN_CMD[@]}"; printf '\n'
printf 'Proposed GEMSWEEP command:'; printf ' %q' "${GEMSWEEP_CMD[@]}"; printf '\n'
printf 'Proposed finalizer command:'; printf ' %q' "${FINALIZE_CMD[@]}"; printf '\n'

if [[ "$EXECUTE" != true ]]; then
  echo "Dry run only. Re-run with --execute after reviewing all three commands."
  exit 0
fi

[[ -d "$GEMSWEEP_DIR/.git" ]] || {
  echo "ERROR: Run bin/bootstrap_gemsweep.sh first: $GEMSWEEP_DIR" >&2
  exit 2
}
[[ "$(git -C "$GEMSWEEP_DIR" rev-parse --short=7 HEAD)" == d06885f ]] || {
  echo "ERROR: GEMSWEEP checkout is not pinned to d06885f" >&2
  exit 2
}
[[ -e "$INDEX_PREFIX" || -n "$(find "$(dirname "$INDEX_PREFIX")" -maxdepth 1 -name "$(basename "$INDEX_PREFIX")*" -print -quit)" ]] || {
  echo "ERROR: Themisto index prefix not found: $INDEX_PREFIX" >&2
  exit 2
}
[[ -s "$REF_GROUPS" ]] || { echo "ERROR: Reference groups missing: $REF_GROUPS" >&2; exit 2; }
[[ -s "$CARRIAGE_REFERENCES" ]] || { echo "ERROR: Carriage reference manifest missing: $CARRIAGE_REFERENCES" >&2; exit 2; }

mkdir -p "$RESULTS_ROOT/cohort/current/sweeps"
"$PROJECT_DIR/bin/check_sweep_reference_coverage.py" \
  --carriage-references "$CARRIAGE_REFERENCES" \
  --ref-groups "$REF_GROUPS" \
  --output "$RESULTS_ROOT/cohort/current/sweeps/frozen_reference_coverage.tsv"

"${MAIN_CMD[@]}"
if [[ "$(wc -l < "$GEMSWEEP_MANIFEST" | tr -d ' ')" -le 1 ]]; then
  echo "No sweep samples are ready for GEMSWEEP; main workflow report is final."
  exit 0
fi
"${GEMSWEEP_CMD[@]}"
"${FINALIZE_CMD[@]}"
