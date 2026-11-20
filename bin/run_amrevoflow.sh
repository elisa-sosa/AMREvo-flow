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

usage() {
  cat <<'EOF'
Usage: bin/run_amrevoflow.sh [options]
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

Without --execute the launcher prints every Nextflow command and changes no
workflow state. It never calls sbatch directly.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
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
