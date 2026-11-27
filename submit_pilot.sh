#!/usr/bin/env bash
# Review and submit the tracked AMREvoFlow pilot controller job.
set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
EXECUTE=false
RUN_ID="pilot_$(date +%Y%m%d_%H%M%S)"

usage() {
  cat <<'EOF'
Usage: ./submit_pilot.sh [--run-id ID] [--execute]

Without --execute, validate prerequisites and print the exact sbatch command.
With --execute, submit the tracked submit/pilot.sbatch controller script.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id) RUN_ID=$2; shift 2 ;;
    --execute) EXECUTE=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

cd "$PROJECT_DIR"

for file in \
  inputs/pilot.isolates.samples.local.tsv \
  inputs/pilot.isolates.read_units.local.tsv \
  inputs/pilot.sweeps.samples.local.tsv \
  inputs/pilot.sweeps.read_units.local.tsv \
  containers/images/ktypr_1.0.1_b2e492530d53.sif \
  containers/images/reporting_multiqc-1.33_matplotlib-3.11.2_procps.sif; do
  [[ -s "$file" ]] || { echo "ERROR: Required pilot file missing or empty: $file" >&2; exit 2; }
done

[[ -d vendor/gemsweep/checkout/.git ]] || {
  echo "ERROR: GEMSWEEP checkout missing. Run bin/bootstrap_gemsweep.sh first." >&2
  exit 2
}
[[ "$(git -C vendor/gemsweep/checkout rev-parse --short=7 HEAD)" == d06885f ]] || {
  echo "ERROR: GEMSWEEP checkout is not pinned to d06885f" >&2
  exit 2
}

mapfile -t isolate_ids < <(tail -n +2 inputs/pilot.isolates.samples.local.tsv | cut -f1 | sort)
[[ "${isolate_ids[*]}" == "P001-COS-12-10 P006-MCA-07-03" ]] || {
  echo "ERROR: Pilot isolate manifest does not contain the reviewed two IDs" >&2
  exit 2
}
mapfile -t sweep_ids < <(tail -n +2 inputs/pilot.sweeps.samples.local.tsv | cut -f1 | sort)
[[ "${sweep_ids[*]}" == "P001-MCA-19" ]] || {
  echo "ERROR: Pilot sweep manifest does not contain the reviewed ID" >&2
  exit 2
}

mkdir -p logs
SBATCH_CMD=(sbatch
  --export="ALL,AMREVO_RUN_ID=$RUN_ID,AMREVOFLOW_ROOT=$PROJECT_DIR"
  --output="$PROJECT_DIR/logs/pilot-controller-%j.out"
  "$PROJECT_DIR/submit/pilot.sbatch")

printf 'Proposed pilot submission:'; printf ' %q' "${SBATCH_CMD[@]}"; printf '\n'
if [[ "$EXECUTE" != true ]]; then
  echo "Dry run only. Re-run with --execute after reviewing the command."
  exit 0
fi

"${SBATCH_CMD[@]}"
