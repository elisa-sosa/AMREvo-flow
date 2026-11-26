#!/usr/bin/env bash
# Apply the AMREvoFlow facility-led Illumina QC policy without Python.
set -euo pipefail

sample_id="" stage="" trim_metrics="" facility_status="accepted"
output="" ready_marker="" override=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) sample_id=$2; shift 2 ;;
    --stage) stage=$2; shift 2 ;;
    --trim-metrics) trim_metrics=$2; shift 2 ;;
    --facility-status) facility_status=$2; shift 2 ;;
    --override) override=true; shift ;;
    --output) output=$2; shift 2 ;;
    --ready-marker) ready_marker=$2; shift 2 ;;
    *) echo "ERROR: Unknown Illumina-QC option: $1" >&2; exit 2 ;;
  esac
done
for required in sample_id stage trim_metrics output ready_marker; do
  [[ -n ${!required} ]] || { echo "ERROR: Missing --${required//_/-}" >&2; exit 2; }
done

facility_status=$(printf '%s' "${facility_status:-accepted}" | tr '[:upper:]' '[:lower:]')
percentage=""
if [[ -s $trim_metrics ]]; then
  percentage=$(sed -nE \
    's/.*Both Surviving:[[:space:]]*[0-9,]+[[:space:]]*\(([0-9.]+)%\).*/\1/ip' \
    "$trim_metrics" | head -n 1)
fi
if [[ -z $percentage && -s $trim_metrics ]]; then
  percentage=$(sed -nE \
    's/.*both[_ ]surviving([^0-9]*)([0-9.]+).*/\2/ip' \
    "$trim_metrics" | head -n 1)
fi

metric_flag=NOT_AVAILABLE
if [[ -n $percentage ]]; then
  percentage=$(awk -v value="$percentage" 'BEGIN { printf "%.6g", value }')
  metric_flag=$(awk -v value="$percentage" 'BEGIN {
    if (value >= 95) print "GOOD";
    else if (value >= 90) print "REVIEW";
    else print "STRONG_REVIEW";
  }')
fi

accepted=false override_applied=false reason=""
if [[ $facility_status == accepted ]]; then
  accepted=true
elif [[ $override == true ]]; then
  accepted=true
  override_applied=true
  reason="override:facility_qc_status=$facility_status"
else
  reason="facility_qc_status=$facility_status"
fi
qc_status=FAIL
if [[ $accepted == true ]]; then qc_status=PASS; fi
printf 'sample_id\tstage\tstatus\tfacility_qc_status\tpaired_survival_pct\tmetric_flag\toverride_applied\treason\n' > "$output"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$sample_id" "$stage" "$qc_status" "$facility_status" "$percentage" \
  "$metric_flag" "$override_applied" "$reason" >> "$output"
if [[ $accepted == true ]]; then printf 'PASS\n' > "$ready_marker"; fi
