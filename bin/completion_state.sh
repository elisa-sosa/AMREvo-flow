#!/usr/bin/env bash
# Create checksummed AMREvoFlow stage records without requiring Python.
set -euo pipefail

[[ ${1:-} == create ]] || {
  echo "ERROR: completion_state.sh currently supports only: create" >&2
  exit 2
}
shift

record="" sample_id="" sample_class="" stage="" status=""
parameter_hash="" software_hash="" run_id=""
inputs=() outputs=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --record) record=$2; shift 2 ;;
    --sample-id) sample_id=$2; shift 2 ;;
    --sample-class) sample_class=$2; shift 2 ;;
    --stage) stage=$2; shift 2 ;;
    --status) status=$2; shift 2 ;;
    --input) inputs+=("$2"); shift 2 ;;
    --output) outputs+=("$2"); shift 2 ;;
    --parameter-hash) parameter_hash=$2; shift 2 ;;
    --software-hash) software_hash=$2; shift 2 ;;
    --run-id) run_id=$2; shift 2 ;;
    *) echo "ERROR: Unknown completion-state option: $1" >&2; exit 2 ;;
  esac
done

for required in record sample_id sample_class stage status parameter_hash software_hash run_id; do
  [[ -n ${!required} ]] || { echo "ERROR: Missing --${required//_/-}" >&2; exit 2; }
done
[[ ${#outputs[@]} -gt 0 ]] || { echo "ERROR: At least one --output is required" >&2; exit 2; }

sha256_stream() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum | awk '{print $1}'
  else
    shasum -a 256 | awk '{print $1}'
  fi
}

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

canonical_path() {
  if command -v realpath >/dev/null 2>&1; then
    realpath "$1"
  else
    readlink -f "$1"
  fi
}

file_size() {
  wc -c < "$1" | awk '{print $1}'
}

json_escape() {
  local value=$1
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  value=${value//$'\t'/\\t}
  value=${value//$'\r'/\\r}
  value=${value//$'\n'/\\n}
  printf '%s' "$value"
}

snapshot() {
  local destination=$1
  shift
  local candidates=() original resolved file first=true
  for original in "$@"; do
    [[ -e $original ]] || { echo "ERROR: Required path is missing: $original" >&2; exit 2; }
    resolved=$(canonical_path "$original")
    if [[ -f $resolved ]]; then
      candidates+=("$resolved")
    elif [[ -d $resolved ]]; then
      while IFS= read -r file; do candidates+=("$file"); done < <(find "$resolved" -type f -print | LC_ALL=C sort)
    fi
  done
  printf '[' > "$destination"
  if [[ ${#candidates[@]} -gt 0 ]]; then
    while IFS= read -r file; do
      [[ -n $file ]] || continue
      if [[ $first == true ]]; then first=false; else printf ',' >> "$destination"; fi
      printf '{"path":"%s","sha256":"%s","size":%s}' \
        "$(json_escape "$file")" "$(sha256_file "$file")" "$(file_size "$file")" \
        >> "$destination"
    done < <(printf '%s\n' "${candidates[@]}" | LC_ALL=C sort -u)
  fi
  printf ']' >> "$destination"
}

record_dir=$(dirname "$record")
mkdir -p "$record_dir"
temporary_dir=$(mktemp -d "${TMPDIR:-/tmp}/amrevo-completion.XXXXXX")
trap 'rm -r "$temporary_dir"' EXIT
snapshot "$temporary_dir/inputs.json" "${inputs[@]}"
snapshot "$temporary_dir/outputs.json" "${outputs[@]}"
[[ $(cat "$temporary_dir/outputs.json") != '[]' ]] || {
  echo "ERROR: Completion record has no output files" >&2
  exit 2
}
input_fingerprint=$(sha256_stream < "$temporary_dir/inputs.json")
output_fingerprint=$(sha256_stream < "$temporary_dir/outputs.json")
completed_at=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
temporary_record="$record_dir/.${record##*/}.$$.tmp"
printf '{\n' > "$temporary_record"
printf '  "schema_version": 1,\n' >> "$temporary_record"
printf '  "sample_id": "%s",\n' "$(json_escape "$sample_id")" >> "$temporary_record"
printf '  "sample_class": "%s",\n' "$(json_escape "$sample_class")" >> "$temporary_record"
printf '  "stage": "%s",\n' "$(json_escape "$stage")" >> "$temporary_record"
printf '  "status": "%s",\n' "$(json_escape "$status")" >> "$temporary_record"
printf '  "run_id": "%s",\n' "$(json_escape "$run_id")" >> "$temporary_record"
printf '  "completed_at": "%s",\n' "$completed_at" >> "$temporary_record"
printf '  "parameter_hash": "%s",\n' "$(json_escape "$parameter_hash")" >> "$temporary_record"
printf '  "software_hash": "%s",\n' "$(json_escape "$software_hash")" >> "$temporary_record"
printf '  "input_fingerprint": "%s",\n' "$input_fingerprint" >> "$temporary_record"
printf '  "output_fingerprint": "%s",\n' "$output_fingerprint" >> "$temporary_record"
printf '  "inputs": %s,\n' "$(cat "$temporary_dir/inputs.json")" >> "$temporary_record"
printf '  "outputs": %s,\n' "$(cat "$temporary_dir/outputs.json")" >> "$temporary_record"
printf '  "metadata": {}\n}\n' >> "$temporary_record"
mv "$temporary_record" "$record"
