/* Manifest validation and approved read-unit preparation. */

process VALIDATE_MANIFESTS {
    tag "${sample_class}"
    label 'validation'

    publishDir "${params.outdir}/cohort/current/provenance/validation/${sample_class}",
        mode: 'copy', overwrite: true

    input:
    tuple val(sample_class), path(samples), path(read_units)
    path qc_overrides
    val allowed_roots_json
    val run_id
    val parameter_hash
    val software_hash
    val results_root
    val fastq_validation
    val gemsweep_index
    val gemsweep_ref_groups

    output:
    tuple val(sample_class), path("${sample_class}.resolved_samples.tsv"), emit: resolved
    tuple val(sample_class), path("${sample_class}.input_checksums.tsv"), emit: checksums
    tuple val(sample_class), path("${sample_class}.validation_summary.json"), emit: summary

    shell:
    '''
    set -euo pipefail
    mkdir validation
    validator_args=(
      --sample-class "!{sample_class}"
      --samples "!{samples}"
      --read-units "!{read_units}"
      --qc-overrides "!{qc_overrides}"
      --allowed-roots-json '!{allowed_roots_json}'
      --results-root "!{results_root}"
      --run-id "!{run_id}"
      --parameter-hash "!{parameter_hash}"
      --software-hash "!{software_hash}"
      --fastq-validation "!{fastq_validation}"
      --output-dir validation
    )
    if [[ "!{sample_class}" == sweep ]]; then
      if [[ -n "!{gemsweep_index}" ]]; then
        validator_args+=(--dependency-prefix "!{gemsweep_index}")
      fi
      if [[ -n "!{gemsweep_ref_groups}" ]]; then
        validator_args+=(--dependency "!{gemsweep_ref_groups}")
      fi
    fi
    validate_manifests.py "${validator_args[@]}"
    mv validation/resolved_samples.tsv "!{sample_class}.resolved_samples.tsv"
    mv validation/input_checksums.tsv "!{sample_class}.input_checksums.tsv"
    mv validation/validation_summary.json "!{sample_class}.validation_summary.json"
    '''
}


process PREPARE_SAMPLE_STATE {
    tag "${sample_class}"
    label 'validation'

    publishDir "${params.outdir}/cohort/current/provenance/validation/${sample_class}",
        mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.archive_report.tsv') ? name : null }

    input:
    tuple val(sample_class), path(resolved)
    val results_root
    val run_id

    output:
    tuple val(sample_class), path("${sample_class}.prepared_samples.tsv"), emit: resolved
    tuple val(sample_class), path("${sample_class}.archive_report.tsv"), emit: archive_report

    shell:
    '''
    archive_invalidated.py \
      --sample-class "!{sample_class}" \
      --resolved "!{resolved}" \
      --results-root "!{results_root}" \
      --run-id "!{run_id}" \
      --output "!{sample_class}.prepared_samples.tsv" \
      --report "!{sample_class}.archive_report.tsv"
    '''
}


process MERGE_ILLUMINA_UNITS {
    tag "${meta.sample_id} (${r1_files.size()} unit${r1_files.size() == 1 ? '' : 's'})"
    label 'validation'

    input:
    tuple val(meta), path(r1_files, stageAs: 'r1??/*'), path(r2_files, stageAs: 'r2??/*'), path(trim_files, stageAs: 'tm??/*')

    output:
    tuple val(meta), path("${meta.sample_id}.R1.fastq.gz"), path("${meta.sample_id}.R2.fastq.gz"),
        path("${meta.sample_id}.quality_read_trimm_info"), emit: reads

    shell:
    '''
    set -euo pipefail
    shopt -s nullglob
    r1=(r1??/*)
    r2=(r2??/*)
    metrics=(tm??/*)
    test "${#r1[@]}" -eq "${#r2[@]}"
    test "${#r1[@]}" -gt 0

    combine_fastq() {
      local output="$1"
      shift
      local all_gzip=true
      local input
      for input in "$@"; do
        if [[ "$input" != *.gz ]]; then all_gzip=false; fi
      done
      if [[ "$all_gzip" == true ]]; then
        cat "$@" > "$output"
      else
        python3 - "$output" "$@" <<'PY'
    import gzip
    import shutil
    import sys
    from pathlib import Path
    output, *inputs = sys.argv[1:]
    with gzip.open(output, "wb") as destination:
        for raw in inputs:
            path = Path(raw)
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rb") as source:
                shutil.copyfileobj(source, destination, length=8 * 1024 * 1024)
    PY
      fi
    }

    combine_fastq "!{meta.sample_id}.R1.fastq.gz" "${r1[@]}"
    combine_fastq "!{meta.sample_id}.R2.fastq.gz" "${r2[@]}"
    : > "!{meta.sample_id}.quality_read_trimm_info"
    for metric in "${metrics[@]}"; do
      if [[ -s "$metric" ]]; then
        printf '### %s\n' "$metric" >> "!{meta.sample_id}.quality_read_trimm_info"
        cat "$metric" >> "!{meta.sample_id}.quality_read_trimm_info"
        printf '\n' >> "!{meta.sample_id}.quality_read_trimm_info"
      fi
    done
    '''
}


process MERGE_ONT_UNITS {
    tag "${meta.sample_id} (${ont_files.size()} unit${ont_files.size() == 1 ? '' : 's'})"
    label 'validation'

    input:
    tuple val(meta), path(ont_files, stageAs: 'ont??/*')

    output:
    tuple val(meta), path("${meta.sample_id}.ONT.fastq.gz"), emit: reads

    shell:
    '''
    set -euo pipefail
    shopt -s nullglob
    reads=(ont??/*)
    test "${#reads[@]}" -gt 0
    all_gzip=true
    for input in "${reads[@]}"; do
      if [[ "$input" != *.gz ]]; then all_gzip=false; fi
    done
    if [[ "$all_gzip" == true ]]; then
      cat "${reads[@]}" > "!{meta.sample_id}.ONT.fastq.gz"
    else
      python3 - "!{meta.sample_id}.ONT.fastq.gz" "${reads[@]}" <<'PY'
    import gzip
    import shutil
    import sys
    from pathlib import Path
    output, *inputs = sys.argv[1:]
    with gzip.open(output, "wb") as destination:
        for raw in inputs:
            path = Path(raw)
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rb") as source:
                shutil.copyfileobj(source, destination, length=8 * 1024 * 1024)
    PY
    fi
    '''
}
