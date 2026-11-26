/* Provisional sweep QC and GEMSWEEP hand-off processes. */

process SWEEP_ILLUMINA_QC {
    tag "${meta.sample_id}"
    label 'qc_light'
    container params.fastqc_image

    publishDir { "${params.outdir}/samples/sweeps/${meta.sample_id}/01_read_qc" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(r1), path(r2), path(trim_metrics)

    output:
    tuple val(meta), path("${meta.sample_id}.sweep_illumina_qc.tsv"), emit: status
    tuple val(meta), path("${meta.sample_id}.SWEEP_READY"), emit: ready, optional: true
    tuple val(meta), path("${meta.sample_id}.sweep_fastqc"), emit: reports
    path "${meta.sample_id}.sweep.software_versions.txt", emit: versions
    path "${meta.sample_id}.SWEEP_ILLUMINA_QC_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    facility_status="!{meta.facility_qc_status}"
    override="!{meta.illumina_qc_override}"
    mkdir -p "${sample_id}.sweep_fastqc"
    fastqc --threads "!{task.cpus}" --outdir "${sample_id}.sweep_fastqc" "!{r1}" "!{r2}"
    override_flag=""
    if [[ "${override}" == "true" ]]; then override_flag="--override"; fi
    assess_illumina_qc.sh \
      --sample-id "${sample_id}" --stage SWEEP_ILLUMINA_QC \
      --trim-metrics "!{trim_metrics}" --facility-status "${facility_status}" \
      ${override_flag} --output "${sample_id}.sweep_illumina_qc.tsv" \
      --ready-marker "${sample_id}.SWEEP_READY"
    {
      printf 'fastqc\t'
      fastqc --version 2>&1 | head -n 1
      printf 'classification\tPROVISIONAL_GEMSWEEP\n'
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.sweep.software_versions.txt"
    stage_status=$(tail -n 1 "${sample_id}.sweep_illumina_qc.tsv" | cut -f3)
    completion_state.sh create \
      --record "${sample_id}.SWEEP_ILLUMINA_QC_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class sweep --stage SWEEP_ILLUMINA_QC \
      --status "$stage_status" --input "!{r1}" --input "!{r2}" --input "!{trim_metrics}" \
      --output "${sample_id}.sweep_illumina_qc.tsv" --output "${sample_id}.sweep_fastqc" \
      --output "${sample_id}.sweep.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' \
      --run-id '!{meta.run_id}'
    '''
}


process GEMSWEEP_MANIFEST_ROW {
    tag "${meta.sample_id}"
    label 'validation'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path "${meta.sample_id}.gemsweep.csv", emit: row

    shell:
    '''
    set -euo pipefail
    r1_path="$(readlink -f "!{r1}")"
    r2_path="$(readlink -f "!{r2}")"
    printf '%s,%s,%s\n' '!{meta.sample_id}' "$r1_path" "$r2_path" \
      > "!{meta.sample_id}.gemsweep.csv"
    '''
}


process COLLECT_GEMSWEEP_MANIFEST {
    tag 'provisional-gemsweep-manifest'
    label 'validation'

    publishDir "${params.outdir}/cohort/current/sweeps", mode: 'copy', overwrite: true

    input:
    path rows

    output:
    path 'gemsweep_manifest.csv', emit: manifest
    path 'gemsweep_manifest.sha256', emit: checksum

    shell:
    '''
    set -euo pipefail
    printf 'ID,R1,R2\n' > gemsweep_manifest.csv
    for row in $(printf '%s\n' !{rows} | tr ' ' '\n' | LC_ALL=C sort); do
      cat "$row" >> gemsweep_manifest.csv
    done
    sha256sum gemsweep_manifest.csv > gemsweep_manifest.sha256
    '''
}
