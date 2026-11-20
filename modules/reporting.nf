process MULTIQC_COHORT {
    tag 'cohort-read-qc'
    label 'reporting'
    container params.multiqc_image
    errorStrategy 'finish'

    publishDir "${params.outdir}/cohort/current/technical_qc", mode: 'copy', overwrite: true

    input:
    path fastqc_reports, stageAs: 'fastqc??/*'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data', emit: data

    shell:
    '''
    set -euo pipefail
    if ! find . -type f -name '*_fastqc.zip' -print -quit | grep -q .; then
      mkdir multiqc_data
      printf '<!doctype html><html><body><h1>AMREvoFlow MultiQC</h1><p>No FastQC reports were produced in this run.</p></body></html>\n' \
        > multiqc_report.html
      printf 'No FastQC reports were produced.\n' > multiqc_data/README.txt
      exit 0
    fi
    multiqc . --force --filename multiqc_report.html --outdir multiqc_output
    mv multiqc_output/multiqc_report.html .
    mv multiqc_output/multiqc_data .
    '''
}


process BUILD_REPORTS {
    tag "${run_id}"
    label 'reporting'
    container params.reporting_image
    errorStrategy 'finish'

    input:
    path resolved_files, stageAs: 'resolved??/*'
    path artifacts, stageAs: 'artifact??/*'
    val results_root
    val run_id
    val configuration_b64

    output:
    path 'report.done', emit: done
    path 'run_summary.txt', emit: summary
    path 'sample_status.tsv', emit: status

    shell:
    '''
    set -euo pipefail
    resolved_args=()
    for resolved in resolved??/*; do
      resolved_args+=(--resolved "$resolved")
    done
    plot_flag=""
    if [[ "!{params.render_plots}" != true ]]; then plot_flag="--no-plots"; fi
    build_reports.py \
      "${resolved_args[@]}" \
      --artifact-dir . \
      --results-root "!{results_root}" \
      --output-dir "!{results_root}/cohort/current" \
      --snapshot-dir "!{results_root}/cohort/snapshots/!{run_id}" \
      --run-id "!{run_id}" \
      --configuration-base64 '!{configuration_b64}' \
      ${plot_flag}
    cp "!{results_root}/cohort/current/run_summary.txt" run_summary.txt
    cp "!{results_root}/cohort/current/sample_status.tsv" sample_status.tsv
    touch report.done
    '''
}
