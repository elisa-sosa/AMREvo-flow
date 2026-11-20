nextflow.enable.dsl = 2


def requiredParameter(String name, value) {
    if (value == null || value.toString().trim() == '') {
        error "Missing required parameter: --${name}"
    }
}

process FINALIZE_GEMSWEEP {
    tag params.run_id
    label 'reporting'
    container params.reporting_image

    input:
    path manifest
    path resolved
    val gemsweep_output
    val results_root
    val index_prefix
    path ref_groups
    val run_id

    output:
    path 'gemsweep_finalize.done'
    path 'run_summary.txt'

    shell:
    '''
    finalize_gemsweep.py \
      --manifest "!{manifest}" \
      --resolved "!{resolved}" \
      --gemsweep-output "!{gemsweep_output}" \
      --results-root "!{results_root}" \
      --index-prefix "!{index_prefix}" \
      --ref-groups "!{ref_groups}" \
      --commit d06885f \
      --run-id "!{run_id}"
    cp "!{results_root}/cohort/current/run_summary.txt" run_summary.txt
    touch gemsweep_finalize.done
    '''
}

workflow {
    requiredParameter('manifest', params.manifest)
    requiredParameter('resolved', params.resolved)
    requiredParameter('gemsweep_output', params.gemsweep_output)
    requiredParameter('index_prefix', params.index_prefix)
    requiredParameter('ref_groups', params.ref_groups)
    requiredParameter('run_id', params.run_id)
    def resultsRoot = file(params.outdir).toAbsolutePath().toString()
    FINALIZE_GEMSWEEP(
        file(params.manifest, checkIfExists: true),
        file(params.resolved, checkIfExists: true),
        file(params.gemsweep_output, checkIfExists: true).toString(),
        resultsRoot,
        file(params.index_prefix).toAbsolutePath().toString(),
        file(params.ref_groups, checkIfExists: true),
        params.run_id,
    )
    FINALIZE_GEMSWEEP.out.run_summary.view { summary -> summary.text }
}
