include { MULTIQC_COHORT; BUILD_REPORTS } from '../modules/reporting'


workflow REPORTING {
    take:
    resolved_files
    artifacts
    fastqc_reports
    results_root
    run_id
    configuration_b64

    main:
    MULTIQC_COHORT(fastqc_reports.collect())
    BUILD_REPORTS(resolved_files.collect(), artifacts.collect(), results_root, run_id, configuration_b64)

    emit:
    done = BUILD_REPORTS.out.done
    summary = BUILD_REPORTS.out.summary
    status = BUILD_REPORTS.out.status
    multiqc = MULTIQC_COHORT.out.report
}
