include {
    SWEEP_ILLUMINA_QC;
    GEMSWEEP_MANIFEST_ROW;
    COLLECT_GEMSWEEP_MANIFEST
} from '../modules/sweep'


workflow SWEEPS {
    take:
    illumina_reads

    main:
    SWEEP_ILLUMINA_QC(illumina_reads)

    sources = illumina_reads.map { meta, r1, r2, trim_metrics ->
        tuple(meta.sample_id as String, meta, r1, r2)
    }
    ready = SWEEP_ILLUMINA_QC.out.ready.map { meta, marker ->
        tuple(meta.sample_id as String, marker)
    }
    usable = sources
        .join(ready, failOnDuplicate: true, failOnMismatch: false)
        .map { sample_id, meta, r1, r2, marker -> tuple(meta, r1, r2) }

    GEMSWEEP_MANIFEST_ROW(usable)
    manifest_rows = GEMSWEEP_MANIFEST_ROW.out.row
        .ifEmpty(file("${projectDir}/assets/no_gemsweep_rows.csv"))
    COLLECT_GEMSWEEP_MANIFEST(manifest_rows.collect())

    emit:
    status = SWEEP_ILLUMINA_QC.out.status
    versions = SWEEP_ILLUMINA_QC.out.versions
    manifest = COLLECT_GEMSWEEP_MANIFEST.out.manifest
    checksum = COLLECT_GEMSWEEP_MANIFEST.out.checksum
    reports = SWEEP_ILLUMINA_QC.out.reports
}
