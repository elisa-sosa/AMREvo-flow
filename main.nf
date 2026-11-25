nextflow.enable.dsl = 2

include { VALIDATE_MANIFESTS; PREPARE_SAMPLE_STATE; MERGE_ILLUMINA_UNITS; MERGE_ONT_UNITS } from './modules/validate'
include { ISOLATES } from './subworkflows/isolates'
include { SWEEPS } from './subworkflows/sweeps'
include { REPORTING } from './subworkflows/reporting'
include { PHYLOGENY } from './subworkflows/phylogeny'


def requiredParameter(String name, value) {
    if (value == null || value.toString().trim() == '') {
        error "Missing required parameter: --${name}"
    }
}


def normaliseRoots(value) {
    if (value == null) {
        return []
    }
    if (value instanceof Collection) {
        return value.collect { it.toString() }.findAll { it }
    }
    def text = value.toString().trim()
    if (!text) {
        return []
    }
    if (text.startsWith('[')) {
        return new groovy.json.JsonSlurper().parseText(text).collect { it.toString() }
    }
    return text.split(',').collect { it.trim() }.findAll { it }
}


def sha256Object(value) {
    def json = groovy.json.JsonOutput.toJson(value)
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest(json.getBytes('UTF-8'))
    return digest.encodeHex().toString()
}


def metadataFor(row, String branch) {
    def branchFingerprint = branch == 'illumina' ? row.illumina_fingerprint : row.ont_fingerprint
    def retryingFailure = row.status == 'RETRIED_PREVIOUS_FAILURE'
    def executionToken = retryingFailure ? row.retry_token : branchFingerprint
    return [
        sample_id: row.sample_id,
        patient_id: row.patient_id,
        timepoint: row.timepoint,
        material: row.material,
        agar_type: row.agar_type,
        facility_qc_status: row.facility_qc_status ?: 'accepted',
        retry_token: executionToken,
        parameter_hash: row.parameter_hash,
        software_hash: row.software_hash,
        // Stable for unchanged inputs; a previous failure uses the new run ID.
        run_id: retryingFailure ? row.run_id : "input_${executionToken.take(16)}",
        illumina_qc_override: row.illumina_qc_override,
        ont_qc_override: row.ont_qc_override,
        assembly_qc_override: row.assembly_qc_override,
    ]
}


workflow {
    if (params.help) {
        log.info '''
AMREvoFlow

Required for each enabled sample class:
  --isolate_samples      inputs/isolates.samples.local.tsv
  --isolate_read_units   inputs/isolates.read_units.local.tsv
  --sweep_samples        inputs/sweeps.samples.local.tsv
  --sweep_read_units     inputs/sweeps.read_units.local.tsv
  --allowed_roots        JSON array or comma-separated read-only data roots

Example validation-only run:
  nextflow run main.nf -profile slurm \
    --isolate_samples inputs/isolates.samples.local.tsv \
    --isolate_read_units inputs/isolates.read_units.local.tsv \
    --allowed_roots '["/shares/amr.imm.uzh/data/projects/NGS00401-Elisa-AMREvo/05_sequencing"]' \
    --enable_sweeps false --validate_only true
'''.stripIndent()
        return
    }

    def roots = normaliseRoots(params.allowed_roots)
    if (!roots) {
        error 'At least one explicit --allowed_roots path is required'
    }
    def runId = params.run_id ?: new Date().format('yyyyMMdd_HHmmss', TimeZone.getTimeZone('Europe/Zurich'))
    def resultsRoot = file(params.outdir).toAbsolutePath().toString()
    def parameterConfig = [
        ont_min_bases: params.ont_min_bases,
        ont_min_n50: params.ont_min_n50,
        ont_min_mean_q: params.ont_min_mean_q,
        assembly_pass_completeness: params.assembly_pass_completeness,
        assembly_min_completeness: params.assembly_min_completeness,
        assembly_pass_contamination: params.assembly_pass_contamination,
        assembly_max_contamination: params.assembly_max_contamination,
        fastq_validation: params.fastq_validation,
        ktypr_mode: 'flanking_30000bp',
        mlst_scheme: 'ecoli_achtman_4',
        amrfinder_organism: 'Escherichia',
    ]
    def softwareConfig = [
        fastqc: params.fastqc_image,
        multiqc: params.multiqc_image,
        nanoplot: params.nanoplot_image,
        unicycler: params.unicycler_image,
        medaka: params.medaka_image,
        bwa_mem2: params.bwa_mem2_image,
        polypolish: params.polypolish_image,
        pypolca: params.pypolca_image,
        quast: params.quast_image,
        checkm2: params.checkm2_image,
        checkm2_db: params.checkm2_db,
        gtdbtk: params.gtdbtk_image,
        gtdbtk_db: params.gtdbtk_db,
        mlst: params.mlst_image,
        amrfinder: params.amrfinder_image,
        ktypr: params.ktypr_image,
        reporting: params.reporting_image,
        gemsweep_commit: params.gemsweep_commit,
        gemsweep_index: params.gemsweep_index,
        gemsweep_ref_groups: params.gemsweep_ref_groups,
        carriage_references: params.carriage_references,
    ]
    def parameterHash = sha256Object(parameterConfig)
    def softwareHash = sha256Object(softwareConfig)
    def configurationB64 = groovy.json.JsonOutput.toJson([
        parameters: parameterConfig,
        software_and_databases: softwareConfig,
    ]).getBytes('UTF-8').encodeBase64().toString()

    def validationEntries = []
    if (params.enable_isolates) {
        requiredParameter('isolate_samples', params.isolate_samples)
        requiredParameter('isolate_read_units', params.isolate_read_units)
        validationEntries << tuple(
            'isolate',
            file(params.isolate_samples, checkIfExists: true),
            file(params.isolate_read_units, checkIfExists: true),
        )
    }
    if (params.enable_sweeps) {
        requiredParameter('sweep_samples', params.sweep_samples)
        requiredParameter('sweep_read_units', params.sweep_read_units)
        requiredParameter('gemsweep_index', params.gemsweep_index)
        requiredParameter('gemsweep_ref_groups', params.gemsweep_ref_groups)
        validationEntries << tuple(
            'sweep',
            file(params.sweep_samples, checkIfExists: true),
            file(params.sweep_read_units, checkIfExists: true),
        )
    }
    if (!validationEntries) {
        error 'At least one of --enable_isolates or --enable_sweeps must be true'
    }

    validationInputs = channel.fromList(validationEntries)
    overrideFile = params.qc_overrides \
        ? file(params.qc_overrides, checkIfExists: true) \
        : file("${projectDir}/assets/no_overrides.tsv", checkIfExists: true)

    VALIDATE_MANIFESTS(
        validationInputs,
        channel.value(overrideFile),
        groovy.json.JsonOutput.toJson(roots),
        runId,
        parameterHash,
        softwareHash,
        resultsRoot,
        params.fastq_validation,
        params.gemsweep_index ?: '',
        params.gemsweep_ref_groups ?: '',
    )
    VALIDATE_MANIFESTS.out.summary.view { sampleClass, summary ->
        "Validated ${sampleClass} manifest:\n${summary.text}"
    }

    if (!params.validate_only) {
        /* Validation-only stops before any existing result can be archived. */
        PREPARE_SAMPLE_STATE(VALIDATE_MANIFESTS.out.resolved, resultsRoot, runId)
        preparedFiles = PREPARE_SAMPLE_STATE.out.resolved
        reportResolved = preparedFiles.map { sampleClass, resolved -> resolved }

        isolateRows = preparedFiles
            .filter { sampleClass, resolved -> sampleClass == 'isolate' }
            .map { sampleClass, resolved -> resolved }
            .splitCsv(header: true, sep: '\t', quote: '"')
        sweepRows = preparedFiles
            .filter { sampleClass, resolved -> sampleClass == 'sweep' }
            .map { sampleClass, resolved -> resolved }
            .splitCsv(header: true, sep: '\t', quote: '"')

        isolateIlluminaUnits = isolateRows
            .filter { row -> row.active == 'true' && row.skip_complete != 'true' && new groovy.json.JsonSlurper().parseText(row.illumina_r1_json) }
            .map { row ->
                def parser = new groovy.json.JsonSlurper()
                def r1 = parser.parseText(row.illumina_r1_json).collect { file(it, checkIfExists: true) }
                def r2 = parser.parseText(row.illumina_r2_json).collect { file(it, checkIfExists: true) }
                def trims = parser.parseText(row.trim_metrics_json).collect { value ->
                    value ? file(value, checkIfExists: true) : file("${projectDir}/assets/no_trim_metrics.txt", checkIfExists: true)
                }
                tuple(metadataFor(row, 'illumina'), r1, r2, trims)
            }
        sweepIlluminaUnits = sweepRows
            .filter { row -> row.active == 'true' && row.skip_complete != 'true' && new groovy.json.JsonSlurper().parseText(row.illumina_r1_json) }
            .map { row ->
                def parser = new groovy.json.JsonSlurper()
                def r1 = parser.parseText(row.illumina_r1_json).collect { file(it, checkIfExists: true) }
                def r2 = parser.parseText(row.illumina_r2_json).collect { file(it, checkIfExists: true) }
                def trims = parser.parseText(row.trim_metrics_json).collect { value ->
                    value ? file(value, checkIfExists: true) : file("${projectDir}/assets/no_trim_metrics.txt", checkIfExists: true)
                }
                tuple(metadataFor(row, 'illumina'), r1, r2, trims)
            }
        isolateOntUnits = isolateRows
            .filter { row -> row.active == 'true' && row.skip_complete != 'true' && new groovy.json.JsonSlurper().parseText(row.ont_fastq_json) }
            .map { row ->
                def paths = new groovy.json.JsonSlurper().parseText(row.ont_fastq_json).collect { file(it, checkIfExists: true) }
                tuple(metadataFor(row, 'ont'), paths)
            }

        MERGE_ILLUMINA_UNITS(isolateIlluminaUnits.mix(sweepIlluminaUnits))
        MERGE_ONT_UNITS(isolateOntUnits)

        isolateIllumina = MERGE_ILLUMINA_UNITS.out.reads
            .filter { meta, r1, r2, trim -> meta.material == 'isolate' }
        sweepIllumina = MERGE_ILLUMINA_UNITS.out.reads
            .filter { meta, r1, r2, trim -> meta.material == 'sweep' }
        isolateOnt = MERGE_ONT_UNITS.out.reads

        analysisArtifacts = channel.of(file("${projectDir}/assets/no_stage_artifacts.tsv", checkIfExists: true))
            .mix(VALIDATE_MANIFESTS.out.checksums.map { _sampleClass, path -> path })
            .mix(VALIDATE_MANIFESTS.out.summary.map { _sampleClass, path -> path })
            .mix(PREPARE_SAMPLE_STATE.out.archive_report.map { _sampleClass, path -> path })
        fastqcReports = channel.of(file("${projectDir}/assets/no_fastqc_reports.txt", checkIfExists: true))

        if (params.enable_isolates) {
            ISOLATES(isolateIllumina, isolateOnt)
            analysisArtifacts = analysisArtifacts
                .mix(ISOLATES.out.status.map { meta, path -> path })
                .mix(ISOLATES.out.versions)
                .mix(ISOLATES.out.polished.map { meta, path -> path })
                .mix(ISOLATES.out.quast.map { meta, path -> path })
                .mix(ISOLATES.out.checkm2.map { meta, path -> path })
                .mix(ISOLATES.out.gtdbtk.map { meta, path -> path })
                .mix(ISOLATES.out.mlst.map { meta, path -> path })
                .mix(ISOLATES.out.amrfinderplus.map { meta, path -> path })
                .mix(ISOLATES.out.ktypr.map { meta, path -> path })
            fastqcReports = fastqcReports.mix(ISOLATES.out.fastqc_reports.map { meta, path -> path })
        }
        if (params.enable_sweeps) {
            SWEEPS(sweepIllumina)
            analysisArtifacts = analysisArtifacts
                .mix(SWEEPS.out.status.map { meta, path -> path })
                .mix(SWEEPS.out.versions)
                .mix(SWEEPS.out.manifest)
            fastqcReports = fastqcReports.mix(SWEEPS.out.reports.map { meta, path -> path })
        }

        REPORTING(reportResolved, analysisArtifacts, fastqcReports, resultsRoot, runId, configurationB64)
        REPORTING.out.summary.view { summary -> "\n${summary.text}" }
        if (params.enable_isolates) {
            PHYLOGENY(REPORTING.out.status, resultsRoot, runId)
        }
    }
}
