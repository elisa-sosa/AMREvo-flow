nextflow.enable.dsl = 2

include { VALIDATE_MANIFEST } from './modules/local/validate_manifest'

params.help = false

def requireParameter(name, value) {
    if (!value) {
        error "Missing required parameter: --${name}"
    }
}

workflow {
    if (params.help) {
        log.info '''
Usage:
  nextflow run main.nf -profile slurm \\
    --samples samples.tsv \\
    --illumina_root /path/to/illumina_reads \\
    --ont_root /path/to/nanopore_fastq \\
    --outdir /path/to/results
'''.stripIndent()
        return
    }

    requireParameter('samples', params.samples)
    requireParameter('illumina_root', params.illumina_root)
    requireParameter('ont_root', params.ont_root)

    samples_file = file(params.samples, checkIfExists: true)
    validated = VALIDATE_MANIFEST(samples_file, params.illumina_root, params.ont_root)

    validated.view { "Validated manifest: ${it}" }

    /*
     * Planned modules:
     * input QC -> Unicycler -> polishing -> assembly QC -> MLST -> metadata merge.
     * Each module will consume one row from the validated manifest.
     */
}
