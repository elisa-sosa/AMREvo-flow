/*
 * AMREvoFlow isolate subworkflow.
 *
 * Inputs
 * ------
 * illumina_reads:
 *   tuple(meta, read_1, read_2, trim_metrics)
 * ont_reads:
 *   tuple(meta, ont_fastq)
 *
 * `trim_metrics` is always a staged path.  The top-level input builder uses
 * assets/no_trim_metrics.txt when no provider metric was delivered; an empty
 * sentinel is reported as NOT_AVAILABLE and does not block a facility-accepted
 * Illumina pair.
 *
 * The streams are intentionally independent through platform QC.  An isolate
 * with one modality therefore completes that QC branch and remains available
 * to the cohort status layer, while only samples with both READY markers enter
 * the hybrid chain.  Because `meta.retry_token` is a `val` input to every
 * analytical process, changing it for a previously failed sample invalidates
 * the corresponding Nextflow cache key without disturbing successful samples.
 */

include {
    ILLUMINA_QC;
    ONT_QC;
    UNICYCLER_HYBRID;
    MEDAKA_POLISH;
    BWA_MEM2_ALIGN;
    POLYPOLISH;
    PYPOLCA;
    QUAST_QC;
    CHECKM2_QC;
    GTDBTK_CLASSIFY;
    MLST_ACHTMAN;
    AMRFINDERPLUS;
    KTYPR
} from '../modules/isolate'


workflow ISOLATES {
    take:
    illumina_reads
    ont_reads

    main:
    /* Each platform QC branch is valid and publishable on its own. */
    ILLUMINA_QC(illumina_reads)
    ONT_QC(ont_reads)

    illumina_source = illumina_reads.map { meta, r1, r2, _trim_metrics ->
        tuple(meta.sample_id as String, meta, r1, r2)
    }
    illumina_ready = ILLUMINA_QC.out.ready.map { meta, marker ->
        tuple(meta.sample_id as String, marker)
    }
    usable_illumina = illumina_source
        .join(illumina_ready, failOnDuplicate: true, failOnMismatch: false)
        .map { sample_id, meta, r1, r2, _marker ->
            tuple(sample_id, meta, r1, r2)
        }

    ont_source = ont_reads.map { meta, ont_fastq ->
        tuple(meta.sample_id as String, meta, ont_fastq)
    }
    ont_ready = ONT_QC.out.ready.map { meta, marker ->
        tuple(meta.sample_id as String, marker)
    }
    usable_ont = ont_source
        .join(ont_ready, failOnDuplicate: true, failOnMismatch: false)
        .map { sample_id, meta, ont_fastq, _marker ->
            tuple(sample_id, meta, ont_fastq)
        }

    /* Inner join is the hybrid eligibility gate: both QC branches must pass. */
    hybrid_sources = usable_illumina
        .join(usable_ont, failOnDuplicate: true, failOnMismatch: false)
        .map { sample_id, illumina_meta, r1, r2, ont_meta, ont_fastq ->
            if (illumina_meta.sample_id != ont_meta.sample_id) {
                error "Internal isolate join mismatch for ${sample_id}"
            }
            tuple(sample_id, illumina_meta, r1, r2, ont_fastq)
        }

    UNICYCLER_HYBRID(
        hybrid_sources.map { _sample_id, meta, r1, r2, ont_fastq ->
            tuple(meta, r1, r2, ont_fastq)
        }
    )

    unicycler_assemblies = UNICYCLER_HYBRID.out.assembly.map { meta, assembly ->
        tuple(meta.sample_id as String, meta, assembly)
    }
    medaka_inputs = hybrid_sources
        .join(unicycler_assemblies, failOnDuplicate: true)
        .map { _sample_id, source_meta, _r1, _r2, ont_fastq, _result_meta, assembly ->
            tuple(source_meta, assembly, ont_fastq)
        }
    MEDAKA_POLISH(medaka_inputs)

    medaka_assemblies = MEDAKA_POLISH.out.assembly.map { meta, assembly ->
        tuple(meta.sample_id as String, meta, assembly)
    }
    bwa_inputs = hybrid_sources
        .join(medaka_assemblies, failOnDuplicate: true)
        .map { _sample_id, source_meta, r1, r2, _ont_fastq, _result_meta, assembly ->
            tuple(source_meta, assembly, r1, r2)
        }
    BWA_MEM2_ALIGN(bwa_inputs)
    POLYPOLISH(BWA_MEM2_ALIGN.out.alignments)

    polypolish_assemblies = POLYPOLISH.out.assembly.map { meta, assembly ->
        tuple(meta.sample_id as String, meta, assembly)
    }
    pypolca_inputs = hybrid_sources
        .join(polypolish_assemblies, failOnDuplicate: true)
        .map { _sample_id, source_meta, r1, r2, _ont_fastq, _result_meta, assembly ->
            tuple(source_meta, assembly, r1, r2)
        }
    PYPOLCA(pypolca_inputs)

    /* Assembly QC fans out, but GTDB-Tk waits for a PASS/WARN eligibility marker. */
    polished_assemblies = PYPOLCA.out.assembly
    QUAST_QC(polished_assemblies)
    CHECKM2_QC(polished_assemblies)

    polished_by_sample = polished_assemblies.map { meta, assembly ->
        tuple(meta.sample_id as String, meta, assembly)
    }
    assembly_ready = CHECKM2_QC.out.ready.map { meta, marker ->
        tuple(meta.sample_id as String, marker)
    }
    gtdbtk_inputs = polished_by_sample
        .join(assembly_ready, failOnDuplicate: true, failOnMismatch: false)
        .map { _sample_id, meta, assembly, _marker ->
            tuple(meta, assembly)
        }
    GTDBTK_CLASSIFY(gtdbtk_inputs)

    /* Exact GTDB species label is the sole E. coli-specific analysis gate. */
    ecoli_ready = GTDBTK_CLASSIFY.out.ecoli.map { meta, marker ->
        tuple(meta.sample_id as String, marker)
    }
    ecoli_assemblies = polished_by_sample
        .join(ecoli_ready, failOnDuplicate: true, failOnMismatch: false)
        .map { _sample_id, meta, assembly, _marker ->
            tuple(meta, assembly)
        }
    MLST_ACHTMAN(ecoli_assemblies)
    AMRFINDERPLUS(ecoli_assemblies)
    KTYPR(ecoli_assemblies)

    status_records = ILLUMINA_QC.out.status
        .mix(ONT_QC.out.status)
        .mix(UNICYCLER_HYBRID.out.status)
        .mix(MEDAKA_POLISH.out.status)
        .mix(BWA_MEM2_ALIGN.out.status)
        .mix(POLYPOLISH.out.status)
        .mix(PYPOLCA.out.status)
        .mix(QUAST_QC.out.status)
        .mix(CHECKM2_QC.out.status)
        .mix(GTDBTK_CLASSIFY.out.status)
        .mix(MLST_ACHTMAN.out.status)
        .mix(AMRFINDERPLUS.out.status)
        .mix(KTYPR.out.status)

    software_versions = ILLUMINA_QC.out.versions
        .mix(ONT_QC.out.versions)
        .mix(UNICYCLER_HYBRID.out.versions)
        .mix(MEDAKA_POLISH.out.versions)
        .mix(BWA_MEM2_ALIGN.out.versions)
        .mix(POLYPOLISH.out.versions)
        .mix(PYPOLCA.out.versions)
        .mix(QUAST_QC.out.versions)
        .mix(CHECKM2_QC.out.versions)
        .mix(GTDBTK_CLASSIFY.out.versions)
        .mix(MLST_ACHTMAN.out.versions)
        .mix(AMRFINDERPLUS.out.versions)
        .mix(KTYPR.out.versions)

    emit:
    status = status_records
    versions = software_versions
    illumina_qc = ILLUMINA_QC.out.status
    ont_qc = ONT_QC.out.status
    polished = polished_assemblies
    quast = QUAST_QC.out.report
    checkm2 = CHECKM2_QC.out.status
    gtdbtk = GTDBTK_CLASSIFY.out.status
    ecoli = ecoli_assemblies
    mlst = MLST_ACHTMAN.out.result
    amrfinderplus = AMRFINDERPLUS.out.result
    ktypr = KTYPR.out.result
    fastqc_reports = ILLUMINA_QC.out.reports
}
