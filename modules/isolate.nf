/*
 * Atomic isolate processes for AMREvoFlow.
 *
 * Channel contract
 * ----------------
 * `meta` is a map containing at least `sample_id`, `patient_id`, `timepoint`,
 * `facility_qc_status`, `retry_token`, and `illumina_qc_override` /
 * `ont_qc_override` / `assembly_qc_override` booleans.  The retry token is
 * changed by the manifest/state layer only for a previously failed sample, so
 * a failed analytical result is not silently reused with `-resume`.
 *
 * These processes never submit Slurm jobs themselves.  Resource selection and
 * execution are entirely delegated to Nextflow profiles and process labels.
 */

process ILLUMINA_QC {
    tag "${meta.sample_id}"
    label 'qc_light'
    container params.fastqc_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/01_read_qc/illumina" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(r1), path(r2), path(trim_metrics)

    output:
    tuple val(meta), path("${meta.sample_id}.illumina_qc.tsv"), emit: status
    tuple val(meta), path("${meta.sample_id}.ILLUMINA_READY"), emit: ready, optional: true
    tuple val(meta), path("${meta.sample_id}.fastqc"), emit: reports
    path "${meta.sample_id}.illumina.software_versions.txt", emit: versions
    path "${meta.sample_id}.ILLUMINA_QC_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail

    sample_id="!{meta.sample_id}"
    facility_status="!{meta.facility_qc_status}"
    override="!{meta.illumina_qc_override}"

    mkdir -p "${sample_id}.fastqc"
    fastqc --threads "!{task.cpus}" --outdir "${sample_id}.fastqc" "!{r1}" "!{r2}"

    override_flag=""
    if [[ "${override}" == "true" ]]; then override_flag="--override"; fi
    assess_illumina_qc.py \
      --sample-id "${sample_id}" --stage ILLUMINA_QC \
      --trim-metrics "!{trim_metrics}" --facility-status "${facility_status}" \
      ${override_flag} --output "${sample_id}.illumina_qc.tsv" \
      --ready-marker "${sample_id}.ILLUMINA_READY"

    {
      printf 'fastqc\t'
      fastqc --version 2>&1 | head -n 1
      printf 'facility_qc_status\t%s\n' "${facility_status:-accepted}"
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.illumina.software_versions.txt"
    stage_status=$(tail -n 1 "${sample_id}.illumina_qc.tsv" | cut -f3)
    completion_state.py create \
      --record "${sample_id}.ILLUMINA_QC_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage ILLUMINA_QC \
      --status "$stage_status" --input "!{r1}" --input "!{r2}" --input "!{trim_metrics}" \
      --output "${sample_id}.illumina_qc.tsv" --output "${sample_id}.fastqc" \
      --output "${sample_id}.illumina.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' \
      --run-id '!{meta.run_id}'
    '''
}


process ONT_QC {
    tag "${meta.sample_id}"
    label 'qc_medium'
    container params.nanoplot_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/01_read_qc/ont" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(ont_reads)

    output:
    tuple val(meta), path("${meta.sample_id}.ont_qc.tsv"), emit: status
    tuple val(meta), path("${meta.sample_id}.ONT_READY"), emit: ready, optional: true
    tuple val(meta), path("${meta.sample_id}.nanoplot"), emit: reports
    path "${meta.sample_id}.nanoplot.software_versions.txt", emit: versions
    path "${meta.sample_id}.ONT_QC_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail

    sample_id="!{meta.sample_id}"
    override="!{meta.ont_qc_override}"
    mkdir -p "${sample_id}.nanoplot"

    NanoPlot --fastq "!{ont_reads}" --outdir "${sample_id}.nanoplot" \
      --threads "!{task.cpus}"

    parse_ont_qc.py \
      --stats "${sample_id}.nanoplot/NanoStats.txt" \
      --output "${sample_id}.ont_qc.tsv" \
      --sample-id "${sample_id}" \
      --override "${override}" \
      --min-bases "!{params.ont_min_bases}" \
      --min-n50 "!{params.ont_min_n50}" \
      --min-mean-qual "!{params.ont_min_mean_q}"

    {
      printf 'nanoplot\t'
      NanoPlot --version 2>&1 | head -n 1
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.nanoplot.software_versions.txt"
    stage_status=$(tail -n 1 "${sample_id}.ont_qc.tsv" | cut -f3)
    completion_state.py create \
      --record "${sample_id}.ONT_QC_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage ONT_QC \
      --status "$stage_status" --input "!{ont_reads}" \
      --output "${sample_id}.ont_qc.tsv" --output "${sample_id}.nanoplot" \
      --output "${sample_id}.nanoplot.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' \
      --run-id '!{meta.run_id}'
    '''
}


process UNICYCLER_HYBRID {
    tag "${meta.sample_id}"
    label 'assembly_heavy'
    container params.unicycler_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/02_hybrid_assembly/01_unicycler" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(r1), path(r2), path(ont_reads)

    output:
    tuple val(meta), path("${meta.sample_id}.unicycler.fasta"), emit: assembly
    tuple val(meta), path("${meta.sample_id}.unicycler.stage.tsv"), emit: status
    path "${meta.sample_id}.unicycler.software_versions.txt", emit: versions
    path "${meta.sample_id}.UNICYCLER_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"

    unicycler -1 "!{r1}" -2 "!{r2}" -l "!{ont_reads}" \
      -o unicycler --threads "!{task.cpus}" --mode normal --keep 1
    test -s unicycler/assembly.fasta
    cp unicycler/assembly.fasta "${sample_id}.unicycler.fasta"
    printf 'sample_id\tstage\tstatus\treason\n%s\tUNICYCLER\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.unicycler.stage.tsv"
    {
      printf 'unicycler\t'
      unicycler --version 2>&1 | head -n 1
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.unicycler.software_versions.txt"
    completion_state.py create --record "${sample_id}.UNICYCLER_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage UNICYCLER --status PASS \
      --input "!{r1}" --input "!{r2}" --input "!{ont_reads}" \
      --output "${sample_id}.unicycler.fasta" --output "${sample_id}.unicycler.stage.tsv" \
      --output "${sample_id}.unicycler.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process MEDAKA_POLISH {
    tag "${meta.sample_id}"
    label 'polish_heavy'
    container params.medaka_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/02_hybrid_assembly/02_medaka" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly), path(ont_reads)

    output:
    tuple val(meta), path("${meta.sample_id}.medaka.fasta"), emit: assembly
    tuple val(meta), path("${meta.sample_id}.medaka.stage.tsv"), emit: status
    path "${meta.sample_id}.medaka.software_versions.txt", emit: versions
    path "${meta.sample_id}.MEDAKA_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    medaka_consensus -i "!{ont_reads}" -d "!{assembly}" -o medaka \
      -t "!{task.cpus}" --bacteria
    test -s medaka/consensus.fasta
    cp medaka/consensus.fasta "${sample_id}.medaka.fasta"
    printf 'sample_id\tstage\tstatus\treason\n%s\tMEDAKA\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.medaka.stage.tsv"
    {
      printf 'medaka\t'
      medaka --version 2>&1 | head -n 1 || true
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.medaka.software_versions.txt"
    completion_state.py create --record "${sample_id}.MEDAKA_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage MEDAKA --status PASS \
      --input "!{assembly}" --input "!{ont_reads}" \
      --output "${sample_id}.medaka.fasta" --output "${sample_id}.medaka.stage.tsv" \
      --output "${sample_id}.medaka.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process BWA_MEM2_ALIGN {
    tag "${meta.sample_id}"
    label 'polish_heavy'
    container params.bwa_mem2_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/02_hybrid_assembly/03_polypolish" },
        mode: 'copy', overwrite: true,
        saveAs: { name ->
            (name.endsWith('.sam') || name.endsWith('.polypolish_input.fasta')) ? null : name
        }

    input:
    tuple val(meta), path(assembly), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample_id}.polypolish_input.fasta"),
        path("${meta.sample_id}.r1.sam"), path("${meta.sample_id}.r2.sam"), emit: alignments
    tuple val(meta), path("${meta.sample_id}.bwa_mem2.stage.tsv"), emit: status
    path "${meta.sample_id}.bwa_mem2.software_versions.txt", emit: versions
    path "${meta.sample_id}.BWA_MEM2_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    cp "!{assembly}" "${sample_id}.polypolish_input.fasta"
    bwa-mem2 index "${sample_id}.polypolish_input.fasta"
    bwa-mem2 mem -t "!{task.cpus}" -a "${sample_id}.polypolish_input.fasta" \
      "!{r1}" > "${sample_id}.r1.sam"
    bwa-mem2 mem -t "!{task.cpus}" -a "${sample_id}.polypolish_input.fasta" \
      "!{r2}" > "${sample_id}.r2.sam"
    test -s "${sample_id}.r1.sam"
    test -s "${sample_id}.r2.sam"
    printf 'sample_id\tstage\tstatus\treason\n%s\tBWA_MEM2\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.bwa_mem2.stage.tsv"
    {
      printf 'bwa-mem2\t'
      bwa-mem2 version 2>&1 | head -n 1 || true
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.bwa_mem2.software_versions.txt"
    completion_state.py create --record "${sample_id}.BWA_MEM2_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage BWA_MEM2 --status PASS \
      --input "!{assembly}" --input "!{r1}" --input "!{r2}" \
      --output "${sample_id}.r1.sam" --output "${sample_id}.r2.sam" \
      --output "${sample_id}.bwa_mem2.stage.tsv" --output "${sample_id}.bwa_mem2.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process POLYPOLISH {
    tag "${meta.sample_id}"
    label 'polish_heavy'
    container params.polypolish_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/02_hybrid_assembly/03_polypolish" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly), path(r1_sam), path(r2_sam)

    output:
    tuple val(meta), path("${meta.sample_id}.polypolish.fasta"), emit: assembly
    tuple val(meta), path("${meta.sample_id}.polypolish.stage.tsv"), emit: status
    path "${meta.sample_id}.polypolish.software_versions.txt", emit: versions
    path "${meta.sample_id}.POLYPOLISH_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"

    polypolish filter --in1 "!{r1_sam}" --in2 "!{r2_sam}" \
      --out1 r1.filtered.sam --out2 r2.filtered.sam
    polypolish polish "!{assembly}" r1.filtered.sam r2.filtered.sam \
      > "${sample_id}.polypolish.fasta"
    test -s "${sample_id}.polypolish.fasta"
    printf 'sample_id\tstage\tstatus\treason\n%s\tPOLYPOLISH\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.polypolish.stage.tsv"
    {
      printf 'polypolish\t'
      polypolish --version 2>&1 | head -n 1 || true
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.polypolish.software_versions.txt"
    completion_state.py create --record "${sample_id}.POLYPOLISH_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage POLYPOLISH --status PASS \
      --input "!{assembly}" --input "!{r1_sam}" --input "!{r2_sam}" \
      --output "${sample_id}.polypolish.fasta" --output "${sample_id}.polypolish.stage.tsv" \
      --output "${sample_id}.polypolish.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process PYPOLCA {
    tag "${meta.sample_id}"
    label 'polish_heavy'
    container params.pypolca_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/02_hybrid_assembly/04_pypolca" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.sample_id}.polished.fasta"), emit: assembly
    tuple val(meta), path("${meta.sample_id}.pypolca.stage.tsv"), emit: status
    path "${meta.sample_id}.pypolca.software_versions.txt", emit: versions
    path "${meta.sample_id}.PYPOLCA_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    pypolca run -a "!{assembly}" -1 "!{r1}" -2 "!{r2}" \
      -t "!{task.cpus}" -o pypolca -f -p "${sample_id}_pypolca" --careful
    corrected="pypolca/${sample_id}_pypolca_corrected.fasta"
    test -s "${corrected}"
    cp "${corrected}" "${sample_id}.polished.fasta"
    printf 'sample_id\tstage\tstatus\treason\n%s\tPYPOLCA\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.pypolca.stage.tsv"
    {
      printf 'pypolca\t'
      pypolca --version 2>&1 | head -n 1 || true
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.pypolca.software_versions.txt"
    completion_state.py create --record "${sample_id}.PYPOLCA_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage PYPOLCA --status PASS \
      --input "!{assembly}" --input "!{r1}" --input "!{r2}" \
      --output "${sample_id}.polished.fasta" --output "${sample_id}.pypolca.stage.tsv" \
      --output "${sample_id}.pypolca.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process QUAST_QC {
    tag "${meta.sample_id}"
    label 'qc_medium'
    container params.quast_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/03_assembly_qc/quast" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.sample_id}.quast.tsv"), emit: report
    tuple val(meta), path("${meta.sample_id}.quast.stage.tsv"), emit: status
    path "${meta.sample_id}.quast.software_versions.txt", emit: versions
    path "${meta.sample_id}.QUAST_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    quast.py "!{assembly}" -o quast --threads "!{task.cpus}"
    test -s quast/report.tsv
    cp quast/report.tsv "${sample_id}.quast.tsv"
    printf 'sample_id\tstage\tstatus\treason\n%s\tQUAST\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.quast.stage.tsv"
    {
      printf 'quast\t'
      quast.py --version 2>&1 | head -n 1
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.quast.software_versions.txt"
    completion_state.py create --record "${sample_id}.QUAST_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage QUAST --status PASS \
      --input "!{assembly}" --output "${sample_id}.quast.tsv" \
      --output "${sample_id}.quast.stage.tsv" --output "${sample_id}.quast.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process CHECKM2_QC {
    tag "${meta.sample_id}"
    label 'checkm2_heavy'
    container params.checkm2_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/03_assembly_qc/checkm2" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.sample_id}.checkm2.tsv"), emit: status
    tuple val(meta), path("${meta.sample_id}.ASSEMBLY_READY"), emit: ready, optional: true
    path "${meta.sample_id}.checkm2.software_versions.txt", emit: versions
    path "${meta.sample_id}.ASSEMBLY_QC_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    override="!{meta.assembly_qc_override}"

    mkdir -p checkm2_input
    ln -s "!{assembly}" "checkm2_input/${sample_id}.fasta"
    checkm2 predict --threads "!{task.cpus}" --input checkm2_input \
      --output-directory checkm2 --database_path "!{params.checkm2_db}" --force
    test -s checkm2/quality_report.tsv

    parse_checkm2_qc.py \
      --report checkm2/quality_report.tsv \
      --output "${sample_id}.checkm2.tsv" \
      --sample-id "${sample_id}" \
      --override "${override}" \
      --pass-completeness "!{params.assembly_pass_completeness}" \
      --min-completeness "!{params.assembly_min_completeness}" \
      --pass-contamination "!{params.assembly_pass_contamination}" \
      --max-contamination "!{params.assembly_max_contamination}"

    {
      printf 'checkm2\t'
      checkm2 --version 2>&1 | head -n 1 || true
      printf 'database\t%s\n' '!{params.checkm2_db}'
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.checkm2.software_versions.txt"
    stage_status=$(tail -n 1 "${sample_id}.checkm2.tsv" | cut -f3)
    completion_state.py create --record "${sample_id}.ASSEMBLY_QC_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage ASSEMBLY_QC --status "$stage_status" \
      --input "!{assembly}" --output "${sample_id}.checkm2.tsv" \
      --output "${sample_id}.checkm2.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process GTDBTK_CLASSIFY {
    tag "${meta.sample_id}"
    label 'gtdbtk_heavy'
    container params.gtdbtk_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/03_assembly_qc/gtdbtk" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.sample_id}.gtdbtk.tsv"), emit: status
    tuple val(meta), path("${meta.sample_id}.ECOLI_READY"), emit: ecoli, optional: true
    path "${meta.sample_id}.gtdbtk.software_versions.txt", emit: versions
    path "${meta.sample_id}.GTDBTK_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    mkdir -p genomes
    ln -s "!{assembly}" "genomes/${sample_id}.fasta"
    export GTDBTK_DATA_PATH="!{params.gtdbtk_db}"

    gtdbtk classify_wf --genome_dir genomes --extension fasta --out_dir gtdbtk \
      --prefix amrevo --cpus "!{task.cpus}" --pplacer_cpus 1
    summary="$(find gtdbtk -type f -name '*.bac120.summary.tsv' -print -quit)"
    test -s "${summary}"

    parse_gtdbtk.py \
      --summary "${summary}" \
      --output "${sample_id}.gtdbtk.tsv" \
      --sample-id "${sample_id}"

    {
      printf 'gtdbtk\t'
      gtdbtk --version 2>&1 | head -n 1
      printf 'database\t%s\n' '!{params.gtdbtk_db}'
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.gtdbtk.software_versions.txt"
    stage_status=$(tail -n 1 "${sample_id}.gtdbtk.tsv" | cut -f3)
    completion_state.py create --record "${sample_id}.GTDBTK_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage GTDBTK --status "$stage_status" \
      --input "!{assembly}" --output "${sample_id}.gtdbtk.tsv" \
      --output "${sample_id}.gtdbtk.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process MLST_ACHTMAN {
    tag "${meta.sample_id}"
    label 'typing_light'
    container params.mlst_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/04_typing/mlst" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.sample_id}.mlst.tsv"), emit: result
    tuple val(meta), path("${meta.sample_id}.mlst.stage.tsv"), emit: status
    path "${meta.sample_id}.mlst.software_versions.txt", emit: versions
    path "${meta.sample_id}.MLST_ACHTMAN_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    mlst --scheme ecoli_achtman_4 --nopath "!{assembly}" \
      > "${sample_id}.mlst.tsv"
    test -s "${sample_id}.mlst.tsv"
    printf 'sample_id\tstage\tstatus\treason\n%s\tMLST_ACHTMAN\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.mlst.stage.tsv"
    {
      printf 'mlst\t'
      mlst --version 2>&1 | head -n 1
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.mlst.software_versions.txt"
    completion_state.py create --record "${sample_id}.MLST_ACHTMAN_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage MLST_ACHTMAN --status PASS \
      --input "!{assembly}" --output "${sample_id}.mlst.tsv" --output "${sample_id}.mlst.stage.tsv" \
      --output "${sample_id}.mlst.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process AMRFINDERPLUS {
    tag "${meta.sample_id}"
    label 'typing_medium'
    container params.amrfinder_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/05_amr/amrfinderplus" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.sample_id}.amrfinderplus.tsv"), emit: result
    tuple val(meta), path("${meta.sample_id}.amrfinderplus.stage.tsv"), emit: status
    path "${meta.sample_id}.amrfinderplus.software_versions.txt", emit: versions
    path "${meta.sample_id}.AMRFINDERPLUS_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    amrfinder --nucleotide "!{assembly}" --organism Escherichia \
      --threads "!{task.cpus}" > "${sample_id}.amrfinderplus.tsv"
    test -s "${sample_id}.amrfinderplus.tsv"
    printf 'sample_id\tstage\tstatus\treason\n%s\tAMRFINDERPLUS\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.amrfinderplus.stage.tsv"
    {
      printf 'amrfinderplus\t'
      amrfinder --version 2>&1 | head -n 1
      amrfinder --database_version 2>&1 | sed 's/^/database\t/'
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
    } > "${sample_id}.amrfinderplus.software_versions.txt"
    completion_state.py create --record "${sample_id}.AMRFINDERPLUS_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage AMRFINDERPLUS --status PASS \
      --input "!{assembly}" --output "${sample_id}.amrfinderplus.tsv" \
      --output "${sample_id}.amrfinderplus.stage.tsv" --output "${sample_id}.amrfinderplus.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}


process KTYPR {
    tag "${meta.sample_id}"
    label 'typing_medium'
    container params.ktypr_image

    publishDir { "${params.outdir}/samples/isolates/${meta.sample_id}/06_capsule/ktypr" },
        mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.sample_id}.ktypr.tsv"), emit: result
    tuple val(meta), path("${meta.sample_id}.ktypr"), emit: details
    tuple val(meta), path("${meta.sample_id}.ktypr.stage.tsv"), emit: status
    path "${meta.sample_id}.ktypr.software_versions.txt", emit: versions
    path "${meta.sample_id}.KTYPR_COMPLETE.json", emit: completion

    shell:
    '''
    set -euo pipefail
    sample_id="!{meta.sample_id}"
    mkdir -p "${sample_id}.ktypr"
    ktypr -i "!{assembly}" -o "${sample_id}.ktypr" -m 0 -f 30000 \
      -n "!{task.cpus}" -p "${sample_id}_" -ko 1
    result="$(find "${sample_id}.ktypr" -type f -name '*_ktypr.tsv' -print \
      | LC_ALL=C sort | head -n 1)"
    test -n "${result}"
    test -s "${result}"
    cp "${result}" "${sample_id}.ktypr.tsv"
    printf 'sample_id\tstage\tstatus\treason\n%s\tKTYPR\tPASS\t\n' \
      "${sample_id}" > "${sample_id}.ktypr.stage.tsv"
    {
      printf 'ktypr\t1.0.1\n'
      printf 'mode\tflanking_30000bp\n'
      printf 'retry_token\t%s\n' '!{meta.retry_token}'
      if test -f /opt/ktypr_source_commit.txt; then
        printf 'source_commit\t'
        cat /opt/ktypr_source_commit.txt
      fi
    } > "${sample_id}.ktypr.software_versions.txt"
    completion_state.py create --record "${sample_id}.KTYPR_COMPLETE.json" \
      --sample-id "${sample_id}" --sample-class isolate --stage KTYPR --status PASS \
      --input "!{assembly}" --output "${sample_id}.ktypr.tsv" --output "${sample_id}.ktypr" \
      --output "${sample_id}.ktypr.stage.tsv" --output "${sample_id}.ktypr.software_versions.txt" \
      --parameter-hash '!{meta.parameter_hash}' --software-hash '!{meta.software_hash}' --run-id '!{meta.run_id}'
    '''
}
