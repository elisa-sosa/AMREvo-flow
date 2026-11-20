process PREPARE_PHYLOGENY {
    tag 'current-ecoli-cohort'
    label 'reporting'
    container params.reporting_image
    errorStrategy 'finish'

    publishDir "${params.outdir}/cohort/current/phylogeny", mode: 'copy', overwrite: true

    input:
    path sample_status
    val results_root

    output:
    path 'phylogeny_jobs.tsv', emit: jobs
    path 'reference_manifest.tsv', emit: references
    path 'phylogeny_group_status.tsv', emit: group_status

    shell:
    '''
    prepare_phylogeny.py \
      --status "!{sample_status}" \
      --results-root "!{results_root}" \
      --output-dir . \
      --minimum-samples 3
    '''
}


process SNIPPY_CONTIGS {
    tag "${group_id}:${sample_id}"
    label 'phylogeny'
    container params.snippy_image
    errorStrategy 'finish'

    input:
    tuple val(group_id), val(sample_id), path(reference), path(assembly)

    output:
    tuple val(group_id), path("${group_id}__${sample_id}.snippy"), emit: result

    shell:
    '''
    set -euo pipefail
    snippy --cpus "!{task.cpus}" --outdir "!{group_id}__!{sample_id}.snippy" \
      --prefix "!{sample_id}" --ref "!{reference}" --ctgs "!{assembly}"
    test -s "!{group_id}__!{sample_id}.snippy/!{sample_id}.aligned.fa"
    '''
}


process SNIPPY_CORE_GROUP {
    tag "${group_id}"
    label 'phylogeny'
    container params.snippy_image
    errorStrategy 'finish'

    input:
    tuple val(group_id), path(snippy_dirs, stageAs: 'snippy??/*')

    output:
    tuple val(group_id), path("${group_id}.core.full.aln"), emit: alignment

    shell:
    '''
    set -euo pipefail
    snippy-core --prefix core snippy??/*
    test -s core.full.aln
    mv core.full.aln "!{group_id}.core.full.aln"
    '''
}


process GUBBINS_GROUP {
    tag "${group_id}"
    label 'phylogeny'
    container params.gubbins_image
    errorStrategy 'finish'

    publishDir { "${params.outdir}/cohort/current/phylogeny/groups/${group_id}" },
        mode: 'copy', overwrite: true

    input:
    tuple val(group_id), path(alignment)

    output:
    tuple val(group_id), path("${group_id}.phylogeny_status.tsv"), emit: status
    tuple val(group_id), path("${group_id}.final_tree.tre"), emit: tree, optional: true

    shell:
    '''
    set -euo pipefail
    variable_sites=$(python3 - "!{alignment}" <<'PY'
    import sys
    sequences = []
    current = []
    for line in open(sys.argv[1]):
        line = line.strip()
        if line.startswith('>'):
            if current: sequences.append(''.join(current))
            current = []
        else:
            current.append(line.upper())
    if current: sequences.append(''.join(current))
    count = 0
    if sequences:
        for column in zip(*sequences):
            states = {base for base in column if base not in {'N', '-', '?'}}
            if len(states) > 1: count += 1
    print(count)
    PY
    )
    if [[ "$variable_sites" -eq 0 ]]; then
      printf 'group_id\tstage\tstatus\tvariable_sites\treason\n%s\tGUBBINS\tNO_VARIABLE_SITES\t0\tcore alignment has no variable sites\n' \
        '!{group_id}' > "!{group_id}.phylogeny_status.tsv"
      exit 0
    fi
    run_gubbins.py --prefix "!{group_id}" --threads "!{task.cpus}" "!{alignment}"
    tree=$(find . -maxdepth 1 -type f -name '*final_tree.tre' -print -quit)
    test -s "$tree"
    cp "$tree" "!{group_id}.final_tree.tre"
    printf 'group_id\tstage\tstatus\tvariable_sites\treason\n%s\tGUBBINS\tCOMPLETED\t%s\t\n' \
      '!{group_id}' "$variable_sites" > "!{group_id}.phylogeny_status.tsv"
    '''
}


process PLOT_PHYLOGENY {
    tag "${group_id}"
    label 'reporting'
    container params.reporting_image
    errorStrategy 'finish'

    publishDir { "${params.outdir}/cohort/current/plots/trees" },
        mode: 'copy', overwrite: true
    publishDir { "${params.outdir}/cohort/snapshots/${run_id}/plots/trees" },
        mode: 'copy', overwrite: true

    input:
    tuple val(group_id), path(tree)
    val run_id

    output:
    path "${group_id}.tree.png", emit: plot

    shell:
    '''
    plot_tree.py --tree "!{tree}" --output "!{group_id}.tree.png" --title "AMREvo !{group_id}"
    '''
}
