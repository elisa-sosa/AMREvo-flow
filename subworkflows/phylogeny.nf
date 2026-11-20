include { PREPARE_PHYLOGENY; SNIPPY_CONTIGS; SNIPPY_CORE_GROUP; GUBBINS_GROUP; PLOT_PHYLOGENY } from '../modules/phylogeny'


workflow PHYLOGENY {
    take:
    sample_status
    results_root
    run_id

    main:
    PREPARE_PHYLOGENY(sample_status, results_root)
    jobs = PREPARE_PHYLOGENY.out.jobs
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.group_id,
                row.sample_id,
                file(row.reference, checkIfExists: true),
                file(row.assembly, checkIfExists: true),
            )
        }
    SNIPPY_CONTIGS(jobs)
    grouped = SNIPPY_CONTIGS.out.result
        .map { group_id, directory -> tuple(group_id, directory) }
        .groupTuple()
    SNIPPY_CORE_GROUP(grouped)
    GUBBINS_GROUP(SNIPPY_CORE_GROUP.out.alignment)
    PLOT_PHYLOGENY(GUBBINS_GROUP.out.tree, run_id)

    emit:
    references = PREPARE_PHYLOGENY.out.references
    initial_status = PREPARE_PHYLOGENY.out.group_status
    final_status = GUBBINS_GROUP.out.status
    plots = PLOT_PHYLOGENY.out.plot
}
