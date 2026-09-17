# Development Roadmap

1. Validate the manifest and write resolved absolute input paths.
2. Add FastQC, NanoPlot and trimming-summary modules.
3. Add per-sample Unicycler and polishing modules.
4. Add QUAST or MetaQUAST, BUSCO and CheckM2 modules.
5. Add MLST and a final metadata/audit merge.
6. Add a small synthetic test dataset and continuous testing.
7. Freeze containers, tag a release, add a licence and archive the release.

The pipeline must run one sample row per task. An input ambiguity is an error,
not a reason to choose the latest matching directory.
