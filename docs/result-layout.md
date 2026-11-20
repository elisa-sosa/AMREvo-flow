# Result layout

Source code is kept shallow, but results need limited nesting so approximately
1,614 isolates and the sweep cohort do not write thousands of ambiguous files
into one directory. The hierarchy answers three questions in order: sample or
cohort, isolate or sweep, and analytical stage.

```text
results/
├── samples/
│   ├── isolates/<sample_id>/
│   │   ├── 01_read_qc/
│   │   ├── 02_hybrid_assembly/
│   │   ├── 03_assembly_qc/
│   │   ├── 04_typing/
│   │   ├── 05_amr/
│   │   ├── 06_capsule/
│   │   ├── sample_status.tsv
│   │   ├── software_versions.txt
│   │   ├── FINAL_COMPLETE.json
│   │   └── history/<old-input-fingerprint>-<run_id>/
│   └── sweeps/<sample_id>/
│       ├── 01_read_qc/
│       ├── 02_gemsweep/
│       ├── sample_status.tsv
│       ├── software_versions.txt
│       └── FINAL_COMPLETE.json
└── cohort/
    ├── current/
    │   ├── sample_status.tsv
    │   ├── run_summary.txt
    │   ├── software_versions.tsv
    │   ├── plots/
    │   ├── technical_qc/
    │   ├── phylogeny/
    │   ├── sweeps/
    │   └── provenance/
    └── snapshots/<run_id>/
        ├── sample_status.tsv
        ├── run_summary.txt
        ├── software_versions.tsv
        ├── manifests/
        ├── provenance/
        └── plots/
```

The numbered sample folders mirror the order of processing and make recovery
and inspection predictable. `current/` is the durable view for day-to-day use;
dated snapshots retain the normalized manifests, checksums, configuration,
status, provenance, software table, and PNGs used for that report.

Raw reads are never copied into this tree. If approved inputs change, the old
derived sample result is moved under that sample's checksum-derived `history/`
location before replacement results are published. Git ignores the entire
`results/` directory.
