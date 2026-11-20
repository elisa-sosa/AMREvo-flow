# Local inputs

This directory is intentionally ignored except for this README. Keep the real
AMREvo manifests and any project-local raw-data staging below this directory;
none of it belongs in Git.

Expected manifest names are:

```text
isolates.samples.local.tsv
isolates.read_units.local.tsv
sweeps.samples.local.tsv
sweeps.read_units.local.tsv
qc_overrides.local.tsv        # optional
```

Paths in read-unit manifests may point to the legacy AMREvo sequencing area or
to future deliveries stored below `inputs/raw/`. The workflow reads those files
without editing, renaming, moving, or deleting them.
