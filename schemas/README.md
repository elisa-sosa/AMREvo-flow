# Manifest schemas

The JSON files describe individual TSV rows. The command-line validator also
enforces relationships that JSON Schema cannot express: sample/read-unit
membership, unique approved merge orders, platform-specific paths, path-root
containment, file integrity, mate identifiers, duplicate read identifiers and
completion-state consistency.

Real manifests use the ignored `.local.tsv` suffix and must not be committed.
