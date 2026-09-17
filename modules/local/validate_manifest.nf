process VALIDATE_MANIFEST {
    tag 'validate-manifest'

    input:
    path samples
    val illumina_root
    val ont_root

    output:
    path 'validated_samples.tsv', emit: manifest

    script:
    """
    validate_manifest.py \\
      --samples "$samples" \\
      --illumina-root "$illumina_root" \\
      --ont-root "$ont_root" \\
      --output validated_samples.tsv
    """
}
