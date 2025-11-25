process ASSEMBLY_STATS {
    tag { meta.sample_id }
    label 'process_single'
    container 'quay.io/biocontainers/assembly-stats:1.0.1--h9948957_10'

    publishDir "${params.outdir}/assembly_stats", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(report)

    script:
    report="${meta.sample_id}.${meta.type}.assembly_stats.txt"

    """
    assembly-stats -t "${fasta}" > ${report}
    """
}
