process ASSEMBLY_STATS {
    tag { meta.sample_id }
    label 'process_single'

    def containerImage = 'quay.io/biocontainers/assembly-stats:1.0.1--h9948957_10'
    def resolvedContainer = (workflow.profile ?: '')
        .tokenize(',')
        .contains('bmrc') ? "${params.singularity_dir}/${containerImage.replace('/', '_').replace(':', '_')}.sif" : containerImage
    container resolvedContainer

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
