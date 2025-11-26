process SIMULATE_ART {
    tag { meta.sample_id }

    def containerImage = 'quay.io/biocontainers/art:3.11.14--h2d50403_1'
    def resolvedContainer = (workflow.profile ?: '')
        .tokenize(',')
        .contains('bmrc') ? "${params.singularity_dir}/${containerImage.replace('/', '_').replace(':', '_')}.sif" : containerImage
    container resolvedContainer

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(reads1), path(reads2)

    script:
    prefix = "${meta.sample_id}"
    """
    art_illumina \\
        -ss HS25 \\
        -i ${fasta} \\
        -l 150 \\
        -f 40 \\
        -m 350 \\
        -s 50 \\
        -o ${prefix}

    mv ${prefix}1.fq ${prefix}_R1.fq
    mv ${prefix}2.fq ${prefix}_R2.fq
    gzip ${prefix}_R1.fq ${prefix}_R2.fq
    """
}
