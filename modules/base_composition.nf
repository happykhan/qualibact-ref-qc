process BASE_COMPOSITION {
    tag { meta.sample_id }

    label 'process_mini'

    publishDir "${params.outdir}/base_composition", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(report)

    script:
    report="${meta.sample_id}.${meta.type}.base_comp.tsv"

    """
    ${projectDir}/bin/count_bases.sh "${fasta}" > ${report}
    """
}
