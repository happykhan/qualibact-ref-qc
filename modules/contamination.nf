process CONTAMINATION_CHECKM {
    tag { meta.sample_id }
    label 'process_medium'
    container 'happykhan/checkm2:0.1.0'
    
    publishDir "${params.outdir}/checkm_summary", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(report)

    script:
    outdir="checkm_out"
    report="${meta.sample_id}.${meta.type}.tsv"

    """
    checkm2 predict -i . -o ${outdir} -x fasta --force -t $task.cpus 
    mv ${outdir}/quality_report.tsv ${report}
    """
}
