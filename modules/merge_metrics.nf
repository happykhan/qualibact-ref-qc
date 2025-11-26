process MERGE_METRICS {
    label 'python_container'
    label 'pixi'
    publishDir "${params.outdir}/metrics_summary", mode: 'copy'

    input:
    path assembly_stats_files
    path checkm_files
    path base_comp_files
    val metadata_json
    val species_name
    val samplesheet_path

    output:
    path "merged_metrics.tsv"
    path "${species_name}_metrics.csv"
    path "summary.csv"
    path "genome_size_histogram.csv"
    path "cds_vs_genome_size.csv"
    path "genome_size_histogram.png", optional: true
    path "cds_vs_genome_size.png", optional: true

    script:
    """
    python ${projectDir}/bin/merge_metrics.py \\
      --assembly-stats ${assembly_stats_files} \\
      --base-comp ${base_comp_files} \\
      --checkm ${checkm_files} \\
      --species "${species_name}" \\
      --samplesheet "${samplesheet_path}" \\
      ${metadata_json ? "--metadata-json "+metadata_json : ""} \\
      --outdir .
    """
}
