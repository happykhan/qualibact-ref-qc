#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { SIMULATE_ART         } from './modules/simulate_reads'
include { ASSEMBLY_SHOVILL     } from './modules/short_read_assembly'
include { ASSEMBLY_STATS       } from './modules/assembly_stats'
include { CONTAMINATION_CHECKM } from './modules/contamination'
include { BASE_COMPOSITION     } from './modules/base_composition'
include { MERGE_METRICS        } from './modules/merge_metrics'
include { PREPARE_FASTA        } from './modules/prepare_fasta'

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // --- validate samplesheet early with a clear error ---
    def samplesheetPath = file(params.samplesheet ?: '')
    if (!samplesheetPath || !samplesheetPath.exists()) {
        error "Samplesheet not found or unreadable: ${params.samplesheet}\n" +
              "Working dir: ${new File('.').canonicalPath}\n" +
              "Please check the path and permissions."
    }
    log.info "Using samplesheet: ${samplesheetPath.toString()}"
    // Build channel from CSV, fail fast if rows missing or sample_id missing
    def assemblies = Channel
        .fromPath(samplesheetPath.toString(), checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sampleId = (row.sample_id ?: row.sample ?: row.id ?: "").toString().trim()
            if (!sampleId) error "Samplesheet row missing sample_id: ${row}"
            def fastaFile = file(row.fasta)
            if (!fastaFile.exists()) error "FASTA for sample ${sampleId} not found: ${row.fasta}"
            def meta = [sample_id: sampleId, type: 'assembly']
            tuple(meta, fastaFile)
        }
    assemblies
        | PREPARE_FASTA
        | set { prepared_assemblies }
    def metrics_input
    if (params.simulate) {
        log.info "Simulating reads with ART and assembling with Shovill"
        SIMULATE_ART(prepared_assemblies)
        def simulated_reads = SIMULATE_ART.out
            .map { meta, r1, r2 -> tuple(meta, r1, r2) }
        ASSEMBLY_SHOVILL(simulated_reads, params.min_contig_length)
        metrics_input = ASSEMBLY_SHOVILL.out
            .map { meta, fasta -> tuple(meta, fasta) }
    } else {
        log.info "Using provided assemblies for metrics only"
        metrics_input = prepared_assemblies
    }

    ASSEMBLY_STATS(metrics_input)
    CONTAMINATION_CHECKM(metrics_input)
    BASE_COMPOSITION(metrics_input)

    def stats_by_id = ASSEMBLY_STATS.out.map { meta, stats -> tuple(meta.sample_id, [meta, stats]) }
    def contam_by_id = CONTAMINATION_CHECKM.out.map { meta, contam -> tuple(meta.sample_id, contam) }
    def base_by_id = BASE_COMPOSITION.out.map { meta, base -> tuple(meta.sample_id, base) }

    stats_by_id
        .join(contam_by_id)
        .map { id, statsPair, contamFile -> tuple(id, [statsPair[0], statsPair[1], contamFile]) }
        .join(base_by_id)
        .map { id, merged, baseFile ->
            def meta = merged[0]
            def statsFile = merged[1]
            def contamFile = merged[2]
            tuple(meta, statsFile, contamFile, baseFile)
        }
        .set { merge_input }

    def asm_files = merge_input.map { meta, stats, contam, base -> stats }
    def checkm_files = merge_input.map { meta, stats, contam, base -> contam }
    def base_files = merge_input.map { meta, stats, contam, base -> base }

    asm_files.collect().set { asm_collected }
    checkm_files.collect().set { checkm_collected }
    base_files.collect().set { base_collected }

    def speciesName = file(params.samplesheet).getParent()?.getName() ?: "species"
    def metadataJson = file(params.samplesheet).getParent() ? file(params.samplesheet).getParent().resolve("assembly_data_report.json") : null

    // MERGE_METRICS(asm_collected, checkm_collected, base_collected, metadataJson, speciesName, file(params.samplesheet))
    // Need to work on this a bit more. Let's just get the inputs done for now.
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/
 workflow.onComplete {
        workDir = new File("${workflow.workDir}")

        println """
        Qualibact-ref Execution Summary
        ---------------------------
        Pipeline Version : ${workflow.manifest.version}
        Nextflow Version : ${nextflow.version}
        Command Line     : ${workflow.commandLine}
        Resumed          : ${workflow.resume}
        Completed At     : ${workflow.complete}
        Duration         : ${workflow.duration}
        Success          : ${workflow.success}
        Exit Code        : ${workflow.exitStatus}
        Error Report     : ${workflow.errorReport ?: '-'}
        Launch Dir       : ${workflow.launchDir}
        """
    }
