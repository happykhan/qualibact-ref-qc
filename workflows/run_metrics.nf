// Run metrics directly on existing assemblies listed in a samplesheet.
// Expected samplesheet columns: sample (sample id) and fasta (path to assembly fasta).

include { CONTAMINATION_CHECKM } from '../modules/contamination'
include { ASSEMBLY_STATS       } from '../modules/assembly_stats'
include { BASE_COMPOSITION     } from '../modules/base_composition'
include { MERGE_METRICS        } from '../modules/merge_metrics'

workflow {

    // Build a channel of (meta, fasta) tuples from the samplesheet
    Channel
        .fromPath(params.samplesheet)
        .ifEmpty { error "Samplesheet not found: ${params.samplesheet}" }
        .splitCsv(header: true)
        .map { row ->
            def sampleId = (row.sample ?: row.sample_id ?: row.id ?: row.name ?: "").toString().trim()
            if (!sampleId) error "Samplesheet row missing sample id: ${row}"
            def fastaFile = file(row.fasta)
            if (!fastaFile.exists()) error "FASTA for sample ${sampleId} not found: ${row.fasta}"
            def meta = [sample_id: sampleId, type: 'assembly']
            tuple(meta, fastaFile)
        }
        .set { assemblies }

    // Metrics on existing assemblies
    ASSEMBLY_STATS(assemblies)
    CONTAMINATION_CHECKM(assemblies)
    BASE_COMPOSITION(assemblies)

    // Key each output by sample_id so results stay aligned
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

    MERGE_METRICS(merge_input)
}
