//import modules for the short read only assembly workflow


include { ASSEMBLY_SHOVILL               } from '../modules/short_read_assembly'
include { CONTAMINATION_CHECKM           } from '../modules/contamination'
include { ASSEMBLY_STATS               } from '../modules/assembly_stats'
include { BASE_COMPOSITION              } from '../modules/base_composition'
include { MERGE_METRICS                 } from '../modules/merge_metrics'

workflow RUN_ASSEMBLY{

    take:

    //take the short reads read channel from the main
    srt_reads

    //main workflow for short read assembly
    main:
    // Simulate reads using ART 

    // Do assembly using shovill
    ASSEMBLY_SHOVILL(srt_reads, params.min_contig_length)
    
    // Get metrics from assembly_stats
    ASSEMBLY_STATS(ASSEMBLY_SHOVILL.out)
    
    // Get contamination check checkm
    CONTAMINATION_CHECKM(ASSEMBLY_SHOVILL.out)

    // Get base composition with bin/count_bases.sh 
    BASE_COMPOSITION(ASSEMBLY_SHOVILL.out)

    // I think we need a final script to merge this all into a dataframe (that we then pass to qualibact)
    MERGE_METRICS(ASSEMBLY_STATS.out, CONTAMINATION_CHECKM.out, BASE_COMPOSITION.out)


}
