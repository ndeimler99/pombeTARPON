/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { IDENTIFY_TELO_CORDINATES } from "../bin/process.nf"
include { GET_STATS_FROM_BAM } from "../bin/process.nf"
include { CONVERT_BAM_TO_FASTA_AND_TRIM } from "../bin/process.nf"
include { BLAST } from "../bin/process.nf"
include { MESHCLUST } from "../bin/process.nf"
include { PLOT_CLUSTERS } from "../bin/process.nf"
include { CONSENSUS_SEQ } from "../bin/process.nf"
include { GENERATE_R_PLOTS } from "../bin/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/ 

workflow TELOMERE_STATS {

    take:
        telo_reads

    main:
        Pinguscript.ping_start(nextflow, workflow, params)
        print("Telomere Analysis")

        telo_coordinates = IDENTIFY_TELO_CORDINATES(telo_reads)
        telo_stats = GET_STATS_FROM_BAM(telo_coordinates.coordinates_identified)

        telo_fasta = CONVERT_BAM_TO_FASTA_AND_TRIM(telo_stats.telo_bam, telo_stats.telo_stats)

        blast_results = BLAST(telo_fasta.telo_fasta, file(params.rDNA_TAS_file))

        meshclust_results = MESHCLUST(telo_fasta.trimmed_telo_fasta)

        telo_stats.telo_stats.mix(telo_fasta.telo_fasta, blast_results.blast_results, meshclust_results.clusters).
            groupTuple().
            map { prefix, files -> tuple (prefix, files.sort{it.name})}.
            set{ merged_fh_ch }

        PLOT_CLUSTERS(merged_fh_ch)

        if (params.consensus == true){
            CONSENSUS_SEQ(merged_fh_ch, file(params.rDNA_TAS_file)) //- create via python script, blast, and plot
        }  

        GENERATE_R_PLOTS(PLOT_CLUSTERS.out.stats)
        // create plots and figures
            // bulk telomere length histogram
            // stranded bulk telo length histogram
            // by cluster telomere length histogram
            // stranded by cluster telomere length histogram
            // bulk telomere length barplot
            // stranded telomere length barplot
            // by cluster telomere length barplot
            // stranded by cluster telomere length barplot
            // C vs G strand telo quality boxplot
            // C vs G strand read quality boxplot
            // read len vs telo length
            // read qual vs telo quality
            // read ln vs read qual
            // telo len vs telo qual


            // telomere read compostion (repeat percentage GGTTAC vs GGTTACA etc)
            // Slippage of telomeric repeats


            // sample comparison plots
                // bulk telomere length
                // read counts (stranded) as different colors

    //emit:
        //telo_bam = telo_stats.coordinates_identified
        //no_telo_coordinates = telo_stats.no_coordinates_identified
        //telo_stats = telo_stats_file_ch.telo_stats
}