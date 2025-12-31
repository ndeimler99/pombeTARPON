#!/usr/bin/env nextflow

/* The following pipeline is intended for research purposes only */
nextflow.enable.dsl=2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

println """\
    pombeTARPON - Telomere Analysis Pipeline on Nanopore Sequencing Data
    ================================================
    v0.0.1
    """.stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validate_parameters } from "./subworkflows/parameter_validation.nf"
include { PREPROCESS_FILES } from "./subworkflows/file_preprocessing.nf"
//include { isolate_telomeres } from "./subworkflows/isolate_telomeres.nf"
//include { analyze_telomeres,cluster_telomeres } from "./subworkflows/telomere_analysis.nf"
include { paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { ALIGNMENT } from "./bin/process.nf"
include { TELOMERE_STATS } from "./subworkflows/telomere_stats.nf"
include { GENERATE_HTML_REPORT } from "./bin/process.nf"
include { getParams; getVersions; getManifest } from "./bin/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//WorkflowMain.initialise(workflow, params, log)

workflow {

    //parameter validation
        // input file is provided (or directory)
    valid_params = validate_parameters()

    if (valid_params.passed.value == false){
        exit 1, "Parameter Validation Failed"
    }

    parameters = getParams()
    versions = getVersions()
    manifest = getManifest()

    // process input files and isolate telomeric sequences
        // returns channel filled with sample - putative read pairs
    preprocess_out = PREPROCESS_FILES()

    // align to reference genome for html report
    alignment = ALIGNMENT(preprocess_out.input, file(params.pombe_genome))

    // take putative reads and identify telo start, telo end, and filter low quality telo sequences
    telo_results = TELOMERE_STATS(preprocess_out.reverse_complemented_reads)

    // generate html report
    //alignment.alignment.collect().view()
    //telo_results.telo_stats.collect().view()
    GENERATE_HTML_REPORT(parameters.params, versions.versions, manifest.manifest, \
                        alignment.alignment.collect(), \
                        telo_results.telo_stats.collect())







    // for each sample
        // isolate telomeric sequences
    
        // identify end of telomere
            // one python script with multiple functions depending on input paramters sample files etc
        
        // identify start of telomere
        // filter telomeres by quality/composition
        // cluster using meshclust
        // blast against TAS/rDNA or custom parameter file

        // create all relevant plots in R/python
            // bulk telomere length histograms
            // bulk telomere length barcharts similar to that of TARPON
            // chromosome cluster specific telomere length boxplot
            // pycairo cluster charts
    // generate html report
    




    // to do:
        // add strand comparison options

}

// When workflow finishes return basic description of finished or not and if it works remove the work directory if specified during run
workflow.onComplete {
    println "Analysis Complete at: $workflow.complete"
    println "Execution Status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Open the Following Report in your Browser ${ params.outdir }/report.html"

    if (workflow.success){
        if (params.remove_wd) {
            "rm -rf ${baseDir}/work".execute()
        }
       
    }
    Pinguscript.ping_complete(nextflow, workflow, params)
}

workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
