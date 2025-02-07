#!/usr/bin/env nextflow

/* The following pipeline is intended for research purposes only */
nextflow.enable.dsl=2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

println """\
    TArPON - Telomere Analysis Pipeline on Nanopore Sequencing Data
    ================================================
    v0.0.1
    """.stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validate_parameters,modify_file_type } from "./subworkflows/parameter_validation.nf"
include { isolate_telomeres } from "./subworkflows/isolate_telomeres.nf"
include { analyze_telomeres,cluster_telomeres } from "./subworkflows/telomere_analysis.nf"
include { paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowMain.initialise(workflow, params, log)

workflow {

    //parameter validation
        // input file is provided (or directory)
    valid_params = validate_parameters()

    if (valid_params.passed.value == false){
        exit 1, "Parameter Validation Failed"
    }

    // if input is directory and --demux is true
        // do not combine files
        // convert each file to bam if fastq
        // pass each file into telomere pipeline

    // if input is single bam file and --nanopore_barcodes is true
        // demux using nanopore dorado demux
    
    // if input is single bam file and barcodes file is provided and adaptor_sequence is none
        // convert file to bam if fastq
        // pass file into telomere pipeline

    // if input is single bam file and barcodes file is provided and adaptor_sequence is provided
        // convert file to bam if fastq
        // pass single bam file into telomere pipeline

    
    // if demux is true or dorado demux has already been run (--nanopore_barcodes is true)
        // align to reference genome
        // isolate telomeric sequences
        // if adaptor is not none or sample file is not none
            // identify end of telomere
        // else
            // identify end of teloemre using regex

    // else
        // align bulk sequencing to reference genome
        // isolate telomeric sequences
        // if barcodes file is provided
            // if adaptor sequence is provided
                // identify adaptor sequence and then demux
            // else
                // demux by barcodes at end of telomere
    
    // identify start of telomeres
    // filter telomeres by composition
    //cluster using meshclust
    // blast TAS/rDNA elements
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
