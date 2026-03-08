/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_abrnaseq_pipeline'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'
include { SALMON_QUANT           } from '../modules/nf-core/salmon/quant/main'
include { DUPRADAR               } from '../modules/nf-core/dupradar/main'
include { QUALIMAP_RNASEQ        } from '../modules/nf-core/qualimap/rnaseq/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ABRNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:
    /*
    * MODULE SPECIFICATION
    ! FASTQC -> in: sample sheet ||| out: FASTQC.out.zip (QC plots)
    ! TRIMAGOLRE -> in: sample sheet ||| out: TRIMAGALORE.out.reads (Trimmed reads ), TRIMAGALORE.out.log (Log file)
    ! STAR_ALIGN -> in: Trimmed reads, STAR index, GTF, ignore GTF (false) ||| out: STAR_ALIGN.out.log_final (Log file), STAR_ALIGN.out.bam (Alignment BAM files)
    ! SALMON_QUANT -> in: Trimmed reads, GTF, Transcriptome, Alignment mode (false), Override library type (false) ||| out: SALMON_QUANT.out.results (Quantification results)
    ? DUPRADAR -> in: Alignment BAM files, GTF ||| out: DUPRADAR.out.multiqc (MultiQC files)
    ? QUALIMAP_RNASEQ -> in: Alignment BAM files, GTF ||| out: QUALIMAP_RNASEQ.out.results (Output data)
    */
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // * MODULE 1: FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{_meta, zip -> zip})

    //
    // * MODULE 2: TRIMGALORE
    //
    TRIMGALORE (
        ch_samplesheet
    )

    //
    // * MODULE 3: STAR_ALIGN
    // * Note: set params.star_index, params.gtf in nextflow.config
    //
    STAR_ALIGN (
        TRIMGALORE.out.reads,
        file(params.star_index),
        file(params.gtf),
        false
    )

    //
    // * MODULE 4: SALMAON_QUANT
    //
    SALMON_QUANT ( 
        TRIMGALORE.out.reads, 
        params.salmon_index,
        params.gtf, 
        params.transcriptome, 
        false, 
        false 
    )

    //
    // * MODULE 6: DUPRADAR
    //
    DUPRADAR ( STAR_ALIGN.out.bam, params.gtf )
    ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc.collect{ _meta, mqc -> mqc })


    //
    // * MODULE 7: QUALIMAP_RNASEQ
    //
    QUALIMAP_RNASEQ ( STAR_ALIGN.out.bam, params.gtf )
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results.collect{ _meta, res -> res })

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'abrnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
