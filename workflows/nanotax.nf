/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { MERGE_AND_GROUP_SAMPLES } from '../modules/local/mergeandgroupsamples'
include { PLOT_CORE               } from '../modules/local/plotcore'
include { PLOT_TAXONOMY           } from '../modules/local/plottaxonomy'
include { DIVERSITY               } from '../modules/local/diversity'
include { SEQKIT                  } from '../modules/local/seqkit'
include { SEQKIT_GREP             } from '../modules/nf-core/seqkit/grep'
include { SEQKIT_SEQ              } from '../modules/nf-core/seqkit/seq'
include { SEQKIT_FQ2FA            } from '../modules/nf-core/seqkit/fq2fa'
include { OBTAIN_IDS              } from '../modules/local/obtainids'
include { MERGE_PICRUST_OUT       } from '../modules/local/mergepicrustout'
include { LEFSE                   } from '../modules/local/lefse'
include { PLOT_LEFSE              } from '../modules/local/plotlefse'

include { BASECALLING             } from '../subworkflows/local/basecalling/main'
include { PREPARE_DATABASES       } from '../subworkflows/local/prepare_databases/main'
include { QUALITY_CONTROL         } from '../subworkflows/local/quality_control/main'
include { TAXONOMIC_ASSIGNMENT    } from '../subworkflows/local/taxonomic_assignment/main'

include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_nanotax_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOTAX {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    /*
     * Database preparation
     */
    PREPARE_DATABASES(
        params.emu_db_name,
        params.mmseqs2_db_name,
        params.skip_emu,
        params.skip_mmseqs2,
    )

    /*
     * Basecalling and demultiplexing using Dorado
     */
    if (!params.skip_basecalling) {
        ch_pod5_dir = channel.value(
            [[id: 'reads'], file(params.dorado_pod5_dir, checkIfExists: true, type: 'dir')]
        )

        BASECALLING(ch_pod5_dir, ch_samplesheet, params.dorado_barcoding_kit)

        ch_samples = BASECALLING.out.samples
    }
    else {
        ch_samples = ch_samplesheet
    }

    /*
     * Quality control and filtering
     */
    if (!params.skip_qc) {
        QUALITY_CONTROL(ch_samples, params.filtlong_sampling)
    }

    /*
     * Taxonomic assignment
     */
    TAXONOMIC_ASSIGNMENT(
        QUALITY_CONTROL.out.reads,
        PREPARE_DATABASES.out.emu_database,
        PREPARE_DATABASES.out.mmseqs2_database,
        params.skip_emu,
        params.skip_mmseqs2,
        params.mmseqs2_min_aln,
        params.mmseqs2_min_identity,
    )

    // // Plots for Taxonomic assignment
    // ch_groups = ch_samplesheet.map{meta,path-> "${meta.id}:${meta.group}"}.collect()
    // PLOT_TAXONOMY((MERGE_AND_GROUP_SAMPLES.out.csv_sample.mix(MERGE_AND_GROUP_SAMPLES.out.csv_group)).flatten()) //csv_group
    // PLOT_CORE(MERGE_AND_GROUP_SAMPLES.out.csv_core,SUMMARY_MMSEQS.out.taxlineage.collect(),ch_groups)

    // // Diversity
    // // ToDo: Solo si hay grupos
    // if(!params.skip_diversity){
    //     ch_groups_info_all = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> "${meta.id}:${meta.group}:${meta.subgroup}:${meta.subsubgroup}"}.collect()
    //     DIVERSITY(MERGE_AND_GROUP_SAMPLES.out.csv_div_nreads,ch_groups_info_all)//ch_groups)
    //     ch_versions = ch_versions.mix(DIVERSITY.out.versions.first())
    // }
    // // Functional prediction
    // if(params.skip_functional_prediction){
    //     ch_input_picrust = (SUMMARY_MMSEQS.out.abundance_picrust.join(ch_input_tax)).map{meta,tsv,fastq -> [tsv,fastq]}
    //     SEQKIT(ch_input_picrust) //ch_input_tax.map{meta, path -> path}.collect(),SUMMARY_MMSEQS.out.abundance_picrust.collect())
    //     (SEQKIT.out.abundance, SEQKIT.out.fasta)
    //     ch_versions = ch_versions.mix(PICRUST2.out.versions.first())

    //     MERGE_PICRUST_OUT(PICRUST2.out.dir.collect(), ch_groups)
    //     LEFSE(MERGE_PICRUST_OUT.out.lefse_input.flatten())
    //     PLOT_LEFSE(LEFSE.out.lefse_output.flatten())
    //     // ToDo:LEFSE SOLO SI HAY GRUPOS
    // }
    //     OBTAIN_IDS(SUMMARY_MMSEQS.out.abundance_picrust)
    //     ch_input_seqkit =  OBTAIN_IDS.out.abundance.join(ch_input_tax)
    //                         .join(OBTAIN_IDS.out.ids)
    //                         .multiMap{meta,tsv,fastq,ids ->
    //                             sequences: [meta, fastq]
    //                             ids: [ids]
    //     }
    //     SEQKIT_GREP(ch_input_seqkit.sequences,ch_input_seqkit.ids)
    //     SEQKIT_FQ2FA(SEQKIT_GREP.out.filter)
    //     SEQKIT_SEQ(SEQKIT_FQ2FA.out.fasta)
    //     ch_input_picrust = OBTAIN_IDS.out.abundance.join( SEQKIT_SEQ.out.fastx)
    //                         .multiMap{meta,abundance_table,fasta->
    //                         abundance_table: [meta, abundance_table]
    //                         fasta: [fasta]
    //                         }
    //     PICRUST2(ch_input_picrust.abundance_table,ch_input_picrust.fasta)
    //     ch_versions = ch_versions.mix(PICRUST2.out.versions.first())
    //     MERGE_PICRUST_OUT(PICRUST2.out.dir.collect(), ch_groups)
    //     LEFSE(MERGE_PICRUST_OUT.out.lefse_input.flatten())

    //     // ToDo:LEFSE: SOLO SI HAY GRUPOS ; version
    // }

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
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nanotax_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'nanotax'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
