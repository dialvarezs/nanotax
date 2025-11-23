include { EMU_ABUNDANCE     } from '../../../modules/nf-core/emu/abundance/main'
include { MMSEQS_EASYSEARCH } from '../../../modules/nf-core/mmseqs/easysearch/main'
include { MMSEQS_SUMMARISE  } from '../../../modules/local/mmseqs_summarise/main'

workflow TAXONOMIC_ASSIGNMENT {
    take:
    ch_reads // channel: [ val(meta), path(fastq) ]
    ch_emu_db
    ch_mmseqs2_db
    val_skip_emu
    val_skip_mmseqs2
    val_min_alignment_length
    val_min_identity

    main:
    ch_versions = channel.empty()

    /*
     * EMU
     */
    if (!val_skip_emu) {
        EMU_ABUNDANCE(ch_reads, ch_emu_db.map { _meta, file -> file })
        ch_versions = ch_versions.mix(EMU_ABUNDANCE.out.versions)
    }

    /*
     * MMseqs2
     */
    if (!val_skip_mmseqs2) {
        MMSEQS_EASYSEARCH(ch_reads, ch_mmseqs2_db)
        ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions)

        MMSEQS_SUMMARISE(MMSEQS_EASYSEARCH.out.tsv, val_min_identity, val_min_alignment_length)
        ch_versions = ch_versions.mix(MMSEQS_SUMMARISE.out.versions)
    }

    // // Taxonomic assignment
    // MMSEQS_EASYSEARCH(ch_input_tax,MMSEQS_CREATE16SDB.out.path_db)
    // ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions.first())

    // ch_mmseqs_output = MMSEQS_EASYSEARCH.out.tsv//.map{meta,tsv -> tsv}.collect()
    // //ch_first_group = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> meta}.first()
    // ch_groups_info = MMSEQS_EASYSEARCH.out.tsv.map{meta,tsv -> "${meta.id}:${meta.group}"}.collect()
    // ch_samples = ch_samplesheet.map{meta,path-> "${meta.id}"}.collect()
    // SUMMARY_MMSEQS(ch_mmseqs_output,ch_samples)
    // MERGE_AND_GROUP_SAMPLES(SUMMARY_MMSEQS.out.summary_csv.collect())//, SUMMARY_MMSEQS.out.abundance_picrust.collect())

    emit:
    versions = ch_versions
}
