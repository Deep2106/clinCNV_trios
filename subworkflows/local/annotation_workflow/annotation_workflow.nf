include { VEP_ANNOTATE } from '../../../modules/local/vep_annotate/vep_annotate'

workflow ANNOTATION_WORKFLOW {

    take:
    ch_vcf            // channel: [ meta, vcf ]
    ch_vep_cache      // channel: value(path)
    ch_vep_plugins    // channel: value(path)
    ch_sv_overlap     // channel: value(vcf.gz)
    ch_sv_overlap_tbi // channel: value(vcf.gz.tbi)
    ch_gene_panel     // channel: value(path) or NO_FILE  unused, kept for interface compatibility
    ch_rscript        // channel: value(path) or NO_FILE  unused, kept for interface compatibility

    main:
    ch_versions = Channel.empty()
    ch_vep_vcf  = Channel.empty()

    if (params.run_vep) {
        VEP_ANNOTATE(
            ch_vcf,
            ch_vep_cache,
            ch_vep_plugins,
            ch_sv_overlap,
            ch_sv_overlap_tbi
        )
        ch_vep_vcf  = VEP_ANNOTATE.out.vcf
        ch_versions = ch_versions.mix(VEP_ANNOTATE.out.versions)
    }

    emit:
    vep_vcf  = ch_vep_vcf
    versions = ch_versions
}
