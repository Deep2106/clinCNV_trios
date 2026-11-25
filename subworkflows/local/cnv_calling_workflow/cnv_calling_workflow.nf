include { CLINCNV  } from '../../../modules/local/clincnv/clincnv'
include { MAKE_VCF } from '../../../modules/local/make_vcf/make_vcf'
include { SORT_VCF } from '../../../modules/local/sort_vcf/sort_vcf'

workflow CNV_CALLING_WORKFLOW {

    take:
    ch_merged_cov
    ch_annotated_bed
    ch_sample_meta
    ch_vcf_header
    ch_ped_file
    ch_fasta_fai

    main:
    ch_versions = Channel.empty()

    CLINCNV(
        ch_merged_cov,
        ch_annotated_bed,
        ch_ped_file
    )
    ch_versions = ch_versions.mix(CLINCNV.out.versions)

    // Trio mode only: cnv_output/trios/CHILD-MOTHER-FATHER/FAMILY_CNVs.txt
    // Proband = first element of directory name before first '-'
    // Singletons/duos: run through clinCNV_singleton pipeline separately
    ch_tsv = CLINCNV.out.cnv_output_dir
        .flatMap { dir ->
            def tsv_files = []
            def trios_dir = new File("${dir}/trios")
            if (trios_dir.exists()) {
                trios_dir.eachDir { family_dir ->
                    family_dir.eachFileMatch(~/.*_CNVs\.txt$/) { f ->
                        def proband_id = family_dir.name.split('-')[0]
                        tsv_files << [proband_id, 'trio', file(f.absolutePath)]
                    }
                }
            }
            return tsv_files
        }

    // Only probands (phenotype=2) match ch_sample_meta — parents silently excluded
    ch_cnv_matched = ch_tsv
        .join(
            ch_sample_meta.map { meta -> [meta.id, meta] },
            by: 0
        )
        .map { sid, mode, f, meta -> [meta, mode, f] }

    MAKE_VCF(
        ch_cnv_matched,
        ch_vcf_header,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(MAKE_VCF.out.versions)

    SORT_VCF(
        MAKE_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(SORT_VCF.out.versions)

    emit:
    vcf      = SORT_VCF.out.vcf
    versions = ch_versions
}
