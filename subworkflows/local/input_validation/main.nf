/*
========================================================================================
    INPUT_VALIDATION SUBWORKFLOW
========================================================================================
    Detects cohort type from samplesheet headers and auto-generates PED file
    Supported formats:
        Singleton : sample_id,bam,bai,sex
        Mixed     : family_id,sample_id,paternal_id,maternal_id,sex,phenotype,bam,bai
----------------------------------------------------------------------------------------
*/

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheet_check/main'

workflow INPUT_VALIDATION {

    take:
    samplesheet // file: path to input CSV

    main:

    ch_versions = Channel.empty()

    SAMPLESHEET_CHECK ( samplesheet )
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    // - Parse validated CSV into sample channel -------------─
    ch_samples = SAMPLESHEET_CHECK.out.csv
        .splitCsv(header: true)
        .map { row -> _parse_row(row) }

    // - Emit BAM channel -------------------------
    ch_bam = ch_samples
        .map { meta, bam, bai -> [ meta, file(bam), file(bai) ] }

    // - PED as value channel — shared safely across all downstream processes
    ch_ped = SAMPLESHEET_CHECK.out.ped.first()

    // - Cohort type: read file content as string -------------
    ch_cohort_type = SAMPLESHEET_CHECK.out.cohort_type
        .map { f -> f.text.trim() }

    emit:
    bam         = ch_bam          // [ meta, bam, bai ]
    ped         = ch_ped          // path: cohort.ped (value channel)
    cohort_type = ch_cohort_type  // string: "singleton" or "mixed"
    versions    = ch_versions
}

// - Private helper: build meta map from CSV row -------------─
def _parse_row(row) {
    def meta = [:]
    meta.id         = row.sample_id
    meta.sex        = row.sex ? row.sex.toInteger() : 0
    meta.family_id  = row.family_id  ?: row.sample_id
    meta.paternal   = row.paternal_id ?: '0'
    meta.maternal   = row.maternal_id ?: '0'
    meta.phenotype  = row.phenotype  ? row.phenotype.toInteger() : 1
    return [ meta, row.bam, row.bai ]
}

