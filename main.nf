#!/usr/bin/env nextflow

/*
========================================================================================
    CNV CALLING PIPELINE (ClinCNV WES) v2.0 — Mixed Cohort
========================================================================================
    Production pipeline:
    - Auto-detects singleton vs mixed cohort from samplesheet
    - Auto-generates PED file for trios/duos
    - GC-content BED annotation (once)
    - Per-sample coverage calculation (scatter)
    - Coverage merge (gather)
    - ClinCNV cohort-level CNV calling (trio-aware if mixed cohort)
    - TSV → VCF conversion (scatter)
    - bcftools sort
    - Optional: VEP annotation
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

import groovy.transform.Field

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================================================
     CNV Calling Pipeline v${workflow.manifest.version}
    =========================================================================

    Usage:
        nextflow run main.nf -profile slurm --input samplesheet.csv [options]

    Required:
        --input             Sample CSV (see formats below)
        --fasta             Reference genome FASTA
        --fasta_fai         Reference genome FASTA index
        --cnv_bed           CNV target BED with gene names in col 4
                            Do NOT use SNV padded BED

    Optional:
        --outdir            Output directory (default: ./results)
        --previous_cov_dir  Path to previous batch .cov files (incremental mode)

    VEP Annotation (--run_vep true):
        --vep_cache         VEP cache directory (cache version auto-detected)
        --sv_overlap_vcf    gnomAD-SV v4.1 DEL/DUP VCF built by prepare_gnomad_sv.sh
                            (gnomad.v4.1.sv.sites.DELDUP.noME.vcf.gz, .tbi alongside)
        --vep_sv_plugin     Patched StructuralVariantOverlap.pm
                            (default: assets/vep_plugins/StructuralVariantOverlap.pm)
        --sv_overlap_pct    Min % of CNV covered by one gnomAD SV (default: 80)

    ClinCNV tuning:
        --clincnv_score_threshold   Minimum score (default: 60)
        --clincnv_min_cluster       Minimum cluster size (default: 10)

    Samplesheet formats:

      Singleton cohort:
        sample_id,bam,bai,sex
        LH0001,/path/to/LH0001.bam,/path/to/LH0001.bai,2

      Mixed cohort (trios/duos/singletons):
        family_id,sample_id,paternal_id,maternal_id,sex,phenotype,bam,bai
        FAM001,FAM001_PROBAND,FAM001_FATHER,FAM001_MOTHER,1,2,/path/bam,/path/bai
        FAM001,FAM001_FATHER,0,0,1,1,/path/bam,/path/bai
        FAM001,FAM001_MOTHER,0,0,2,1,/path/bam,/path/bai

        sex:       1=male, 2=female
        phenotype: 0=unknown, 1=unaffected, 2=affected
        paternal_id/maternal_id: 0 = parent not in cohort

    Outputs:
        results/1_clincnv_output/   Always produced
        results/2_vep_annotated/    If --run_vep true
    =========================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    INCLUDE SUBWORKFLOWS
========================================================================================
*/

include { INPUT_VALIDATION       } from './subworkflows/local/input_validation/main'
include { COVERAGE_WORKFLOW      } from './subworkflows/local/coverage_workflow/coverage_workflow'
include { CNV_CALLING_WORKFLOW   } from './subworkflows/local/cnv_calling_workflow/cnv_calling_workflow'
include { ANNOTATION_WORKFLOW    } from './subworkflows/local/annotation_workflow/annotation_workflow'

/*
========================================================================================
    VALIDATE PARAMETERS
========================================================================================
*/

def validateParameters() {
    def errors = []

    if (!params.input)
        errors << "ERROR: --input samplesheet required"
    else if (!file(params.input).exists())
        errors << "ERROR: Samplesheet not found: ${params.input}"

    if (!params.fasta)
        errors << "ERROR: --fasta required"

    if (!params.fasta_fai)
        errors << "ERROR: --fasta_fai required"

    if (!params.cnv_bed)
        errors << "ERROR: --cnv_bed required (CNV target BED with comma-separated gene names in col 4)"
    else if (!file(params.cnv_bed).exists())
        errors << "ERROR: cnv_bed not found: ${params.cnv_bed}"

    // vcf_header has a config default (assets/header.vcf) but assets/ only gets
    // NO_FILE created at runtime — if header.vcf is not committed, the failure
    // surfaces inside MAKE_VCF instead of here.
    if (!params.vcf_header)
        errors << "ERROR: vcf_header not set"
    else if (!file(params.vcf_header).exists())
        errors << "ERROR: vcf_header not found: ${params.vcf_header}"

    if (params.run_vep) {
        if (!params.vep_cache)      errors << "ERROR: --vep_cache required when --run_vep true"
        if (!params.sv_overlap_vcf) errors << "ERROR: --sv_overlap_vcf required when --run_vep true"
        if (!file(params.vep_sv_plugin).exists())
            errors << "ERROR: patched VEP plugin not found: ${params.vep_sv_plugin}"
        else if (!file(params.vep_sv_plugin).text.contains('PATCH: avoid chrchr'))
            errors << "ERROR: ${params.vep_sv_plugin} is not the patched StructuralVariantOverlap.pm"
        if (params.sv_overlap_vcf && !file("${params.sv_overlap_vcf}.tbi").exists())
            errors << "ERROR: index not found: ${params.sv_overlap_vcf}.tbi"
    }

    if (errors) {
        log.error "\nPARAMETER VALIDATION FAILED\n" + errors.join("\n")
        exit 1
    }
    log.info "Parameters validated successfully"
}

validateParameters()

/*
========================================================================================
    PRINT STATIC BANNER (before cohort type is known)
========================================================================================
*/

log.info """
╔═══════════════════════════════════════════════════════════════════╗
║       ClinCNV WES CNV Calling Pipeline v${workflow.manifest.version}
╠═══════════════════════════════════════════════════════════════════╣
║  Input samplesheet  : ${params.input}
║  Reference FASTA    : ${params.fasta ?: 'NOT SET'}
║  CNV target BED     : ${params.cnv_bed ?: 'NOT SET'}
║  ClinCNV score >=   : ${params.clincnv_score_threshold}
║  ClinCNV min cluster: ${params.clincnv_min_cluster}
║  Run VEP            : ${params.run_vep}
║  gnomAD-SV VCF      : ${params.run_vep ? params.sv_overlap_vcf : 'N/A'}
║  Output directory   : ${params.outdir}
║  Cohort type        : auto-detecting from samplesheet...
╚═══════════════════════════════════════════════════════════════════╝
""".stripIndent()

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    ch_versions = Channel.empty()

    // Ensure assets directory exists
    file("${projectDir}/assets").mkdirs()
    if (!file("${projectDir}/assets/NO_FILE").exists()) {
        file("${projectDir}/assets/NO_FILE").text = ""
    }

    /*
    ============================================================================
        STAGE 1: INPUT VALIDATION
        Auto-detects singleton vs mixed cohort
        Auto-generates PED file for mixed cohort
    ============================================================================
    */

    INPUT_VALIDATION( file(params.input) )

    ch_bam         = INPUT_VALIDATION.out.bam
    ch_ped         = INPUT_VALIDATION.out.ped
    ch_cohort_type = INPUT_VALIDATION.out.cohort_type
    ch_versions    = ch_versions.mix(INPUT_VALIDATION.out.versions)

    // Log cohort type once it is resolved from the process output
    ch_cohort_type.subscribe { type ->
        log.info "Cohort type detected: ${type.toUpperCase()}"
    }

    /*
    ============================================================================
        STAGE 2: COVERAGE (GC annotate → per-sample BedCoverage → merge)
    ============================================================================
    */

    ch_fasta      = Channel.value(file(params.fasta,     checkIfExists: true))
    ch_fasta_fai  = Channel.value(file(params.fasta_fai, checkIfExists: true))
    ch_vcf_header = Channel.value(file(params.vcf_header))
    ch_cnv_bed    = Channel.value(file(params.cnv_bed,   checkIfExists: true))

    COVERAGE_WORKFLOW(
        ch_bam,
        ch_cnv_bed,
        ch_fasta,
        ch_fasta_fai
    )
    ch_annotated_bed = COVERAGE_WORKFLOW.out.annotated_bed
    ch_merged_cov    = COVERAGE_WORKFLOW.out.merged_cov
    ch_versions      = ch_versions.mix(COVERAGE_WORKFLOW.out.versions)

    /*
    ============================================================================
        STAGE 3: CNV CALLING (ClinCNV → TSV→VCF → sort)
        PED auto-passed from INPUT_VALIDATION (empty file for singletons)
    ============================================================================
    */

    ch_sample_meta = COVERAGE_WORKFLOW.out.bam_with_sex.map { meta, bam, bai -> meta }

    CNV_CALLING_WORKFLOW(
        ch_merged_cov,
        ch_annotated_bed,
        ch_sample_meta,
        ch_vcf_header,
        ch_ped,
        ch_fasta_fai
    )

    ch_sorted_vcf = CNV_CALLING_WORKFLOW.out.vcf
    ch_versions   = ch_versions.mix(CNV_CALLING_WORKFLOW.out.versions)

    /*
    ============================================================================
        STAGE 4: VEP ANNOTATION (optional)
    ============================================================================
    */

    ch_vep_cache      = params.run_vep ? Channel.value(file(params.vep_cache,       checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_vep_sv_plugin  = params.run_vep ? Channel.value(file(params.vep_sv_plugin,   checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_sv_overlap     = params.run_vep ? Channel.value(file(params.sv_overlap_vcf,  checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_sv_overlap_tbi = params.run_vep ? Channel.value(file("${params.sv_overlap_vcf}.tbi", checkIfExists: true)) : Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_gene_panel     = Channel.value(file("${projectDir}/assets/NO_FILE"))
    ch_rscript        = Channel.value(file("${projectDir}/assets/NO_FILE"))

    ANNOTATION_WORKFLOW(
        ch_sorted_vcf,
        ch_vep_cache,
        ch_vep_sv_plugin,
        ch_sv_overlap,
        ch_sv_overlap_tbi,
        ch_gene_panel,
        ch_rscript
    )

    ch_versions = ch_versions.mix(ANNOTATION_WORKFLOW.out.versions)

    // Cohort table of gnomAD-SV overlap rates (sample, total, with_overlap, fraction)
    ANNOTATION_WORKFLOW.out.sv_qc
        .map { meta, f -> f.text }
        .collectFile(name: 'sv_overlap_qc.tsv', storeDir: "${params.outdir}/pipeline_info",
                     seed: "sample_id\ttotal_cnvs\tcnvs_with_overlap\tfraction\n", sort: true)

    /*
    ============================================================================
        STAGE 5: SOFTWARE VERSIONS
    ============================================================================
    */

    ch_versions
        .flatten()
        .filter  { it != null }
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}

/*
========================================================================================
    COMPLETION HANDLERS
========================================================================================
*/

workflow.onComplete {
    log.info """
╔═══════════════════════════════════════════════════════════════════╗
║                    PIPELINE COMPLETED                             ║
╠═══════════════════════════════════════════════════════════════════╣
║  Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
║  Duration : ${workflow.duration}
║  Output   : ${params.outdir}
╚═══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

workflow.onError {
    log.error """
╔═══════════════════════════════════════════════════════════════════╗
║                     PIPELINE FAILED                               ║
║  Error: ${workflow.errorMessage}
╚═══════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
