#!/bin/bash
#SBATCH --job-name=CNDcnv_trios
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --partition=ProdQ
#SBATCH --output=nextflow_%j.out
#SBATCH --error=nextflow_%j.err

# ============================================================================
# ClinCNV WES CNV Calling Pipeline v2.0  Mixed Cohort
# ============================================================================
# Cohort type is AUTO-DETECTED from samplesheet column headers:
#   Singleton : sample_id,bam,bai,sex
#   Mixed     : family_id,sample_id,paternal_id,maternal_id,sex,phenotype,bam,bai
# PED file is AUTO-GENERATED  no manual PED preparation needed.
#
# Usage:
#   FIRST RUN:  sbatch run_pipeline.slurm
#   NEXT BATCH: Edit BATCH_NAME, INPUT_CSV, PREVIOUS_COV_DIR
# ============================================================================

# ============================================================================
# USER CONFIGURATION  EDIT THESE FOR EACH RUN
# ============================================================================

BATCH_NAME="CNV_Legacy_batch1"          # log label only  not passed to Nextflow
INPUT_CSV="sample_trios.csv"

RUN_VEP="true"

# ClinCNV tuning  set explicitly so each batch records the thresholds used
CLINCNV_SCORE="50"
CLINCNV_MIN_CLUSTER="10"

# Incremental coverage store
# First run:       leave empty  batch is archived to ${OUTDIR}/coverage_archive
# Subsequent runs: point at that archive  NOTE it lives under OUTDIR, not PROJECT_DIR:
#                  ${PROJECT_DIR}/results_trios/coverage_archive
# A wrong path is now a hard error, not a silent empty cohort.
PREVIOUS_COV_DIR=""

# ============================================================================
# PATHS  do not edit unless directory structure changes
# ============================================================================

PROJECT_DIR=/path/to
PIPELINE_DIR=${PROJECT_DIR}/clinCNV_trios

export NXF_WORK=${PROJECT_DIR}/work_trios
mkdir -p "${NXF_WORK}"

OUTDIR=${PROJECT_DIR}/CNV_results/results_trios_legacy

# ============================================================================
# ENVIRONMENT
# ============================================================================

module purge
module --ignore-cache load Java/17

export NXF_SINGULARITY_LIBRARYDIR=/path/to/Docker_Images
export NXF_SINGULARITY_CACHEDIR=${PROJECT_DIR}/.singularity_cache
mkdir -p "${NXF_SINGULARITY_CACHEDIR}"

export NXF_ANSI_LOG=false
export _JAVA_OPTIONS="-Djava.io.tmpdir=${TMPDIR:-/tmp}"

# ============================================================================
# PRE-RUN CHECKS
# ============================================================================

cd "${PIPELINE_DIR}" || { echo "ERROR: Cannot cd to ${PIPELINE_DIR}"; exit 1; }

if [ ! -f "main.nf" ]; then
    echo "ERROR: main.nf not found in ${PIPELINE_DIR}"
    exit 1
fi

if [ ! -f "${PROJECT_DIR}/${INPUT_CSV}" ]; then
    echo "ERROR: Samplesheet not found: ${PROJECT_DIR}/${INPUT_CSV}"
    exit 1
fi

echo "=============================================================================="
echo "ClinCNV WES CNV Calling Pipeline v2.0  trios Cohort"
echo "=============================================================================="
echo "Date:              $(date)"
echo "Node:              $(hostname)"
echo "SLURM Job ID:      ${SLURM_JOB_ID}"
echo "Batch:             ${BATCH_NAME}"
echo "Input CSV:         ${INPUT_CSV}"
echo "Previous cov dir:  ${PREVIOUS_COV_DIR:-NOT SET (first run)}"
echo "Run VEP:           ${RUN_VEP}"
echo "ClinCNV score:     ${CLINCNV_SCORE}"
echo "ClinCNV cluster:   ${CLINCNV_MIN_CLUSTER}"
echo "Output:            ${OUTDIR}"
echo "NXF_WORK:          ${NXF_WORK}"
echo "Submit dir:        ${SLURM_SUBMIT_DIR:-$(pwd)}   (SLURM .out/.err land here)"
echo "NOTE: Cohort type (singleton+duos/trios) auto-detected from samplesheet headers"
echo "=============================================================================="

# ============================================================================
# BUILD OPTIONAL ARGUMENTS
# ============================================================================

OPTIONAL_ARGS=()

if [ -n "${PREVIOUS_COV_DIR}" ]; then
    OPTIONAL_ARGS+=(--previous_cov_dir "${PREVIOUS_COV_DIR}")
fi

# ============================================================================
# RUN PIPELINE
# ============================================================================

nextflow run main.nf \
    -profile slurm \
    -w "${NXF_WORK}" \
    --input                     "${PROJECT_DIR}/${INPUT_CSV}" \
    --fasta                     /path/to/reference/hg38.fasta \
    --fasta_fai                 /path/to/reference/hg38.fasta.fai \
    --cnv_bed        /path/to/reference/EXOME_Xtra_EPIFUN/clincnv_beds/legacy_consensus.min20.annotated.bed \
    --vep_cache                 /path/to/resources/vep_cache \
    --vep_plugins               /path/to/resources/VEP_plugins \
    --sv_overlap_vcf /path/to/resources/sv_cache/gnomad.v4.1.sv.sites.DELDUP.noME.vcf.gz \
    --outdir                    "${OUTDIR}" \
    --run_vep                   "${RUN_VEP}" \
    --clincnv_score_threshold   "${CLINCNV_SCORE}" \
    --clincnv_min_cluster       "${CLINCNV_MIN_CLUSTER}" \
    "${OPTIONAL_ARGS[@]}" \
    -resume

EXIT_CODE=$?

echo "=============================================================================="
echo "Pipeline exit code: ${EXIT_CODE}"
echo "Date: $(date)"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "SUCCESS"
    echo "  1_clincnv_output:  ${OUTDIR}/1_clincnv_output/"
    echo "  2_vep_annotated:   ${OUTDIR}/2_vep_annotated/  (if run_vep=true)"
else
    echo "FAILED  check:"
    echo "  ${SLURM_SUBMIT_DIR:-.}/nextflow_${SLURM_JOB_ID}.err"
    echo "  ${PIPELINE_DIR}/.nextflow.log"
fi
echo "=============================================================================="

exit ${EXIT_CODE}
