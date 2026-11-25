# ClinCNV WES CNV Calling Pipeline  Trios

A production-grade DSL2 Nextflow pipeline for germline CNV calling from WES data using [ClinCNV](https://github.com/imgag/ClinCNV) in trio-aware mode (`GERMLINE_TRIOS`), with optional VEP annotation.

> **For singletons and duos without family structure, use the [clinCNV_singletons](https://github.com/Deep2106/clinCNV_singletons) repository instead.**

---

## Pipeline Overview

```
samplesheet.csv (family_id, sample_id, paternal_id, maternal_id, sex, phenotype, bam, bai)
    │
    ├── SAMPLESHEET_CHECK   (validate CSV, auto-detect cohort type, auto-generate PED)
    ├── GC_ANNOTATE_BED     (ngs-bits BedAnnotateGC + BedSort  once per run)
    ├── INFER_SEX           (ngs-bits SampleGender -method xy  per sample)
    ├── BED_COVERAGE        (ngs-bits BedCoverage  per sample, scatter)
    ├── ARCHIVE_COVERAGE    (persist .cov files for incremental runs)
    ├── MERGE_COVERAGE      (ClinCNV mergeFilesFromFolder  full cohort gather)
    │
    ├── CLINCNV             (cohort-level CNV calling  GERMLINE_TRIOS)
    │                       (all samples used for normalisation; only trios get output)
    │
    ├── MAKE_VCF            (ClinCNV _CNVs.txt → VCF  proband only, per trio)
    ├── SORT_VCF            (bcftools sort)
    │
    └── [optional] VEP_ANNOTATE    (StructuralVariantOverlap + MANE Select)
```

### Outputs

| Directory | Contents | When |
|---|---|---|
| `results/1_clincnv_output/` | Coverage files, ClinCNV calls (`trios/CHILD-MOTHER-FATHER/`), sorted VCFs | Always |
| `results/2_vep_annotated/` | VEP-annotated VCFs (proband only) | `--run_vep true` |

> **Note:** Only proband (phenotype=2) VCFs are produced. Parents are used for normalisation and de novo scoring only.

> **Note:** Singletons and duos present in the samplesheet will produce a WARNING and must be run through `clinCNV_singletons` separately.

---

## Requirements

### Software
- Nextflow ≥ 25.04.0
- Singularity ≥ 3.8
- Java 17+

### Containers (Singularity `.sif`)

| Tool | Container |
|---|---|
| Python 3.14 | `python_3.14.2.sif` |
| ngs-bits 2025_12 | `ngsbits.sif` |
| ClinCNV (master) + R 4.2 | `clincnv.sif` |
| BCFtools 1.23 | `bcftool_1.23.sif` |
| VEP 115 | `vep.sif` |

---

## Quick Start

### 1. Prepare samplesheet

```csv
family_id,sample_id,paternal_id,maternal_id,sex,phenotype,bam,bai
FAM001,FAM001_PROBAND,FAM001_FATHER,FAM001_MOTHER,1,2,/data/FAM001_PROBAND.bam,/data/FAM001_PROBAND.bai
FAM001,FAM001_FATHER,0,0,1,1,/data/FAM001_FATHER.bam,/data/FAM001_FATHER.bai
FAM001,FAM001_MOTHER,0,0,2,1,/data/FAM001_MOTHER.bam,/data/FAM001_MOTHER.bai
FAM002,FAM002_PROBAND,FAM002_FATHER,FAM002_MOTHER,2,2,/data/FAM002_PROBAND.bam,/data/FAM002_PROBAND.bai
FAM002,FAM002_FATHER,0,0,1,1,/data/FAM002_FATHER.bam,/data/FAM002_FATHER.bai
FAM002,FAM002_MOTHER,0,0,2,0,/data/FAM002_MOTHER.bam,/data/FAM002_MOTHER.bai
```

**Sex codes:** `1`=male, `2`=female
**Phenotype codes:** `0`=unknown, `1`=unaffected (parent), `2`=affected (proband)
**Parent IDs:** use `0` when parent is not in the cohort

The PED file is **auto-generated** from the samplesheet  no manual preparation needed.

### 2. Edit and submit
```bash
# Edit paths in run_pipeline.slurm, then:
sbatch run_pipeline.slurm
```

---

## Samplesheet Format

| Column | Required | Description |
|---|---|---|
| `family_id` | Yes | Family identifier |
| `sample_id` | Yes | Unique sample identifier |
| `paternal_id` | Yes | Father's sample_id, or `0` if absent |
| `maternal_id` | Yes | Mother's sample_id, or `0` if absent |
| `sex` | Yes | 1=male, 2=female |
| `phenotype` | Yes | 0=unknown, 1=unaffected, 2=affected |
| `bam` | Yes | Absolute path to BAM file |
| `bai` | Yes | Absolute path to BAM index |

---

## Cohort Type Auto-Detection

`SAMPLESHEET_CHECK` detects cohort type from column headers automatically:

| Samplesheet columns | Detected mode | Pipeline to use |
|---|---|---|
| `sample_id,bam,bai,sex` | Singleton+duos | `clinCNV_singletons` |
| `family_id,sample_id,paternal_id,...` | trios | `clinCNV_trios` (this pipeline) |

For mixed samplesheets, the pipeline:
- Auto-generates PED from samplesheet
- Builds comma-separated `trios.txt` for ClinCNV from complete trios only
- Warns about singletons/duos and excludes them from trio calling
- Exits with error if no complete trios are found

---

## Parameters

### Required
```bash
--input          # samplesheet CSV
--fasta          # hg38 reference FASTA
--fasta_fai      # hg38 FASTA index
--cnv_bed        # CNV target BED (not SNV padded BED)
--outdir         # output directory
```

### Optional  VEP Annotation
```bash
--run_vep        # true/false (default: false)
--vep_cache      # path to VEP cache (homo_sapiens_vep_115_GRCh38)
--vep_plugins    # path to VEP plugins directory
--sv_overlap_vcf # gnomAD-SV v4.1 VCF (tabix-indexed)
```

### Optional  Incremental Runs
```bash
--previous_cov_dir  # NFS path to archived .cov files from previous batches
```

### ClinCNV Tuning
```bash
--clincnv_score_threshold    # default: 60
--clincnv_min_length         # default: 0 (allow single-exon CNVs)
--clincnv_min_cluster        # default: 10 (raise to 15 at >100 samples)
--clincnv_max_germ_cnvs      # default: 2000
```

---

## CNV Target BED

Use the **CNV-specific** target BED  not the SNV padded BED. Padding dilutes coverage signal and inflates false positive CNV calls.

---

## Trio vs Singleton/duos Output Differences

| Aspect | Singleton | Trio |
|---|---|---|
| ClinCNV mode | `GERMLINE_SINGLE` | `GERMLINE_TRIOS` |
| Output dir | `cnv_output/normal/SAMPLE/` | `cnv_output/trios/CHILD-MOTHER-FATHER/` |
| CNV file | `SAMPLE_cnvs.tsv` | `FAMILY_CNVs.txt` |
| CN column | col4 | col4 (kid_CN_change) |
| Loglikelihood column | col5 | col8 |
| VCF produced for | all samples | proband only |
| De novo detection | no | yes (via family normalisation) |

---

## Incremental Cohort Design

**First run:** leave `PREVIOUS_COV_DIR=""` in `run_pipeline.slurm`  coverage files are archived to `results/coverage_archive/` automatically.

**Subsequent runs:** set `PREVIOUS_COV_DIR="${OUTDIR}/coverage_archive"`  new batch merges with all previous batches for ClinCNV normalisation.

> ClinCNV re-processes the full cohort every batch. Always include complete family units per batch  never split a trio across batches.

---

## ClinCNV Parameters Guide

| Parameter | Default | Rationale |
|---|---|---|
| `--scoreG` | 60 | Clinical grade; lower (50) for discovery |
| `--lengthG` | 0 | Allow single-exon CNVs |
| `--minimumNumOfElemsInCluster` | 10 | Raise to 15 at >100 samples |
| `--maxNumGermCNVs` | 2000 | Permissive for affected cohorts |
| `--noPlot` | on | Prevents write failures in read-only containers |

---

## VEP Annotation

- **Mode:** `--offline --cache`
- **Cache:** homo_sapiens_vep_115_GRCh38
- **Assembly:** GRCh38 (explicit `--assembly GRCh38 --species homo_sapiens`)
- **Plugin:** `StructuralVariantOverlap` with gnomAD-SV v4.1

### Download VEP cache
```bash
singularity exec -B /data vep.sif \
    perl /opt/vep/src/ensembl-vep/INSTALL.pl \
    -a cf -s homo_sapiens -y GRCh38 \
    --CACHE_VERSION 115 \
    -c /data/db/vep_cache
```

### Download gnomAD-SV v4.1
```bash
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
```

---

## Pipeline Structure

```
clinCNV_trios/
├── main.nf
├── nextflow.config
├── run_pipeline.slurm
├── assets/
│   ├── header.vcf
│   └── NO_FILE
├── conf/
│   ├── base.config
│   ├── slurm.config
│   ├── modules.config
│   └── test.config
├── containers/
│   ├── clincnv.def
│   └── ngsbits.def
├── modules/local/
│   ├── samplesheet_check/
│   ├── gc_annotate_bed/
│   ├── bed_coverage/
│   ├── merge_coverage/
│   ├── archive_coverage/
│   ├── clincnv/
│   ├── infer_sex/
│   ├── make_vcf/
│   ├── sort_vcf/
│   └── vep_annotate/
└── subworkflows/local/
    ├── input_validation/
    ├── coverage_workflow/
    ├── cnv_calling_workflow/
    └── annotation_workflow/
```

---

---

## Citation

If you use this pipeline, please cite:

- **ClinCNV:** [Bessonov et al. 2023](https://github.com/imgag/ClinCNV)
- **ngs-bits:** [IMGAG](https://github.com/imgag/ngs-bits)
- **VEP:** [McLaren et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0974-4)
- **gnomAD-SV:** [Collins et al. 2020](https://doi.org/10.1038/s41586-020-2287-8)
- **Nextflow:** [Di Tommaso et al. 2017](https://doi.org/10.1038/nbt.3820)

---


---

## Citing this pipeline

If you use this pipeline in your work, please cite it:


BibTeX:

```bibtex
@software{bharti_clinCNV_trios_pipeline,
  author  = {Bharti, Deepak},
  title   = { ClinCNV CNV Calling Pipeline for trios},
  version = {1.0.0},
  year    = {2026},
  url     = {https://github.com/Deep2106/clinCNV_trios}
}
```

You can also use the **"Cite this repository"** button on the GitHub page. It reads [`CITATION.cff`](CITATION.cff) and gives the citation in APA and BibTeX.

If the pipeline is published in a paper, cite that paper instead; this section will be updated with the reference.




## Author

Deepak Bharti, Clinical Bioinformatician, RCSI Dublin

Contact: [deepakbharti@rcsi.com](mailto:deepakbharti@rcsi.com)

> CNV calling was performed using a DSL2 Nextflow pipeline based on ClinCNV (Bessonov et al. 2023), developed at RCSI (Bharti D, 2026). Trio-aware calling used `--triosFile` for de novo CNV detection. VEP annotation used StructuralVariantOverlap with gnomAD-SV v4.1. Sex was inferred from BAM files using ngs-bits SampleGender.
