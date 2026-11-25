include { GC_ANNOTATE_BED  } from '../../../modules/local/gc_annotate/gc_annotate_bed'
include { INFER_SEX        } from '../../../modules/local/infer_sex/infer_sex'
include { BED_COVERAGE     } from '../../../modules/local/bed_coverage/bed_coverage'
include { MERGE_COVERAGE   } from '../../../modules/local/merge_coverage/merge_coverage'
include { ARCHIVE_COVERAGE } from '../../../modules/local/archive_coverage/archive_coverage'

workflow COVERAGE_WORKFLOW {

    take:
    ch_bam         // channel: [ meta, bam, bai ]
    ch_bed         // channel: value(cnv_bed) — gene names in col 4 (comma-separated)
    ch_fasta       // channel: value(fasta)
    ch_fasta_fai   // channel: value(fasta.fai)

    main:
    ch_versions = Channel.empty()

    // - Incremental mode guard ----------------------
    // A wrong/moved previous_cov_dir must never silently yield zero archived
    // samples: ClinCNV would succeed on a smaller normalisation cohort and
    // emit different calls with nothing in the log to say history was lost.
    def n_previous_cov = 0
    if (params.previous_cov_dir) {
        def prev_dir = file(params.previous_cov_dir)
        if (!prev_dir.exists()) {
            error "--previous_cov_dir does not exist: ${params.previous_cov_dir}"
        }
        n_previous_cov = prev_dir.list().findAll { it.endsWith('.cov') }.size()
        if (n_previous_cov == 0) {
            error "no .cov files found in --previous_cov_dir: ${params.previous_cov_dir}"
        }
        log.info "INFO: incremental mode — ${n_previous_cov} archived .cov file(s) from ${params.previous_cov_dir}"
    } else {
        log.info "INFO: current batch only — no --previous_cov_dir set"
    }

    // Step 1: Annotate BED with GC content — runs once
    GC_ANNOTATE_BED(
        ch_bed,
        ch_fasta,
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(GC_ANNOTATE_BED.out.versions)

    // Step 2: Infer sex only for samples where sex is unknown (0)
    //   - Known sex (1/2): keep samplesheet value, skip SampleGender
    //   - Unknown sex (0): run SampleGender and replace with inferred value
    ch_known_sex   = ch_bam.filter { meta, bam, bai -> meta.sex != null && meta.sex != 0 }
    ch_unknown_sex = ch_bam.filter { meta, bam, bai -> meta.sex == null || meta.sex == 0  }

    INFER_SEX(ch_unknown_sex)
    ch_versions = ch_versions.mix(INFER_SEX.out.versions)

    // Parse SampleGender output and update meta.sex for unknown-sex samples
    ch_inferred = ch_unknown_sex
        .join(INFER_SEX.out.sex, by: 0)
        .map { meta, bam, bai, sex_txt ->
            def inferred_sex = 0
            sex_txt.eachLine { line ->
                if (!line.startsWith('#')) {
                    def cols = line.split('\t')
                    if (cols.size() >= 3) {
                        inferred_sex = cols[2].trim() == 'male'   ? 1 :
                                       cols[2].trim() == 'female' ? 2 : 0
                    }
                }
            }
            if (inferred_sex == 0) {
                log.warn "Sex could not be inferred for '${meta.id}' — sex-chromosome CN calls will be skipped for this sample"
            } else {
                log.info "INFO: inferred sex for '${meta.id}': ${inferred_sex == 1 ? 'male' : 'female'}"
            }
            [ meta + [sex: inferred_sex], bam, bai ]
        }

    // Merge known-sex and inferred-sex channels
    ch_sex_updated = ch_known_sex.mix(ch_inferred)

    // Step 3: Per-sample coverage (scatter) — current batch only
    BED_COVERAGE(
        ch_sex_updated,
        GC_ANNOTATE_BED.out.annotated_bed
    )
    ch_versions = ch_versions.mix(BED_COVERAGE.out.versions)

    // Step 4: Build full cohort coverage channel
    // ARCHIVE_COVERAGE writes this batch into the same directory that
    // previous_cov_dir points at, so on any -resume or re-run of an already
    // archived batch the glob below picks this batch up a second time.
    // Dedupe by basename (BED_COVERAGE emits "${meta.id}.cov"), current wins:
    // the current batch's freshly computed coverage supersedes any archived copy.
    ch_current_cov  = BED_COVERAGE.out.cov
    ch_previous_cov = params.previous_cov_dir
        ? Channel.fromPath("${params.previous_cov_dir}/*.cov", checkIfExists: true)
        : Channel.empty()

    // NB: do NOT use .collect().combine(.collect()) here — combine flattens both
    // lists into one tuple, so the map would bind a single Path instead of a list.
    // Tag each file with its origin and gather once with toList(), which does not flatten.
    ch_merge_input = ch_current_cov.map { f -> [ 'current', f ] }
        .mix( ch_previous_cov.map { f -> [ 'archive', f ] } )
        .toList()
        .map { items ->
            def current  = items.findAll { it[0] == 'current' }.collect { it[1] }
            def names    = current.collect { it.name } as Set
            def archived = items.findAll { it[0] == 'archive' }.collect { it[1] }
            def kept     = archived.findAll { !names.contains(it.name) }
            def dropped  = archived.size() - kept.size()
            if (dropped > 0) {
                log.info "INFO: ${dropped} archived .cov file(s) superseded by the current batch"
            }
            log.info "INFO: cohort for ClinCNV = ${current.size()} current + ${kept.size()} archived = ${current.size() + kept.size()} samples"
            current + kept
        }

    // Step 5: Archive current batch .cov files to persistent store
    ARCHIVE_COVERAGE(
        ch_current_cov.collect()
    )

    // Step 6: Merge full cohort (deduped: current batch + non-overlapping archive)
    MERGE_COVERAGE(
        ch_merge_input
    )
    ch_versions = ch_versions.mix(MERGE_COVERAGE.out.versions)

    emit:
    annotated_bed  = GC_ANNOTATE_BED.out.annotated_bed
    merged_cov     = MERGE_COVERAGE.out.merged_cov
    bam_with_sex   = ch_sex_updated   // updated meta with resolved sex
    versions       = ch_versions
}
