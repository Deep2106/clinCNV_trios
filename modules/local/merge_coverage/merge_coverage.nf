process MERGE_COVERAGE {
    tag "merge_coverage"
    label 'process_low'

    container params.containers.clincnv

    input:
    path cov_files   // collected list from all samples

    output:
    path "merged.cov", emit: merged_cov
    path "versions.yml", emit: versions

    script:
    """
    # - Region-set integrity check --------------------
    # Archived .cov files from earlier batches were computed against whatever
    # BED was current at the time. If --cnv_bed has changed since (new capture
    # kit, re-annotation), the row sets no longer align and mergeFilesFromFolder
    # produces a misaligned matrix that ClinCNV will happily consume.
    # Every .cov entering the merge must describe the same regions.
    REF_SUM=""
    REF_FILE=""
    N_FILES=0
    for f in *.cov; do
        [ -e "\$f" ] || continue
        N_FILES=\$(( N_FILES + 1 ))
        SUM=\$(tail -n +2 "\$f" | cut -f1-3 | md5sum | cut -d' ' -f1)
        if [ -z "\$REF_SUM" ]; then
            REF_SUM="\$SUM"
            REF_FILE="\$f"
            REF_N=\$(tail -n +2 "\$f" | wc -l)
        elif [ "\$SUM" != "\$REF_SUM" ]; then
            THIS_N=\$(tail -n +2 "\$f" | wc -l)
            echo "ERROR: region set mismatch between coverage files." >&2
            echo "  \${REF_FILE}: \${REF_N} regions" >&2
            echo "  \${f}: \${THIS_N} regions" >&2
            echo "  These were computed against different BED files." >&2
            echo "  Either rerun the whole cohort with one --cnv_bed, or clear" >&2
            echo "  the coverage archive and recompute." >&2
            exit 1
        fi
    done
    echo "INFO: region-set check passed — \${N_FILES} coverage file(s), \${REF_N} regions each"

    Rscript \${EBROOTCLINCNV}/mergeFilesFromFolder.R \\
        -i . \\
        -o merged.cov

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clincnv: \$(Rscript -e "packageVersion('ClinCNV')" 2>/dev/null || echo "1.18.3")
    END_VERSIONS
    """
}
