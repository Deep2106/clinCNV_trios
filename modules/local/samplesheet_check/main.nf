process SAMPLESHEET_CHECK {
    tag "samplesheet_check"
    label 'process_single'

    // Pure bash — use the smallest available container
    container params.containers.python

    input:
    path samplesheet

    output:
    path "validated.csv",  emit: csv
    path "cohort.ped",     emit: ped
    path "cohort_type.txt",emit: cohort_type
    path "versions.yml",   emit: versions

    script:
    """
    python3 << 'EOF'
import csv
import sys

SINGLETON_COLS = {'sample_id', 'bam', 'bai', 'sex'}
MIXED_COLS     = {'family_id', 'sample_id', 'paternal_id', 'maternal_id',
                  'sex', 'phenotype', 'bam', 'bai'}

rows = []
with open("${samplesheet}") as fh:
    # restval='' so short rows yield '' rather than None (None.strip() crashes)
    reader = csv.DictReader(fh, restval='')
    headers = set(reader.fieldnames or [])
    for row in reader:
        rows.append(row)

if not rows:
    sys.exit("ERROR: Samplesheet is empty")

# - Detect cohort type --------------------------
is_mixed = MIXED_COLS.issubset(headers)
is_singleton = SINGLETON_COLS.issubset(headers) and not is_mixed

if not is_mixed and not is_singleton:
    sys.exit(
        "ERROR: Unrecognised samplesheet format.\\n"
        "  Singleton columns required : " + str(sorted(SINGLETON_COLS)) + "\\n"
        "  Mixed columns required     : " + str(sorted(MIXED_COLS))
    )

cohort_type = "mixed" if is_mixed else "singleton"
print(f"INFO: Cohort type detected: {cohort_type.upper()}")

# - Validate rows ----------------------------─
sample_ids = set()
errors     = []
warnings   = []
infer_sex  = []

def cell(row, key, default=''):
    v = row.get(key, default)
    return (v if v is not None else default).strip()

for i, row in enumerate(rows, start=2):
    sid = cell(row, 'sample_id')
    if not sid:
        errors.append(f"Row {i}: sample_id is empty")
        continue
    if sid in sample_ids:
        errors.append(f"Row {i}: duplicate sample_id '{sid}'")
    sample_ids.add(sid)

    if not cell(row, 'bam'): errors.append(f"Row {i} ({sid}): bam path is empty")
    if not cell(row, 'bai'): errors.append(f"Row {i} ({sid}): bai path is empty")

    # sex 0 = unknown; inferred downstream by INFER_SEX from the BAM
    sex = cell(row, 'sex')
    if sex not in ('0', '1', '2'):
        errors.append(
            f"Row {i} ({sid}): sex must be 0 (unknown), 1 (male) or 2 (female), got '{sex}'"
        )
    elif sex == '0':
        infer_sex.append(sid)

    if is_mixed:
        if not cell(row, 'family_id'):
            errors.append(f"Row {i} ({sid}): family_id is empty (required for mixed cohort)")
        phen = cell(row, 'phenotype')
        if phen not in ('0', '1', '2'):
            errors.append(
                f"Row {i} ({sid}): phenotype must be 0 (unknown), 1 (unaffected) "
                f"or 2 (affected), got '{phen}'"
            )

if errors:
    print("\\nERRORS found in samplesheet:")
    for e in errors: print(f"  {e}")
    sys.exit(1)

if infer_sex:
    print(f"INFO: {len(infer_sex)} sample(s) with unknown sex — will be inferred from BAM: "
          + ", ".join(infer_sex))

# - Family structure classification (mixed cohort only) ---------
n_trio = n_duo = n_single = 0

if is_mixed:
    for row in rows:
        sid = cell(row, 'sample_id')
        pat = cell(row, 'paternal_id', '0')
        mat = cell(row, 'maternal_id', '0')
        has_pat = pat != '0' and pat in sample_ids
        has_mat = mat != '0' and mat in sample_ids

        if has_pat and has_mat:
            n_trio += 1
        elif pat != '0' or mat != '0':
            n_duo += 1
            missing = []
            if pat != '0' and not has_pat: missing.append(f"father '{pat}'")
            if mat != '0' and not has_mat: missing.append(f"mother '{mat}'")
            detail = (" — " + " and ".join(missing) + " not in samplesheet") if missing else ""
            warnings.append(
                f"WARNING: '{sid}' is a DUO (incomplete family){detail} — "
                f"excluded from de novo calling, will run as singleton"
            )
        else:
            n_single += 1

for w in warnings:
    print(w)

# - Write validated CSV -------------------------─
out_fields = list(reader.fieldnames)
with open("validated.csv", "w", newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=out_fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)

# - Generate PED file (mixed only) -------------------─
# PED format: family_id  sample_id  paternal_id  maternal_id  sex  phenotype
# Duos are written to the PED but are not eligible for de novo calling.
if is_mixed:
    with open("cohort.ped", "w") as fh:
        for row in rows:
            fam  = cell(row, 'family_id')
            sid  = cell(row, 'sample_id')
            pat  = cell(row, 'paternal_id', '0')
            mat  = cell(row, 'maternal_id', '0')
            sex  = cell(row, 'sex')
            phen_raw = cell(row, 'phenotype', '1')
            phen = '1' if phen_raw == '0' else phen_raw  # 0=unknown -> unaffected in PED
            # Write 0 for parents not present in samplesheet
            pat_out = pat if pat in sample_ids else '0'
            mat_out = mat if mat in sample_ids else '0'
            fh.write(f"{fam}\\t{sid}\\t{pat_out}\\t{mat_out}\\t{sex}\\t{phen}\\n")
    print(f"INFO: PED file generated with {len(rows)} samples")
    print(f"INFO: Probands in complete trios (eligible for de novo): {n_trio}")
    print(f"INFO: Duos (no de novo, run as singleton):               {n_duo}")
    print(f"INFO: Unrelated / parent samples:                        {n_single}")
else:
    # Singleton — emit empty PED (NO_FILE equivalent)
    open("cohort.ped", "w").close()
    print("INFO: No PED file generated (singleton cohort)")

# - Write cohort type --------------------------─
with open("cohort_type.txt", "w") as fh:
    fh.write(cohort_type)

print("INFO: Samplesheet validation complete")
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """
}
