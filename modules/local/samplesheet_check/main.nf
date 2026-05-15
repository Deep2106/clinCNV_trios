process SAMPLESHEET_CHECK {
    tag "samplesheet_check"
    label 'process_single'

    // Pure bash  use the smallest available container
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
import os

SINGLETON_COLS = {'sample_id', 'bam', 'bai', 'sex'}
MIXED_COLS     = {'family_id', 'sample_id', 'paternal_id', 'maternal_id',
                  'sex', 'phenotype', 'bam', 'bai'}

rows = []
with open("${samplesheet}") as fh:
    reader = csv.DictReader(fh)
    headers = set(reader.fieldnames or [])
    for row in reader:
        rows.append(row)

if not rows:
    sys.exit("ERROR: Samplesheet is empty")

# ── Detect cohort type ────────────────────────────────────────────────────
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

# ── Validate rows ─────────────────────────────────────────────────────────
sample_ids = set()
errors     = []
warnings   = []

for i, row in enumerate(rows, start=2):
    sid = row.get('sample_id', '').strip()
    if not sid:
        errors.append(f"Row {i}: sample_id is empty")
        continue
    if sid in sample_ids:
        errors.append(f"Row {i}: duplicate sample_id '{sid}'")
    sample_ids.add(sid)

    bam = row.get('bam', '').strip()
    bai = row.get('bai', '').strip()
    if not bam: errors.append(f"Row {i} ({sid}): bam path is empty")
    if not bai: errors.append(f"Row {i} ({sid}): bai path is empty")

    sex = row.get('sex', '').strip()
    if sex not in ('1', '2'):
        errors.append(f"Row {i} ({sid}): sex must be 1 (male) or 2 (female), got '{sex}'")

    if is_mixed:
        pat = row.get('paternal_id', '0').strip()
        mat = row.get('maternal_id', '0').strip()
        phen = row.get('phenotype', '').strip()
        if phen not in ('0', '1', '2'):
            errors.append(f"Row {i} ({sid}): phenotype must be 0 (unknown), 1 (unaffected) or 2 (affected), got '{phen}'")

if errors:
    print("\\nERRORS found in samplesheet:")
    for e in errors: print(f"  {e}")
    sys.exit(1)

# ── Duo detection (mixed cohort only) ────────────────────────────────────
if is_mixed:
    for row in rows:
        sid = row['sample_id'].strip()
        pat = row.get('paternal_id', '0').strip()
        mat = row.get('maternal_id', '0').strip()
        has_pat = pat != '0' and pat in sample_ids
        has_mat = mat != '0' and mat in sample_ids
        if (pat != '0' or mat != '0') and not (has_pat and has_mat):
            warnings.append(
                f"WARNING: '{sid}' is a DUO (incomplete family)  "
                f"excluded from de novo calling, will run as singleton"
            )

for w in warnings:
    print(w)

# ── Write validated CSV ───────────────────────────────────────────────────
out_fields = list(rows[0].keys())
with open("validated.csv", "w", newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=out_fields)
    writer.writeheader()
    writer.writerows(rows)

# ── Generate PED file (mixed only) ───────────────────────────────────────
# PED format: family_id  sample_id  paternal_id  maternal_id  sex  phenotype
# Duos excluded from trios but still written to PED for completeness
if is_mixed:
    with open("cohort.ped", "w") as fh:
        for row in rows:
            fam  = row.get('family_id',   row['sample_id']).strip()
            sid  = row['sample_id'].strip()
            pat  = row.get('paternal_id', '0').strip()
            mat  = row.get('maternal_id', '0').strip()
            sex  = row['sex'].strip()
            phen_raw = row.get('phenotype', '1').strip()
            phen = '1' if phen_raw == '0' else phen_raw  # 0=unknown → treat as unaffected in PED
            # Write 0 for parents not present in samplesheet
            pat_out = pat if pat in sample_ids else '0'
            mat_out = mat if mat in sample_ids else '0'
            fh.write(f"{fam}\\t{sid}\\t{pat_out}\\t{mat_out}\\t{sex}\\t{phen}\\n")
    print(f"INFO: PED file generated with {len(rows)} samples")
    # Summary
    complete_trios = sum(
        1 for r in rows
        if r.get('paternal_id','0').strip() in sample_ids
        and r.get('maternal_id','0').strip() in sample_ids
    )
    print(f"INFO: Complete trios (eligible for de novo): {complete_trios}")
    print(f"INFO: Duos/singletons (no de novo): {len(rows) - complete_trios}")
else:
    # Singleton  emit empty PED (NO_FILE equivalent)
    open("cohort.ped", "w").close()
    print("INFO: No PED file generated (singleton cohort)")

# ── Write cohort type ─────────────────────────────────────────────────────
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
