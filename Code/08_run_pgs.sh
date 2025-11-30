#!/usr/bin/env bash
set -euo pipefail

# Directory settings (use relative paths)
BASE="."
DATA="$BASE/data"
JOBS="$DATA/jobs.csv"

OUTDIR="$DATA/PGS"
mkdir -p "$OUTDIR"

# Loop through jobs.csv (skip header)
tail -n +2 "$JOBS" | while IFS=',' read -r CODE TRAIT; do
    TRAIT=$(echo "$TRAIT" | xargs)
    [[ -z "$TRAIT" ]] && continue

    # Input files
    GWAS_FILE="$DATA/qc_output/${TRAIT}.QC.Transformed"
    TARGET_PREFIX="$DATA/prePGS/${TRAIT}_target_QC"
    VALID_SNP="$DATA/prePGS/${TRAIT}.valid.snp"

    OUT_PREFIX="$OUTDIR/${TRAIT}_PGS"

    # Required file checks
    [[ -f "$GWAS_FILE" ]] || continue
    [[ -f "${TARGET_PREFIX}.bim" ]] || continue
    [[ -f "$VALID_SNP" ]] || continue

    # PGS scoring
    plink \
      --bfile "$TARGET_PREFIX" \
      --extract "$VALID_SNP" \
      --score "$GWAS_FILE" 3 5 6 header \
      --out "$OUT_PREFIX"
done

