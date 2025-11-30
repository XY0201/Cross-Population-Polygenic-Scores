#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <trait>"
    exit 1
fi

TRAIT="$1"

# Paths (relative, GitHub-safe)
BASE="data"
TARGET="data/1000G_QC/target_1000G_finalQC"

GWAS_GZ="${BASE}/qc_output/${TRAIT}.QC.gz"
GWAS_TRANS="${BASE}/qc_output/${TRAIT}.QC.Transformed"

A1="${BASE}/mismatch/${TRAIT}.a1"
MIS="${BASE}/mismatch/${TRAIT}.mismatch"

OUTDIR="${BASE}/prePGS"
mkdir -p "$OUTDIR"


# 1. Transform: uncompress QC.gz → QC.Transformed
gunzip -c "$GWAS_GZ" > "$GWAS_TRANS"


# 2. Create trait-specific target dataset
plink \
  --bfile "$TARGET" \
  --exclude "$MIS" \
  --a1-allele "$A1" 1 2 \
  --make-bed \
  --out "$OUTDIR/${TRAIT}_target_QC"


# 3. Clumping
plink \
  --bfile "$OUTDIR/${TRAIT}_target_QC" \
  --clump "$GWAS_TRANS" \
  --clump-snp-field SNP \
  --clump-field P \
  --clump-p1 1e-4 \
  --clump-p2 1e-4 \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --out "$OUTDIR/${TRAIT}_clump"


# 4. Extract clumped SNP list
awk 'NR>1 {print $3}' \
  "$OUTDIR/${TRAIT}_clump.clumped" \
  > "$OUTDIR/${TRAIT}.valid.snp"

