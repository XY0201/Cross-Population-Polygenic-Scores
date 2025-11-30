#!/usr/bin/env bash
set -euo pipefail

PANEL="data/1000G_phase3/integrated_call_samples_v3.20130502.ALL.4col.panel"
TARGET="data/1000G_QC/target_1000G_finalQC"
OUTDIR="data/pop_AF"

mkdir -p "$OUTDIR"

# Create .keep files
for POP in AFR AMR EAS EUR SAS; do
    awk -v p=$POP 'NR>1 && $3==p {print $1, $1}' "$PANEL" \
        > "$OUTDIR/${POP}.keep"
done

# Compute allele frequencies
for POP in AFR AMR EAS EUR SAS; do
    plink \
      --bfile "$TARGET" \
      --keep "$OUTDIR/${POP}.keep" \
      --freq \
      --out "$OUTDIR/${POP}_freq"
done

