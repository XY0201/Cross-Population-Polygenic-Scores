#!/usr/bin/env bash

# =====================================================================
# Script: run_full_pipeline.sh
# Description:
#   Full processing pipeline for 1000 Genomes Phase 3:
#     Step 1: Clean per-chromosome VCF (chr1–22)
#     Step 2: Convert VCF to PLINK format for each chromosome
#     Step 3: Merge chromosomes into a genome-wide dataset
#     Step 4: Perform standard QC and create final QC dataset
#
# Requirements:
#   bcftools, plink, tabix, Rscript (for heterozygosity filter)
#
# Directory structure expected:
#   data/
#     ├── 1000G_phase3/     (raw VCF)
#     └── 1000G_QC/         (outputs)
#
# Usage:
#   bash run_full_pipeline.sh
# =====================================================================

set -euo pipefail

RAW_DIR="data/1000G_phase3"
QC_DIR="data/1000G_QC"
mkdir -p "${QC_DIR}"

echo "============================================================"
echo " 1000 Genomes Processing Pipeline"
echo "============================================================"


# ============================================================
# Step 1: Clean chromosomes 1–22
# ============================================================
echo "[STEP 1] Cleaning per-chromosome VCFs..."

for chr in {1..22}; do
    echo "[INFO] chr${chr} ..."

    IN_VCF="${RAW_DIR}/ALL.chr${chr}.vcf.gz"
    CLEAN_VCF="${QC_DIR}/chr${chr}.clean.vcf.gz"
    OUT_PREFIX="${QC_DIR}/chr${chr}"

    # Normalize and filter SNPs
    bcftools norm -m -any -O z "${IN_VCF}" \
    | bcftools view \
        --types snps \
        --min-alleles 2 --max-alleles 2 \
        --exclude 'strlen(REF)!=1 || strlen(ALT)!=1' \
        -O z -o "${CLEAN_VCF}"

    tabix -p vcf "${CLEAN_VCF}"

    # Assign variant IDs + convert to PLINK
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' "${CLEAN_VCF}" -O z \
    | plink \
        --vcf /dev/stdin \
        --double-id \
        --make-bed \
        --out "${OUT_PREFIX}"

done

echo "[STEP 1 DONE] All per-chromosome files cleaned."


# ============================================================
# Step 2: Merge chr1–22 PLINK files
# ============================================================
echo "[STEP 2] Merging chromosomes..."

OUT_PREFIX="${QC_DIR}/target_1000G_raw"

# Base: chr1
cp "${QC_DIR}/chr1.bed" "${OUT_PREFIX}.bed"
cp "${QC_DIR}/chr1.bim" "${OUT_PREFIX}.bim"
cp "${QC_DIR}/chr1.fam" "${OUT_PREFIX}.fam"

# Merge list
MERGE_LIST="${QC_DIR}/merge_list.txt"
rm -f "${MERGE_LIST}"

for chr in {2..22}; do
    echo "${QC_DIR}/chr${chr}" >> "${MERGE_LIST}"
done

plink \
  --bfile "${OUT_PREFIX}" \
  --merge-list "${MERGE_LIST}" \
  --make-bed \
  --out "${OUT_PREFIX}"

echo "[STEP 2 DONE] Genome-wide merged PLINK dataset created."


# ============================================================
# Step 3: Standard QC
# ============================================================
echo "[STEP 3] Running QC..."

RAW="${QC_DIR}/target_1000G_raw"
PREFIX="${QC_DIR}/target_1000G"

# 3.1 SNP & sample QC
plink \
  --bfile "${RAW}" \
  --maf 0.01 \
  --hwe 1e-6 \
  --geno 0.01 \
  --mind 0.01 \
  --write-snplist \
  --make-just-fam \
  --out "${PREFIX}.QC"

# 3.2 LD pruning
plink \
  --bfile "${RAW}" \
  --keep "${PREFIX}.QC.fam" \
  --extract "${PREFIX}.QC.snplist" \
  --indep-pairwise 200 50 0.25 \
  --out "${PREFIX}.QC"

# 3.3 Heterozygosity (F statistic)
plink \
  --bfile "${RAW}" \
  --extract "${PREFIX}.QC.prune.in" \
  --keep "${PREFIX}.QC.fam" \
  --het \
  --out "${PREFIX}.QC"

# 3.4 Run user-defined heterozygosity filter
echo "[INFO] Running R heterozygosity filter: f_filter.R"
Rscript f_filter.R


# 3.5 Find duplicated SNPs
awk 'NR==FNR{keep[$1]=1; next} ($2 in keep){print $2}' \
  "${PREFIX}.QC.snplist" "${RAW}.bim" | sort | uniq -d \
  > "${PREFIX}.dup"

# 3.6 Final QC dataset
plink \
  --bfile "${RAW}" \
  --keep "${PREFIX}.valid.sample" \
  --extract "${PREFIX}.QC.snplist" \
  --exclude "${PREFIX}.dup" \
  --make-bed \
  --out "${PREFIX}_finalQC"

echo "[STEP 3 DONE] Final QC dataset generated:"
echo "  ${PREFIX}_finalQC.bed"
echo "  ${PREFIX}_finalQC.bim"
echo "  ${PREFIX}_finalQC.fam"

echo "============================================================"
echo " Pipeline completed successfully!"
echo "============================================================"

