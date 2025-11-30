#!/usr/bin/env bash

# ============================================================
# Script: gwas_qc.sh
# Purpose:
#   Extract a single-sample GWAS-VCF, apply QC filters, and
#   standardize variant IDs into CHR:BP:REF:ALT format.
#
# Output:
#   <trait>.QC.gz  (QC'd GWAS summary statistics)
#
# Usage:
#   bash gwas_qc.sh <path/to/file.vcf.gz> <trait_name>
# ============================================================

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <path/to/file.vcf.gz> <trait_name>"
  exit 1
fi

VCF="$1"
TRAIT="$2"
OUT_RAW="${TRAIT}.raw.tsv"
OUT_QC="${TRAIT}.QC.gz"

# ============================================================
# Step 1 — Extract fields from GWAS-VCF
# ============================================================

SAMPLE=$(bcftools query -l "$VCF" | head -1)
echo "Extracting using sample: $SAMPLE"

HEADER="CHR\tBP\trsID\tREF\tALT\tES\tSE\tLP\tP\tAF\tSS\tEZ"

{
  echo -e "$HEADER"
  bcftools query -s "$SAMPLE" \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%ES]\t[%SE]\t[%LP]\t[%AF]\t[%SS]\t[%EZ]\n' \
    "$VCF"
} > "$OUT_RAW"

echo "Extracted fields -> $OUT_RAW"


# ============================================================
# Step 2 — QC + Add standardized SNP ID (CHR:BP:REF:ALT)
# ============================================================

echo "[INFO] Applying QC filters and generating normalized SNP IDs..."

awk 'BEGIN { FS=OFS="\t" }
NR==1 {
  print "CHR","BP","SNP","REF","ALT","ES","SE","LP","P","AF","SS","EZ","rsID";
  next;
}

{
  len_ref = length($4);
  len_alt = length($5);

  # Autosomal single-base SNPs only
  if($1 < 1 || $1 > 22 || len_ref != 1 || len_alt != 1)
    next;

  # Remove A/T and C/G ambiguous SNPs
  if( ($4=="A" && $5=="T") || ($4=="T" && $5=="A") ||
      ($4=="C" && $5=="G") || ($4=="G" && $5=="C") )
    next;

  # Convert LP to P (p = 10^-LP)
  pval = ($8=="." || $8=="" ? "." : 10^(-$8));

  # Normalized variant ID
  norm_snp = $1 ":" $2 ":" $4 ":" $5;

  # Deduplicate
  if(!seen[norm_snp]++)
    print $1, $2, norm_snp, $4, $5, $6, $7, $8, pval, $9, $10, $11, $3;
}' "$OUT_RAW" | gzip > "$OUT_QC"

echo "Final QC output -> $OUT_QC"

