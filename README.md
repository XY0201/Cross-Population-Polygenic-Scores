# Cross-Population-Polygenic-Scores
Hi

| Script                       | Description                                                                                               |
| ---------------------------- | --------------------------------------------------------------------------------------------------------- |
| **Quality Control**          |                                                                                                           |
| `01_target_prepare_1000G.sh` | Cleans raw 1000G VCFs, converts to PLINK, merges chr1–22, and performs standard target QC.                |
| `02_gwas_qc.sh`              | Extracts fields from GWAS-VCF, filters SNPs, assigns `CHR:BP:REF:ALT`, outputs cleaned `<trait>.QC.gz`.   |
| `03_gwas_mismatch.R`         | Harmonizes GWAS alleles with target data; outputs `.a1` (correct effect allele) and `.mismatch` SNP list. |
| **SNP selection**            |                                                                                                           |
| `04_clumping.sh`             | Performs LD clumping using two p-value thresholds (`1e-4` and `1e-5`); outputs `<trait>.valid.snp`.       |
| **Allele Frequency**         |                                                                                                           |
| `05_compute_pop_AF.sh`       | Computes allele frequency (AF) for AFR/AMR/EAS/EUR/SAS using 1000G.                                       |
| `06_make_all_AF.R`           | Builds a master AF table containing variant AF for all populations.                                       |
| `07_plot_AF_distributions.R` | Extracts AF for clumped SNPs and plots AF distribution histograms for each trait.                         |
| **PGS computation**          |                                                                                                           |
| `08_run_pgs.sh`              | Calculates PGS for each trait using valid SNPs and harmonized alleles.                                    |
| `09_compute_PGS_mean.R`      | Computes mean PGS for each (population × sex), merges subpopulations, outputs global z-scored means.      |
| **Regression analyses**      |                                                                                                           |
| `regression-Height.R`        | Regresses mean Height on Height PGS.                                                                      |
| `regression-BMI.R`           | Regresses BMI prevalence on BMI PGS.                                                                      |
| `regression-HDL.R`           | Regresses mean HDL on HDL PGS.                                                                            |
| `regression-SBP.R`           | Regresses mean SBP on SBP PGS.                                                                            |
