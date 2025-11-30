# Cross-Population-Polygenic-Scores

## **Abstract**

This project compares polygenic scores (PGS) derived from standard GWAS and sibling-based GWAS to evaluate their cross-population predictive correlations. Using 1000 Genomes data, the analysis relates population-level mean PGS to corresponding traits, including height, BMI, SBP, and HDL.

---

## **Code**

### Environment

This project was run using R 4.5 and PLINK 1.9.

### Overview of all scripts used in this project

<table>
  <tr><th>Script</th><th>Description</th></tr>

  <!-- Quality Control -->
  <tr>
    <td colspan="2" align="center"><b>Quality Control</b></td>
  </tr>
  <tr><td>01_target_prepare_1000G.sh</td><td>Cleans raw 1000G VCFs, converts to PLINK, merges chr1–22, and performs standard target QC.</td></tr>
  <tr><td>02_gwas_qc.sh</td><td>Extracts fields from GWAS-VCF, filters SNPs, assigns CHR:BP:REF:ALT.</td></tr>
  <tr><td>03_gwas_mismatch.R</td><td>Harmonizes GWAS alleles with target data and identifies allele mismatches.</td></tr>

  <!-- SNP selection -->
  <tr>
    <td colspan="2" align="center"><b>SNP selection</b></td>
  </tr>
  <tr><td>04_clumping.sh</td><td>Performs LD clumping to select independent SNPs for each trait using p-value threshold = 1e-4.</td></tr>

  <!-- Allele Frequency -->
  <tr>
    <td colspan="2" align="center"><b>Allele Frequency</b></td>
  </tr>
  <tr><td>05_compute_pop_AF.sh</td><td>Computes allele frequency for AFR/AMR/EAS/EUR/SAS populations using 1000G.</td></tr>
  <tr><td>06_make_all_AF.R</td><td>Builds a master AF table containing variant AF across all populations.</td></tr>
  <tr><td>07_plot_AF_distributions.R</td><td>Extracts AF for clumped SNPs and creates AF distribution histograms for each trait.</td></tr>

  <!-- PGS computation -->
  <tr>
    <td colspan="2" align="center"><b>PGS computation</b></td>
  </tr>
  <tr><td>08_run_pgs.sh</td><td>Calculates PGS for each trait based on valid SNPs and harmonized alleles.</td></tr>
  <tr><td>09_compute_PGS_mean.R</td><td>Computes mean PGS for each (population × sex), merges subpopulations, and generates global z-scored means.</td></tr>

  <!-- Regression analyses -->
  <tr>
    <td colspan="2" align="center"><b>Regression analyses</b></td>
  </tr>
  <tr><td>regression-Height.R</td><td>Regresses mean Height on Height PGS.</td></tr>
  <tr><td>regression-BMI.R</td><td>Regresses BMI prevalence on BMI PGS.</td></tr>
  <tr><td>regression-HDL.R</td><td>Regresses mean HDL on HDL PGS.</td></tr>
  <tr><td>regression-SBP.R</td><td>Regresses mean SBP on SBP PGS.</td></tr>

</table>


