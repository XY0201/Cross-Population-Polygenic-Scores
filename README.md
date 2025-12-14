# Cross-Population-Polygenic-Scores

## 📄 **Abstract**

This project compares polygenic scores (PGS) derived from standard GWAS and sibling-based GWAS to evaluate their cross-population predictive performance. Using 1000 Genomes data, the analysis relates population-level mean PGS to corresponding traits, including height, BMI, SBP, and HDL.

---

## 📄 **Data**

### 1. Target data
1000 Genomes Phase 3 genotype data used as the target dataset.  
**Source:** https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

### 2. Base GWAS data
GWAS summary statistics used as base datasets for computing polygenic scores.  
**Source:** https://opengwas.io/datasets/

| GWAS ID     | Trait     |
|-------------|-----------|
| ieu-b-4813  | Height    |
| ieu-b-4814  | Height_sd |
| ieu-b-4815  | BMI       |
| ieu-b-4816  | BMI_sd    |
| ieu-b-4817  | SBP       |
| ieu-b-4818  | SBP_sd    |
| ieu-b-4843  | HDL       |
| ieu-b-4844  | HDL_sd    |

### 3. Phenotype data
Population-level trait values (height, BMI, SBP, HDL) used for regression analyses.  
**Source:** https://ncdrisc.org/index.html

### 4. Uploaded Data
- **panel_with_pop_filter.csv** — Processed 1000G sample panel with merged sub-population labels.  
- **region.xlsx** — Mapping between 1000G populations and their corresponding countries/regions for trait assignment.  
- **regression/** — Final sex-stratified population-level tables linking phenotypes with PGS results.

---

## 📄 **Code**

### Environment

This project was run using R 4.5 and PLINK 1.9.

### Overview of all scripts used in this project

<table>
  <tr><th>Script</th><th>Description</th></tr>

  <!-- Quality Control -->
  <tr>
    <td colspan="2" align="center"><b>Quality Control</b></td>
  </tr>
  <tr><td><code>01_target_prepare_1000G.sh</code></td><td>Cleans raw 1000G VCFs, converts to PLINK, merges chr1–22, and performs standard target QC.</td></tr>
  <tr><td><code>02_gwas_qc.sh</code></td><td>Extracts fields from GWAS-VCF, filters SNPs, assigns CHR:BP:REF:ALT.</td></tr>
  <tr><td><code>03_gwas_mismatch.R</code></td><td>Harmonizes GWAS alleles with target data and identifies allele mismatches.</td></tr>

  <!-- SNP selection -->
  <tr>
    <td colspan="2" align="center"><b>SNP selection</b></td>
  </tr>
  <tr><td><code>04_clumping.sh</code></td><td>Performs LD clumping to select independent SNPs for each trait using p-value threshold = 1e-4.</td></tr>

  <!-- Allele Frequency -->
  <tr>
    <td colspan="2" align="center"><b>Allele Frequency</b></td>
  </tr>
  <tr><td><code>05_compute_pop_AF.sh</code></td><td>Computes allele frequency for AFR/AMR/EAS/EUR/SAS populations using 1000G.</td></tr>
  <tr><td><code>06_make_all_AF.R</code></td><td>Builds a master AF table containing variant AF across all populations.</td></tr>
  <tr><td><code>07_plot_AF_distributions.R</code></td><td>Extracts AF for clumped SNPs and creates AF distribution histograms for each trait.</td></tr>

  <!-- PGS computation -->
  <tr>
    <td colspan="2" align="center"><b>PGS computation</b></td>
  </tr>
  <tr><td><code>08_run_pgs.sh</code></td><td>Calculates PGS for each trait based on valid SNPs and harmonized alleles.</td></tr>
  <tr><td><code>09_compute_PGS_mean.R</code></td><td>Computes mean PGS for each (population × sex), merges subpopulations, and generates global z-scored means.</td></tr>
  <tr><td><code>10_PGS_correlation.R</code></td><td>Analyzes and visualizes correlations between population-mean PGS from standard and sibling GWAS.</td></tr>
  <tr><td><code>11_Random_AF.R</code></td><td>Generates EUR allele-frequency–matched random SNP sets and analyzes their allele frequency distributions as a control for GWAS SNPs.</td></tr>

  <!-- Regression analyses -->
  <tr>
    <td colspan="2" align="center"><b>Regression analyses</b></td>
  </tr>
  <tr><td><code>regression-Height.R</code></td><td>Regresses mean Height on Height PGS.</td></tr>
  <tr><td><code>regression-BMI.R</code></td><td>Regresses BMI prevalence on BMI PGS.</td></tr>
  <tr><td><code>regression-HDL.R</code></td><td>Regresses mean HDL on HDL PGS.</td></tr>
  <tr><td><code>regression-SBP.R</code></td><td>Regresses mean SBP on SBP PGS.</td></tr>

</table>

---

## 📄 **Results**

- **af_distributions/** — Allele-frequency distribution plots for each trait.

- **regression/** — Regression summary tables (`*_reg.csv`, `*_sd_reg.csv`) and sex-stratified scatterplots comparing population-level PGS with phenotype means.

- **PGS correlation/** — Sex-stratified scatterplots comparing population-mean PGS from standard and sibling-based GWAS across super-populations.
