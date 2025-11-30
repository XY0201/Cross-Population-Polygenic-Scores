#!/usr/bin/env Rscript

# Batch extract allele frequencies for clumped SNPs
# and generate distribution plots for each trait.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(scales)
})

# Paths (modify as needed)
BASE <- "data"
AF_FILE   <- file.path(BASE, "all_snp_AF.tsv")
VALID_DIR <- file.path(BASE, "prePGS")
OUT_AF    <- file.path(BASE, "base_with_AF_clumped")
OUT_FIG   <- file.path(BASE, "figures", "af_distributions")

dir.create(OUT_AF, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# Read master AF table
af_all <- read_tsv(AF_FILE, col_types = cols())

# Find trait SNP files
valid_files <- list.files(VALID_DIR, pattern = "\\.valid\\.snp$", full.names = TRUE)
if (length(valid_files) == 0) stop("No .valid.snp files found in prePGS.")

# Plot function
plot_AF <- function(df, trait, pop_cols) {
  pops <- intersect(c("AFR", "AMR", "EAS", "EUR", "SAS"), pop_cols)
  if (length(pops) == 0) stop(paste("No population AF columns found for", trait))
  
  df <- df %>% rename_with(~ paste0("AF_", .), all_of(pops))
  df_long <- df %>%
    select(starts_with("AF_")) %>%
    pivot_longer(everything(), names_to = "Population", values_to = "AF") %>%
    mutate(
      Population = sub("AF_", "", Population),
      AF = as.numeric(AF)
    ) %>%
    filter(!is.na(AF), AF >= 0, AF < 0.22)
  
  breaks_seq <- seq(0, 0.22, by = 0.02)
  bin_labels <- paste0(
    round(head(breaks_seq, -1) * 100), "-",
    round(tail(breaks_seq, -1) * 100), "%"
  )
  
  df_long$bin <- cut(
    df_long$AF,
    breaks = breaks_seq,
    include.lowest = TRUE,
    right = FALSE,
    labels = bin_labels
  )
  
  df_sum <- df_long %>%
    group_by(Population, bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Population) %>%
    mutate(Proportion = count / sum(count)) %>%
    ungroup()
  
  p <- ggplot(df_sum, aes(x = bin, y = Proportion, fill = Population)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = paste0("Allele Frequency Distribution - ", trait),
      x = "Effect Allele Frequency (EAF)",
      y = "Proportion"
    ) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  outfile <- file.path(OUT_FIG, paste0(trait, ".pdf"))
  ggsave(outfile, p, width = 12, height = 5)
}

# Main loop
for (vf in valid_files) {
  trait <- sub("\\.valid\\.snp$", "", basename(vf))
  
  valid <- read_tsv(vf, col_names = FALSE, col_types = cols())
  colnames(valid) <- "SNP"
  
  df_trait <- af_all %>% semi_join(valid, by = "SNP")
  out_tsv <- file.path(OUT_AF, paste0(trait, "_AF.tsv"))
  write_tsv(df_trait, out_tsv)
  
  pop_cols <- intersect(names(df_trait), c("AFR", "AMR", "EAS", "EUR", "SAS"))
  plot_AF(df_trait, trait, pop_cols)
}
