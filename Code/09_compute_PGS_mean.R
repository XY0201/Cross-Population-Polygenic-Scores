#!/usr/bin/env Rscript

library(dplyr)
library(readr)

# Read trait name
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: compute_pgs_mean.R <trait>")
}
TRAIT <- args[1]

# File paths
BASE <- "data"
PGS_FILE   <- file.path(BASE, "PGS", paste0(TRAIT, "_PGS.profile"))
PANEL_FILE <- file.path(BASE, "panel.csv")
OUTDIR     <- file.path(BASE, "PGS_mean")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Load data
pgs <- read.table(PGS_FILE, header = TRUE)
panel <- read.csv(PANEL_FILE)

# Merge CHN / IND / YEN
panel2 <- panel %>%
  mutate(pop_filter = case_when(
    pop %in% c("CHB", "CHS", "CDX") ~ "CHN",
    pop %in% c("GIH", "ITU") ~ "IND",
    pop %in% c("YRI", "ESN") ~ "YEN",
    TRUE ~ pop
  ))

# Merge PGS + population info
merged <- pgs %>%
  rename(ID = IID) %>%
  left_join(panel2, by = c("ID" = "sample"))

# Step 1: raw means per population × gender
subpop_raw <- merged %>%
  filter(!is.na(gender)) %>%
  group_by(pop_filter, gender) %>%
  summarize(
    PGS_mean_raw = mean(SCORE, na.rm = TRUE),
    PGS_sd_raw   = sd(SCORE,   na.rm = TRUE),
    n            = n(),
    .groups = "drop"
  )

# Step 2: global z-score
subpop_raw <- subpop_raw %>%
  mutate(PGS_global_z = as.numeric(scale(PGS_mean_raw)))

# Output
OUTFILE <- file.path(OUTDIR, paste0(TRAIT, "_PGS_global_mean_by_sex.csv"))
write.csv(subpop_raw, OUTFILE, row.names = FALSE)

print(subpop_raw)
