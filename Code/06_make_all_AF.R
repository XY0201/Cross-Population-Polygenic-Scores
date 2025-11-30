#!/usr/bin/env Rscript

library(dplyr)
library(readr)

base <- "data/pop_AF"

afr <- read_table(file.path(base, "AFR_freq.frq")) %>% select(SNP, AFR = MAF)
amr <- read_table(file.path(base, "AMR_freq.frq")) %>% select(SNP, AMR = MAF)
eas <- read_table(file.path(base, "EAS_freq.frq")) %>% select(SNP, EAS = MAF)
eur <- read_table(file.path(base, "EUR_freq.frq")) %>% select(SNP, EUR = MAF)
sas <- read_table(file.path(base, "SAS_freq.frq")) %>% select(SNP, SAS = MAF)

all <- afr %>%
  full_join(amr, by="SNP") %>%
  full_join(eas, by="SNP") %>%
  full_join(eur, by="SNP") %>%
  full_join(sas, by="SNP")

write_tsv(all, "data/all_snp_AF.tsv")
