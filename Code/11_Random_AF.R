#library(dplyr)
#library(readr)
#library(ggplot2)
#library(tidyr)
#library(scales)

TRAIT <- "BMI"   
BASE <- "data/"

CLUMP_FILE <- file.path(BASE, "prePGS", paste0(TRAIT, "_clump.clumped"))
AF_FILE <- file.path(BASE, "all_snp_AF.tsv")

OUT_SNPLIST <- file.path(BASE, "Random SNPs", paste0(TRAIT, "_clumped.snplist"))
OUT_RANDOM <- file.path(BASE, "Random SNPs", paste0(TRAIT, "_random_EURmatched.snplist"))
OUT_PLOT <- file.path(BASE, "figures/af_distributions", paste0(TRAIT, "_Random.pdf"))
pops <- c("AFR","AMR","EAS","EUR","SAS")

clump <- read_table(CLUMP_FILE, comment = "")

gwas <- clump %>%
  filter(!is.na(SNP), SNP != "") %>%   
  select(SNP)

write.table(
  gwas,
  file = OUT_SNPLIST,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

#af <- read_tsv(AF_FILE)

breaks_seq <- seq(0, 1, by = 0.02)

af2 <- af %>%
  mutate(
    bin = cut(
      EUR,                      
      breaks = breaks_seq,
      include.lowest = TRUE,
      right = FALSE
    ),
    in_gwas = SNP %in% gwas$SNP  
  )

set.seed(2002) 

af2_counts <- af2 %>%
  filter(!in_gwas, !is.na(bin)) %>%
  left_join(
    af2 %>%
      group_by(bin) %>%
      summarise(n_gwas = sum(in_gwas), .groups = "drop"),
    by = "bin"
  )

random_df <- data.frame()

bins <- unique(af2_counts$bin)

set.seed(2002)

for (b in bins) {
  
  df_bin <- af2_counts[af2_counts$bin == b, ]
  
  k <- unique(df_bin$n_gwas)
  
  sampled <- df_bin[sample(seq_len(nrow(df_bin)), size = k, replace = FALSE), ]
  
  random_df <- rbind(random_df, sampled)
}

write.table(
  random_df$SNP,
  file = OUT_RANDOM,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

random_af <- dplyr::semi_join(af, random_df, by = "SNP")
df_af <- random_af[, c("AFR", "AMR", "EAS", "EUR", "SAS")]

df_long <- df_af %>%
  pivot_longer(
    cols = everything(),
    names_to = "Population",
    values_to = "AF"
  ) %>%
  filter(!is.na(AF), AF >= 0, AF < 0.22)

breaks_seq_plot <- seq(0, 0.22, by = 0.02)

bin_labels <- paste0(
  round(head(breaks_seq_plot, -1) * 100), "-",
  round(tail(breaks_seq_plot, -1) * 100), "%"
)

df_long <- df_long %>%
  mutate(
    bin = cut(
      AF,
      breaks = breaks_seq_plot,
      include.lowest = TRUE,
      right = FALSE,
      labels = bin_labels
    )
  )

df_long$bin <- factor(df_long$bin, levels = bin_labels, ordered = TRUE)

df_sum <- df_long %>%
  group_by(Population, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Population) %>%
  mutate(Proportion = count / sum(count)) %>%
  ungroup()

total_n <- nrow(random_df)

legend_labels <- setNames(
  paste0(pops, " (n=", format(total_n, big.mark=","), ")"),
  pops
)

p <- ggplot(df_sum, aes(x = bin, y = Proportion, fill = Population)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_brewer(palette="Set1", labels=legend_labels) +
  labs(
    title = paste0("Distribution of Effect Allele Frequencies - ", TRAIT, " (Random SNPs)"),
    x = "Effect Allele Frequency (EAF)",
    y = "Proportion"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

ggsave(OUT_PLOT, p, width = 12, height = 5)
