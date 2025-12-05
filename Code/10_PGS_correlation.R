library(dplyr)
library(ggplot2)
library(readr)
library(showtext)

font_add_google("Noto Sans", "noto")
showtext_auto()

sib_file <- "data/trait_data/BMI_pop_male_2013_final.csv"

std_file <- sub("BMI_", "BMI_sd_", sib_file)

cat("Sibling file:  ", sib_file, "\n")
cat("Standard file: ", std_file, "\n")

sib <- read_csv(sib_file) %>%
  select(pop_filter, Country, SuperPop, PGS_sib = PGS_global_z)

std <- read_csv(std_file) %>%
  select(pop_filter, Country, SuperPop, PGS_std = PGS_global_z)

df <- inner_join(std, sib, by = c("pop_filter", "Country", "SuperPop"))

cor_pearson  <- cor(df$PGS_std, df$PGS_sib, method = "pearson")
cor_spearman <- cor(df$PGS_std, df$PGS_sib, method = "spearman")

cat("Pearson correlation =", round(cor_pearson, 3), "\n")
cat("Spearman correlation =", round(cor_spearman, 3), "\n")

text_label <- paste0(
  "Pearson r = ", round(cor_pearson, 3), "\n",
  "Spearman ρ = ", round(cor_spearman, 3)
)

p <- ggplot(df, aes(x = PGS_std, y = PGS_sib, color = SuperPop)) +
  geom_point(size = 3) +
  geom_text(aes(label = pop_filter),
            vjust = -0.6, size = 3,show.legend = FALSE
  ) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  annotate("text",
           x = max(df$PGS_std),
           y = min(df$PGS_sib),
           label = text_label,
           hjust = 1, vjust = 0, size = 4) +
  labs(
    x = "Mean PGS (Standard GWAS, global z-score)",
    y = "Mean PGS (Sibling GWAS, global z-score)",
    title = "Standard vs Sibling Mean PGS (BMI, Male)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

print(p)


pdf("BMI_male_standard_vs_sibling_PGS.pdf", width = 6.5, height = 5)
print(p)
dev.off()
