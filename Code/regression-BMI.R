############################################################
# Process BMI data to 1000G subpopulation-level traits
# Year = 2013; Men/Women separately
############################################################

library(dplyr)
library(readr)
library(ggplot2)
library(showtext)

font_add_google("Noto Sans", "noto")
showtext_auto()

# 1. Read raw BMI data
BMI <- read_csv("data/trait_data/NCD_RisC_Lancet_2024_BMI_child_adolescent_country_ageStd.csv")
colnames(BMI)[3] <- 'Country'

# 2. Filter Year = 2013
BMI_2013_country_mean <- BMI %>%
  filter(Year == 2013) %>%
  group_by(Country, Sex) %>%
  summarise(mean_BMI = mean(`Prevalence of BMI > 2SD (obesity)`, na.rm = TRUE)) %>%
  arrange(Country, Sex)

# 3. 1000G pop → country map
pop_map <- tribble(
  ~pop, ~country,
  "ACB", "Barbados",
  "BEB", "Bangladesh",
  "CHN", "China",
  "CEU", "United States of America",
  "CLM", "Colombia",
  "FIN", "Finland",
  "GBR", "United Kingdom",
  "IND", "India",
  "GWD", "Gambia",
  "IBS", "Spain",
  "JPT", "Japan",
  "KHV", "Viet Nam",
  "LWK", "Kenya",
  "MSL", "Sierra Leone",
  "MXL", "Mexico",
  "PEL", "Peru",
  "PJL", "Pakistan",
  "PUR", "Puerto Rico",
  "STU", "Sri Lanka",
  "TSI", "Italy",
  "YEN", "Nigeria"
)

# 4. Merge BMI with 1000G mapping
merged <- BMI_2013_country_mean %>%
  inner_join(pop_map, by = c("Country" = "country"))

# Split by sex
male_BMI    <- merged %>% filter(Sex == "Boys")
female_BMI  <- merged %>% filter(Sex == "Girls")

############################################################
# ↓↓↓ Below is your PGS merging logic (renamed to BMI) ↓↓↓
############################################################

#sibling
BMI_PGS <- read_csv("Desktop/GWAS/data/PGS_mean/BMI_PGS_global_mean_by_sex.csv")
BMI_PGS$pop <- BMI_PGS$pop_filter
#standard
BMI_PGS <- read_csv("Desktop/GWAS/data/PGS_mean/BMI_sd_PGS_global_mean_by_sex.csv")
BMI_PGS$pop <- BMI_PGS$pop_filter

# Map superpop
superpop_map <- tribble(
  ~pop,  ~SuperPop,
  "BEB", "SAS",
  "ACB", "AFR",
  "CEU", "EUR",
  "CHN", "EAS",
  "CLM", "AMR",
  "FIN", "EUR",
  "GWD", "AFR",
  "IND", "SAS",
  "TSI", "EUR",
  "JPT", "EAS",
  "LWK", "AFR",
  "MXL", "AMR",
  "YEN", "AFR",
  "PJL", "SAS",
  "PEL", "AMR",
  "PUR", "AMR",
  "MSL", "AFR",
  "IBS", "EUR",
  "STU", "SAS",
  "GBR", "EUR",
  "KHV", "EAS"
)

# Merge phenotype + PGS
male_dat <- male_BMI %>%
  inner_join(BMI_PGS %>% filter(gender == "male"),
             by = c("pop" = "pop")) %>%
  inner_join(superpop_map, by = "pop")

female_dat <- female_BMI %>%
  inner_join(BMI_PGS %>% filter(gender == "female"),
             by = c("pop" = "pop")) %>%
  inner_join(superpop_map, by = "pop")

#sibling
write_csv(male_dat, "data/trait_data/BMI_pop_male_2013_final.csv")
write_csv(female_dat, "data/trait_data/BMI_pop_female_2013_final.csv")
#standard
write_csv(male_dat, "data/trait_data/BMI_sd_pop_male_2013_final.csv")
write_csv(female_dat, "data/trait_data/BMI_sd_pop_female_2013_final.csv")

############################################################
# Regression models
############################################################

model_male   <- lm(mean_BMI ~ PGS_global_z, data = male_dat)
model_female <- lm(mean_BMI ~ PGS_global_z, data = female_dat)

cor_male   <- cor.test(male_dat$PGS_global_z, male_dat$mean_BMI)
cor_female <- cor.test(female_dat$PGS_global_z, female_dat$mean_BMI)

# Create summary
male_summary <- data.frame(
  group = "Men",
  beta  = coef(model_male)[2],
  r     = cor_male$estimate,
  p     = summary(model_male)$coefficients[2,4],
  R2    = summary(model_male)$r.squared
)

female_summary <- data.frame(
  group = "Women",
  beta  = coef(model_female)[2],
  r     = cor_female$estimate,
  p     = summary(model_female)$coefficients[2,4],
  R2    = summary(model_female)$r.squared
)

summary_table <- rbind(male_summary, female_summary)
summary_table

#sibling
write.csv(summary_table,'data/figures/regression/BMI_reg.csv')
#standard
write.csv(summary_table,'data/figures/regression/BMI_sd_reg.csv')

############################################################
# Plot function (HEIGHT → BMI version)
############################################################

plot_regression <- function(data, model, cor_obj, title, outfile) {
  
  txt <- paste(
    sprintf("β₁ = %.3f", coef(model)[2]),
    sprintf("r = %.3f", cor_obj$estimate),
    sprintf("p = %s", format.pval(summary(model)$coefficients[2,4], digits = 2, eps = 0.001)),
    sep = "\n"
  )
  
  slope <- coef(model)[2]
  
  if (slope > 0) {
    ann_x <- max(data$PGS_global_z)
    ann_y <- min(data$mean_BMI)
    ann_hjust <- 1
    ann_vjust <- 0
  } else {
    ann_x <- max(data$PGS_global_z)
    ann_y <- max(data$mean_BMI)
    ann_hjust <- 1
    ann_vjust <- 1
  }
  
  p <- ggplot(data, aes(x = PGS_global_z, y = mean_BMI, color = SuperPop)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_smooth(method = "lm",
                color = "black",
                fill = "grey80",
                linewidth = 0.8,
                alpha = 0.4) +
    geom_text(
      aes(label = pop),
      vjust = -0.6,
      size = 3,
      show.legend = FALSE
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = title,
      x = "PGS (global z-score)",
      y = "Mean obesity prevalence (BMI > 2SD)",
      color = "SuperPop"
    ) +
    annotate(
      "text",
      x = ann_x,
      y = ann_y,
      label = txt,
      hjust = ann_hjust,
      vjust = ann_vjust,
      size = 4
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  pdf(outfile, height = 5, width = 6.5)
  print(p)
  dev.off()
}

############################################################
# Draw plots
############################################################

# 1. Male – sibling GWAS
plot_regression(
  data = male_dat,
  model = model_male,
  cor_obj = cor_male,
  title = "BMI vs PGS (Sibling GWAS) - Male",
  outfile = "data/figures/regression/BMI_sibling_male.pdf"
)

# 2. Female – sibling GWAS
plot_regression(
  data = female_dat,
  model = model_female,
  cor_obj = cor_female,
  title = "BMI vs PGS (Sibling GWAS) - Female",
  outfile = "data/figures/regression/BMI_sibling_female.pdf"
)

# 3. Male – standard GWAS
plot_regression(
  data = male_dat,
  model = model_male,
  cor_obj = cor_male,
  title = "BMI vs PGS (Standard GWAS) - Male",
  outfile = "data/figures/regression/BMI_standard_male.pdf"
)

# 4. Female – standard GWAS
plot_regression(
  data = female_dat,
  model = model_female,
  cor_obj = cor_female,
  title = "BMI vs PGS (Standard GWAS) - Female",
  outfile = "data/figures/regression/BMI_standard_female.pdf"
)

