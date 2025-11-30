############################################################
# Process height data to 1000G subpopulation-level traits
# Adult = age 18–19; Year = 2013; male and female separately
############################################################

library(dplyr)
library(readr)
library(ggplot2)
library(showtext)

font_add_google("Noto Sans", "noto")
showtext_auto()

# 1. Read raw height data
height <- read_csv("data/trait_data/NCD_RisC_Lancet_2020_height_child_adolescent_country.csv")

# 2. Filter Year = 2013 and Age group = 18 or 19
height_2013_country_mean <- height %>%
  # Keep 2013 and age 18–19
  filter(Year == 2013, `Age group` %in% c(18, 19)) %>%
  # Group by Country and Sex
  group_by(Country, Sex) %>%
  # Directly compute mean height for 18–19
  summarise(mean_height = mean(`Mean height`, na.rm = TRUE)) %>%
  arrange(Country, Sex)

# 3. Prepare 1000G population → country mapping
pop_map <- tribble(
  ~pop, ~country,
  "ACB", "Barbados",
  #"ASW", "United States of America",
  "BEB", "Bangladesh",
  #"CDX", "China",
  "CHN", "China",
  #"CHB", "China",
  #"CHS", "China",
  "CEU", "United States of America",
  "CLM", "Colombia",
  #"ESN", "Nigeria",
  "FIN", "Finland",
  "GBR", "United Kingdom",
  "IND", "India",
  #"GIH", "India",
  #"ITU", "India",
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
  #"YRI", "Nigeria",
  "YEN", "Nigeria"
)

# 4. Merge height data with population mapping
merged <- height_2013_country_mean %>%
  inner_join(pop_map, by = c("Country" = "country"))

# Boys
male_height <- merged %>%
  filter(Sex == "Boys")

# Girls
female_height <- merged %>%
  filter(Sex == "Girls")

# Save
#write_csv(male_height, "data/trait_data/height_pop_male_2013_18_19_raw.csv")
#write_csv(female_height, "data/trait_data/height_pop_female_2013_18_19_raw.csv")

### regression
#sibling
Height_subpop_mean_filter <- read_csv("Desktop/GWAS/data/PGS_mean/Height_PGS_global_mean_by_sex.csv")
Height_subpop_mean_filter$pop <- Height_subpop_mean_filter$pop_filter
#standard
Height_subpop_mean_filter <- read_csv("Desktop/GWAS/data/PGS_mean/Height_sd_PGS_global_mean_by_sex.csv")
Height_subpop_mean_filter$pop <- Height_subpop_mean_filter$pop_filter
#EUR
Height_subpop_mean_filter <- subset(Height_subpop_mean_filter, pop %in% c('FIN','GBR','CEU','TSI','IBS'))

superpop_map <- tribble(
  ~pop,  ~SuperPop,
  "BEB", "SAS",
  "ACB", "AFR",
  "CEU", "EUR",
  
  # China merged → CHN
  "CHN", "EAS",
  
  "CLM", "AMR",
  "FIN", "EUR",
  "GWD", "AFR",
  
  # India merged → IND
  "IND", "SAS",
  
  "TSI", "EUR",
  "JPT", "EAS",
  "LWK", "AFR",
  "MXL", "AMR",
  
  # ESN + YRI merged → YEN
  "YEN", "AFR",
  
  "PJL", "SAS",
  "PEL", "AMR",
  "PUR", "AMR",
  "MSL", "AFR",
  "IBS", "EUR",
  "STU", "SAS",
  "GBR", "EUR",
  #"ASW", "AFR",
  "KHV", "EAS"
)

# Merge files
male_dat <- male_height %>%
  inner_join(Height_subpop_mean_filter %>% filter(gender == "male"),
             by = c("pop" = "pop")) %>%
  inner_join(superpop_map, by="pop")

female_dat <- female_height %>%
  inner_join(Height_subpop_mean_filter %>% filter(gender == "female"),
             by = c("pop" = "pop")) %>%
  inner_join(superpop_map, by="pop")

#sibling
write_csv(male_dat, "data/trait_data/height_pop_male_2013_18_19_final.csv")
write_csv(female_dat, "data/trait_data/height_pop_female_2013_18_19_final.csv")
#standard
write_csv(male_dat, "data/trait_data/height_sd_pop_male_2013_18_19_final.csv")
write_csv(female_dat, "data/trait_data/height_sd_pop_female_2013_18_19_final.csv")
#EUR_sibling
write_csv(male_dat, "data/trait_data/height_EUR_pop_male_2013_18_19_final.csv")
write_csv(female_dat, "data/trait_data/height_EUR_pop_female_2013_18_19_final.csv")
#EUR_standard
write_csv(male_dat, "data/trait_data/height_EUR_sd_pop_male_2013_18_19_final.csv")
write_csv(female_dat, "data/trait_data/height_EUR_sd_pop_female_2013_18_19_final.csv")

# regression
model_male <- lm(mean_height ~ PGS_global_z, data = male_dat)
#summary(model_male)

model_female <- lm(mean_height ~ PGS_global_z, data = female_dat)
#summary(model_female)

# Boys correlation
cor_male <- cor.test(male_dat$PGS_global_z, male_dat$mean_height)

# Girls correlation
cor_female <- cor.test(female_dat$PGS_global_z, female_dat$mean_height)

# Boys summary
male_summary <- data.frame(
  group = "Boys",
  beta = coef(model_male)[2],
  r = cor_male$estimate,
  p = summary(model_male)$coefficients[2,4],
  R2 = summary(model_male)$r.squared
)

# Girls summary
female_summary <- data.frame(
  group = "Girls",
  beta = coef(model_female)[2],
  r = cor_female$estimate,
  p = summary(model_female)$coefficients[2,4],
  R2 = summary(model_female)$r.squared
)

summary_table <- rbind(male_summary, female_summary)
summary_table
#sibling
write.csv(summary_table,'data/figures/regression/Height_reg.csv')
#standard
write.csv(summary_table,'data/figures/regression/Height_sd_reg.csv')
#EUR_sibling
write.csv(summary_table,'data/figures/regression/Height_EUR_reg.csv')
#EUR_standard
write.csv(summary_table,'data/figures/regression/Height_EUR_sd_reg.csv')

lm_m <- summary(model_male)
#txt_m <- paste0("β1 = ", round(coef(model_male)[2], 2),"\nr = ", round(cor_male$estimate, 2), "\np = ", signif(lm_m$coefficients[2,4], 3))


###
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
    ann_y <- min(data$mean_height)
    ann_hjust <- 1
    ann_vjust <- 0
  } else {
    ann_x <- max(data$PGS_global_z)
    ann_y <- max(data$mean_height)
    ann_hjust <- 1
    ann_vjust <- 1
  }
  
  p <- ggplot(data, aes(x = PGS_global_z, y = mean_height, color = SuperPop)) +
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
      y = "Mean height (cm)",
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


# 1. Male – sibling GWAS
plot_regression(
  data = male_dat,
  model = model_male,
  cor_obj = cor_male,
  title = "Height vs PGS (Sibling GWAS) - Male",
  outfile = "data/figures/regression/Height_sibling_male.pdf"
)

# 2. Female – sibling GWAS
plot_regression(
  data = female_dat,
  model = model_female,
  cor_obj = cor_female,
  title = "Height vs PGS (Sibling GWAS) - Female",
  outfile = "data/figures/regression/Height_sibling_female.pdf"
)

# 3. Male – standard GWAS
plot_regression(
  data = male_dat,
  model = model_male,
  cor_obj = cor_male,
  title = "Height vs PGS (Standard GWAS) - Male",
  outfile = "data/figures/regression/Height_standard_male.pdf"
)

# 4. Female – standard GWAS
plot_regression(
  data = female_dat,
  model = model_female,
  cor_obj = cor_female,
  title = "Height vs PGS (Standard GWAS) - Female",
  outfile = "data/figures/regression/Height_standard_female.pdf"
)

# 5. Male – EUR sibling GWAS
plot_regression(
  data = male_dat,
  model = model_male,
  cor_obj = cor_male,
  title = "Height vs PGS (Sibling GWAS) - Male",
  outfile = "data/figures/regression/Height_EUR_male.pdf"
)

# 6. Female – EUR sibling GWAS
plot_regression(
  data = female_dat,
  model = model_female,
  cor_obj = cor_female,
  title = "Height vs PGS (Sibling GWAS) - Female",
  outfile = "data/figures/regression/Height_EUR_female.pdf"
)

# 7. Male – EUR sibling GWAS
plot_regression(
  data = male_dat,
  model = model_male,
  cor_obj = cor_male,
  title = "Height vs PGS (Standard GWAS) - Male",
  outfile = "data/figures/regression/Height_EUR_sd_male.pdf"
)

# 8. Female – EUR sibling GWAS
plot_regression(
  data = female_dat,
  model = model_female,
  cor_obj = cor_female,
  title = "Height vs PGS (Standard GWAS) - Female",
  outfile = "data/figures/regression/Height_EUR_sd_female.pdf"
)
