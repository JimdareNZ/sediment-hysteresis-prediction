# ══════════════════════════════════════════════════════════════════════════════
# 05_pca_analysis.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Principal component analysis of event-scale hydrological metrics and
#   spatial predictor indices, with ANOVA tests for hysteresis direction
#   and sub-catchment effects on PC1-PC4.
#
# Outputs:
#   - output/pca_results.rds (PCA object, scores, variance explained)
#   - Printed ANOVA summary for Table S6
#
# Inputs:
#   - data/Final_Hysteresis_Dataset_Spatial_Ante.csv
#
# Dependencies:
#   tidyverse
#
# Reference:
#   Dare et al. (2026). Event-scale prediction of sediment hysteresis
#   regimes using hydrological and spatial indices. Journal of Hydrology.
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)

# =========================================================================
# 1. Load and prepare data
# =========================================================================
df <- read_csv("data/Final_Hysteresis_Dataset_Spatial_Ante.csv")

df <- df %>%
  mutate(
    Site_Code = case_when(
      Site == "Tuangi" ~ "PKO-D",  Site == "Inanga" ~ "PKO-UW",
      Site == "Kokopu" ~ "PKO-UE", Site == "Kokako" ~ "MGT",
      Site == "Pipi"   ~ "PUN",    Site == "Kotuku" ~ "WHR-D",
      Site == "Matata" ~ "PNG-M",  Site == "Tuna"   ~ "PNG-TR1"
    ),
    Sub_Catchment = case_when(
      Site %in% c("Tuangi", "Inanga", "Kokopu") ~ "Pokopoko",
      Site == "Kokako"                           ~ "Mangatoetoe",
      Site == "Pipi"                             ~ "Puanene",
      Site == "Kotuku"                           ~ "Wharere",
      Site %in% c("Matata", "Tuna")              ~ "Pongakawa"
    ),
    Direction = factor(
      ifelse(Zuec_class %in% c(1, 2), "Clockwise", "Anticlockwise"),
      levels = c("Clockwise", "Anticlockwise")
    ),
    Class_Label = factor(
      case_when(
        Zuec_class == 1 ~ "I",  Zuec_class == 2 ~ "II",
        Zuec_class == 3 ~ "III", Zuec_class == 4 ~ "IV"
      ),
      levels = c("I", "II", "III", "IV")
    )
  ) %>%
  filter(Zuec_class %in% 1:4)

# =========================================================================
# 2. PCA
# =========================================================================
pca_vars <- c(
  "RB_flashiness", "rising_time_h", "baseflow_index", "recession_k_per_h",
  "symmetry_index", "Zuec_h", "Yield_kg_ha", "SEI", "CCI", "LUSP"
)

df_pca <- df %>%
  select(Site_Code, Sub_Catchment, Direction, Class_Label, all_of(pca_vars)) %>%
  filter(complete.cases(.))

pca_result <- prcomp(df_pca %>% select(all_of(pca_vars)),
                     center = TRUE, scale. = TRUE)

var_exp <- summary(pca_result)$importance[2, ] * 100

for (i in 1:4) {
  df_pca[[paste0("PC", i)]] <- pca_result$x[, i]
}

cat("Variance explained:\n")
for (i in 1:4) cat(paste0("  PC", i, ": ", round(var_exp[i], 1), "%\n"))
cat(paste0("  Cumulative PC1-PC4: ", round(sum(var_exp[1:4]), 1), "%\n"))

# =========================================================================
# 3. ANOVA tests
# =========================================================================
cat("\n── ANOVA: PC ~ Direction ──\n")
direction_anova <- tibble(PC = paste0("PC", 1:4), Variance = round(var_exp[1:4], 1))

for (i in 1:4) {
  s <- summary(aov(as.formula(paste0("PC", i, " ~ Direction")), data = df_pca))[[1]]
  direction_anova$F_value[i] <- round(s$`F value`[1], 1)
  direction_anova$p_value[i] <- signif(s$`Pr(>F)`[1], 3)
  direction_anova$eta_sq[i]  <- round(s$`Sum Sq`[1] / sum(s$`Sum Sq`), 3)
}
print(direction_anova)

cat("\n── ANOVA: PC ~ Sub_Catchment ──\n")
sc_anova <- tibble(PC = paste0("PC", 1:4), Variance = round(var_exp[1:4], 1))

for (i in 1:4) {
  s <- summary(aov(as.formula(paste0("PC", i, " ~ Sub_Catchment")), data = df_pca))[[1]]
  sc_anova$F_value[i] <- round(s$`F value`[1], 1)
  sc_anova$p_value[i] <- signif(s$`Pr(>F)`[1], 3)
  sc_anova$eta_sq[i]  <- round(s$`Sum Sq`[1] / sum(s$`Sum Sq`), 3)
}
print(sc_anova)

cat("\n── Summary table for Table S6 ──\n")
summary_table <- tibble(
  PC                     = paste0("PC", 1:4),
  `Variance (%)`         = round(var_exp[1:4], 1),
  `F (Direction)`        = direction_anova$F_value,
  `p (Direction)`        = direction_anova$p_value,
  `eta2 (Direction)`     = direction_anova$eta_sq,
  `F (Sub-catchment)`    = sc_anova$F_value,
  `p (Sub-catchment)`    = sc_anova$p_value,
  `eta2 (Sub-catchment)` = sc_anova$eta_sq
)
print(summary_table, width = Inf)

# =========================================================================
# 4. Save results for figure generation
# =========================================================================
dir.create("output", showWarnings = FALSE)

saveRDS(
  list(pca_result = pca_result, df_pca = df_pca, var_exp = var_exp),
  file = "output/pca_results.rds"
)

cat("\nPCA results saved to output/pca_results.rds\n")
cat("PCA analysis complete.\n")
