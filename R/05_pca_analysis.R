# ══════════════════════════════════════════════════════════════════════════════
# 05_pca_analysis.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Principal component analysis of event-scale hydrological metrics and
#   spatial predictor indices, with ANOVA tests for hysteresis direction
#   and sub-catchment effects on PC1-PC4.
#
# Outputs:
#   Figures:
#     - figures/Fig5_PCA_Combined.png (main: PC1 vs PC2 + PC3 vs PC4)
#     - figures/FigS1_PCA_PC1_PC3.png (supplementary)
#     - figures/FigS2_PCA_PC1_PC4.png (supplementary)
#   Tables:
#     - Printed ANOVA summary for Table S6
#
# Inputs:
#   - data/Final_Hysteresis_Dataset_Spatial_Ante.csv
#
# Dependencies:
#   tidyverse, ggrepel, cowplot
#
# Reference:
#   Dare et al. (2026). Event-scale prediction of sediment hysteresis
#   regimes using hydrological and spatial indices. Journal of Hydrology.
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggrepel)
library(cowplot)

# =========================================================================
# 1. Load and prepare data
# =========================================================================
df <- read_csv("data/Final_Hysteresis_Dataset_Spatial_Ante.csv")

# Site name and sub-catchment mapping
df <- df %>%
  mutate(
    Site_Code = case_when(
      Site == "Tuangi" ~ "PKO-D",
      Site == "Inanga" ~ "PKO-UW",
      Site == "Kokopu" ~ "PKO-UE",
      Site == "Kokako" ~ "MGT",
      Site == "Pipi"   ~ "PUN",
      Site == "Kotuku" ~ "WHR-D",
      Site == "Matata" ~ "PNG-M",
      Site == "Tuna"   ~ "PNG-TR1"
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
        Zuec_class == 1 ~ "I",
        Zuec_class == 2 ~ "II",
        Zuec_class == 3 ~ "III",
        Zuec_class == 4 ~ "IV"
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

# Add scores to dataframe
for (i in 1:4) {
  df_pca[[paste0("PC", i)]] <- pca_result$x[, i]
}

cat("Variance explained:\n")
for (i in 1:4) cat(paste0("  PC", i, ": ", round(var_exp[i], 1), "%\n"))
cat(paste0("  Cumulative PC1-PC4: ", round(sum(var_exp[1:4]), 1), "%\n"))

# =========================================================================
# 3. ANOVA tests
# =========================================================================

# --- Direction ---
cat("\n── ANOVA: PC ~ Direction ──\n")
direction_anova <- tibble(PC = paste0("PC", 1:4), Variance = round(var_exp[1:4], 1))

for (i in 1:4) {
  s <- summary(aov(as.formula(paste0("PC", i, " ~ Direction")), data = df_pca))[[1]]
  direction_anova$F_value[i] <- round(s$`F value`[1], 1)
  direction_anova$p_value[i] <- signif(s$`Pr(>F)`[1], 3)
  direction_anova$eta_sq[i]  <- round(s$`Sum Sq`[1] / sum(s$`Sum Sq`), 3)
}
print(direction_anova)

# --- Sub-catchment ---
cat("\n── ANOVA: PC ~ Sub_Catchment ──\n")
sc_anova <- tibble(PC = paste0("PC", 1:4), Variance = round(var_exp[1:4], 1))

for (i in 1:4) {
  s <- summary(aov(as.formula(paste0("PC", i, " ~ Sub_Catchment")), data = df_pca))[[1]]
  sc_anova$F_value[i] <- round(s$`F value`[1], 1)
  sc_anova$p_value[i] <- signif(s$`Pr(>F)`[1], 3)
  sc_anova$eta_sq[i]  <- round(s$`Sum Sq`[1] / sum(s$`Sum Sq`), 3)
}
print(sc_anova)

# --- Summary table for manuscript (Table S6) ---
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
# 4. Variable labels and colours
# =========================================================================
var_labels <- c(
  "RB_flashiness"     = "Flashiness",
  "rising_time_h"     = "Rising time",
  "baseflow_index"    = "Baseflow index",
  "recession_k_per_h" = "Recession rate",
  "symmetry_index"    = "Hydrograph symmetry",
  "Zuec_h"            = "Hysteresis index (h)",
  "Yield_kg_ha"       = "Sediment yield",
  "SEI"               = "SEI",
  "CCI"               = "CCI",
  "LUSP"              = "LUSP"
)

sc_colours <- c(
  "Pokopoko"    = "#f4a582",
  "Mangatoetoe" = "#ca0020",
  "Puanene"     = "#9970ab",
  "Wharere"     = "#404040",
  "Pongakawa"   = "#92c5de"
)

# =========================================================================
# 5. Reusable biplot function
# =========================================================================
make_pca_biplot <- function(pca_result, df_pca, pc_x = 1, pc_y = 2,
                            var_labels, sc_colours, var_exp,
                            nudge_overrides = NULL) {

  pc_x_name <- paste0("PC", pc_x)
  pc_y_name <- paste0("PC", pc_y)

  scores <- as_tibble(pca_result$x) %>%
    bind_cols(df_pca %>% select(Site_Code, Sub_Catchment)) %>%
    rename(PCx = !!pc_x_name, PCy = !!pc_y_name)

  loadings <- as_tibble(pca_result$rotation, rownames = "Variable") %>%
    select(Variable, PCx = !!pc_x_name, PCy = !!pc_y_name)

  scale_factor <- min(
    max(abs(scores$PCx)) / max(abs(loadings$PCx)),
    max(abs(scores$PCy)) / max(abs(loadings$PCy))
  ) * 0.55

  loadings <- loadings %>%
    mutate(
      PCx_scaled = PCx * scale_factor,
      PCy_scaled = PCy * scale_factor,
      Label = var_labels[Variable]
    )

  site_centroids <- scores %>%
    group_by(Site_Code, Sub_Catchment) %>%
    summarise(PCx = mean(PCx), PCy = mean(PCy), .groups = "drop")

  # Split labels into auto-repelled and manually nudged
  if (!is.null(nudge_overrides)) {
    nudge_labels   <- names(nudge_overrides)
    loadings_clean <- loadings %>% filter(!Label %in% nudge_labels)
    loadings_nudge <- loadings %>%
      filter(Label %in% nudge_labels) %>%
      mutate(
        nudge_x = map_dbl(Label, ~ nudge_overrides[[.x]][1]),
        nudge_y = map_dbl(Label, ~ nudge_overrides[[.x]][2])
      )
  } else {
    loadings_clean <- loadings
    loadings_nudge <- loadings %>% filter(FALSE)
  }

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey80", linewidth = 0.3) +
    stat_ellipse(
      data = scores,
      aes(x = PCx, y = PCy, colour = Sub_Catchment, fill = Sub_Catchment),
      type = "t", level = 0.95, geom = "polygon",
      alpha = 0.12, linewidth = 0.6
    ) +
    geom_point(
      data = scores,
      aes(x = PCx, y = PCy, fill = Sub_Catchment),
      shape = 21, size = 1.8, alpha = 0.6, colour = "white", stroke = 0.3
    ) +
    geom_segment(
      data = loadings,
      aes(x = 0, y = 0, xend = PCx_scaled, yend = PCy_scaled),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      colour = "grey20", linewidth = 0.6
    ) +
    geom_point(
      data = site_centroids,
      aes(x = PCx, y = PCy),
      shape = 23, size = 3.5, fill = "white", colour = "grey20", stroke = 0.8
    ) +
    geom_label_repel(
      data = site_centroids,
      aes(x = PCx, y = PCy, label = Site_Code, fill = Sub_Catchment),
      size = 3, label.padding = unit(0.12, "lines"),
      box.padding = unit(0.3, "lines"), min.segment.length = 0,
      seed = 42, colour = "#fdf6e3", segment.color = "grey20",
      fontface = "bold", label.size = 0.4, show.legend = FALSE
    ) +
    geom_label_repel(
      data = loadings_clean,
      aes(x = PCx_scaled, y = PCy_scaled, label = Label),
      size = 2.8, label.padding = unit(0.1, "lines"),
      box.padding = unit(0.5, "lines"), min.segment.length = 0,
      force = 3, seed = 42, colour = "grey20", fill = "white", alpha = 0.85
    )

  if (nrow(loadings_nudge) > 0) {
    p <- p +
      geom_label(
        data = loadings_nudge,
        aes(x = PCx_scaled + nudge_x, y = PCy_scaled + nudge_y, label = Label),
        size = 2.8, label.padding = unit(0.1, "lines"),
        colour = "grey20", fill = "white", alpha = 0.85
      )
  }

  p + scale_fill_manual(values = sc_colours, name = "Sub-catchment") +
    scale_colour_manual(values = sc_colours, name = "Sub-catchment") +
    labs(
      x = paste0("PC", pc_x, " (", round(var_exp[pc_x], 1), "% variance explained)"),
      y = paste0("PC", pc_y, " (", round(var_exp[pc_y], 1), "% variance explained)")
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey95", linewidth = 0.3),
      axis.text        = element_text(size = 9),
      axis.title       = element_text(size = 10),
      legend.position  = "right",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.5, "cm"),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(5, 5, 5, 5, "mm")
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.8)), colour = "none")
}

# =========================================================================
# 6. Generate figures
# =========================================================================
dir.create("figures", showWarnings = FALSE)

# Main figure panels
fig5a <- make_pca_biplot(
  pca_result, df_pca, pc_x = 1, pc_y = 2,
  var_labels, sc_colours, var_exp,
  nudge_overrides = list(
    "Hydrograph symmetry" = c(-0.8, 0.3),
    "LUSP"                = c(0.5, -0.4)
  )
)

fig5b <- make_pca_biplot(
  pca_result, df_pca, pc_x = 3, pc_y = 4,
  var_labels, sc_colours, var_exp,
  nudge_overrides = NULL
)

# Combined Figure 5
fig5 <- plot_grid(
  fig5a + theme(legend.position = "none"),
  fig5b + theme(legend.position = "none"),
  labels = c("(a)", "(b)"),
  ncol = 2, align = "h", axis = "tb"
)

legend     <- get_legend(fig5a)
fig5_final <- plot_grid(fig5, legend, ncol = 2, rel_widths = c(1, 0.15))

ggsave("figures/Fig5_PCA_Combined.png", fig5_final,
       dpi = 300, width = 14, height = 7, units = "in", bg = "white")

# Supplementary figures
fig_s1 <- make_pca_biplot(
  pca_result, df_pca, pc_x = 1, pc_y = 3,
  var_labels, sc_colours, var_exp
)
ggsave("figures/FigS1_PCA_PC1_PC3.png", fig_s1,
       dpi = 300, width = 10, height = 7, units = "in", bg = "white")

fig_s2 <- make_pca_biplot(
  pca_result, df_pca, pc_x = 1, pc_y = 4,
  var_labels, sc_colours, var_exp
)
ggsave("figures/FigS2_PCA_PC1_PC4.png", fig_s2,
       dpi = 300, width = 10, height = 7, units = "in", bg = "white")

cat("All PCA figures saved to figures/\n")
cat("PCA analysis complete.\n")
