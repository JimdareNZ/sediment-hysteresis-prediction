# ══════════════════════════════════════════════════════════════════════════════
# 06_fig5_pca_biplots.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 5: PCA biplots (PC1 vs PC2, PC3 vs PC4) — main text
# Figures S1-S2: Supplementary PCA biplots (PC1 vs PC3, PC1 vs PC4)
#
# Input:  output/pca_results.rds (from 05_pca_analysis.R)
# Output:
#   figures/Fig5_PCA_Combined.png
#   figures/FigS1_PCA_PC1_PC3.png
#   figures/FigS2_PCA_PC1_PC4.png
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(ggrepel)
library(cowplot)

# =========================================================================
# 1. Load PCA results
# =========================================================================
pca_data   <- readRDS("output/pca_results.rds")
pca_result <- pca_data$pca_result
df_pca     <- pca_data$df_pca
var_exp    <- pca_data$var_exp

# =========================================================================
# 2. Variable labels and colours
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
# 3. Reusable biplot function
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
# 4. Generate figures
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
cat("Figure 5 saved\n")

# Supplementary figures
fig_s1 <- make_pca_biplot(
  pca_result, df_pca, pc_x = 1, pc_y = 3,
  var_labels, sc_colours, var_exp
)
ggsave("figures/FigS1_PCA_PC1_PC3.png", fig_s1,
       dpi = 300, width = 10, height = 7, units = "in", bg = "white")
cat("Figure S1 saved\n")

fig_s2 <- make_pca_biplot(
  pca_result, df_pca, pc_x = 1, pc_y = 4,
  var_labels, sc_colours, var_exp
)
ggsave("figures/FigS2_PCA_PC1_PC4.png", fig_s2,
       dpi = 300, width = 10, height = 7, units = "in", bg = "white")
cat("Figure S2 saved\n")

cat("All PCA figures complete.\n")
