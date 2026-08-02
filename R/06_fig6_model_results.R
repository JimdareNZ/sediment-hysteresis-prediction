# ══════════════════════════════════════════════════════════════════════════════
# 06_fig6_model_results.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 6: XGBoost model performance (ROC, confusion matrix, variable
# importance) for the pure unidirectional model (I vs IV).
# Also generates supplementary Figures S3 (all classes) and S4 (figure-eight).
#
# Inputs: output/models/*_figure_data.rds (from 04_xgboost_modelling.R)
# Outputs:
#   figures/Fig6_Model_Results.png
#   figures/FigS3_AllClasses_Model_Results.png
#   figures/FigS4_FigureEight_Model_Results.png
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(tidyverse)
library(cowplot)

# ── Load model outputs ────────────────────────────────────────────────────────
full_data <- readRDS("output/models/all_classes_figure_data.rds")
pure_data <- readRDS("output/models/pure_I_vs_IV_figure_data.rds")
fig8_data <- readRDS("output/models/fig_eight_II_vs_III_figure_data.rds")

# ── Feature groups and labels ─────────────────────────────────────────────────
hydrological <- c(
  "Q_peak", "rising_time_h", "recession_time_h", "symmetry_index",
  "RB_flashiness", "peak_sharpness", "baseflow_index", "recession_k_per_h",
  "n_peaks", "VQ_m3", "vQ_mm", "duration_h", "n_obs", "mean_dt_s",
  "Q_fixed_min", "Q_fixed_max"
)
antecedent <- c(
  "dt_prev_event_h", "Q_peak_prev", "VQ_m3_prev",
  "n_events_prev_30d", "n_events_prev_14d", "n_events_prev_7d"
)
spatial <- c(
  "QCI", "Slope", "DDI", "CCI", "SEI", "LU_high", "LU_mid", "LU_low",
  "LUSP", "DI", "SI", "Sed_Risk_mean", "Sed_Risk_weighted_mean", "Sed_Risk_tail_2_5"
)

var_labels <- c(
  "QCI" = "Quick Connectivity Index", "vQ_mm" = "Specific runoff",
  "baseflow_index" = "Baseflow index", "Slope" = "Slope",
  "Q_peak" = "Peak discharge", "Q_peak_prev" = "Peak discharge (prev.)",
  "symmetry_index" = "Hydrograph symmetry", "duration_h" = "Event duration",
  "recession_time_h" = "Recession time", "rising_time_h" = "Rising time",
  "peak_sharpness" = "Peak sharpness", "VQ_m3" = "Event volume",
  "VQ_m3_prev" = "Event volume (prev.)", "recession_k_per_h" = "Recession rate",
  "n_obs" = "No. observations", "RB_flashiness" = "Flashiness",
  "SEI" = "Soil Erodibility Index", "CCI" = "Composite Connectivity Index",
  "LUSP" = "Land Use Sediment Potential", "DDI" = "Deep Drainage Index",
  "dt_prev_event_h" = "Time since prev. event", "mean_dt_s" = "Mean time step",
  "Sed_Risk_tail_2_5" = "Sediment risk (tail)", "Sed_Risk_mean" = "Sediment risk (mean)",
  "n_events_prev_30d" = "Events prev. 30 days", "n_events_prev_14d" = "Events prev. 14 days",
  "n_events_prev_7d" = "Events prev. 7 days"
)

group_colours <- c("Hydrological" = "#4393c3", "Spatial" = "#d73027", "Antecedent" = "#74c476")

# ── Helper: smooth ROC ────────────────────────────────────────────────────────
interp_roc <- function(roc_df, cv_name) {
  d <- roc_df %>% filter(cv_type == cv_name) %>% arrange(1 - specificity)
  x_new <- seq(min(1 - d$specificity), max(1 - d$specificity), length.out = 500)
  y_new <- approx(1 - d$specificity, d$sensitivity, xout = x_new)$y
  tibble(specificity = 1 - x_new, sensitivity = y_new, cv_type = cv_name)
}

# ── Helper: confusion matrix panel ────────────────────────────────────────────
make_conf_panel <- function(conf_df, title, recode_levels = NULL) {
  if (!is.null(recode_levels)) {
    conf_df <- conf_df %>%
      mutate(Truth = recode(Truth, !!!recode_levels),
             Prediction = recode(Prediction, !!!recode_levels))
  }
  conf_df <- conf_df %>%
    mutate(Truth = factor(Truth, levels = c("Clockwise", "Anticlockwise")),
           Prediction = factor(Prediction, levels = c("Clockwise", "Anticlockwise"))) %>%
    group_by(Truth) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()

  ggplot(conf_df, aes(x = Truth, y = Prediction, fill = pct)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = paste0(n, "\n(", round(pct, 1), "%)")), size = 3.5, colour = "grey20") +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 100), name = "% of true class") +
    labs(x = "Observed", y = "Predicted", title = title) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 9),
          axis.title = element_text(size = 9), legend.position = "right",
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          legend.key.size = unit(0.5, "cm"), plot.title = element_text(size = 9, face = "bold"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(3, 3, 3, 3, "mm"))
}

# ── Helper: variable importance panel ─────────────────────────────────────────
make_imp_panel <- function(imp_df, title) {
  imp <- imp_df %>%
    slice_max(Importance, n = 15) %>%
    mutate(
      Label = ifelse(Variable %in% names(var_labels), var_labels[Variable], Variable),
      Label = factor(Label, levels = Label[order(Importance)]),
      Group = case_when(Variable %in% hydrological ~ "Hydrological",
                        Variable %in% antecedent ~ "Antecedent",
                        Variable %in% spatial ~ "Spatial")
    )
  ggplot(imp, aes(x = Label, y = Importance, fill = Group)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = group_colours, name = "Predictor group") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(imp$Importance) * 1.15)) +
    coord_flip() +
    labs(x = NULL, y = "Importance (gain)", title = title) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
          axis.text = element_text(size = 8), axis.title = element_text(size = 9),
          legend.position = "bottom", legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 7), legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(size = 9, face = "bold"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(3, 10, 3, 3, "mm"))
}

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6 — Pure model (I vs IV): ROC + confusion + importance
# ══════════════════════════════════════════════════════════════════════════════
dir.create("figures", showWarnings = FALSE)

roc_all <- bind_rows(
  interp_roc(full_data$roc_combined, "LOSO CV") %>% mutate(Model = "All classes"),
  interp_roc(pure_data$roc_combined, "LOSO CV") %>% mutate(Model = "Pure (I vs IV)"),
  interp_roc(fig8_data$roc_combined, "LOSO CV") %>% mutate(Model = "Figure-eight (II vs III)")
)

panel_a <- ggplot(roc_all, aes(x = 1 - specificity, y = sensitivity,
                               colour = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  annotate("text", x = 0.02, y = 0.98,
           label = paste0("Pure (I vs IV): AUC = ", round(pure_data$auc_loso, 3)),
           hjust = 0, size = 2.8, colour = "#2166ac", fontface = "bold") +
  annotate("text", x = 0.02, y = 0.91,
           label = paste0("All classes: AUC = ", round(full_data$auc_loso, 3)),
           hjust = 0, size = 2.8, colour = "#636363", fontface = "bold") +
  annotate("text", x = 0.02, y = 0.84,
           label = paste0("Figure-eight (II vs III): AUC = ", round(fig8_data$auc_loso, 3)),
           hjust = 0, size = 2.8, colour = "#d73027", fontface = "bold") +
  scale_colour_manual(values = c("All classes" = "#636363", "Pure (I vs IV)" = "#2166ac",
                                 "Figure-eight (II vs III)" = "#d73027")) +
  scale_linetype_manual(values = c("All classes" = "dotted", "Pure (I vs IV)" = "solid",
                                   "Figure-eight (II vs III)" = "dashed")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "(a) LOSO ROC curves by class subset") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        axis.text = element_text(size = 8), axis.title = element_text(size = 9),
        plot.title = element_text(size = 9, face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA))

panel_b <- make_conf_panel(pure_data$conf_df, "(b) Confusion matrix, I vs IV (LOSO)")
panel_c <- make_imp_panel(pure_data$importance_df, "(c) Variable importance, I vs IV")

top_row <- plot_grid(panel_a, panel_b, ncol = 2, align = "hv", rel_widths = c(1, 1))
fig6 <- plot_grid(top_row, panel_c, ncol = 1, rel_heights = c(1, 1))

ggsave("figures/Fig6_Model_Results.png", fig6,
       dpi = 300, width = 12, height = 10, units = "in", bg = "white")
cat("Figure 6 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S3 — All classes (I+II vs III+IV)
# ══════════════════════════════════════════════════════════════════════════════
fig_s3 <- plot_grid(
  make_conf_panel(full_data$conf_df, "(a) Confusion matrix, all classes (LOSO)"),
  make_imp_panel(full_data$importance_df, "(b) Variable importance, all classes"),
  ncol = 1, rel_heights = c(0.8, 1)
)

ggsave("figures/FigS3_AllClasses_Model_Results.png", fig_s3,
       dpi = 300, width = 10, height = 10, units = "in", bg = "white")
cat("Figure S3 saved\n")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S4 — Figure-eight (II vs III)
# ══════════════════════════════════════════════════════════════════════════════
fig8_recode <- c("Clockwise_Fig_8" = "Clockwise", "Anticlockwise_Fig_8" = "Anticlockwise")

fig_s4 <- plot_grid(
  make_conf_panel(fig8_data$conf_df, "(a) Confusion matrix, II vs III (LOSO)", recode_levels = fig8_recode),
  make_imp_panel(fig8_data$importance_df, "(b) Variable importance, II vs III"),
  ncol = 1, rel_heights = c(0.8, 1)
)

ggsave("figures/FigS4_FigureEight_Model_Results.png", fig_s4,
       dpi = 300, width = 10, height = 10, units = "in", bg = "white")
cat("Figure S4 saved\n")
