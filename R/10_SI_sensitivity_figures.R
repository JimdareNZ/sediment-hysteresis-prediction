# ══════════════════════════════════════════════════════════════════════════════
# 10_SI_sensitivity_figures.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure S5: Label perturbation sensitivity (Reviewer 2, Comment 1)
# Figure S6: CCI/LUSP weighting sensitivity (Reviewer 2, Comment 3)
#
# Input:  results/sensitivity_results.csv            (08_sensitivity_label_perturbation.R)
#         results/weighting_sensitivity_results.csv  (09_sensitivity_weighting.R)
# Output: figures/FigS5_sensitivity.png
#         figures/FigS6_weighting_sensitivity.png
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)

dir.create("figures", showWarnings = FALSE)

subset_levels <- c("pure", "full", "fig8")
subset_labels <- c(
  pure = "Pure unidirectional (I vs IV)",
  full = "Full binary (I+II vs III+IV)",
  fig8 = "Figure-eight (II vs III)"
)

# ── Figure S5: Label perturbation sensitivity ────────────────────────────────
results_label <- read.csv("results/sensitivity_results.csv", stringsAsFactors = FALSE)

baseline <- results_label |>
  filter(flip_rate == 0) |>
  select(subset, baseline_auc = roc_auc) |>
  mutate(subset = factor(subset, levels = subset_levels))

fig_s5 <- results_label |>
  filter(flip_rate > 0) |>
  mutate(
    flip_pct = factor(paste0(flip_rate * 100, "%"), levels = c("10%", "20%", "33%")),
    subset   = factor(subset, levels = subset_levels)
  ) |>
  ggplot(aes(x = flip_pct, y = roc_auc)) +
  geom_boxplot(fill = "grey80", outlier.size = 0.5) +
  geom_hline(data = baseline, aes(yintercept = baseline_auc),
             linetype = "dashed", colour = "black") +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "red") +
  facet_wrap(~ subset, labeller = labeller(subset = subset_labels)) +
  labs(x = "Label flip rate", y = "LOSO ROC AUC",
       caption = "Dashed line = unperturbed baseline; dotted red line = random chance (0.5)") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))

ggsave("figures/FigS5_sensitivity.png", fig_s5,
       width = 8, height = 3.5, dpi = 300, bg = "white")

# ── Figure S6: Weighting sensitivity ─────────────────────────────────────────
results_weight <- read.csv("results/weighting_sensitivity_results.csv", stringsAsFactors = FALSE)

scenario_order <- c("Original", "CCI equal", "CCI reversed",
                    "LUSP equal", "LUSP steep", "Both equal", "Both adversarial")

fig_s6 <- results_weight |>
  mutate(
    scenario = factor(scenario, levels = rev(scenario_order)),
    subset   = factor(subset, levels = subset_levels)
  ) |>
  ggplot(aes(x = roc_auc, y = scenario)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "red") +
  facet_wrap(~ subset, labeller = labeller(subset = subset_labels)) +
  labs(x = "LOSO ROC AUC", y = NULL,
       caption = "Dotted red line = random chance (0.5)") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9), axis.text.y = element_text(size = 8))

ggsave("figures/FigS6_weighting_sensitivity.png", fig_s6,
       width = 9, height = 4, dpi = 300, bg = "white")