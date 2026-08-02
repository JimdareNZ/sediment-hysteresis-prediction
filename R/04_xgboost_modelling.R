# ══════════════════════════════════════════════════════════════════════════════
# 04_xgboost_modelling.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Train an XGBoost classifier to predict hysteresis direction from
#   hydrological and spatial predictors. Evaluates performance using
#   5-fold stratified CV and leave-one-site-out (LOSO) CV.
#
# Class subsets:
#   Run this script three times, changing the 'class_subset' parameter:
#     "all"       — classes I+II (CW) vs III+IV (ACW), n = 678
#     "pure"      — class I (CW) vs class IV (ACW), n = 227
#     "fig_eight" — class II (CW) vs class III (ACW), n = 451
#
# Inputs:
#   - data/Final_Hysteresis_Dataset_Spatial_Ante.csv
#
# Outputs:
#   - output/models/<subset>_figure_data.rds (ROC, AUC, confusion, importance)
#   - output/models/<subset>_final_fit.rds
#   - output/models/<subset>_workflow.rds
#
# Notes:
#   - Load_kg, Yield_kg_ha, and Load_kg_prev are excluded from the predictor
#     set because they require continuous SSC monitoring.
#   - Hyperparameters are tuned via space-filling grid search (60 candidates).
#
# Dependencies:
#   tidymodels, xgboost, vip, dials
#
# Reference:
#   Dare et al. (2026). Event-scale prediction of sediment hysteresis
#   regimes using hydrological and spatial indices. Journal of Hydrology.
# ══════════════════════════════════════════════════════════════════════════════

library(tidymodels)
library(xgboost)
library(vip)
library(dials)

tidymodels_prefer()

# =========================================================================
# Configuration — change this to run each subset
# =========================================================================
class_subset <- "all"  # Options: "all", "pure", "fig_eight"

output_dir <- "output/models"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =========================================================================
# 1. Load and prepare data
# =========================================================================
raw <- read.csv("data/Final_Hysteresis_Dataset_Spatial_Ante.csv",
                stringsAsFactors = FALSE)

# Predictor groups (Load_kg, Yield_kg_ha, Load_kg_prev excluded)
hydrological <- c(
  "Q_peak", "rising_time_h", "recession_time_h", "symmetry_index",
  "RB_flashiness", "peak_sharpness", "baseflow_index", "recession_k_per_h",
  "n_peaks", "VQ_m3", "vQ_mm", "duration_h",
  "n_obs", "mean_dt_s", "Q_fixed_min", "Q_fixed_max"
)

antecedent <- c(
  "dt_prev_event_h", "Q_peak_prev", "VQ_m3_prev",
  "n_events_prev_30d", "n_events_prev_14d", "n_events_prev_7d"
)

spatial <- c(
  "QCI", "Slope", "DDI", "CCI", "SEI",
  "LU_high", "LU_mid", "LU_low", "LUSP",
  "DI", "SI", "Sed_Risk_mean", "Sed_Risk_weighted_mean", "Sed_Risk_tail_2_5"
)

all_features <- c(hydrological, antecedent, spatial)

# Filter and collapse classes based on subset
if (class_subset == "all") {
  dat <- raw |>
    dplyr::filter(Zuec_class %in% 1:4) |>
    dplyr::mutate(
      direction = factor(
        ifelse(Zuec_class %in% c(1, 2), "Clockwise", "Anticlockwise"),
        levels = c("Clockwise", "Anticlockwise")
      )
    )
  subset_label <- "all_classes"

} else if (class_subset == "pure") {
  dat <- raw |>
    dplyr::filter(Zuec_class %in% c(1, 4)) |>
    dplyr::mutate(
      direction = factor(
        ifelse(Zuec_class == 1, "Clockwise", "Anticlockwise"),
        levels = c("Clockwise", "Anticlockwise")
      )
    )
  subset_label <- "pure_I_vs_IV"

} else if (class_subset == "fig_eight") {
  dat <- raw |>
    dplyr::filter(Zuec_class %in% c(2, 3)) |>
    dplyr::mutate(
      direction = factor(
        ifelse(Zuec_class == 2, "Clockwise_Fig_8", "Anticlockwise_Fig_8"),
        levels = c("Clockwise_Fig_8", "Anticlockwise_Fig_8")
      )
    )
  subset_label <- "fig_eight_II_vs_III"

} else {
  stop("class_subset must be 'all', 'pure', or 'fig_eight'")
}

dat <- dat |>
  dplyr::select(all_of(c("Site", all_features, "direction"))) |>
  tidyr::drop_na()

cat("Subset:", subset_label, "\n")
cat("N =", nrow(dat), "\n")
cat("Class distribution:\n")
print(table(dat$direction))

# =========================================================================
# 2. Recipe
# =========================================================================
xgb_recipe <- recipe(direction ~ ., data = dat) |>
  update_role(Site, new_role = "ID") |>
  step_zv(all_predictors()) |>
  step_nzv(all_predictors())

# =========================================================================
# 3. CV folds
# =========================================================================
set.seed(42)
folds_5 <- vfold_cv(dat, v = 5, strata = direction)

tune_metrics <- metric_set(roc_auc, accuracy, f_meas)

# =========================================================================
# 4. Tunable model specification
# =========================================================================
xgb_spec_tune <- boost_tree(
  trees          = tune(),
  tree_depth     = tune(),
  learn_rate     = tune(),
  sample_size    = tune(),
  mtry           = 0.7,
  min_n          = tune(),
  loss_reduction = tune()
) |>
  set_engine(
    "xgboost",
    counts      = FALSE,
    alpha       = 0.01,
    lambda      = 1.0,
    eval_metric = "logloss",
    nthread     = parallel::detectCores() - 1
  ) |>
  set_mode("classification")

# =========================================================================
# 5. Workflow
# =========================================================================
xgb_wf_tune <- workflow() |>
  add_recipe(xgb_recipe) |>
  add_model(xgb_spec_tune)

# =========================================================================
# 6. Hyperparameter grid
# =========================================================================
xgb_params <- parameters(
  trees(range          = c(100L, 1000L)),
  tree_depth(range     = c(2L, 6L)),
  learn_rate(range     = c(-3, -1), trans = log10_trans()),
  dials::sample_prop(range = c(0.5, 1.0)),
  min_n(range          = c(1L, 20L)),
  loss_reduction(range = c(-5, 0), trans = log10_trans())
)

set.seed(42)
xgb_grid <- grid_space_filling(xgb_params, size = 60)

# =========================================================================
# 7. Hyperparameter tuning
# =========================================================================
cat("\nTuning hyperparameters (60 candidates x 5 folds)...\n")

tune_results <- tune_grid(
  xgb_wf_tune,
  resamples = folds_5,
  grid      = xgb_grid,
  metrics   = tune_metrics,
  control   = control_grid(save_pred = TRUE, verbose = TRUE, allow_par = FALSE)
)

best_params <- select_best(tune_results, metric = "roc_auc")
cat("\nBest hyperparameters:\n")
print(best_params)

# =========================================================================
# 8. Finalise and evaluate
# =========================================================================
xgb_wf_final <- finalize_workflow(xgb_wf_tune, best_params)

# 5-fold CV with tuned model
cat("\nFitting tuned model with 5-fold CV...\n")
set.seed(42)
cv_results <- fit_resamples(
  xgb_wf_final,
  resamples = folds_5,
  metrics   = tune_metrics,
  control   = control_resamples(save_pred = TRUE)
)

cat("\n5-fold CV results:\n")
print(collect_metrics(cv_results))

# LOSO CV
cat("\nFitting LOSO CV...\n")
set.seed(42)
loso_splits <- group_vfold_cv(dat, group = Site)

loso_results <- fit_resamples(
  xgb_wf_final,
  resamples = loso_splits,
  metrics   = tune_metrics,
  control   = control_resamples(save_pred = TRUE)
)

cat("\nLOSO results:\n")
print(collect_metrics(loso_results))

# =========================================================================
# 9. Final model fit (all data)
# =========================================================================
cat("\nFitting final model on all data...\n")
final_fit <- fit(xgb_wf_final, data = dat)

# Variable importance
importance_df <- final_fit |>
  extract_fit_parsnip() |>
  vip::vi()

cat("\nTop 15 features:\n")
print(head(importance_df, 15))

# =========================================================================
# 10. Compile and save outputs
# =========================================================================

# ROC data
cv_preds   <- collect_predictions(cv_results)
loso_preds <- collect_predictions(loso_results)

# Determine probability column name (varies by factor levels)
prob_col <- names(cv_preds)[grepl("^\\.pred_", names(cv_preds)) & !grepl("class", names(cv_preds))][1]

roc_5fold_data <- cv_preds |>
  roc_curve(truth = direction, !!sym(prob_col)) |>
  mutate(cv_type = "5-fold CV")

roc_loso_data <- loso_preds |>
  roc_curve(truth = direction, !!sym(prob_col)) |>
  mutate(cv_type = "LOSO CV")

roc_combined <- bind_rows(roc_5fold_data, roc_loso_data)

auc_5fold <- cv_preds |>
  roc_auc(truth = direction, !!sym(prob_col)) |>
  pull(.estimate)

auc_loso <- loso_preds |>
  roc_auc(truth = direction, !!sym(prob_col)) |>
  pull(.estimate)

# Confusion matrix
conf_mat_loso <- loso_preds |>
  conf_mat(truth = direction, estimate = .pred_class)

lvls <- levels(dat$direction)
conf_df <- tibble(
  Truth      = rep(lvls, each = 2),
  Prediction = rep(lvls, times = 2),
  n          = as.vector(conf_mat_loso$table)
)

cat("\n5-fold vs LOSO comparison:\n")
print(bind_rows(
  collect_metrics(cv_results)   |> mutate(cv_type = "5-fold"),
  collect_metrics(loso_results) |> mutate(cv_type = "LOSO")
) |> select(cv_type, .metric, mean, std_err))

# Save
saveRDS(
  list(
    roc_combined  = roc_combined,
    auc_5fold     = auc_5fold,
    auc_loso      = auc_loso,
    conf_df       = conf_df,
    importance_df = importance_df
  ),
  file = file.path(output_dir, paste0(subset_label, "_figure_data.rds"))
)

saveRDS(final_fit,    file.path(output_dir, paste0(subset_label, "_final_fit.rds")))
saveRDS(xgb_wf_final, file.path(output_dir, paste0(subset_label, "_workflow.rds")))

cat("\nOutputs saved to", output_dir, "with prefix:", subset_label, "\n")
cat("Modelling complete for", subset_label, "\n")
