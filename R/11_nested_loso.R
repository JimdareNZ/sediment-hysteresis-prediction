# ══════════════════════════════════════════════════════════════════════════════
# 11_nested_loso.R
# ══════════════════════════════════════════════════════════════════════════════
# Nested leave-one-site-out (LOSO) cross-validation (Reviewer 2, round 2,
# Comment 1).
#
# In 04_xgboost_modelling.R, hyperparameters were tuned by 5-fold CV on the
# full dataset before LOSO evaluation, so each held-out site contributed to
# model selection. Here, tuning is repeated inside every LOSO training fold:
# the held-out site is never seen during tuning or fitting.
#
# Everything else matches 04 exactly: data prep, recipe, feature set, engine
# settings, parameter ranges, 60-candidate space-filling grid (seed 653),
# inner 5-fold CV stratified by direction (seed 653), selection by ROC AUC,
# Clockwise as the positive class.
#
# Input:  data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv
# Output: results/nested_loso_results.csv         (per site x subset)
#         results/nested_loso_summary.csv         (per subset)
#         results/nested_loso_hyperparameters.csv (selected per outer fold)
#         tables/TableS11_nested_loso.tex
#
# Runtime: ~2 hours (3 subsets x 8 outer folds x 60 candidates x 5 inner folds).
# Completed outer folds are cached in results/nested_loso_cache/, so an
# interrupted run resumes where it stopped.
# ══════════════════════════════════════════════════════════════════════════════

library(tidymodels)
library(xgboost)
library(dials)

tidymodels_prefer()

dir.create("results/nested_loso_cache", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# ── 1. Data and features (as in 04) ──────────────────────────────────────────
raw <- read.csv("data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv",
                stringsAsFactors = FALSE)

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

subset_defs <- list(
  full = list(classes = 1:4,
              label_fn = function(zc) ifelse(zc %in% c(1, 2), "Clockwise", "Anticlockwise")),
  pure = list(classes = c(1, 4),
              label_fn = function(zc) ifelse(zc == 1, "Clockwise", "Anticlockwise")),
  fig8 = list(classes = c(2, 3),
              label_fn = function(zc) ifelse(zc == 2, "Clockwise", "Anticlockwise"))
)

prepare_subset <- function(raw_data, subset_name) {
  sdef <- subset_defs[[subset_name]]
  raw_data |>
    dplyr::filter(Zuec_class %in% sdef$classes) |>
    dplyr::mutate(direction = factor(sdef$label_fn(Zuec_class),
                                     levels = c("Clockwise", "Anticlockwise"))) |>
    dplyr::select(all_of(c("Site", all_features, "direction"))) |>
    tidyr::drop_na()
}

# ── 2. Tunable workflow, parameter ranges and grid (as in 04) ────────────────
build_tune_wf <- function(train) {
  rec <- recipe(direction ~ ., data = train) |>
    update_role(Site, new_role = "ID") |>
    step_zv(all_predictors()) |>
    step_nzv(all_predictors())

  spec <- boost_tree(
    trees = tune(), tree_depth = tune(), learn_rate = tune(),
    sample_size = tune(), mtry = 0.7, min_n = tune(), loss_reduction = tune()
  ) |>
    set_engine("xgboost", counts = FALSE, alpha = 0.01, lambda = 1.0,
               eval_metric = "logloss",
               nthread = parallel::detectCores() - 1) |>
    set_mode("classification")

  workflow() |> add_recipe(rec) |> add_model(spec)
}

xgb_params <- parameters(
  trees(range          = c(100L, 1000L)),
  tree_depth(range     = c(2L, 6L)),
  learn_rate(range     = c(-3, -1), trans = log10_trans()),
  dials::sample_prop(range = c(0.5, 1.0)),
  min_n(range          = c(1L, 20L)),
  loss_reduction(range = c(-5, 0), trans = log10_trans())
)

set.seed(653)
xgb_grid <- grid_space_filling(xgb_params, size = 60)

tune_metrics <- metric_set(roc_auc, accuracy, f_meas)

# ── 3. One outer fold: tune on 7 sites, predict the held-out site ────────────
run_outer_fold <- function(dat, held_out, subset_name) {
  cache_file <- file.path("results/nested_loso_cache",
                          paste0(subset_name, "_", held_out, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))

  train <- dplyr::filter(dat, Site != held_out)
  test  <- dplyr::filter(dat, Site == held_out)

  set.seed(653)
  inner_folds <- vfold_cv(train, v = 5, strata = direction)

  wf <- build_tune_wf(train)
  tuned <- tune_grid(wf, resamples = inner_folds, grid = xgb_grid,
                     metrics = tune_metrics,
                     control = control_grid(save_pred = FALSE, allow_par = FALSE))

  best    <- select_best(tuned, metric = "roc_auc")
  fit_wf  <- finalize_workflow(wf, best) |> fit(data = train)

  preds <- predict(fit_wf, test, type = "prob") |>
    dplyr::bind_cols(predict(fit_wf, test)) |>
    dplyr::mutate(direction = test$direction, Site = held_out)

  out <- list(
    preds  = preds,
    params = best |> dplyr::select(-.config) |>
      dplyr::mutate(subset = subset_name, held_out = held_out)
  )
  saveRDS(out, cache_file)
  out
}

# Site-level metric that returns NA when undefined (e.g. single-class site)
safe_metric <- function(fn, data, ...) {
  tryCatch(suppressWarnings(fn(data, truth = direction, ...)$.estimate),
           error = function(e) NA_real_)
}

# ── 4. Main loop ─────────────────────────────────────────────────────────────
set.seed(2026)
all_preds  <- list()
all_params <- list()
t_start    <- Sys.time()

for (sub in names(subset_defs)) {
  dat   <- prepare_subset(raw, sub)
  sites <- sort(unique(dat$Site))
  cat(sprintf("── %s (N = %d, %d sites) ──\n", sub, nrow(dat), length(sites)))

  for (s in sites) {
    res <- run_outer_fold(dat, s, sub)
    all_preds[[paste(sub, s)]]  <- dplyr::mutate(res$preds, subset = sub)
    all_params[[paste(sub, s)]] <- res$params
    cat(sprintf("  held out %-7s done [%.1f min]\n", s,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}

preds_df  <- dplyr::bind_rows(all_preds)
params_df <- dplyr::bind_rows(all_params)

# ── 5. Metrics: per-site (Table 5 convention) and pooled ─────────────────────
per_site <- preds_df |>
  dplyr::group_by(subset, Site) |>
  dplyr::group_modify(~ tibble(
    n        = nrow(.x),
    roc_auc  = safe_metric(roc_auc,  .x, .pred_Clockwise),
    accuracy = safe_metric(accuracy, .x, estimate = .pred_class),
    f1       = safe_metric(f_meas,   .x, estimate = .pred_class)
  )) |>
  dplyr::ungroup()

pooled <- preds_df |>
  dplyr::group_by(subset) |>
  dplyr::group_modify(~ tibble(
    pooled_roc_auc  = safe_metric(roc_auc,  .x, .pred_Clockwise),
    pooled_accuracy = safe_metric(accuracy, .x, estimate = .pred_class),
    pooled_f1       = safe_metric(f_meas,   .x, estimate = .pred_class)
  )) |>
  dplyr::ungroup()

# Original LOSO values from Table 5 (tuning on the full dataset)
original <- tribble(
  ~subset, ~orig_roc_auc, ~orig_accuracy, ~orig_f1,
  "full",  0.619,         0.561,          0.337,
  "pure",  0.792,         0.710,          0.578,
  "fig8",  0.525,         0.553,          0.273
)

summary_df <- per_site |>
  dplyr::group_by(subset) |>
  dplyr::summarise(
    nested_roc_auc  = mean(roc_auc,  na.rm = TRUE),
    nested_accuracy = mean(accuracy, na.rm = TRUE),
    nested_f1       = mean(f1,       na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(pooled,   by = "subset") |>
  dplyr::left_join(original, by = "subset") |>
  dplyr::mutate(subset = factor(subset, levels = c("full", "pure", "fig8"))) |>
  dplyr::arrange(subset)

cat("\n══ Nested vs original LOSO ══\n")
print(summary_df, width = Inf)

write.csv(per_site,   "results/nested_loso_results.csv",         row.names = FALSE)
write.csv(summary_df, "results/nested_loso_summary.csv",         row.names = FALSE)
write.csv(params_df,  "results/nested_loso_hyperparameters.csv", row.names = FALSE)

# ── 6. Table S11 ─────────────────────────────────────────────────────────────
f3 <- function(x) formatC(x, format = "f", digits = 3)
subset_lab <- c(full = "All classes (I+II vs III+IV)",
                pure = "Pure (I vs IV)",
                fig8 = "Figure-eight (II vs III)")

rows <- with(summary_df, paste0(
  subset_lab[as.character(subset)], " & ",
  f3(orig_roc_auc), " & ", f3(orig_accuracy), " & ", f3(orig_f1), " & ",
  f3(nested_roc_auc), " & ", f3(nested_accuracy), " & ", f3(nested_f1), " \\\\"
))

latex <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{Leave-one-site-out (LOSO) performance with hyperparameters tuned ",
         "on the full dataset (original, as in Table 5 of the main text) and with ",
         "hyperparameters re-tuned within each LOSO training fold (nested). In the ",
         "nested scheme, the held-out site was excluded from both tuning and fitting. ",
         "Values are means across held-out sites; F1 is calculated with clockwise as ",
         "the positive class.}"),
  "\\label{tab:nested_loso}", "\\small",
  "\\begin{tabular}{lrrrrrr}", "\\toprule",
  " & \\multicolumn{3}{c}{Original LOSO} & \\multicolumn{3}{c}{Nested LOSO} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
  "Class subset & ROC AUC & Accuracy & F1 & ROC AUC & Accuracy & F1 \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"
)

writeLines(latex, "tables/TableS11_nested_loso.tex")
cat("\nWritten: tables/TableS11_nested_loso.tex\n")