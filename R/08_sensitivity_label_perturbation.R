# ══════════════════════════════════════════════════════════════════════════════
# 08_sensitivity_label_perturbation.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Label perturbation sensitivity analysis (Reviewer 2, Comment 1).
#   Wraps the existing tidymodels XGBoost pipeline in a perturbation loop.
#   For each flip rate, randomly flips CW/ACW direction labels and reruns
#   the LOSO evaluation using the original tuned hyperparameters.
#
# Inputs:
#   - data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv
#
# Outputs:
#   - results/sensitivity_results.csv
#   - results/sensitivity_summary.csv
#   Figure S5 is produced from these outputs by 10_SI_sensitivity_figures.R.
#
# Notes:
#   - Runtime: ~8 hours.
#   - The reported results were generated before a seed was fixed; the
#     archived CSVs in results/ reproduce the reported values exactly. Reruns
#     with the seed below are reproducible but will differ slightly from the
#     reported medians. Rerunning overwrites the archived CSVs.
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Packages ──────────────────────────────────────────────────────────────
library(tidymodels)
library(xgboost)
library(scales)

tidymodels_prefer()

set.seed(2026)
dir.create("results", showWarnings = FALSE)

# ── 1. Load & prepare data ───────────────────────────────────────────────────
raw <- read.csv(
  "data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv",
  stringsAsFactors = FALSE
)

# Feature groups — matches original pipeline exactly
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

# ── 2. Class subset definitions ──────────────────────────────────────────────
subset_defs <- list(
  full = list(
    classes = c(1, 2, 3, 4),
    label_fn = function(zc) ifelse(zc %in% c(1, 2), "Clockwise", "Anticlockwise")
  ),
  pure = list(
    classes = c(1, 4),
    label_fn = function(zc) ifelse(zc == 1, "Clockwise", "Anticlockwise")
  ),
  fig8 = list(
    classes = c(2, 3),
    label_fn = function(zc) ifelse(zc == 2, "Clockwise", "Anticlockwise")
  )
)

# ── 3. Tuned hyperparameters from Table S5 ───────────────────────────────────
tuned_params <- list(
  full = list(
    trees          = 954,
    tree_depth     = 4,
    min_n          = 4,
    learn_rate     = 0.00148,
    loss_reduction = 8.90e-04,
    mtry           = 0.70,
    sample_size    = 0.737
  ),
  pure = list(
    trees          = 771,
    tree_depth     = 5,
    min_n          = 1,
    learn_rate     = 0.0310,
    loss_reduction = 1.60e-03,
    mtry           = 0.70,
    sample_size    = 0.610
  ),
  fig8 = list(
    trees          = 252,
    tree_depth     = 4,
    min_n          = 18,
    learn_rate     = 0.00255,
    loss_reduction = 2.18e-05,
    mtry           = 0.70,
    sample_size    = 0.746
  )
)

# ── 4. Build finalised workflow for a given subset ───────────────────────────
build_workflow <- function(dat, hp) {
  rec <- recipe(direction ~ ., data = dat) |>
    update_role(Site, new_role = "ID") |>
    step_zv(all_predictors()) |>
    step_nzv(all_predictors())

  spec <- boost_tree(
    trees          = hp$trees,
    tree_depth     = hp$tree_depth,
    learn_rate     = hp$learn_rate,
    sample_size    = hp$sample_size,
    mtry           = hp$mtry,
    min_n          = hp$min_n,
    loss_reduction = hp$loss_reduction
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

  workflow() |>
    add_recipe(rec) |>
    add_model(spec)
}

# ── 5. Prepare data for a given subset ───────────────────────────────────────
prepare_subset <- function(raw_data, subset_name) {
  sdef <- subset_defs[[subset_name]]
  raw_data |>
    dplyr::filter(Zuec_class %in% sdef$classes) |>
    dplyr::mutate(
      direction = factor(
        sdef$label_fn(Zuec_class),
        levels = c("Clockwise", "Anticlockwise")
      )
    ) |>
    dplyr::select(all_of(c("Site", all_features, "direction"))) |>
    tidyr::drop_na()
}

# ── 6. Label perturbation function ───────────────────────────────────────────
perturb_labels <- function(data, flip_rate) {
  if (flip_rate == 0) return(data)
  n <- nrow(data)
  n_flip <- round(n * flip_rate)
  flip_idx <- sample(seq_len(n), n_flip, replace = FALSE)
  data$direction[flip_idx] <- factor(
    ifelse(data$direction[flip_idx] == "Clockwise",
           "Anticlockwise", "Clockwise"),
    levels = c("Clockwise", "Anticlockwise")
  )
  data
}

# ── 7. Run one LOSO evaluation ───────────────────────────────────────────────
run_one_loso <- function(dat, wf) {
  loso_splits <- group_vfold_cv(dat, group = Site)
  loso_results <- fit_resamples(
    wf,
    resamples = loso_splits,
    metrics   = metric_set(roc_auc, accuracy, f_meas),
    control   = control_resamples(save_pred = FALSE)
  )
  collect_metrics(loso_results) |>
    dplyr::select(.metric, mean) |>
    tidyr::pivot_wider(names_from = .metric, values_from = mean)
}

# ── 8. Main sensitivity loop ────────────────────────────────────────────────
flip_rates   <- c(0, 0.10, 0.20, 0.33)
n_iter       <- 200
subset_names <- c("full", "pure", "fig8")

results <- list()
idx <- 0

total_runs <- length(subset_names) * (1 + sum(flip_rates > 0) * n_iter)
cat(sprintf("Starting sensitivity analysis: %d LOSO evaluations\n", total_runs))
cat(sprintf("(3 subsets x [1 baseline + %d perturbation rates x %d iterations])\n\n",
            sum(flip_rates > 0), n_iter))

t_start <- Sys.time()

for (sub in subset_names) {

  cat(sprintf("── Subset: %s ──\n", sub))

  dat_base <- prepare_subset(raw, sub)
  hp <- tuned_params[[sub]]

  cat(sprintf("  N = %d events\n", nrow(dat_base)))

  for (fr in flip_rates) {

    n_runs <- ifelse(fr == 0, 1, n_iter)
    cat(sprintf("  Flip rate: %.0f%% (%d runs)\n", fr * 100, n_runs))

    for (iter in seq_len(n_runs)) {

      idx <- idx + 1

      dat_perturbed <- perturb_labels(dat_base, fr)
      wf_iter <- build_workflow(dat_perturbed, hp)

      metrics <- tryCatch(
        run_one_loso(dat_perturbed, wf_iter),
        error = function(e) {
          warning(sprintf("Error in %s, rate=%.2f, iter=%d: %s",
                          sub, fr, iter, e$message))
          tibble(roc_auc = NA_real_, accuracy = NA_real_, f_meas = NA_real_)
        }
      )

      results[[idx]] <- tibble(
        subset    = sub,
        flip_rate = fr,
        iteration = iter,
        n_events  = nrow(dat_perturbed),
        roc_auc   = metrics$roc_auc,
        accuracy  = metrics$accuracy,
        f1        = metrics$f_meas
      )

      if (idx %% 25 == 0) {
        elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
        cat(sprintf("    [%5.1f min] %d / %d runs complete\n",
                    elapsed, idx, total_runs))
      }
    }
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nCompleted in %.1f minutes\n\n", elapsed_total))

# ── 9. Combine and summarise ────────────────────────────────────────────────
results_df <- bind_rows(results)
write.csv(results_df, "results/sensitivity_results.csv", row.names = FALSE)

summary_df <- results_df |>
  group_by(subset, flip_rate) |>
  summarise(
    n_iter     = n(),
    auc_median = median(roc_auc, na.rm = TRUE),
    auc_q05    = quantile(roc_auc, 0.05, na.rm = TRUE),
    auc_q95    = quantile(roc_auc, 0.95, na.rm = TRUE),
    acc_median = median(accuracy, na.rm = TRUE),
    f1_median  = median(f1, na.rm = TRUE),
    .groups    = "drop"
  )

cat("══ Summary ══\n")
print(summary_df, width = Inf)

write.csv(summary_df, "results/sensitivity_summary.csv", row.names = FALSE)