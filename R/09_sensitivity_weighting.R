# ══════════════════════════════════════════════════════════════════════════════
# 09_sensitivity_weighting.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Weighting sensitivity analysis (Reviewer 2, Comment 3). Tests whether
#   XGBoost LOSO results are robust to alternative weighting schemes for CCI
#   (Composite Connectivity Index) and LUSP (Land Use Sediment Potential).
#
# Approach:
#   - Back-derive overland flow (O) and lateral flow (L) components
#     from stored QCI, DDI, CCI values
#   - Recompute CCI and/or LUSP with alternative weights
#   - Rerun the XGBoost LOSO pipeline for each scenario
#   - Compare ROC AUC and variable importance across scenarios
#
# Inputs:
#   - data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv
#
# Outputs:
#   - results/weighting_sensitivity_results.csv
#   - results/weighting_sensitivity_importance.csv
#   Figure S6 is produced from these outputs by 10_SI_sensitivity_figures.R.
#
# Notes:
#   - Runtime: ~7 minutes (7 scenarios x 3 subsets).
#   - LUSP is recomputed and z-scored across the eight monitored sites in
#     every scenario, including "Original". The stored LUSP in the dataset is
#     z-scored across ten sub-catchments, so the "Original" scenario here is
#     the internal reference for the weighting comparison and differs slightly
#     from the main-text model.
#   - The reported results were generated before a seed was fixed; the
#     archived CSVs in results/ reproduce the reported values exactly.
#     Rerunning overwrites the archived CSVs.
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Packages ──────────────────────────────────────────────────────────────
library(tidymodels)
library(xgboost)
library(vip)
library(scales)

tidymodels_prefer()

set.seed(2026)
dir.create("results", showWarnings = FALSE)

# ── 1. Load data ─────────────────────────────────────────────────────────────
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

# ── 2. Back-derive CCI flow pathway components ──────────────────────────────
# Original formula (polygon-level, preserved through area-weighted mean):
#   CCI = 0.5*O + 0.3*L + 0.2*(1 - DDI)
#   QCI = O + L
#   DDI = D
# Solve for L and O:
#   L = (0.5*QCI + 0.2 - 0.2*DDI - CCI) / 0.2
#   O = QCI - L

site_components <- raw |>
  dplyr::distinct(Site, QCI, DDI, CCI) |>
  dplyr::mutate(
    L_component = (0.5 * QCI + 0.2 - 0.2 * DDI - CCI) / 0.2,
    O_component = QCI - L_component
  )

# Verify back-derivation
site_components <- site_components |>
  dplyr::mutate(
    CCI_check = 0.5 * O_component + 0.3 * L_component + 0.2 * (1 - DDI)
  )

cat("── Back-derivation verification ──\n")
cat(sprintf("Max CCI error: %.2e\n",
            max(abs(site_components$CCI - site_components$CCI_check))))
print(site_components |> dplyr::select(Site, O_component, L_component, DDI, CCI, CCI_check))

# ── 3. Define weighting scenarios ────────────────────────────────────────────
# CCI weights: overland / lateral / deep drainage
cci_schemes <- list(
  original = c(O = 0.5, L = 0.3, D = 0.2),
  equal    = c(O = 1/3, L = 1/3, D = 1/3),
  reversed = c(O = 0.2, L = 0.3, D = 0.5)
)

# LUSP weights: high / moderate / low sediment production
lusp_schemes <- list(
  original = c(high = 3, mid = 2, low = 1),
  equal    = c(high = 1, mid = 1, low = 1),
  steep    = c(high = 5, mid = 2, low = 1)
)

# Scenario definitions: which CCI and LUSP scheme to use
scenarios <- tribble(
  ~scenario,               ~cci_scheme,  ~lusp_scheme,  ~description,
  "Original",              "original",   "original",    "Original weights (baseline)",
  "CCI equal",             "equal",      "original",    "CCI: equal pathway weights",
  "CCI reversed",          "reversed",   "original",    "CCI: deep drainage dominant",
  "LUSP equal",            "original",   "equal",       "LUSP: no LU differentiation",
  "LUSP steep",            "original",   "steep",       "LUSP: 5/2/1 high contrast",
  "Both equal",            "equal",      "equal",       "Both indices with equal weights",
  "Both adversarial",      "reversed",   "equal",       "CCI reversed + LUSP equal (worst case)"
)

# ── 4. Functions to recompute indices ────────────────────────────────────────
recompute_cci <- function(site_comps, weights) {
  # weights = c(O, L, D)
  # CCI = w_O * O + w_L * L + w_D * (1 - DDI)
  site_comps |>
    dplyr::mutate(
      CCI_new = weights["O"] * O_component +
                weights["L"] * L_component +
                weights["D"] * (1 - DDI)
    ) |>
    dplyr::select(Site, CCI_new)
}

recompute_lusp <- function(data, weights) {
  # weights = c(high, mid, low)
  # Compute raw weighted sum per site, then z-score across sites
  site_lu <- data |>
    dplyr::distinct(Site, LU_high, LU_mid, LU_low) |>
    dplyr::mutate(
      LUSP_raw = weights["high"] * LU_high +
                 weights["mid"]  * LU_mid  +
                 weights["low"]  * LU_low,
      LUSP_new = as.numeric(scale(LUSP_raw))
    ) |>
    dplyr::select(Site, LUSP_new)
  site_lu
}

apply_scenario <- function(raw_data, site_comps, cci_weights, lusp_weights) {
  dat <- raw_data

  # Recompute CCI
  cci_new <- recompute_cci(site_comps, cci_weights)
  dat <- dat |>
    dplyr::left_join(cci_new, by = "Site") |>
    dplyr::mutate(CCI = CCI_new) |>
    dplyr::select(-CCI_new)

  # Recompute LUSP
  lusp_new <- recompute_lusp(raw_data, lusp_weights)
  dat <- dat |>
    dplyr::left_join(lusp_new, by = "Site") |>
    dplyr::mutate(LUSP = LUSP_new) |>
    dplyr::select(-LUSP_new)

  dat
}

# ── 5. Class subset definitions ──────────────────────────────────────────────
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

# ── 6. Tuned hyperparameters from Table S5 ───────────────────────────────────
tuned_params <- list(
  full = list(
    trees = 954, tree_depth = 4, min_n = 4,
    learn_rate = 0.00148, loss_reduction = 8.90e-04,
    mtry = 0.70, sample_size = 0.737
  ),
  pure = list(
    trees = 771, tree_depth = 5, min_n = 1,
    learn_rate = 0.0310, loss_reduction = 1.60e-03,
    mtry = 0.70, sample_size = 0.610
  ),
  fig8 = list(
    trees = 252, tree_depth = 4, min_n = 18,
    learn_rate = 0.00255, loss_reduction = 2.18e-05,
    mtry = 0.70, sample_size = 0.746
  )
)

# ── 7. Workflow and evaluation functions ─────────────────────────────────────
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

prepare_subset <- function(data, subset_name) {
  sdef <- subset_defs[[subset_name]]
  data |>
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

run_loso_with_importance <- function(dat, wf) {
  # LOSO evaluation
  loso_splits <- group_vfold_cv(dat, group = Site)
  loso_results <- fit_resamples(
    wf,
    resamples = loso_splits,
    metrics   = metric_set(roc_auc, accuracy, f_meas),
    control   = control_resamples(save_pred = FALSE)
  )
  metrics <- collect_metrics(loso_results) |>
    dplyr::select(.metric, mean) |>
    tidyr::pivot_wider(names_from = .metric, values_from = mean)

  # Variable importance from full-data fit
  full_fit <- fit(wf, data = dat)
  imp <- full_fit |>
    extract_fit_parsnip() |>
    vip::vi() |>
    dplyr::rename(feature = Variable, importance = Importance)

  list(metrics = metrics, importance = imp)
}

# ── 8. Main loop ─────────────────────────────────────────────────────────────
subset_names <- c("full", "pure", "fig8")

results_metrics <- list()
results_importance <- list()
idx_m <- 0
idx_i <- 0

total_runs <- nrow(scenarios) * length(subset_names)
cat(sprintf("Starting weighting sensitivity: %d LOSO evaluations\n",
            total_runs))
cat(sprintf("(%d scenarios x %d subsets)\n\n",
            nrow(scenarios), length(subset_names)))

t_start <- Sys.time()

for (s in seq_len(nrow(scenarios))) {

  sc <- scenarios[s, ]
  cat(sprintf("── Scenario: %s ──\n", sc$scenario))

  # Apply alternative weights to the raw dataset
  dat_reweighted <- apply_scenario(
    raw,
    site_components,
    cci_schemes[[sc$cci_scheme]],
    lusp_schemes[[sc$lusp_scheme]]
  )

  for (sub in subset_names) {

    idx_m <- idx_m + 1
    dat_sub <- prepare_subset(dat_reweighted, sub)
    hp <- tuned_params[[sub]]
    wf <- build_workflow(dat_sub, hp)

    cat(sprintf("  %s (N = %d)...", sub, nrow(dat_sub)))

    res <- tryCatch(
      run_loso_with_importance(dat_sub, wf),
      error = function(e) {
        warning(sprintf("Error in %s/%s: %s", sc$scenario, sub, e$message))
        list(
          metrics = tibble(roc_auc = NA_real_, accuracy = NA_real_,
                           f_meas = NA_real_),
          importance = tibble(feature = character(), importance = numeric())
        )
      }
    )

    results_metrics[[idx_m]] <- tibble(
      scenario    = sc$scenario,
      cci_scheme  = sc$cci_scheme,
      lusp_scheme = sc$lusp_scheme,
      subset      = sub,
      n_events    = nrow(dat_sub),
      roc_auc     = res$metrics$roc_auc,
      accuracy    = res$metrics$accuracy,
      f1          = res$metrics$f_meas
    )

    # Store top 10 features for this scenario/subset
    top_imp <- res$importance |>
      dplyr::slice_head(n = 10) |>
      dplyr::mutate(
        scenario = sc$scenario,
        subset   = sub,
        rank     = dplyr::row_number()
      )
    results_importance <- c(results_importance, list(top_imp))

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    cat(sprintf(" done [%.1f min]\n", elapsed))
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cat(sprintf("\nCompleted in %.1f minutes\n\n", elapsed_total))

# ── 9. Combine and summarise ────────────────────────────────────────────────
metrics_df <- bind_rows(results_metrics)
importance_df <- bind_rows(results_importance)

cat("══ LOSO ROC AUC by scenario and subset ══\n")
metrics_wide <- metrics_df |>
  dplyr::select(scenario, subset, roc_auc) |>
  tidyr::pivot_wider(names_from = subset, values_from = roc_auc)
print(metrics_wide)

write.csv(metrics_df, "results/weighting_sensitivity_results.csv", row.names = FALSE)
write.csv(importance_df, "results/weighting_sensitivity_importance.csv", row.names = FALSE)