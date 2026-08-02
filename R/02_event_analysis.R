# ══════════════════════════════════════════════════════════════════════════════
# 02_event_analysis.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Read individual event CSVs produced by 01_event_definition.R, calculate
#   hydrograph metrics, hysteresis indices, and antecedent condition indices
#   to produce the final analytical dataset.
#
# Workflow:
#   1. Read and combine event CSVs across all sites
#   2. Calculate event-scale hydrological metrics (rising time, recession,
#      symmetry, flashiness, baseflow index, recession constant, volume,
#      load, yield)
#   3. Calculate hysteresis indices (Lloyd HI, Flushing Index, Zuecco
#      classification)
#   4. Calculate antecedent condition indices (time since previous event,
#      preceding event magnitude, event frequency within 7/14/30-day windows)
#   5. Merge and export the final dataset
#
# Inputs:
#   - Event CSVs from 01_event_definition.R (one per event per site)
#   - functions.R (shared utility functions including Zuecco classification)
#
# Outputs:
#   - data/Final_Hysteresis_Dataset.csv (without antecedent indices)
#   - data/Final_Hysteresis_Dataset_Ante.csv (with antecedent indices)
#
# Dependencies:
#   tidyverse, lubridate, purrr
#
# Reference:
#   Dare et al. (2026). Event-scale prediction of sediment hysteresis
#   regimes using hydrological and spatial indices. Journal of Hydrology.
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(lubridate)
library(purrr)

source("R/functions.R")

# =========================================================================
# 1. Read and combine events across all sites
# =========================================================================
# Site configuration: name, event folder, and catchment area (ha)

site_config <- tibble::tribble(
  ~site_name, ~event_folder,       ~area_ha,
  "Tuangi",   "output/events/Tuangi_Events",   9461.728,
  "Kokako",   "output/events/Kokako_Events",   1132.848,
  "Kokopu",   "output/events/Kokopu_Events",   5378.934,
  "Kotuku",   "output/events/Kotuku_Events",   3246.450,
  "Matata",   "output/events/Matata_Events",   9976.359,
  "Pipi",     "output/events/Pipi_Events",     1319.197,
  "Tuna",     "output/events/Tuna_Events",     4182.178,
  "Inanga",   "output/events/Inanga_Events",   3317.304
)

cat("Reading events from all sites...\n")

All_Events <- pmap_dfr(site_config, function(site_name, event_folder, area_ha) {
  df <- Event_Reader(name = event_folder)
  df$Site    <- site_name
  df$Area_ha <- area_ha
  df
})

All_Events <- All_Events %>%
  mutate(Time = parse_date_time(Time, orders = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d"), tz = "Etc/GMT-12"))

cat("Total events loaded:", length(unique(paste(All_Events$Site, All_Events$Final_Event))), "\n")

# =========================================================================
# 2. Calculate event-scale hydrological metrics
# =========================================================================
cat("Calculating hydrograph metrics...\n")

# Parameters for peak counting
min_peak_frac <- 0.80   # ignore peaks below this fraction of Q_peak
min_peak_sep  <- 100     # minimum separation between peaks (minutes)
recession_fit_hours <- 3 # hours of falling limb used to fit recession k

event_metrics <- All_Events %>%
  group_by(Site, Final_Event) %>%
  mutate(
    dt_s      = as.numeric(difftime(lead(Time), Time, units = "secs")),
    dV        = DISC * dt_s,
    dLoad_kg  = .pred * DISC * dt_s * 1e-3
  ) %>%
  summarise(
    Source_File       = unique(Source_File),
    t0                = min(Time, na.rm = TRUE),
    t1                = max(Time, na.rm = TRUE),
    i_peak            = which.max(DISC),
    t_peak            = Time[i_peak],
    Q_peak            = DISC[i_peak],
    rising_time_h     = as.numeric(difftime(t_peak, t0, units = "hours")),
    recession_time_h  = as.numeric(difftime(t1, t_peak, units = "hours")),
    symmetry_index    = rising_time_h / (rising_time_h + recession_time_h),

    RB_flashiness = {
      q   <- DISC
      num <- sum(abs(q - dplyr::lag(q)), na.rm = TRUE)
      den <- sum(q, na.rm = TRUE)
      ifelse(den > 0, num / den, NA_real_)
    },

    peak_sharpness = {
      q         <- DISC
      mean_rise <- mean(q[Time <= t_peak], na.rm = TRUE)
      ifelse(is.finite(mean_rise) && mean_rise > 0, Q_peak / mean_rise, NA_real_)
    },

    baseflow_index = {
      num <- sum(baseflow, na.rm = TRUE)
      den <- sum(DISC, na.rm = TRUE)
      ifelse(den > 0, num / den, NA_real_)
    },

    recession_k_per_h = {
      fall_idx <- which(Time >= t_peak & Time <= (t_peak + hours(recession_fit_hours)))
      q_fall   <- DISC[fall_idx]
      t_fall   <- Time[fall_idx]
      ok       <- is.finite(q_fall) & (q_fall > 0) & is.finite(t_fall)
      q_fall   <- q_fall[ok]
      t_fall   <- t_fall[ok]

      if (length(q_fall) >= 3) {
        tt  <- as.numeric(difftime(t_fall, t_fall[1], units = "hours"))
        fit <- try(lm(log(q_fall) ~ tt), silent = TRUE)
        if (inherits(fit, "try-error")) NA_real_ else as.numeric(-coef(fit)[["tt"]])
      } else {
        NA_real_
      }
    },

    n_peaks    = count_peaks(Time, DISC, min_frac = min_peak_frac, min_sep_min = min_peak_sep),
    VQ_m3      = sum(dV, na.rm = TRUE),
    vQ_mm      = VQ_m3 / (unique(Area_ha) * 10),
    Load_kg    = sum(dLoad_kg, na.rm = TRUE),
    Yield_kg_ha = Load_kg / unique(Area_ha),
    duration_h = as.numeric(difftime(t1, t0, units = "hours")),
    n_obs      = dplyr::n(),
    mean_dt_s  = mean(dt_s, na.rm = TRUE),

    .groups = "drop"
  )

cat("Hydrograph metrics calculated for", nrow(event_metrics), "events\n")

# =========================================================================
# 3. Calculate hysteresis indices
# =========================================================================
cat("Calculating hysteresis indices...\n")

hysteresis_output <- Calc_Hysteresis_Indices(data = All_Events)

# Merge with event metrics
event_metrics <- event_metrics %>% mutate(Unique = paste0(Site, Final_Event))
hysteresis_output <- hysteresis_output %>% mutate(Unique = paste0(Site, Final_Event))

Final_Dataset <- merge(event_metrics, hysteresis_output, by = "Unique") %>%
  select(-Unique, -Site.y, -Final_Event.y) %>%
  rename(
    Site        = Site.x,
    Final_Event = Final_Event.x,
    Q_fixed_min = Q_fixed_evt_min,
    Q_fixed_max = Q_fixed_evt_max
  )

write.csv(Final_Dataset, file = "data/Final_Hysteresis_Dataset.csv", row.names = FALSE)
cat("Saved: data/Final_Hysteresis_Dataset.csv\n")

# =========================================================================
# 4. Calculate antecedent condition indices
# =========================================================================
cat("Calculating antecedent indices...\n")

Hys_dat <- Final_Dataset %>%
  mutate(
    t0 = if_else(
      str_detect(as.character(t0), "^\\d{4}-\\d{2}-\\d{2}$"),
      paste0(t0, " 00:00:00"),
      as.character(t0)
    ),
    t1 = if_else(
      str_detect(as.character(t1), "^\\d{4}-\\d{2}-\\d{2}$"),
      paste0(t1, " 00:00:00"),
      as.character(t1)
    ),
    t0 = ymd_hms(t0, tz = "Etc/GMT-12"),
    t1 = ymd_hms(t1, tz = "Etc/GMT-12")
  )

Hys_dat_ante <- Hys_dat %>%
  group_by(Site) %>%
  arrange(t0, .by_group = TRUE) %>%
  mutate(
    # Time since previous event
    t1_prev          = lag(t1),
    dt_prev_event_h  = as.numeric(difftime(t0, t1_prev, units = "hours")),
    dt_prev_event_d  = dt_prev_event_h / 24,

    # Previous event magnitude
    Q_peak_prev      = lag(Q_peak),
    Load_kg_prev     = lag(Load_kg),
    VQ_m3_prev       = lag(VQ_m3),

    # Overlap flag
    prev_overlap_flag = dt_prev_event_h < 0
  ) %>%
  mutate(
    # Event counts within antecedent windows
    n_events_prev_30d = map_int(
      row_number(),
      ~ sum(t1 < t0[.x] & t1 >= (t0[.x] - days(30)), na.rm = TRUE)
    ),
    n_events_prev_14d = map_int(
      row_number(),
      ~ sum(t1 < t0[.x] & t1 >= (t0[.x] - days(14)), na.rm = TRUE)
    ),
    n_events_prev_7d = map_int(
      row_number(),
      ~ sum(t1 < t0[.x] & t1 >= (t0[.x] - days(7)), na.rm = TRUE)
    )
  ) %>%
  ungroup()

# Diagnostic checks
cat("Ordering check (should all be FALSE):\n")
print(Hys_dat_ante %>% group_by(Site) %>% summarise(any_unsorted = any(diff(t0) < 0)))

cat("Overlap count:", sum(Hys_dat_ante$prev_overlap_flag, na.rm = TRUE), "\n")

write.csv(Hys_dat_ante, file = "data/Final_Hysteresis_Dataset_Ante.csv", row.names = FALSE)
cat("Saved: data/Final_Hysteresis_Dataset_Ante.csv\n")
cat("Event analysis complete.\n")
