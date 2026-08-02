# ══════════════════════════════════════════════════════════════════════════════
# 01_event_definition.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Identify high-flow events from ANN-derived TSS and discharge time series
#   using the hydroEvents package, followed by interactive manual refinement.
#
# Workflow:
#   1. Load ANN-derived TSS and discharge data for a given site
#   2. Pad to regular 15-minute intervals and interpolate short gaps
#   3. Separate baseflow using Lyne-Hollick filter (alpha = 0.99)
#   4. Identify events using eventMinima on the quickflow component
#   5. Interactively review and refine event boundaries
#   6. Export individual event CSVs for subsequent analysis
#
# Inputs:
#   - ANN-derived TSS CSV file with columns: Time, DISC, TURB, .pred
#   - Site-specific eventMinima parameters (delta.y, delta.x, threshold)
#
# Outputs:
#   - One CSV per event in the site output folder
#   - One PNG per event showing the hydrograph with baseflow
#
# Notes:
#   - This script is run once per site. Change the site_name variable below
#     to process each site in turn.
#   - Event boundaries are refined interactively using Event_Modifier(),
#     which prompts the user to accept, adjust, or delete each event.
#   - The eventMinima parameters were tuned on a site-specific basis and
#     refined through visual inspection (see Methods, Section 3.5).
#
# Dependencies:
#   tidyverse, hydroEvents, zoo, padr
#
# Reference:
#   Dare et al. (2026). Event-scale prediction of sediment hysteresis
#   regimes using hydrological and spatial indices. Journal of Hydrology.
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(hydroEvents)
library(zoo)
library(padr)

source("R/functions.R")

# -------------------------
# Configuration
# -------------------------
# Set the site to process. Run this script once per site, changing the
# site_name and eventMinima parameters as needed.

site_name <- "Matata"  # Options: Matata, Tuangi, Inanga, Kokopu, Kokako, Pipi, Kotuku, Tuna

# Site-specific eventMinima parameters
# These control event separation sensitivity and were tuned per site
event_params <- list(
  Matata = list(delta.y = 0.2, delta.x = 300, threshold = 0.2),
  Tuangi = list(delta.y = 0.2, delta.x = 300, threshold = 0.2),
  Inanga = list(delta.y = 0.2, delta.x = 300, threshold = 0.2),
  Kokopu = list(delta.y = 0.2, delta.x = 300, threshold = 0.2),
  Kokako = list(delta.y = 0.2, delta.x = 300, threshold = 0.2),
  Pipi   = list(delta.y = 0.2, delta.x = 300, threshold = 0.2),
  Kotuku = list(delta.y = 0.2, delta.x = 300, threshold = 0.2),
  Tuna   = list(delta.y = 0.2, delta.x = 300, threshold = 0.2)
)

# -------------------------
# Paths
# -------------------------
input_file  <- file.path("data", "raw", paste0(site_name, "_TSS_Final.csv"))
output_dir  <- file.path("output", "events", paste0(site_name, "_Events"))

# -------------------------
# 1. Load and prepare data
# -------------------------
cat("Loading data for", site_name, "\n")

site_data <- read.csv(input_file) %>%
  mutate(
    Time = as.POSIXct(Time, format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT+12"),
    DISC = as.numeric(DISC)
  )

# Pad to regular 15-minute intervals and interpolate short gaps (max 4 hours)
site_data <- pad(site_data)
site_data$TURB  <- na.spline(site_data$TURB,  maxgap = 16)
site_data$DISC  <- na.spline(site_data$DISC,  maxgap = 16)
site_data$.pred <- na.spline(site_data$.pred, maxgap = 16)

# -------------------------
# 2. Baseflow separation
# -------------------------
# Lyne-Hollick recursive digital filter with alpha = 0.99
# Baseflow is subtracted to isolate the quickflow component for event detection

flow_data <- site_data %>%
  filter(!is.na(DISC)) %>%
  mutate(baseflow = baseflowA(DISC, alpha = 0.99)$bf) %>%
  select(Time, DISC, baseflow, .pred, TURB)

# -------------------------
# 3. Event identification
# -------------------------
# eventMinima identifies events based on local minima in the quickflow signal
# Parameters are site-specific (see event_params above)

params <- event_params[[site_name]]

event_results <- eventMinima(
  flow_data$DISC - flow_data$baseflow,
  delta.y   = params$delta.y,
  delta.x   = params$delta.x,
  threshold = params$threshold
)

# Visual check of identified events
plotEvents(
  data    = flow_data$DISC,
  events  = event_results,
  xlab    = "Index",
  ylab    = "Discharge (m3/s)",
  colpnt  = "#E41A1C",
  colline = "#377EB8",
  main    = paste0(site_name, " — eventMinima results")
)

# -------------------------
# 4. Tag events in dataframe
# -------------------------
flow_data$Event_id <- NA

for (i in 1:nrow(event_results)) {
  start_idx <- event_results$srt[i]
  end_idx   <- event_results$end[i]
  flow_data$Event_id[start_idx:end_idx] <- i
}

# -------------------------
# 5. Interactive event refinement
# -------------------------
# Event_Modifier displays each event and prompts the user to:
#   - Accept (Y) — saves the event as-is
#   - Adjust (N) — enter start/end offsets in hours to revise boundaries
#   - Delete — skip the event entirely
#
# Refined events are saved as individual CSVs in the output folder

cat("Starting interactive event refinement for", site_name, "\n")
cat("Output directory:", output_dir, "\n")

Event_Modifier(
  data         = flow_data,
  start_event  = 1,
  folder_name  = output_dir
)

cat("Event definition complete for", site_name, "\n")
