# ══════════════════════════════════════════════════════════════════════════════
# functions.R
# ══════════════════════════════════════════════════════════════════════════════
# Shared utility functions for the hysteresis prediction pipeline.
#
# Contents:
#   Event definition:
#     Event_Modifier()           — Interactive event boundary review
#     Event_Reader()             — Read and combine event CSVs
#
#   Hydrograph metrics:
#     count_peaks()              — Count peaks in discharge time series
#
#   Hysteresis classification:
#     Zuecco_Classification()    — Zuecco et al. (2016) 8-class scheme
#     qfixed_for_event()         — Determine valid Q_fixed grid
#     Hyst_Indicies()            — Lloyd HI and flushing index
#     Calc_Hysteresis_Indices()  — Wrapper to classify all events
#
# References:
#   Zuecco, G., Penna, D., Borga, M., & van Meerveld, H.J. (2016).
#     A versatile index to characterize hysteresis between hydrological
#     variables at the runoff event timescale. Hydrological Processes,
#     30, 1449-1466. https://doi.org/10.1002/hyp.10681
#
#   Lloyd, C.E.M., Freer, J.E., Johnes, P.J., & Collins, A.L. (2016).
#     Testing an improved index for analysing storm discharge-concentration
#     hysteresis. Hydrology and Earth System Sciences, 20, 625-632.
#     https://doi.org/10.5194/hess-20-625-2016
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(dplyr)
library(lubridate)
library(purrr)


# ══════════════════════════════════════════════════════════════════════════════
# EVENT DEFINITION
# ══════════════════════════════════════════════════════════════════════════════

#' Interactive event boundary review and adjustment
#'
#' Displays each identified event as a hydrograph plot and prompts the user
#' to accept, adjust, or delete the event. Accepted and adjusted events are
#' saved as individual CSV and PNG files.
#'
#' @param data A dataframe containing columns: Time, DISC, baseflow, Event_id
#' @param start_event Integer. The first event ID to review (default 1).
#' @param folder_name Character. Output directory for event CSVs and PNGs.
#' @return NULL (side effects: writes CSV and PNG files to folder_name)

Event_Modifier <- function(data, start_event = 1, folder_name = NA) {

  show_plot <- function(p) {
    print(p)
    utils::flush.console()
    try(grDevices::dev.flush(), silent = TRUE)
    Sys.sleep(0.15)
  }

  ask <- function(prompt) {
    tolower(trimws(readline(prompt)))
  }

  if (!is.na(folder_name) && !dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }

  last_event <- max(data$Event_id, na.rm = TRUE)
  tz_used    <- "Etc/GMT+12"

  for (i in seq.int(start_event, last_event)) {

    event_rows  <- dplyr::filter(data, Event_id == i)
    event_start <- min(event_rows$Time)
    event_end   <- max(event_rows$Time)

    initial_df <- data %>%
      dplyr::filter(Time >= event_start, Time <= event_end)

    initial_plot <- ggplot(initial_df) +
      geom_path(aes(x = Time, y = baseflow), lty = 2) +
      geom_point(aes(x = Time, y = DISC), colour = "midnightblue") +
      theme_bw() +
      ggtitle(paste0("Event ", i, " — initial"))

    show_plot(initial_plot)
    happy <- ask("Are you happy with this event? (Y/N or 'delete'): ")

    if (happy %in% c("y", "yes")) {
      message("Event ", i, ": accepted")
      initial_df$Final_Event <- i
      if (!is.na(folder_name)) {
        write.csv(initial_df,
                  file = file.path(folder_name, paste0("Event_", i, ".csv")),
                  row.names = FALSE)
        ggsave(filename = file.path(folder_name, paste0("Event_", i, ".png")),
               plot = initial_plot, device = "png", width = 10, height = 6)
      }

    } else if (happy == "delete") {
      message("Event ", i, ": deleted")
      next

    } else if (happy %in% c("n", "no")) {
      message("Event ", i, ": revising")
      repeat {
        start_offset <- suppressWarnings(
          as.numeric(readline("Enter start offset (hours, e.g. -2): "))
        )
        end_offset <- suppressWarnings(
          as.numeric(readline("Enter end offset (hours, e.g. 3): "))
        )
        if (any(is.na(c(start_offset, end_offset)))) {
          message("Please enter numeric values.")
          next
        }

        revised_start <- event_start + lubridate::hours(start_offset)
        revised_end   <- event_end   + lubridate::hours(end_offset)

        revised_df <- data %>%
          dplyr::filter(Time >= revised_start, Time <= revised_end)

        revised_plot <- ggplot(revised_df) +
          geom_path(aes(x = Time, y = baseflow), lty = 2) +
          geom_point(aes(x = Time, y = DISC), colour = "midnightblue") +
          theme_bw() +
          ggtitle(paste0("Event ", i, " — revised (",
                         start_offset, "h, ", end_offset, "h)"))

        show_plot(revised_plot)
        happy <- ask("Are you happy with the revised plot? (Y/N or 'delete'): ")

        if (happy == "delete") {
          message("Event ", i, ": deleted")
          break
        }
        if (happy %in% c("y", "yes")) {
          revised_df$Final_Event <- i
          if (!is.na(folder_name)) {
            write.csv(revised_df,
                      file = file.path(folder_name, paste0("Event_", i, ".csv")),
                      row.names = FALSE)
            ggsave(filename = file.path(folder_name, paste0("Event_", i, ".png")),
                   plot = revised_plot, device = "png", width = 10, height = 6)
          }
          break
        }
      }
    } else {
      message("Event ", i, ": unrecognised response, skipping")
    }
  }
}


#' Read and combine event CSVs from a site folder
#'
#' Reads all CSV files in the specified folder, orders them chronologically,
#' and assigns sequential Final_Event numbers.
#'
#' @param name Character. Path to the folder containing event CSVs.
#' @return A dataframe with all events combined and a Final_Event column.

Event_Reader <- function(name) {
  files <- list.files(path = paste0(name, "/"), pattern = "\\.csv$", full.names = TRUE)

  Event_List <- purrr::map2(files, basename(files), function(f, fname) {
    df <- read.csv(f, check.names = FALSE)
    df$Source_File <- fname
    df
  })

  get_start_time <- function(df) {
    time_col <- names(df)[grep("time|date", names(df), ignore.case = TRUE)][1]
    if (is.na(time_col)) return(NA)
    suppressWarnings(min(as.POSIXct(df[[time_col]]), na.rm = TRUE))
  }

  start_times <- purrr::map_dbl(Event_List, ~ as.numeric(get_start_time(.x)))
  order_idx   <- order(start_times, na.last = NA)

  Event_List <- Event_List[order_idx]
  Event_List <- purrr::map2(Event_List, seq_along(Event_List), function(df, i) {
    df$Final_Event <- i
    df
  })

  dplyr::bind_rows(Event_List)
}


# ══════════════════════════════════════════════════════════════════════════════
# HYDROGRAPH METRICS
# ══════════════════════════════════════════════════════════════════════════════

#' Count peaks in a discharge time series
#'
#' Identifies local maxima above a minimum prominence threshold with
#' minimum temporal separation between counted peaks.
#'
#' @param time_vec POSIXct vector of timestamps
#' @param q_vec Numeric vector of discharge values
#' @param min_frac Minimum fraction of Q_peak for a peak to be counted
#' @param min_sep_min Minimum separation between peaks in minutes
#' @return Integer count of peaks

count_peaks <- function(time_vec, q_vec, min_frac = 0.2, min_sep_min = 60) {
  n <- length(q_vec)
  if (n < 3) return(0L)

  is_peak <- (q_vec > dplyr::lag(q_vec, default = -Inf)) &
    (q_vec >= dplyr::lead(q_vec, default = -Inf))

  q_peak   <- max(q_vec, na.rm = TRUE)
  cand_idx <- which(is_peak & (q_vec >= min_frac * q_peak))
  if (length(cand_idx) == 0) return(0L)

  kept      <- integer(0)
  last_time <- as.POSIXct(NA)
  for (i in cand_idx) {
    if (is.na(last_time) || as.numeric(difftime(time_vec[i], last_time, units = "mins")) >= min_sep_min) {
      kept      <- c(kept, i)
      last_time <- time_vec[i]
    }
  }
  length(kept)
}


# ══════════════════════════════════════════════════════════════════════════════
# HYSTERESIS CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════════

#' Classify a single event into one of eight Zuecco hysteresis classes
#'
#' Normalises Q and y to [0,1], splits at the discharge peak, interpolates
#' rising and falling limb concentrations at fixed Q intervals, computes
#' trapezoid areas, and assigns a class (1-8) based on the sign pattern.
#'
#' Classes 1-4: enrichment (1=CW, 2=CW fig-8, 3=ACW fig-8, 4=ACW)
#' Classes 5-8: dilution (mirrors 1-4)
#' Class 0: linearity
#'
#' @param Q Numeric vector. Raw discharge time series.
#' @param y Numeric vector. Raw concentration time series.
#' @param Q_fixed Numeric vector. Normalised Q evaluation points.
#' @return List with diff_area, h, hyst_class

Zuecco_Classification <- function(Q, y, Q_fixed) {

  Qrng <- range(Q, na.rm = TRUE)
  Yrng <- range(y, na.rm = TRUE)

  if (!is.finite(Qrng[1]) || !is.finite(Qrng[2]) || diff(Qrng) == 0 ||
      !is.finite(Yrng[1]) || !is.finite(Yrng[2]) || diff(Yrng) == 0) {
    return(list(diff_area = numeric(0), h = 0, hyst_class = 0))
  }

  Qn <- (Q - Qrng[1]) / diff(Qrng)
  Yn <- (y - Yrng[1]) / diff(Yrng)

  i_peak <- which.max(Qn)
  if (i_peak == 1L || i_peak == length(Qn)) {
    return(list(diff_area = numeric(0), h = 0, hyst_class = 0))
  }

  rise <- data.frame(q = Qn[1:i_peak], c = Yn[1:i_peak])
  fall <- data.frame(q = Qn[i_peak:length(Qn)], c = Yn[i_peak:length(Yn)])

  rise <- rise[order(rise$q, rise$c), ]
  rise <- rise[!duplicated(rise$q), ]
  fall <- fall[order(fall$q, fall$c), ]
  fall <- fall[!duplicated(fall$q), ]

  if (nrow(rise) < 2 || nrow(fall) < 2) {
    return(list(diff_area = numeric(0), h = 0, hyst_class = 0))
  }

  q_lo <- max(min(rise$q), min(fall$q))
  q_hi <- min(max(rise$q), max(fall$q))
  if (!(q_hi > q_lo)) {
    return(list(diff_area = numeric(0), h = 0, hyst_class = 0))
  }

  qg <- sort(unique(Q_fixed))
  qg <- qg[qg >= q_lo & qg <= q_hi]
  if (length(qg) < 2) qg <- seq(q_lo, q_hi, length.out = 3L)

  yr <- approx(rise$q, rise$c, xout = qg, rule = 2, ties = "ordered")$y
  yf <- approx(fall$q, fall$c, xout = qg, rule = 2, ties = "ordered")$y

  dx <- diff(qg)
  Ar <- 0.5 * (head(yr, -1) + tail(yr, -1)) * dx
  Af <- 0.5 * (head(yf, -1) + tail(yf, -1)) * dx
  dA <- Ar - Af
  h  <- sum(dA)

  min_dA <- min(dA)
  max_dA <- max(dA)
  min_y_rise <- min(Yn[1:i_peak])
  max_y_rise <- max(Yn[1:i_peak])
  change_max <- abs(max_y_rise - Yn[1])
  change_min <- abs(min_y_rise - Yn[1])

  choose_enrichment <- function(h, min_dA, max_dA) {
    if (min_dA > 0 && max_dA > 0) return(1L)
    if (min_dA < 0 && max_dA < 0) return(4L)
    if (min_dA <= 0 && max_dA > 0 && h >= 0) return(2L)
    if (min_dA < 0 && max_dA >= 0 && h < 0) return(3L)
    return(0L)
  }

  choose_dilution <- function(h, min_dA, max_dA) {
    if (min_dA > 0 && max_dA > 0) return(5L)
    if (min_dA < 0 && max_dA < 0) return(8L)
    if (min_dA <= 0 && max_dA > 0 && h >= 0) return(6L)
    if (min_dA < 0 && max_dA >= 0 && h < 0) return(7L)
    return(0L)
  }

  hyst_class <- if (change_max > change_min) {
    choose_enrichment(h, min_dA, max_dA)
  } else if (change_max < change_min) {
    choose_dilution(h, min_dA, max_dA)
  } else {
    min_y_fall <- min(Yn[i_peak:length(Yn)])
    max_y_fall <- max(Yn[i_peak:length(Yn)])
    change_max_f <- abs(max_y_fall - Yn[1])
    change_min_f <- abs(min_y_fall - Yn[1])
    if (change_max_f > change_min_f) choose_enrichment(h, min_dA, max_dA)
    else if (change_max_f < change_min_f) choose_dilution(h, min_dA, max_dA)
    else 0L
  }

  list(diff_area = dA, h = h, hyst_class = hyst_class)
}


#' Determine the valid Q_fixed grid for a single event
#'
#' @param df Dataframe with Q_Norm and limb columns
#' @param Q_fixed Numeric vector of candidate normalised Q values
#' @return Numeric vector of valid Q_fixed values, or NULL

qfixed_for_event <- function(df, Q_fixed) {
  q_rise <- df$Q_Norm[df$limb == "Rising"]
  q_fall <- df$Q_Norm[df$limb == "Falling"]

  if (length(na.omit(q_rise)) < 2L || length(na.omit(q_fall)) < 2L) return(NULL)

  q_lo <- max(min(q_rise, na.rm = TRUE), min(q_fall, na.rm = TRUE))
  q_hi <- min(max(q_rise, na.rm = TRUE), max(q_fall, na.rm = TRUE))

  if (!is.finite(q_lo) || !is.finite(q_hi) || !(q_hi > q_lo)) return(NULL)

  qg <- sort(unique(Q_fixed))
  qg <- qg[qg >= q_lo & qg <= q_hi]
  if (length(qg) < 2L) qg <- seq(q_lo, q_hi, length.out = 3L)
  qg
}


#' Calculate Lloyd hysteresis index and flushing index
#'
#' @param data Dataframe with columns: Time, Q_Norm, TSS_Norm, limb
#' @param Int Integer. Interval resolution (1, 5, 10, or 25)
#' @return Dataframe with HI_Lloyd_Index and Flushing_Index

Hyst_Indicies <- function(data, Int) {

  require(reshape2)

  names(data) <- c("Time", "Flow", "Contaminant", "Limb")

  Int_Vect <- switch(as.character(Int),
    "1"  = seq(0.01, 0.99, by = 0.01),
    "5"  = seq(0.05, 0.95, by = 0.05),
    "10" = seq(0.10, 0.90, by = 0.10),
    "25" = seq(0.25, 0.75, by = 0.25),
    stop("Unsupported interval. Use 1, 5, 10, or 25.")
  )

  Final_Table <- NULL

  for (x in unique(data$Limb)) {
    Limb_data <- data[data$Limb == x, ]
    for (i in Int_Vect) {
      dist_from_val <- Limb_data$Flow - i
      dist_under <- dist_from_val
      dist_under[dist_under > 0] <- NA
      Under <- which.max(dist_under)
      dist_over <- dist_from_val
      dist_over[dist_over < 0] <- NA
      Over <- which.min(dist_over)

      if (length(Under) == 0 || length(Over) == 0) {
        Y <- NA
      } else {
        m1 <- lm(
          c(Limb_data$Contaminant[Under], Limb_data$Contaminant[Over]) ~
          c(Limb_data$Flow[Under], Limb_data$Flow[Over])
        )
        Y <- (coefficients(m1)[2] * i) + coefficients(m1)[1]
      }

      Final_Table <- rbind(Final_Table, data.frame(
        Limb = x, Interval = i,
        Under = ifelse(is.na(Y), NA, Limb_data$Flow[Under]),
        Over  = ifelse(is.na(Y), NA, Limb_data$Flow[Over]),
        Y     = Y
      ))
    }
  }

  Output <- dcast(Final_Table, Interval ~ Limb, value.var = "Y")
  Output$HI_Lloyd <- Output$Rising - Output$Falling

  HI_Lloyd_Index <- mean(Output$HI_Lloyd, na.rm = TRUE)
  Flushing_Index <- data[which.max(data[data$Limb == "Rising", "Flow"]), "Contaminant"] -
                    data[which.min(data[data$Limb == "Rising", "Flow"]), "Contaminant"]

  data.frame(HI_Lloyd_Index = HI_Lloyd_Index, Flushing_Index = Flushing_Index)
}


#' Calculate hysteresis indices for all events in a dataset
#'
#' @param data Dataframe with columns: Site, Final_Event, Time, DISC, .pred
#' @return Dataframe with one row per event containing hysteresis indices

Calc_Hysteresis_Indices <- function(data) {

  Final_DF <- NULL

  safe_Hyst_Indicies <- possibly(Hyst_Indicies,
    otherwise = tibble(HI_Lloyd_Index = NA_real_, Flushing_Index = NA_real_),
    quiet = TRUE
  )

  safe_Zuc_Indicies <- possibly(Zuecco_Classification,
    otherwise = list(diff_area = NA_real_, h = NA_real_, hyst_class = NA_real_),
    quiet = TRUE
  )

  Q_fixed <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.999)

  for (i in unique(data$Site)) {
    Event_Nos <- data %>% filter(Site == i) %>% select(Final_Event) %>% unique()

    for (x in Event_Nos$Final_Event) {

      Norm_DF <- data %>%
        filter(Site == i, Final_Event == x) %>%
        mutate(
          Q_Norm   = (DISC - min(DISC, na.rm = TRUE)) / (max(DISC, na.rm = TRUE) - min(DISC, na.rm = TRUE)),
          TSS_Norm = (.pred - min(.pred, na.rm = TRUE)) / (max(.pred, na.rm = TRUE) - min(.pred, na.rm = TRUE)),
          t_peak   = Time[which.max(DISC)],
          limb     = if_else(Time <= t_peak, "Rising", "Falling")
        ) %>%
        ungroup() %>%
        select(Time, Q_Norm, TSS_Norm, limb) %>%
        filter(complete.cases(.))

      HysI <- safe_Hyst_Indicies(data = Norm_DF, Int = 1)

      Q_fixed_evt <- tryCatch(qfixed_for_event(Norm_DF, Q_fixed), error = function(e) NULL)

      if (is.null(Q_fixed_evt)) {
        Zuec <- list(diff_area = numeric(0), h = NA, hyst_class = NA)
      } else {
        Zuec <- safe_Zuc_Indicies(Q = Norm_DF$Q_Norm, y = Norm_DF$TSS_Norm, Q_fixed = Q_fixed_evt)
        if (is.null(Zuec$hyst_class) || is.na(Zuec$hyst_class)) {
          Zuec$hyst_class <- NA
          if (is.null(Zuec$h) || is.na(Zuec$h)) Zuec$h <- NA
        }
      }

      Final_DF <- rbind(Final_DF, data.frame(
        Site           = i,
        Final_Event    = x,
        HI_Lloyd_Index = HysI$HI_Lloyd_Index,
        Flushing_Index = HysI$Flushing_Index,
        Zuec_h         = Zuec$h,
        Zuec_class     = Zuec$hyst_class,
        Q_fixed_evt_min = min(Q_fixed_evt, na.rm = TRUE),
        Q_fixed_evt_max = max(Q_fixed_evt, na.rm = TRUE)
      ))
    }
  }

  return(Final_DF)
}
