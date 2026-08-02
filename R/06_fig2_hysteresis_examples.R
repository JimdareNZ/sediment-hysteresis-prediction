# ══════════════════════════════════════════════════════════════════════════════
# 06_fig2_hysteresis_examples.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 2: Representative hysteresis loop examples for four events spanning
# contrasting regimes. Left panels: hydrograph + sedigraph. Right panels:
# normalised Q-TSS hysteresis loop.
#
# Inputs: Individual event CSVs from 01_event_definition.R
# Output: figures/Fig2_Hysteresis_Examples.png
# ══════════════════════════════════════════════════════════════════════════════

library(ggplot2)
library(tidyverse)
library(cowplot)
library(scales)

# -------------------------
# Load representative events
# -------------------------
kokako <- read_csv("output/events/Kokako_Events/Event_84_E.csv") %>%
  mutate(Site = "MGT", Direction = "Clockwise")
pipi <- read_csv("output/events/Pipi_Events/Event_4.csv") %>%
  mutate(Site = "PUN", Direction = "Clockwise") %>%
  slice(-(1:16))
kotuku <- read_csv("output/events/Kotuku_Events/Event_19.csv") %>%
  mutate(Site = "WHR-D", Direction = "Anticlockwise")
inanga <- read_csv("output/events/Inanga_Events/Event_8.csv") %>%
  mutate(Site = "PKO-UW", Direction = "Anticlockwise")

# -------------------------
# Colours and theme
# -------------------------
cw_col  <- "#d73027"
acw_col <- "#4575b4"

theme_hysteresis <- function() {
  theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      axis.title       = element_text(size = 8),
      axis.text        = element_text(size = 7),
      plot.title       = element_text(size = 9, face = "bold", hjust = 0),
      plot.margin      = margin(4, 4, 4, 4, "mm"),
      plot.background  = element_rect(fill = "white", colour = NA),
      legend.position  = "none"
    )
}

# -------------------------
# Helper functions
# -------------------------
normalise <- function(x) (x - min(x, na.rm = TRUE)) /
  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

get_limb <- function(disc) {
  peak_i <- which.max(disc)
  limb <- rep("Falling", length(disc))
  limb[1:peak_i] <- "Rising"
  limb
}

make_loop_panel <- function(df, panel_title, direction) {
  col <- if (direction == "Clockwise") cw_col else acw_col
  df <- df %>%
    mutate(Q_norm = normalise(DISC), TSS_norm = normalise(.pred),
           Limb = get_limb(DISC), Time = as.POSIXct(Time))
  rising  <- df %>% filter(Limb == "Rising")
  falling <- df %>% filter(Limb == "Falling")
  rise_mid  <- rising[round(nrow(rising) / 2), ]
  fall_mid  <- falling[round(nrow(falling) / 2), ]
  rise_next <- rising[round(nrow(rising) / 2) + 1, ]
  fall_next <- falling[round(nrow(falling) / 2) + 1, ]

  ggplot(df, aes(x = Q_norm, y = TSS_norm)) +
    geom_path(colour = col, linewidth = 0.8, alpha = 0.9) +
    geom_segment(data = rise_mid,
      aes(x = Q_norm, y = TSS_norm, xend = rise_next$Q_norm, yend = rise_next$TSS_norm),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = col, linewidth = 0.8) +
    geom_segment(data = fall_mid,
      aes(x = Q_norm, y = TSS_norm, xend = fall_next$Q_norm, yend = fall_next$TSS_norm),
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = col, linewidth = 0.8) +
    geom_point(data = df[1, ], aes(x = Q_norm, y = TSS_norm),
               shape = 21, size = 3, fill = "white", colour = col, stroke = 1) +
    labs(x = "Normalised Q", y = "Normalised TSS", title = panel_title) +
    theme_hysteresis()
}

make_hydro_panel <- function(df, direction) {
  col <- if (direction == "Clockwise") cw_col else acw_col
  df <- df %>% mutate(Time = as.POSIXct(Time))
  q_max <- max(df$DISC, na.rm = TRUE)
  tss_max <- max(df$.pred, na.rm = TRUE)
  scale_factor <- q_max / tss_max

  ggplot(df, aes(x = Time)) +
    geom_area(aes(y = DISC), fill = "#4393c3", alpha = 0.3) +
    geom_line(aes(y = DISC), colour = "#4393c3", linewidth = 0.7) +
    geom_line(aes(y = .pred * scale_factor), colour = col, linewidth = 0.7, lty = 2) +
    scale_y_continuous(
      name = expression(Q ~ (m^3 ~ s^-1)), limits = c(0, NA),
      sec.axis = sec_axis(~ . / scale_factor, name = expression(TSS ~ (mg ~ L^-1)))
    ) +
    scale_x_datetime(date_labels = "%d %b\n%H:%M", date_breaks = "6 hours") +
    labs(x = NULL) +
    theme_hysteresis() +
    theme(axis.title.y.left = element_text(colour = "#4393c3", size = 7),
          axis.title.y.right = element_text(colour = col, size = 7),
          axis.text.x = element_text(size = 6))
}

# -------------------------
# Build and assemble
# -------------------------
final_figure <- plot_grid(
  make_hydro_panel(kokako, "Clockwise"),  make_loop_panel(kokako, "(a) MGT \u2014 Clockwise", "Clockwise"),
  make_hydro_panel(pipi, "Clockwise"),    make_loop_panel(pipi, "(b) PUN \u2014 Clockwise", "Clockwise"),
  make_hydro_panel(kotuku, "Anticlockwise"), make_loop_panel(kotuku, "(c) WHR-D \u2014 Anticlockwise", "Anticlockwise"),
  make_hydro_panel(inanga, "Anticlockwise"), make_loop_panel(inanga, "(d) PKO-UW \u2014 Anticlockwise", "Anticlockwise"),
  ncol = 2, rel_widths = c(1.4, 1), align = "hv"
)

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Fig2_Hysteresis_Examples.png", final_figure,
       dpi = 300, width = 10, height = 14, units = "in", bg = "white")
