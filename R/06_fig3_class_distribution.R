# ══════════════════════════════════════════════════════════════════════════════
# 06_fig3_class_distribution.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 3: Percentage distribution of hysteresis classes across all sites,
# weighted by event count (left panel) and event sediment load (right panel).
# The load-weighted panel was added in response to Reviewer 1, Comment 1.
#
# Input:  data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv
# Output: figures/Fig3_Hysteresis_Class_Distribution.png
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)

# ── 1. Load data ─────────────────────────────────────────────────────────────
raw <- read.csv("data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv",
                stringsAsFactors = FALSE)

# ── 2. Site mapping (field names → manuscript site codes) ────────────────────
site_lookup <- tribble(
  ~Site,     ~site_code, ~sub_catchment,
  "Tuangi",  "PKO-D",    "Pokopoko",
  "Inanga",  "PKO-UW",   "Pokopoko",
  "Kokopu",  "PKO-UE",   "Pokopoko",
  "Kokako",  "MGT",      "Mangatoetoe",
  "Pipi",    "PUN",      "Puanene",
  "Kotuku",  "WHR-D",    "Wharere",
  "Matata",  "PNG-M",    "Pongakawa",
  "Tuna",    "PNG-TR1",  "Pongakawa"
)

site_order <- c("PKO-D", "PKO-UW", "PKO-UE",
                "PNG-M", "PNG-TR1",
                "MGT", "PUN", "WHR-D")

# ── 3. Prepare data ─────────────────────────────────────────────────────────
dat <- raw |>
  left_join(site_lookup, by = "Site") |>
  filter(!is.na(Zuec_class), !is.na(site_code)) |>
  mutate(
    class_label = factor(paste0("Class ", as.roman(Zuec_class)),
                         levels = paste0("Class ", as.roman(1:8))),
    site_code   = factor(site_code, levels = site_order)
  )

# ── 4. Count-weighted and load-weighted proportions ─────────────────────────
count_props <- dat |>
  group_by(site_code, class_label) |>
  summarise(value = n(), .groups = "drop_last") |>
  mutate(total = sum(value), pct = 100 * value / total,
         weight_type = "Event count") |>
  ungroup()

load_props <- dat |>
  group_by(site_code, class_label) |>
  summarise(value = sum(Load_kg, na.rm = TRUE), .groups = "drop_last") |>
  mutate(total = sum(value), pct = 100 * value / total,
         weight_type = "Sediment load") |>
  ungroup()

panel_levels <- c("Event count", "Sediment load")

combined <- bind_rows(count_props, load_props) |>
  mutate(weight_type = factor(weight_type, levels = panel_levels))

# ── 5. Bar labels: n (count panel) and total load in t (load panel) ──────────
bar_labels <- bind_rows(
  count_props |> distinct(site_code, total, weight_type) |>
    mutate(label = paste0("n = ", total)),
  load_props  |> distinct(site_code, total, weight_type) |>
    mutate(label = paste0(round(total / 1000), " t"))
) |>
  mutate(weight_type = factor(weight_type, levels = panel_levels))

# ── 6. Colour scheme: red → blue for enrichment I–IV, greys for dilution ─────
class_colours <- c(
  "Class I"    = "#D73027", "Class II"   = "#FC8D59",
  "Class III"  = "#91BFDB", "Class IV"   = "#4575B4",
  "Class V"    = "#BDBDBD", "Class VI"   = "#969696",
  "Class VII"  = "#737373", "Class VIII" = "#525252"
)

# ── 7. Plot ──────────────────────────────────────────────────────────────────
fig3 <- ggplot(combined, aes(x = site_code, y = pct, fill = class_label)) +
  geom_col(position = "stack", width = 0.8, colour = "white", linewidth = 0.2) +
  geom_text(data = bar_labels,
            aes(x = site_code, y = 103, label = label),
            inherit.aes = FALSE, size = 2.8, vjust = 0) +
  facet_wrap(~ weight_type) +
  scale_fill_manual(values = class_colours, name = "Hysteresis class") +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    legend.position    = "bottom",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8),
    strip.text         = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    plot.background    = element_rect(fill = "white", colour = NA)
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Fig3_Hysteresis_Class_Distribution.png", fig3,
       width = 8, height = 4.5, dpi = 300, bg = "white")

# ── 8. Summary table of load-weighted proportions (reported in Section 4.3) ──
load_props |>
  select(site_code, class_label, pct) |>
  pivot_wider(names_from = class_label, values_from = pct, values_fill = 0) |>
  print(width = Inf)