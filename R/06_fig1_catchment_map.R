# ══════════════════════════════════════════════════════════════════════════════
# 06_fig1_catchment_map.R
# ══════════════════════════════════════════════════════════════════════════════
# Figure 1: Two-column layout. Left: NZ context map + monitoring network.
# Right: Three spatial predictor panels (SEI, CCI, LUSP).
#
# Inputs: Shapefiles (not redistributed; see Data Availability Statement)
# Output: figures/Fig1_Catchment_Map.png
#
# Note: Shapefile paths must be updated to match your local directory.
# ══════════════════════════════════════════════════════════════════════════════

library(sf)
library(ggplot2)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(cowplot)
library(ggrepel)

# ── Load spatial data ─────────────────────────────────────────────────────────
# Update these paths to match your local shapefile locations
shp_base <- "data/shapefiles"

BOP_Region       <- read_sf(file.path(shp_base, "BOPRC/BOP Region - Land Area.shp")) %>% st_transform(4326)
Waihi_Estuary    <- read_sf(file.path(shp_base, "waihi_estuary/Waihi.shp")) %>% st_transform(4326)
monitoring_sites <- read_sf(file.path(shp_base, "waihi_estuary/Monitoring_Sites.shp")) %>% st_transform(4326)
Sub_Catchments   <- read_sf(file.path(shp_base, "waihi_estuary/WaihiEstuaryCatchmentsGenerated_MonitoringSitesSubset_JD_25052022.shp")) %>% st_transform(4326)
REC              <- read_sf(file.path(shp_base, "waihi_estuary/REC.shp")) %>% st_transform(4326) %>% st_intersection(Waihi_Estuary)

# ── Attribute sites and sub-catchments ────────────────────────────────────────
monitoring_sites <- monitoring_sites %>%
  mutate(MapName = case_when(
    SiteName == "Tuna" ~ "PNG-TR1", SiteName == "Kokako" ~ "MGT",
    SiteName == "Kokopu" ~ "PKO-UE", SiteName == "Inanga" ~ "PKO-UW",
    SiteName == "Tuangi" ~ "PKO-D", SiteName == "Matata" ~ "PNG-M",
    SiteName == "Pipi" ~ "PUN", SiteName == "Kotuku" ~ "WHR-D"
  )) %>%
  filter(!is.na(MapName))

Sub_Catchments <- Sub_Catchments %>%
  mutate(Sub_Catchment = case_when(
    Name %in% c("FN956338", "GN086350", "GN083710") ~ "Pokopoko",
    Name == "GN265752" ~ "Mangatoetoe",
    Name %in% c("GN849464", "GN955302", "HN063546") ~ "Pongakawa",
    Name == "GN562808" ~ "Puanene",
    Name %in% c("GN618644", "GN435054") ~ "Wharere",
    TRUE ~ NA_character_
  ))

predictor_indices <- tibble(
  Name = c("FN956338","GN083710","GN086350","GN265752","GN435054","GN562808","GN618644","GN849464","GN955302","HN063546"),
  SEI  = c(0.874, 1.024, 0.786, 1.162, -0.116, 1.729, -0.438, -0.097, -0.247, -1.667),
  CCI  = c(0.199, 0.184, 0.164, 0.150, 0.161, 0.164, 0.147, 0.155, 0.153, 0.161),
  LUSP = c(-0.642, -0.892, -0.706, 1.680, -0.270, 1.832, -0.029, -0.337, -0.920, 0.285)
)

Sub_Catchments <- Sub_Catchments %>%
  left_join(predictor_indices, by = "Name") %>%
  filter(!Name %in% c("HN063546", "GN435054")) %>%
  arrange(Name == "FN956338")

# ── Map settings ──────────────────────────────────────────────────────────────
xlims <- c(176.27, 176.66)
ylims <- c(-38.06, -37.70)

sc_colours <- c("Pokopoko" = "#f4a582", "Mangatoetoe" = "#ca0020",
                "Puanene" = "#9970ab", "Wharere" = "#404040", "Pongakawa" = "#92c5de")

theme_map <- function() {
  theme_bw() +
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, colour = "black"),
          panel.grid.major = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7),
          legend.key.width = unit(1.0, "cm"), legend.key.height = unit(0.3, "cm"),
          plot.title = element_text(size = 9, face = "bold", hjust = 0),
          plot.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = "white", colour = NA))
}

# ── Panel A: NZ context ───────────────────────────────────────────────────────
world <- ne_countries(scale = "medium", returnclass = "sf")

panel_a <- ggplot() +
  geom_sf(data = world, fill = "antiquewhite", linewidth = 0.2) +
  geom_sf(data = BOP_Region, fill = "#d9ead3", linewidth = 0.2, lty = 2) +
  geom_sf(data = Waihi_Estuary, fill = "grey30", linewidth = 0.3) +
  geom_rect(aes(xmin = 176.3, xmax = 176.65, ymin = -38.08, ymax = -37.7),
            fill = NA, colour = "red", linewidth = 1) +
  coord_sf(xlim = c(171.31, 180), ylim = c(-41.5, -34.5), expand = TRUE) +
  annotation_north_arrow(location = "tr", which_north = "true",
    pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
    style = north_arrow_fancy_orienteering, height = unit(0.8, "cm"), width = unit(0.8, "cm")) +
  annotation_scale(location = "br", width_hint = 0.35, text_cex = 0.6) +
  ggtitle("(a) Study location") + theme_map() + theme(legend.position = "none")

# ── Panel B: Monitoring network ───────────────────────────────────────────────
panel_b <- ggplot() +
  geom_sf(data = BOP_Region, fill = "antiquewhite", alpha = 0.6, linewidth = 0.15) +
  geom_sf(data = Waihi_Estuary, fill = "antiquewhite", linewidth = 0.8) +
  geom_sf(data = Sub_Catchments, aes(fill = Sub_Catchment), linewidth = 0.6) +
  scale_fill_manual(values = sc_colours, name = "Sub-catchment", na.value = "antiquewhite") +
  geom_sf(data = REC[REC$StreamOrde > 1, ], colour = "grey60", linewidth = 0.3, alpha = 0.5) +
  geom_sf(data = monitoring_sites, shape = 21, size = 2.5, fill = "black", colour = "white", stroke = 0.8) +
  geom_label_repel(data = monitoring_sites, aes(label = MapName, geometry = geometry),
    stat = "sf_coordinates", size = 2.5, label.padding = unit(0.10, "lines"),
    box.padding = unit(0.30, "lines"), min.segment.length = 0, seed = 42, max.overlaps = 20) +
  annotation_north_arrow(location = "tr", which_north = "true",
    pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
    style = north_arrow_fancy_orienteering, height = unit(0.8, "cm"), width = unit(0.8, "cm")) +
  annotation_scale(location = "br", width_hint = 0.35, text_cex = 0.6) +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) +
  guides(fill = guide_legend(ncol = 5, title.position = "top")) +
  ggtitle("(b) Monitoring network") + theme_map() +
  theme(legend.position = "bottom", legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6), legend.key.size = unit(0.4, "cm"))

# ── Panels C-E: Spatial predictors ────────────────────────────────────────────
make_predictor_panel <- function(fill_var, panel_title, palette, direction = 1, legend_label) {
  ggplot() +
    geom_sf(data = BOP_Region, fill = "antiquewhite", alpha = 0.6, linewidth = 0.15) +
    geom_sf(data = Waihi_Estuary, fill = "antiquewhite", linewidth = 0.8) +
    geom_sf(data = Sub_Catchments, aes(fill = .data[[fill_var]]), linewidth = 0.6) +
    scale_fill_distiller(palette = palette, direction = direction, name = legend_label, na.value = "grey90") +
    geom_sf(data = REC[REC$StreamOrde > 1, ], colour = "grey60", linewidth = 0.3, alpha = 0.5) +
    geom_sf(data = monitoring_sites, shape = 21, size = 2.5, fill = "black", colour = "white", stroke = 0.8) +
    annotation_scale(location = "br", width_hint = 0.35, text_cex = 0.6) +
    annotation_north_arrow(location = "tr", which_north = "true",
      pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
      style = north_arrow_fancy_orienteering, height = unit(0.8, "cm"), width = unit(0.8, "cm")) +
    coord_sf(xlim = xlims, ylim = ylims, expand = FALSE) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                 barwidth = unit(4, "cm"), barheight = unit(0.35, "cm"))) +
    ggtitle(panel_title) + theme_map() +
    theme(legend.position = "bottom", legend.title = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 6), legend.margin = margin(2, 0, 0, 0),
          legend.box.margin = margin(-5, 0, 0, 0))
}

panel_c <- make_predictor_panel("SEI", "(c) Soil Erodibility Index (SEI)", "YlOrRd", legend_label = "SEI")
panel_d <- make_predictor_panel("CCI", "(d) Composite Connectivity Index (CCI)", "Blues", legend_label = "CCI")
panel_e <- make_predictor_panel("LUSP", "(e) Land Use Sediment Potential (LUSP)", "YlGn", legend_label = "LUSP")

# ── Assemble ──────────────────────────────────────────────────────────────────
panel_a <- panel_a + theme(plot.margin = margin(5, 5, 5, 5))
panel_b <- panel_b + theme(plot.margin = margin(5, 5, 5, 5))
panel_c <- panel_c + theme(plot.margin = margin(5, 5, 5, 5))
panel_d <- panel_d + theme(plot.margin = margin(5, 5, 5, 5))
panel_e <- panel_e + theme(plot.margin = margin(5, 5, 5, 5))

left_col <- plot_grid(panel_a, panel_b, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.8, 1))
right_col <- plot_grid(panel_c, panel_d, panel_e, ncol = 1, align = "v", axis = "lr")

left_col  <- left_col  + theme(plot.margin = margin(0, -80, 0, 0, "pt"))
right_col <- right_col + theme(plot.margin = margin(0, -70, 0, 0, "pt"))

final_figure <- plot_grid(left_col, right_col, ncol = 2, rel_widths = c(1, 1)) +
  theme(plot.margin = margin(0, 0, 0, 0))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Fig1_Catchment_Map.png", final_figure,
       dpi = 300, width = 8, height = 11, units = "in", bg = "white")
cat("Figure 1 saved\n")
