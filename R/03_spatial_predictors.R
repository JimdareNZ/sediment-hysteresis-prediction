# ══════════════════════════════════════════════════════════════════════════════
# 03_spatial_predictors.R
# ══════════════════════════════════════════════════════════════════════════════
# Purpose:
#   Derive spatial predictor indices for each sub-catchment from publicly
#   available geospatial datasets, then merge with the event-scale
#   hysteresis dataset.
#
# Indices derived:
#   - Land use: LUSP (Land Use Sediment Potential), DI (Disturbance Index),
#     SI (Stability Index)
#   - Soil: SEI (Soil Erodibility Index)
#   - Flow path: QCI (Quick Connectivity Index), CCI (Composite Connectivity
#     Index), DDI (Deep Drainage Index), Slope
#   - Biophysical: Sed_Risk (mean, weighted mean, tail fraction)
#
# Inputs:
#   - BOPRC Land Use 2021 shapefile
#   - Fundamental Soils Layer (FSL) shapefile
#   - LWS Physiographic Classification shapefile
#   - Biophysical sediment risk shapefile
#   - Sub-catchment boundaries shapefile
#   - data/Final_Hysteresis_Dataset_Ante.csv (from 02_event_analysis.R)
#
# Outputs:
#   - output/spatial/Land_Use2021_Indices_Subcatchment.csv
#   - output/spatial/SEI_Subcatchment.csv
#   - output/spatial/Flow_Path_Subcatchment.csv
#   - output/spatial/Biophysical_sed_risk.csv
#   - data/Final_Hysteresis_Dataset_Spatial_Ante.csv (final merged dataset)
#
# Notes:
#   - Geospatial source datasets are not redistributed in this repository.
#     Users should obtain them from the original sources listed in the
#     manuscript Data Availability Statement.
#   - All indices are aggregated to the sub-catchment scale using an
#     area-weighted mean.
#
# Dependencies:
#   sf, tidyverse
#
# Reference:
#   Dare et al. (2026). Event-scale prediction of sediment hysteresis
#   regimes using hydrological and spatial indices. Journal of Hydrology.
# ══════════════════════════════════════════════════════════════════════════════

library(sf)
library(tidyverse)

# -------------------------
# Paths — update these to match your local data locations
# -------------------------
shp_land_use       <- "data/shapefiles/BOPRC/Land_Use_2021.shp"
shp_fsl            <- "data/shapefiles/FSL/fsl-north-island-v11-all-attributes.shp"
shp_flow_path      <- "data/shapefiles/LWS_Physiographics/FSL_Flow_Path.shp"
shp_biophysical    <- "data/shapefiles/biophysical_risk/waihi_pong_sediment.shp"
shp_sub_catchments <- "data/shapefiles/waihi_estuary/WaihiEstuaryCatchmentsGenerated_MonitoringSitesSubset_JD_25052022.shp"

output_dir <- "output/spatial"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------
# Site name mapping (sub-catchment ID → anonymised site name)
# -------------------------
site_lookup <- c(
  "FN956338" = "Inanga",
  "GN083710" = "Tuangi",
  "GN086350" = "Kokopu",
  "GN265752" = "Kokako",
  "GN435054" = "Koura",
  "GN562808" = "Pipi",
  "GN618644" = "Kotuku",
  "GN849464" = "Matata",
  "GN955302" = "Tuna",
  "HN063546" = "Patiki"
)

add_site_id <- function(df) {
  df %>% mutate(SiteID = site_lookup[Name])
}

# -------------------------
# Load sub-catchments (shared across all sections)
# -------------------------
Sub_Catchments <- read_sf(shp_sub_catchments) %>% st_transform(crs = 2193)

# =========================================================================
# 1. Land use indices (LUSP, DI, SI)
# =========================================================================
cat("Deriving land use indices...\n")

Land_Use <- read_sf(shp_land_use)
Land_Use$Correct_Area <- as.numeric(st_area(Land_Use))

Subcatchment_LU <- st_intersection(st_zm(Land_Use), st_zm(Sub_Catchments)) %>%
  as.data.frame() %>%
  group_by(Name, LU2021) %>%
  summarise(area = sum(Correct_Area), .groups = "drop") %>%
  group_by(Name) %>%
  mutate(Percentage = area / sum(area) * 100) %>%
  select(Name, LU2021, Percentage) %>%
  pivot_wider(names_from = LU2021, id_cols = Name, values_from = Percentage) %>%
  replace(is.na(.), 0) %>%
  add_site_id()

# Land use sediment production categories
Subcatchment_LU_Indices <- Subcatchment_LU %>%
  mutate(
    # High disturbance / high sediment yield (weight 3)
    LU_high = arable + dairy + `sheep/beef` + deer + urban + lifestyle +
              `converted from forest exotic` + `converted from forest native`,
    # Intermediate risk (weight 2)
    LU_mid  = `forest exotic` + scrub + orchard + kiwifruit + avocado,
    # Low risk / protective (weight 1)
    LU_low  = `forest native` + reserve + wetland,
    # Weighted land use sediment potential
    LUSP_raw = 3 * LU_high + 2 * LU_mid + 1 * LU_low
  ) %>%
  mutate(
    LUSP = as.numeric(scale(LUSP_raw)),
    DI   = as.numeric(scale(LU_high)),
    SI   = as.numeric(scale(LU_low))
  ) %>%
  select(Name, SiteID, LU_high, LU_mid, LU_low, LUSP_raw, LUSP, DI, SI)

write.csv(Subcatchment_LU_Indices,
          file = file.path(output_dir, "Land_Use2021_Indices_Subcatchment.csv"),
          row.names = FALSE)

# =========================================================================
# 2. Soil Erodibility Index (SEI)
# =========================================================================
cat("Deriving Soil Erodibility Index...\n")

FSL <- read_sf(shp_fsl) %>% st_transform(crs = 2193)
FSL$Correct_Area <- as.numeric(st_area(FSL))

# Intersect with sub-catchments and select relevant attributes
# DRAIN_CLAS: drainage class (5 = well drained, 1 = very poor)
# DSLO_CLASS: depth to slow permeable layer (1 = shallow, 6 = none)
# ROCK_CLASS: rock fragments (1 = non-rocky, 5 = extremely rocky)
# CARBON_CLA: soil carbon (1 = very high, 5 = very low)
# CEC_CLASS:  cation exchange capacity (1 = very high, 5 = very low)

SEI_data <- st_intersection(st_zm(FSL), st_zm(Sub_Catchments)) %>%
  as.data.frame() %>%
  select(Name, Correct_Area, DSLO_CLASS, DRAIN_CLAS, ROCK_CLASS, CARBON_CLA, CEC_CLASS)

# Reverse-code so higher values = greater erodibility, then standardise
SEI_calc <- SEI_data %>%
  mutate(
    DRAIN_erod  = 6 - as.numeric(as.character(DRAIN_CLAS)),
    DSLO_erod   = 7 - as.numeric(as.character(DSLO_CLASS)),
    ROCK_erod   = 6 - as.numeric(as.character(ROCK_CLASS)),
    CARBON_erod = 6 - as.numeric(as.character(CARBON_CLA)),
    CEC_erod    = 6 - as.numeric(as.character(CEC_CLASS))
  ) %>%
  mutate(across(
    c(DRAIN_erod, DSLO_erod, ROCK_erod, CARBON_erod, CEC_erod),
    ~ scale(.)[, 1],
    .names = "{.col}_std"
  )) %>%
  mutate(
    SEI_polygon = DRAIN_erod_std + DSLO_erod_std + ROCK_erod_std +
                  CARBON_erod_std + CEC_erod_std
  ) %>%
  group_by(Name) %>%
  summarise(
    SEI = weighted.mean(SEI_polygon, Correct_Area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_site_id()

write.csv(SEI_calc,
          file = file.path(output_dir, "SEI_Subcatchment.csv"),
          row.names = FALSE)

# =========================================================================
# 3. Flow path and connectivity indices (QCI, CCI, DDI, Slope)
# =========================================================================
cat("Deriving flow path and connectivity indices...\n")

Flow_Path <- read_sf(shp_flow_path) %>% st_transform(crs = 2193)
Flow_Path$Correct_Area <- as.numeric(st_area(Flow_Path))

# CCI weights: overland flow (0.5), lateral flow (0.3), deep drainage (0.2)
# Assigned based on expert judgement reflecting a conceptual hierarchy
# of hydrological connectivity (see Methods, Section 3.7.2)

Flow_Indices <- st_intersection(st_zm(Flow_Path), st_zm(Sub_Catchments)) %>%
  as.data.frame() %>%
  mutate(
    OverlandFl_s = OverlandFl / 100,
    LateralDra_s = LateralDra / 5,
    DeepDraina_s = DeepDraina / 5,
    QCI   = OverlandFl_s + LateralDra_s,
    Slope = SlopeIndex,
    DDI   = DeepDraina_s,
    CCI   = 0.5 * OverlandFl_s + 0.3 * LateralDra_s + 0.2 * (1 - DeepDraina_s)
  ) %>%
  group_by(Name) %>%
  summarise(
    QCI   = weighted.mean(QCI, Correct_Area, na.rm = TRUE),
    Slope = weighted.mean(Slope, Correct_Area, na.rm = TRUE),
    DDI   = weighted.mean(DDI, Correct_Area, na.rm = TRUE),
    CCI   = weighted.mean(CCI, Correct_Area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_site_id()

write.csv(Flow_Indices,
          file = file.path(output_dir, "Flow_Path_Subcatchment.csv"),
          row.names = FALSE)

# =========================================================================
# 4. Biophysical sediment risk
# =========================================================================
cat("Deriving biophysical sediment risk indices...\n")

Biophysical <- read_sf(shp_biophysical) %>% st_transform(crs = 2193)
Biophysical$Correct_Area <- as.numeric(st_area(Biophysical))

Biophysical_risk <- st_intersection(st_zm(Biophysical), st_zm(Sub_Catchments)) %>%
  as.data.frame() %>%
  group_by(Name) %>%
  summarise(
    Sed_Risk_mean          = mean(value, na.rm = TRUE),
    Sed_Risk_weighted_mean = weighted.mean(value, Correct_Area, na.rm = TRUE),
    Sed_Risk_tail_2_5      = sum(Correct_Area[value > 2.5], na.rm = TRUE) /
                             sum(Correct_Area[!is.na(value)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  add_site_id()

write.csv(Biophysical_risk,
          file = file.path(output_dir, "Biophysical_sed_risk.csv"),
          row.names = FALSE)

# =========================================================================
# 5. Merge all spatial indices with event dataset
# =========================================================================
cat("Merging spatial indices with event dataset...\n")

Combined_Spatial <- Flow_Indices %>%
  left_join(SEI_calc, by = "Name") %>%
  left_join(Subcatchment_LU_Indices, by = "Name") %>%
  left_join(Biophysical_risk, by = "Name") %>%
  select(SiteID.x, QCI, Slope, DDI, CCI, SEI,
         LU_high, LU_mid, LU_low, LUSP, DI, SI,
         Sed_Risk_mean, Sed_Risk_weighted_mean, Sed_Risk_tail_2_5) %>%
  rename(SiteID = SiteID.x)

Hys_dat <- read.csv("data/Final_Hysteresis_Dataset_Ante.csv", check.names = FALSE)
Hys_dat <- merge(Hys_dat, Combined_Spatial, by.x = "Site", by.y = "SiteID")

write.csv(Hys_dat,
          file = "data/Final_Hysteresis_Dataset_Spatial_Ante.csv",
          row.names = FALSE)

cat("Spatial predictor derivation complete.\n")
cat("Final dataset saved: data/Final_Hysteresis_Dataset_Spatial_Ante.csv\n")
