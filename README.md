# Event-scale prediction of sediment hysteresis regimes using hydrological and spatial indices

[![DOI](https://zenodo.org/badge/DOI/PLACEHOLDER.svg)](https://doi.org/PLACEHOLDER)

## Overview

This repository contains the code and data associated with:

> Dare, J.E., Özkundakci, D., McDowell, R.W., & Hussain, E. (2026). Event-scale prediction of sediment hysteresis regimes using hydrological and spatial indices. *Journal of Hydrology* (in review).

The study predicts suspended sediment hysteresis direction from readily available discharge records and static geospatial indices, without requiring continuous suspended sediment concentration monitoring. The analysis covers 678 enrichment events across eight monitoring sites in five sub-catchments of the Waihi Estuary catchment, Bay of Plenty, New Zealand.

## Repository structure

```         
├── README.md
├── LICENSE
├── data/
│   └── Final_Hysteresis_Dataset_Spatial_Ante.csv
├── R/
│   ├── 01_event_definition.R
│   ├── 02_zuecco_classification.R
│   ├── 03_spatial_predictors.R
│   ├── 04_xgboost_modelling.R
│   ├── 05_pca_analysis.R
│   └── 06_figures.R
└── tables/
```

## Data

`Final_Hysteresis_Dataset_Spatial_Ante.csv` contains the processed event-scale dataset used in all analyses. The dataset comprises 716 events across eight monitoring sites, with anonymised site names. Each row represents a single hydrological event and includes:

-   Event-scale hydrological metrics (peak discharge, rising time, recession time, baseflow index, flashiness, etc.)
-   Sediment transport metrics (load, yield)
-   Zuecco hysteresis classification (classes I–VIII)
-   Antecedent condition indices (time since previous event, preceding event magnitude, event frequency)
-   Static spatial predictor indices (SEI, CCI, QCI, DDI, Slope, LUSP, DI, SI)

Raw sensor data and ANN model objects are described in the companion paper: [Dare et al. (2026), *Scientific Reports*](https://doi.org/10.1038/s41598-026-43915-9).

## Scripts

All analyses were carried out in R (R Core Team, 2025). Scripts are numbered in the order they should be run:

| Script | Description |
|----------------------------|--------------------------------------------|
| `01_event_definition.R` | Event identification using hydroEvents package with manual refinement; hydrograph metric extraction |
| `02_zuecco_classification.R` | R implementation of the Zuecco et al. (2016) hysteresis classification framework |
| `03_spatial_predictors.R` | Derivation of spatial predictor indices (SEI, CCI, QCI, DDI, Slope, LUSP, DI, SI) from national geospatial datasets |
| `04_xgboost_modelling.R` | XGBoost classification framework using tidymodels, including 5-fold CV and LOSO CV for three class subsets |
| `05_pca_analysis.R` | Principal component analysis and ANOVA tests on PC1–PC4 |
| `06_figures.R` | All manuscript and supplementary figures |

## Requirements

-   R (≥ 4.3.0)
-   Key packages: `tidymodels`, `xgboost`, `sf`, `ggplot2`, `cowplot`, `hydroEvents`

Install all required packages:

``` r
install.packages(c("tidymodels", "xgboost", "sf", "ggplot2", "cowplot", 
                    "hydroEvents", "ggspatial", "ggrepel", "rnaturalearth",
                    "rnaturalearthdata"))
```

## Geospatial data sources

The spatial predictor indices were derived from the following publicly available datasets:

-   **Fundamental Soils Layer (FSL)**: [Landcare Research (2020)](https://lris.scinfo.org.nz/layer/48104-fundamental-soils-layer-nz/)
-   **Physiographic and Hydrological Classification**: Pearson, L. & Rissmann, C. (2021). Report 2021/25, Land and Water Science, New Zealand.
-   **Land use data**: Bay of Plenty Regional Council

These datasets are not redistributed in this repository. Users wishing to reproduce the spatial predictor derivation should obtain them from the original sources.

## Citation

If you use this code or data, please cite:

```         
Dare, J.E., Özkundakci, D., McDowell, R.W., & Hussain, E. (2026). 
Event-scale prediction of sediment hysteresis regimes using hydrological 
and spatial indices. Journal of Hydrology (in review).
```

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contact

James Dare — jd227\@students.waikato.ac.nz