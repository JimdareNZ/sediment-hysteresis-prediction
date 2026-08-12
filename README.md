# Bimodal predictability of suspended sediment hysteresis direction from landscape and hydrological data enables spatially targeted catchment management

[![DOI](https://zenodo.org/badge/DOI/PLACEHOLDER.svg)](https://doi.org/PLACEHOLDER)

## Overview

This repository contains the code and data associated with:

> Dare, J.E., Özkundakci, D., McDowell, R.W., & Hussain, E. (2026). Bimodal predictability of suspended sediment hysteresis direction from landscape and hydrological data enables spatially targeted catchment management. *Journal of Hydrology* (in review).

The study predicts suspended sediment hysteresis direction from readily available discharge records and static geospatial indices, without requiring continuous suspended sediment concentration monitoring. The analysis covers 678 enrichment events across eight monitoring sites in five sub-catchments of the Waihi Estuary catchment, Bay of Plenty, New Zealand.

## Repository structure

```
├── README.md
├── LICENSE
├── data/
│   └── Final_Hysteresis_Dataset_Spatial_Ante.csv
├── R/
│   ├── functions.R
│   ├── 01_event_definition.R
│   ├── 02_event_analysis.R
│   ├── 03_spatial_predictors.R
│   ├── 04_xgboost_modelling.R
│   ├── 05_pca_analysis.R
│   ├── 06_fig1_catchment_map.R
│   ├── 06_fig2_hysteresis_examples.R
│   ├── 06_fig3_class_distribution.R
│   ├── 06_fig4_hysteresis_index.R
│   ├── 06_fig5_pca_biplots.R
│   ├── 06_fig6_model_results.R
│   └── 06_fig7_conceptual_framework.R
└── figures/
```

## Data

`Final_Hysteresis_Dataset_Spatial_Ante.csv` contains the processed event-scale dataset used in all analyses. The dataset comprises 716 events across eight monitoring sites, with anonymised site names. Each row represents a single hydrological event and includes:

- Event-scale hydrological metrics (peak discharge, rising time, recession time, baseflow index, flashiness, etc.)
- Sediment transport metrics (load, yield)
- Zuecco hysteresis classification (classes I–VIII)
- Antecedent condition indices (time since previous event, preceding event magnitude, event frequency)
- Static spatial predictor indices (SEI, CCI, QCI, DDI, Slope, LUSP, DI, SI)

Raw sensor data and ANN model objects are described in the companion paper: [Dare et al. (2026), *Scientific Reports*](https://doi.org/10.1038/s41598-026-43915-9).

## Scripts

All analyses were carried out in R (R Core Team, 2025). Scripts are numbered in the order they should be run. Figure scripts (prefixed `06_`) can be run independently once the analysis scripts have been executed.

### Analysis pipeline

| Script | Description |
|--------|-------------|
| `functions.R` | Shared utility functions: event definition, hydrograph metrics, and hysteresis classification (Zuecco et al., 2016; Lloyd et al., 2016) |
| `01_event_definition.R` | Event identification using hydroEvents with interactive manual refinement |
| `02_event_analysis.R` | Hydrograph metric extraction, hysteresis index calculation, and antecedent condition indices |
| `03_spatial_predictors.R` | Derivation of spatial predictor indices (SEI, CCI, QCI, DDI, Slope, LUSP, DI, SI) from national geospatial datasets |
| `04_xgboost_modelling.R` | XGBoost classification using tidymodels, with 5-fold CV and LOSO CV for three class subsets (all classes, pure, figure-eight) |
| `05_pca_analysis.R` | Principal component analysis and ANOVA tests; saves PCA results for figure generation |

### Figure scripts

| Script | Figures produced |
|--------|-----------------|
| `06_fig1_catchment_map.R` | Figure 1: Study area and spatial predictor maps |
| `06_fig2_hysteresis_examples.R` | Figure 2: Representative hysteresis loop examples |
| `06_fig3_class_distribution.R` | Figure 3: Hysteresis class distribution across sites |
| `06_fig4_hysteresis_index.R` | Figure 4: Continuous hysteresis index violin plots |
| `06_fig5_pca_biplots.R` | Figure 5: PCA biplots; Figures S1–S2: supplementary PCA biplots |
| `06_fig6_model_results.R` | Figure 6: Model performance (ROC, confusion matrix, variable importance); Figures S3–S4: supplementary model results |
| `06_fig7_conceptual_framework.R` | Figure 7: Conceptual framework diagram |

## Requirements

- R (≥ 4.3.0)
- Key packages: `tidymodels`, `xgboost`, `sf`, `ggplot2`, `cowplot`, `hydroEvents`, `DiagrammeR`

Install all required packages:

```r
install.packages(c("tidymodels", "xgboost", "sf", "ggplot2", "cowplot",
                    "hydroEvents", "ggspatial", "ggrepel", "rnaturalearth",
                    "rnaturalearthdata", "DiagrammeR", "DiagrammeRsvg",
                    "rsvg", "zoo", "padr", "reshape2", "vip", "scales"))
```

## Geospatial data sources

The spatial predictor indices were derived from the following publicly available datasets:

- **Fundamental Soils Layer (FSL)**: [Landcare Research (2020)](https://lris.scinfo.org.nz/layer/48104-fundamental-soils-layer-nz/)
- **Physiographic and Hydrological Classification**: Pearson, L. & Rissmann, C. (2021). Report 2021/25, Land and Water Science, New Zealand.
- **Land use data**: Bay of Plenty Regional Council

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

James Dare — jd227@students.waikato.ac.nz
