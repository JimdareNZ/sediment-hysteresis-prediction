# Bimodal predictability of suspended sediment hysteresis direction from landscape and hydrological data enables spatially targeted catchment management

[![DOI](https://zenodo.org/badge/1318867872.svg)](https://doi.org/10.5281/zenodo.21899318)

Code and data accompanying:

> Dare, J.E., Özkundakci, D., McDowell, R.W., & Hussain, E. (2026). Bimodal predictability of suspended sediment hysteresis direction from landscape and hydrological data enables spatially targeted catchment management. *Journal of Hydrology* (in review).

The study predicts suspended sediment hysteresis direction from discharge records and static geospatial indices, without requiring continuous suspended sediment concentration monitoring. It covers eight monitoring sites in five sub-catchments of the Waihi Estuary catchment, Bay of Plenty, New Zealand.

## Repository structure

```
├── README.md
├── LICENSE
├── data/
│   ├── Final_Hysteresis_Dataset_Spatial_Ante_Final.csv   # analysis dataset (716 events)
│   └── fig2_events/                         # four event time series used in Fig. 2
├── R/
│   ├── functions.R                          # shared functions (sourced by 01, 02)
│   ├── 01_event_definition.R
│   ├── 02_event_analysis.R
│   ├── 03_spatial_predictors.R
│   ├── 04_xgboost_modelling.R
│   ├── 05_pca_analysis.R
│   ├── 06_fig1_catchment_map.R ... 06_fig7_conceptual_framework.R
│   ├── 07_table4_load_partitioning.R
│   ├── 08_sensitivity_label_perturbation.R
│   ├── 09_sensitivity_weighting.R
│   └── 10_SI_sensitivity_figures.R
├── results/        # intermediate outputs; includes archived sensitivity results
├── figures/        # figure outputs
└── tables/         # LaTeX table outputs
```

**All scripts must be run from the repository root** (e.g. `Rscript R/07_table4_load_partitioning.R`). All paths are relative to the root, and each script creates its output folder if it does not exist.

## Data

`data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv` is the processed event-scale dataset used in all analyses, figures and tables. It is the final output of scripts 01–03. Each row is one hydrological event and includes:

- Event-scale hydrological metrics (peak discharge, rising/recession time, baseflow index, flashiness, etc.)
- Sediment transport metrics (load, yield)
- Zuecco hysteresis classification (`Zuec_class`, classes 1–8 = I–VIII) and index (`Zuec_h`), and the Lloyd hysteresis index
- Antecedent condition indices (time since previous event, preceding event magnitude, event frequency)
- Static spatial predictor indices (SEI, CCI, QCI, DDI, Slope, LUSP, DI, SI)

Use only this `_Final` file. An earlier version without the `_Final` suffix had inverted Zuecco classes and is not included.

`data/fig2_events/` contains the 15-minute discharge and ANN-derived TSS time series (`Time`, `DISC`, `.pred`) for the four example events in Figure 2. File names are as written by `01_event_definition.R`:

| File | Site | Class |
|---|---|---|
| `Event_84_E.csv` | Kokako (MGT) | I |
| `Event_4.csv` | Pipi (PUN) | I |
| `Event_19.csv` | Kotuku (WHR-D) | IV |
| `Event_8.csv` | Inanga (PKO-UW) | IV |

### Site codes

Site names in the dataset are anonymised field names. They map to the manuscript site codes as follows:

| Dataset `Site` | Manuscript code | Sub-catchment |
|---|---|---|
| Tuangi | PKO-D | Pokopoko |
| Inanga | PKO-UW | Pokopoko (nested within PKO-D) |
| Kokopu | PKO-UE | Pokopoko (nested within PKO-D) |
| Kokako | MGT | Mangatoetoe |
| Pipi | PUN | Puanene |
| Kotuku | WHR-D | Wharere |
| Matata | PNG-M | Pongakawa |
| Tuna | PNG-TR1 | Pongakawa (nested within PNG-M) |

## Scripts and run order

| Script | Description | Produces |
|---|---|---|
| `functions.R` | Event definition helpers, hydrograph metrics, Zuecco (2016) and Lloyd (2016) hysteresis indices | — |
| `01_event_definition.R` | Event identification with hydroEvents plus interactive manual refinement (run once per site) | `results/events/` |
| `02_event_analysis.R` | Hydrograph metrics, hysteresis classification and antecedent indices | `data/Final_Hysteresis_Dataset_Ante.csv` |
| `03_spatial_predictors.R` | Spatial predictor indices from national geospatial datasets, merged with the event data | `results/spatial/`, `data/Final_Hysteresis_Dataset_Spatial_Ante_Final.csv` |
| `04_xgboost_modelling.R` | XGBoost classification (tidymodels), 5-fold CV and leave-one-site-out CV; run once per class subset (`all`, `pure`, `fig_eight`) | `results/models/` |
| `05_pca_analysis.R` | PCA and ANOVA tests | `results/pca_results.rds`, Table S6 (console) |
| `06_fig1_catchment_map.R` | Study area and spatial predictor maps | Fig. 1 |
| `06_fig2_hysteresis_examples.R` | Representative hysteresis loops | Fig. 2 |
| `06_fig3_class_distribution.R` | Hysteresis class distribution by event count and by load | Fig. 3 |
| `06_fig4_hysteresis_index.R` | Hysteresis index distributions | Fig. 4 |
| `06_fig5_pca_biplots.R` | PCA biplots | Fig. 5, Figs. S1–S2 |
| `06_fig6_model_results.R` | Model performance (ROC, confusion matrix, variable importance) | Fig. 6, Figs. S3–S4 |
| `06_fig7_conceptual_framework.R` | Conceptual framework diagram | Fig. 7 |
| `07_table4_load_partitioning.R` | Sediment load by hysteresis direction and site, with outlet-only network totals | Table 4 |
| `08_sensitivity_label_perturbation.R` | Label-flip sensitivity of LOSO performance | `results/sensitivity_results.csv`, `results/sensitivity_summary.csv` |
| `09_sensitivity_weighting.R` | CCI/LUSP weighting sensitivity | `results/weighting_sensitivity_results.csv`, `results/weighting_sensitivity_importance.csv` |
| `10_SI_sensitivity_figures.R` | SI sensitivity figures from the archived results | Figs. S5–S6 |

Output files:

- Main text: `figures/Fig1_Catchment_Map.png`, `Fig2_Hysteresis_Examples.png`, `Fig3_Hysteresis_Class_Distribution.png`, `Fig4_Hysteresis_Index_Distribution.png`, `Fig5_PCA_Combined.png`, `Fig6_Model_Results.png`, `Fig7_Conceptual_Framework.png`; `tables/Table4_Load_By_Direction.tex`
- SI: `figures/FigS1_PCA_PC1_PC3.png`, `FigS2_PCA_PC1_PC4.png`, `FigS3_AllClasses_Model_Results.png`, `FigS4_FigureEight_Model_Results.png`, `FigS5_sensitivity.png`, `FigS6_weighting_sensitivity.png`

### What can be run from this repository alone

Scripts 04–05, 06 (Figs. 2–7), 07 and 10 run directly from the archived data and results. Scripts 01–03 and Fig. 1 also need raw inputs that are not redistributed:

- **01–02** need the ANN-derived TSS and discharge time series (`data/raw/<site>_TSS_Final.csv`). These are described in the companion paper: [Dare et al. (2026), *Scientific Reports*](https://doi.org/10.1038/s41598-026-43915-9).
- **03 and Fig. 1** need the geospatial shapefiles listed below, placed under `data/shapefiles/`.
- **Table S7** (temporal resolution of validation events) is based on raw discrete laboratory sample data, which is not archived in this repository.

## Runtime and reproducibility

- `08_sensitivity_label_perturbation.R` takes about **8 hours**; `09_sensitivity_weighting.R` takes about **7 minutes**. All other scripts take minutes or less.
- The sensitivity results in `results/` were generated before a random seed was fixed. Scripts 08 and 09 now call `set.seed(2026)`, so reruns are reproducible but will differ slightly from the medians reported in the manuscript. The archived CSVs in `results/` reproduce the reported values exactly, and `10_SI_sensitivity_figures.R` reads them directly. Rerunning 08 or 09 overwrites the archived CSVs, so back them up first.

## Requirements

R (≥ 4.3.0; tested with 4.5.1) and the following packages:

```r
install.packages(c(
  "tidyverse", "tidymodels", "xgboost", "vip", "dials", "scales",
  "sf", "rnaturalearth", "rnaturalearthdata", "ggspatial", "ggrepel", "cowplot",
  "hydroEvents", "zoo", "padr", "reshape2", "lubridate",
  "DiagrammeR", "DiagrammeRsvg", "rsvg"
))
```

(`dplyr`, `tidyr`, `purrr`, `tibble` and `ggplot2` are installed with `tidyverse`.)

## Geospatial data sources

The spatial predictor indices were derived from the following datasets, which are not redistributed here:

- **Fundamental Soils Layer (FSL)**: [Landcare Research (2020)](https://lris.scinfo.org.nz/layer/48104-fundamental-soils-layer-nz/)
- **Physiographic and Hydrological Classification**: Pearson, L. & Rissmann, C. (2021). Report 2021/25, Land and Water Science, New Zealand.
- **Land use data**: Bay of Plenty Regional Council

## Citation

```
Dare, J.E., Özkundakci, D., McDowell, R.W., & Hussain, E. (2026).
Bimodal predictability of suspended sediment hysteresis direction from
landscape and hydrological data enables spatially targeted catchment
management. Journal of Hydrology. Code and data:
https://doi.org/10.5281/zenodo.21899319
```

## License

MIT License — see [LICENSE](LICENSE).

## Contact

James Dare — jd227@students.waikato.ac.nz
