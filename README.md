# Modeling the Impact of Outbreak Response Immunization on Typhoid Fever Using Data from 19 Outbreaks Between 2000–2022

This repository contains the R code and data used in the analysis described in:

> **Modeling the impact of outbreak response immunization on typhoid fever using data from 19 outbreaks between 2000–2022**
> Monica Duong and Jong-Hoon Kim
> Epidemiology, Public Health, Impact, International Vaccine Institute, Seoul 08826, Korea

## Overview

This analysis models the potential impact of typhoid conjugate vaccine (TCV) deployment during typhoid fever outbreaks. Using a systematic review of historical outbreak data, we simulate outbreak response vaccination (ORI) campaigns across a range of timing and coverage scenarios and estimate:

- Cases and deaths averted
- Antimicrobial resistance (AMR) burden averted (including MDR, FQNS, and XDR strains)
- Cost-effectiveness (cost per DALY averted, cost per case averted)

Uncertainty in vaccine direct effectiveness (DVE), campaign duration, and immunological delay is propagated using a Sobol quasi-random sequence (n = 200 samples).

## Repository Structure

```
typhoid_outbreak_modeling/
├── R/
│   ├── 1-data.R          # Data loading, cleaning, and preprocessing
│   ├── 2-functions.R     # Core modeling and CEA helper functions
│   ├── 3-model.R         # Simulation loop (weekly and monthly outbreaks)
│   ├── figures.R         # Code to reproduce all manuscript figures
│   └── tables.R          # Code to reproduce supplementary tables
├── data/
│   ├── typhoid_outbreak_summary_v1.csv   # Outbreak-level summary data
│   ├── typhoid_outbreak_ts_v2.csv        # Outbreak time series data
│   ├── amr-data.csv                      # AMR resistance proportions by outbreak
│   ├── parameters.csv                    # Epidemiological parameters (disability weights, durations, etc.)
│   ├── country_costs.csv                 # Cost-of-illness data by country
│   ├── WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES1.xls  # UN WPP 2024 population by age
│   ├── wpp_life_expectancy_20241022.csv  # UN WPP 2024 life expectancy
│   ├── Workforce_Worldbank.xls           # World Bank labor force participation rates
│   ├── API_NY.GDP.PCAP.xls               # World Bank GDP per capita
│   ├── imf-dm-export-20250216.xls        # IMF CPI inflation data
│   ├── weekly_output_n200_fin.csv        # Model output: weekly outbreaks, base case
│   ├── weekly_output_n200_fin_ive.csv    # Model output: weekly outbreaks, IVE sensitivity
│   ├── monthly_output_n200_fin.csv       # Model output: monthly outbreaks, base case
│   └── monthly_output_n200_fin_ive.csv   # Model output: monthly outbreaks, IVE sensitivity
└── papers/
    └── Typhoid Fever Modeling Manuscript-main_text.pdf
```

## How to Reproduce the Analysis

Scripts should be run in order:

### Step 1 — Data preparation (`R/1-data.R`)

Loads and cleans all input datasets. Key tasks:
- Applies inclusion/exclusion criteria to outbreak records
- Harmonizes time series to weekly or monthly resolution
- Adjusts cost-of-illness (COI) estimates for inflation using IMF CPI data
- Computes COI per patient from severity proportions

**Working directory:** Set your working directory to `data/` before running, or update the file paths in the script.

### Step 2 — Functions (`R/2-functions.R`)

Defines all helper functions used in the simulation:

| Function | Description |
|---|---|
| `compute_yld()` | Years lived with disability per case, by severity |
| `compute_yll()` | Years of life lost from deaths, with discounting |
| `svim_sum_across_week()` | Aggregates weekly simulation results across runs |
| `svim_sum_across_monthly()` | Aggregates monthly simulation results across runs |
| `add_cea_results()` | Computes full CEA outputs (DALYs, costs, ICER) |
| `fin_median()` | Summarizes results to median and 95% uncertainty intervals |

Source this file before running `3-model.R`.

### Step 3 — Simulation (`R/3-model.R`)

Runs the vaccine impact simulation for all outbreaks under specified scenarios:

- **Base case:** Vaccination at week 8 (weekly outbreaks) or month 2 (monthly outbreaks), 80% coverage
- **Sensitivity analysis:** IVE ranging from 0–50% (in 5% increments)

Uncertainty is propagated via 200 Sobol-sampled parameter sets drawn from distributions for:
- Direct vaccine effectiveness (DVE): truncated normal, mean = 0.83, SD = 0.087
- Campaign duration: truncated gamma, mean ~8 days
- Immunological delay: uniform, 14–28 days

Outputs are written to CSV files in `data/`.

### Step 4 — Figures and Tables

Run `R/figures.R` and `R/tables.R` after the model outputs are available. These read the pre-computed CSV outputs and reproduce all manuscript figures and supplementary tables.

## Key Parameters

| Parameter | Value | Source |
|---|---|---|
| Direct vaccine effectiveness (DVE) | 83% (mean) | Literature |
| Indirect vaccine effectiveness (IVE) | 17% (base case) | Khanam et al. |
| Vaccine cost per dose | See `parameters.csv` | — |
| Discount rate | See `parameters.csv` | — |
| AMR disability weight multiplier | 2× | — |
| Cost-effectiveness threshold | 3× GDP per capita | WHO-CHOICE |

## Data Sources

- **Outbreak data:** Systematic review of published typhoid outbreak reports
- **Population:** UN World Population Prospects 2024
- **Life expectancy:** UN World Population Prospects 2024
- **GDP per capita:** World Bank (indicator: NY.GDP.PCAP.CD)
- **Labor force participation:** World Bank / ILO modeled estimates
- **Inflation adjustment:** IMF World Economic Outlook (CPI, 2017–2025)

## Software Requirements

R (≥ 4.2) with the following packages:

```r
install.packages(c(
  "data.table", "tidyverse", "ggplot2", "dplyr", "readxl",
  "TruncExpFam", "randtoolbox", "truncnorm",
  "patchwork", "gridExtra", "ggpattern", "ggtext"
))
```

## Contact

Jong-Hoon Kim
Epidemiology, Public Health, Impact, International Vaccine Institute, Seoul 08826, Korea

For questions about the analysis, please open an issue in this repository.
