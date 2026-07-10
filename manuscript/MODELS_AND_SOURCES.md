# Models & definitive sources — typhoid PPV and ORI

## 1. Latent-class PPV model

**Definitive source:** `latent_class_ppv/stan/model_final.stan` (fit by `latent_class_ppv/run_final.R`).

**Goal.** Without a gold standard, jointly estimate (a) blood-culture sensitivity `Se_BC`
(volume-dependent) and bone-marrow sensitivity `Se_BM`, and (b) the **PPV of the clinical /
suspected case definition** — `P(true typhoid | suspected)` — separately for hospital-referred
(`phi_s`) and community/outbreak (`pi_o`) settings.

**Two data layers sharing transportable test parameters** (see `fig0_latent_class_model.png`):

1. **Paired-culture layer (Hui–Walter conditional-independence LCM).** For each of 10 historic
   studies with paired **blood + bone-marrow** culture, the 2×2 cross-classification
   (BC±, BM±) is `Multinomial` with cell probabilities built from hospital prevalence `phi_s`
   (= the hospital PPV), `Se_BC`, `Se_BM`, and specificities. Identifies the **level** of Se_BC,
   Se_BM, and `phi_s`. Conditional independence of the two cultures given true status is assumed
   (a conditional-dependence variant is the sensitivity check).

2. **Outbreak-positivity layer.** For each outbreak with suspected→tested→confirmed counts, the
   culture-positive count `k_o` among `n_o` tested is
   `k_o ~ Binomial(n_o, pi_o·Se_BC(V_o) + (1-pi_o)(1-Sp))` — true typhoid detected at Se_BC plus
   false positives at `1-Sp`. Identifies the community PPV `pi_o` given the Se_BC curve.

**The bridge — volume-anchored Se_BC.** `logit Se_BC = alpha0 + alpha1·log(volume) + tau·u`.
`alpha1` (volume slope) has an **Antillón-2018-informed** prior `N(0.36, 0.15)`; `alpha0` = level;
`tau` = between-study heterogeneity. Making Se_BC a function of blood volume is what lets the
hospital layer's Se inform the outbreak layer's `pi`.

**Identification & priors.** In the positivity layer Se_BC and `pi` enter only as the product
`pi·Se` (a ridge); the paired layer + the volume prior pin Se_BC and break the ridge.
`Se_BM ~ Beta(12,2)` (favours high — bone marrow is the better reference; also breaks LCM
label-switching). `Sp ~ Beta(200,1)` (≈0.995, pinned near 1).

**Transportable vs local.** *Transportable* (test): `alpha0, alpha1, tau, Se_BM, Sp` — reused
everywhere. *Local* (population): `phi_s` (hospital), `pi_o` (community) — the case-definition PPV
is **not** one transferable constant (spans ~3-fold, 0.23–0.65 in the community set).

**Companion:** `latent_class_ppv/stan/community_ppv.stan` — standalone partially-pooled community
model (`logit pi_o = mu + sigma·z_o`, Se_BC prior-pinned) that yields a **transferable predictive
PPV `pi_new`** for a new outbreak. This is the principled source for the community hyper-
distribution; the ORI integration currently reconstructs `mu/sigma` from `model_final`'s `pi`
draws (documented approximation in `renewal/R/ppv.R`).

**Fit / outputs.** cmdstanr, 4 chains, seed 2026 →
`final_parameters.csv`, `final_phi_hospital.csv`, `final_pi_community.csv`, `final_draws.rds`.

**Headline estimates.** Se_BC 0.52 / 0.62 / 0.68 at 2 / 5 / 10 mL; Se_BM 0.90; Sp 0.99;
hospital φ 0.66–0.99; community π 0.23 (Kabwama) / 0.26 (Aye) / 0.65 (Neil).

### PPV data & literature (definitive)
| Input | File | Provenance |
|---|---|---|
| Paired BC/BM 2×2 + volumes | `latent_class_ppv/data/mogasale2016_paired_bc_bmc_seed.csv` | Mogasale 2016 review (Akoh 1991, Gasem 1995/2002, Wain 2008, …) |
| Outbreak suspected/tested/confirmed | `latent_class_ppv/data/community_surveillance_ppv.csv` | Aye 2004, Neil 2012, Kabwama 2017 |
| Se_BC vs blood volume (slope prior) | `latent_class_ppv/data/antillon2018_volume_sensitivity.csv` | Antillón 2018 |
| Clinical-definition Se/Sp (PPV spectrum) | `latent_class_ppv/data/clinical_definition_accuracy.csv` | Thriemer 2012, Hosoglu 2006, Storey 2015 |

Methods lineage: Hui & Walter 1980 (LCM identification); Mogasale 2016; Antillón 2018; Reitsma
2005 / Rutter–Gatsonis 2001 (DTA pooling for the clinical layer, if triggered).

## 2. Renewal-equation ORI model

**Definitive source:** `renewal/` module — driver `renewal/run_analysis.R`, config `renewal/config.yml`.
- `renewal/R/renewal_core.R` — exact R_t reconstruction (`R_t = I_t / Lambda_t`) + vaccination
  counterfactual with transmission feedback (emergent herd effect).
- `renewal/R/scenario.R` — Sobol scenario grid (coverage × timing × VE draws), both renewal & static.
- `renewal/R/ppv.R` — **PPV (π) posterior propagation**; additive `S=I+B` and multiplicative regimes.
- `renewal/R/gi.R` — typhoid generation interval (mean 14 d, sd 8.4 d, weekly).
- `R/2-functions.R` — **reused** CEA (DALY/cost: `compute_yld/yll`, `add_cea_results`).

### ORI data (definitive)
| Input | File |
|---|---|
| Outbreak time series (39 outbreaks; 13 resolution-eligible) | `data/Typhoid_Outbreak_Time_Series_2000_2022_Timeseries.csv` |
| Outbreak metadata (population, deaths) | `data/Typhoid_Outbreak_Time_Series_2000_2022_Summary.csv` |
| AMR resistance proportions | `data/amr-data.csv` |
| **PPV posterior input** | `latent_class_ppv/outputs/final_draws.rds` (community π draws) |

Outputs: `renewal/tables/*` (`summary_delay_coverage`, `tab_ppv_effect`, …), `renewal/figures/*`.
Methods lineage: Cori 2013 / Wallinga–Teunis renewal equation; Halloran static comparator.
