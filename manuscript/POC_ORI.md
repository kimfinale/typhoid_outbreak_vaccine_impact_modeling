# Proof-of-concept: ORI impact on developing-country outbreaks, with PPV (pi) propagation

Goal: show the end-to-end pipeline — candidate **time-series** data + candidate **PPV** data
→ **ORI (outbreak-response immunization) impact** — works, so it can be re-run when the
datasets grow. (User: "when the analyses are successful, it's not a big deal to redo … when
the dataset is modified.")

## Data readiness

**PPV (pi).** Ready from the `latent_class_ppv` fit — community PPV of the clinical case
definition, itself from developing-country outbreaks:
Aye/Myanmar **0.26**, Neil/Kasese-Uganda **0.65**, Kabwama/Kampala-Uganda **0.23**
(`latent_class_ppv/tables/final_pi_community.csv`). Hospital phi per study also available.

**Time series.** Two streams:
- *Existing* tf dataset (`data/Typhoid_Outbreak_Time_Series_2000_2022`): 39 outbreaks, many
  developing-country. Three overlap exactly with the pi estimates (Aye, Neil, Kabwama).
- *New* 2023+ candidates (`ori_timeseries_candidates.csv`): 32 yes. PDF fetch (OA-first,
  `fetch_pdfs.py`) got **13/32** via Unpaywall; the key new **developing-country outbreak
  curves (Malaysia 2019, Bangladesh 2024) failed** (not OA, Sci-Hub dead) → need institutional
  access. The Netherlands cruise-ship curve (`ts_poc_cruiseship.csv`) was digitized but SET
  ASIDE: high-income, common-source, closed-population — wrong substrate for ORI.

## POC analysis (`poc_ori_kabwama.R`)

Substrate: **Kabwama 2017, Kampala** — large urban **propagated** outbreak, daily SUSPECTED
surveillance (10,152 suspected, Jan–Jun 2015), community pi = 0.23. Reuses the renewal ORI
engine (`renewal/R/renewal_core.R`): reconstruct R_t → vaccination counterfactual with
transmission feedback (emergent herd effect).

Results (coverage 0.6, transmission-VE 0.8, +2 wk immunity):

| campaign protected by | % reduction | suspected averted | TRUE typhoid averted (pi=0.23) |
|---|---|---|---|
| wk 3 | 88% | 8,947 | 2,058 |
| wk 5 | 84% | 8,476 | 1,949 |
| wk 7 | 72% | 7,313 | 1,682 |
| wk 9 | 57% | 5,784 | 1,330 |

pi CI at the earliest campaign: **1,387–3,740** true cases averted (pi 0.155–0.418).

## The PPV integration (the point)

Running ORI on a **suspected**-case curve yields averted counts in *suspected* units.
Under a constant-PPV multiplicative regime, R_t and the **% reduction are pi-invariant**, but
the **absolute TRUE burden averted = pi × suspected-averted** — the CEA/DALY numerator is
~4× smaller than the suspected count implies (pi=0.23). Campaigns are triggered and sized on
suspected counts, yet only ~23% are vaccine-preventable typhoid.

## Honest caveats
- R_t reconstructed from suspected includes early surveillance ramp-up / reporting artifacts
  (peak R_t ~6.8 is inflated); a Bayesian R_t with a reporting model is the rigorous version
  (see `uncertainty/`).
- Constant-PPV (multiplicative) regime assumed; an **additive** non-typhoid background B(t)
  would instead *dilute* the observed % reduction — worth testing when B is estimable.
- Single-outbreak demo. The full analysis runs `renewal/run_analysis.R` across all
  developing-country outbreaks and propagates the pi posterior (not a point value).
