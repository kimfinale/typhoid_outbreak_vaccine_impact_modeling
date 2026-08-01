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

Substrate: **Kabwama 2017, Kampala** — a mixed, common-source-dominant outbreak with daily
suspected-case surveillance (10,152 suspected, Jan–Jun 2015). The updated POC uses two additive
levels: suspected = true typhoid + other febrile illness, followed by true typhoid = source +
propagated incidence inside one renewal recursion. The source fraction is θ = 0.60.

The script regenerates timing-specific true and surveillance-visible reductions for coverage 0.6,
efficacy 0.8, and a two-week immunity delay; generated values should be taken from
`manuscript/data/poc_ori_kabwama_results.csv`, rather than the superseded multiplicative numbers.

## The PPV integration (the point)

The primary observation model is additive: the vaccine changes true typhoid but not the other-
febrile component. Consequently PPV need not cancel from percentage reductions, and the reduction
visible in suspected-case surveillance is smaller than the true-typhoid reduction. At the second
level, source and propagated typhoid are added inside one recursion; θ is not an outcome weight.

## Honest caveats
- R_t reconstructed from suspected includes early surveillance ramp-up / reporting artifacts
  (peak R_t ~6.8 is inflated); a Bayesian R_t with a reporting model is the rigorous version
  (see `uncertainty/`).
- The approximately flat other-febrile allocation and θ = 0.60 source fraction are structural
  assumptions; the multiplicative PPV and historical θ-weighted models are comparators only.
- Single-outbreak demo. The full analysis runs `renewal/run_analysis.R` across all
  developing-country outbreaks and propagates the pi posterior (not a point value).
