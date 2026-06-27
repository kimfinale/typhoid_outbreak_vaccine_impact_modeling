# Evidence memo — parameter ranges for the rigorous R_t analysis

Targeted (not systematic) literature review supporting `params_uncertainty.csv`.
Three uncertainty axes: generation interval (GI), asymptomatic/unreported-but-
transmitting infection, and positive predictive value (PPV) of clinical cases.

## Framing: what actually moves the estimand

The renewal reconstruction `R_t = I_t / Λ_t` and the **percent case reduction**
are **invariant** to any *time-constant* thinning of cases (reporting fraction,
PPV level, asymptomatic fraction) when the unreported infections share the
symptomatic generation interval — the constant cancels in the ratio (same algebra
as the Part D α-invariance). Rigor therefore concentrates on the components that
genuinely change R_t or the impact estimand:

1. **GI shape** (μ_g, CV, and a chronic-carrier long tail) — the dominant lever.
2. **Time-variation in PPV / reporting** — biases the *shape* of R_t (bounded
   sensitivity; constant PPV is tested for invariance).
3. **Asymptomatic transmission** — enters mechanistically through (a) a longer
   effective GI (carrier tail) and (b) the transmission-relevant vaccine efficacy
   ψ_T, since transmission the vaccine cannot block lowers ψ_T relative to ψ.

## 1. Generation interval

No direct empirical typhoid GI/serial-interval estimate exists; it is derived
from natural history (the intrinsic GI is the convolution of the latent period
and the infectiousness profile).

- **Incubation period:** usual 8–14 d; systematic review/meta-analysis subgroup
  means 9.7–21.2 d (BMC Infect Dis 2018; doi:10.1186/s12879-018-3391-3).
- **Latent period:** short — shedding can begin ~1 d post-infection (Marudaiappan 2021, PMC8892528).
- **Infectious/shedding profile:** acute + early convalescent; household study median ~11 d (range 3–61) (Marudaiappan 2021).
- **Endemic-model anchor:** time-series SIR typhoid models use a ~4-week (28 d) generation interval (Phillips 2020, PLoS NTD 14:e0007828); Pitzer 2014 (PLoS NTD 8:e2642).
- **Adopted:** μ_g base 14 d (2 wk), uncertainty 10–28 d; CV base 0.6, range 0.4–0.8.

## 2. Asymptomatic / unreported-but-transmitting infection

- **Infection : clinical-case ratio:** paired serology (STRATAA) median 24.2,
  IQR 11.4–58.9 (medRxiv 2025.03.15.25324021 / PMC12527133); seroincidence ~10×
  blood-culture incidence in Vellore (PLoS NTD; PMC11896026). Part of this gap is
  blood-culture insensitivity (≈0.61), the remainder is genuinely subclinical.
- **Controlled human infection:** ~37% of challenged naïve participants did not
  develop disease (attack rate ~63%), i.e. a substantial subclinical fraction
  (re-challenge CHIM; PMC7598925).
- **Chronic carriers:** 2–4% of cases; shed up to ~1 yr; carrier transmissibility
  estimated **low** but a pivotal uncertainty for indirect/herd effects (Pitzer
  2014; Watson & Edmunds 2015 review, PMC4504000).
- **Adopted:** asymptomatic/subclinical fraction f_a base 0.5 (range 0.25–0.75);
  relative infectiousness c_a base 0.3 (range 0.1–0.5); carrier tail weight
  π_c 0.03 (0.02–0.04) with low relative infectiousness (0–0.1).
- **Derived link:** φ_a = f_a·c_a / (f_a·c_a + (1−f_a)); ψ_T/ψ = 1 − φ_a
  (base ≈0.77; range ≈0.57–0.97 — consistent with, and now mechanistic for, the
  previous 0.6·ψ–ψ bracket).

## 3. Positive predictive value (PPV) of clinical "suspected" cases

- **Blood-culture confirmation among suspected:** ~10% (9.8% Central India,
  PMC5121673; 10.1% Tanzania, PMC6551910).
- **Blood-culture sensitivity:** 0.61 (95% CI 0.52–0.70) (Mogasale et al. 2016,
  Ann Clin Microbiol Antimicrob 15:32).
- **Implied PPV** (true typhoid among suspected) ≈ confirmed/suspected ÷
  sensitivity; in our dataset this is the per-outbreak α (Part D), capped at 1.
  Outbreak settings raise pre-test probability, so PPV is higher than the ~10%
  raw confirmation rate (≈0.2–0.7 across settings).
- **Time-variation:** no strong evidence for a specific within-outbreak PPV trend;
  modelled as a **bounded sensitivity** (relative ±0.3 linear trend over the
  outbreak), with constant-PPV invariance tested explicitly.

## Citations

- BMC Infect Dis 2018, 18:464. doi:10.1186/s12879-018-3391-3 (incubation).
- Marudaiappan 2021, J Infect Dis (PMC8892528) (shedding/latent).
- Phillips MT et al. 2020, PLoS NTD 14:e0007828 (4-wk GI anchor).
- Pitzer VE et al. 2014, PLoS NTD 8:e2642 (carriers; transmission model).
- Watson & Edmunds 2015, Vaccine 33:C42 (PMC4504000) (review; carrier role).
- STRATAA paired serology 2025 (medRxiv 2025.03.15.25324021; PMC12527133) (infection:case ratio).
- Vellore wastewater/seroincidence, PLoS NTD (PMC11896026).
- Typhoidal Salmonella re-challenge CHIM (PMC7598925) (subclinical fraction).
- Mogasale V et al. 2016, Ann Clin Microbiol Antimicrob 15:32 (culture sensitivity).
- Diagnostic confirmation rates: PMC5121673 (Central India), PMC6551910 (Tanzania).
