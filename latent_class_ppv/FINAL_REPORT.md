# Blood-culture sensitivity, bone-marrow sensitivity, and clinical-definition PPV for typhoid — final report

A Bayesian latent-class analysis that estimates, from real data: (a) blood-culture sensitivity
`Se_BC` as a function of blood volume, (b) bone-marrow-culture sensitivity `Se_BM`, and (c) the
positive predictive value (PPV) of the clinical/suspected case definition under **hospital**
(`phi_s`) and **community-surveillance** (`pi_o`) ascertainment. Fit with cmdstanr; converged
cleanly (0 divergences, max R-hat 1.002, min ESS 1826). Reproduce: `Rscript latent_class_ppv/run_final.R`.

---

## 1. Executive summary (headline estimates, median [90% CrI])

| Quantity | Estimate |
|---|---|
| **Se_BC at 2 mL** | **0.52 [0.44, 0.61]** |
| **Se_BC at 5 mL** | **0.62 [0.54, 0.69]** |
| **Se_BC at 10 mL** | **0.68 [0.59, 0.76]** |
| volume slope (per log-mL, logit) | 0.42 [0.18, 0.65] |
| **Se_BM (bone marrow)** | **0.90 [0.87, 0.93]** |
| **Sp_BC (specificity)** | **0.99 [0.98, 1.00]** |
| between-study heterogeneity `tau` (logit) | 0.49 [0.25, 0.85] |
| **Hospital PPV `phi_s`** (10 historic studies) | **0.66 – 0.99** (median 0.83) |
| **Community PPV `pi_o`** (3 outbreaks) | Aye **0.26** [0.09,0.60], Kabwama **0.23** [0.16,0.42], Neil **0.65** [0.44,0.92] |

**Three take-home numbers:** blood culture detects **~60%** of true typhoid at a typical 5 mL draw
(rising with volume); **bone marrow ~90%**; and the clinical case definition's PPV is **high in
hospitals (~0.8) but low in community surveillance (~0.25)** — they are not the same quantity.

---

## 2. Problem and approach

Observed typhoid counts are a layered observation on true infection. For the cultured subset,
`positivity = k/n = pi * Se_BC` (with Sp_BC ~ 1). A single test cannot separate `pi` (PPV, a
population property) from `Se_BC` (a test property) — only their product. The identifying leverage
is **historic studies that cultured the same patients by BOTH blood and bone marrow**: two
conditionally-independent, near-perfectly-specific tests are jointly identified (Hui-Walter), which
pins `Se_BC` and `Se_BM`. That anchored `Se_BC` then lets each outbreak's positivity identify its
`pi`. Transportable parameters (`Se_BC` volume curve, `Se_BM`) are shared; every PPV/prevalence
(`phi_s`, `pi_o`) is local.

## 3. Data consolidated

- **Historic paired blood+bone-marrow (10 studies, 635 patients):** the 2x2 tables (BC+BM+, BC+BM-,
  BC-BM+, BC-BM-) + blood volumes, from Mogasale 2016 (Vallenas both-negative cell corrected to 15;
  Hirsowitz volume imputed at the median). Volumes 2-9 mL.
- **Community outbreak positivity (3 tight syndromic points):** Aye 2004 (3/22), Neil 2012 (25/63),
  Kabwama 2017 (51/364) — culture-confirmed / cultured, from the outbreak Summary sheet. Assumed
  community blood volume 5 mL (swept 3-10 mL).
- **Volume-slope prior:** Antillon 2018 meta-regression (Se 0.51@2mL -> 0.65@10mL) informs `alpha1`.

## 4. Model (final; the severity offset was dropped — see Section 7)

```
Se sub-model:   logit(Se_BC) = alpha0 + alpha1*log(volume) + tau*u
Historic LCM:   (BC+BM+, BC+BM-, BC-BM+, BC-BM-) ~ Multinomial(N, p(phi_s, Se_BC, Se_BM, Sp))
Outbreak:       k_o ~ Binomial(n_o, pi_o*Se_BC(vol_o) + (1-pi_o)*(1-Sp_BC))
```
`Se_BM ~ Beta(12,2)` (favours high, breaks the label-switch symmetry), `Sp ~ Beta(200,1)`,
`alpha1 ~ Normal(0.36, 0.15)` (Antillon), `phi_s, pi_o ~ Beta(1,1)` independent (not pooled).

## 5. Results

### 5a. Blood-culture sensitivity `Se_BC` — volume-dependent
Se_BC rises with blood volume: **0.52 (2 mL) -> 0.62 (5 mL) -> 0.68 (10 mL)**, slope 0.42 [0.18,0.65]
per log-mL (clearly > 0). This reproduces Mogasale (66% vs bone marrow) and Antillon (0.59 pooled;
0.51->0.65 over 2-10 mL) from an independent latent-class fit. Substantial residual heterogeneity
(`tau` 0.49) remains after volume — consistent with unmodelled media, technique, illness week, and
prior antibiotics. See `figures/fig_final_se_volume.png`.

### 5b. Bone-marrow-culture sensitivity `Se_BM`
**0.90 [0.87, 0.93]** — bone marrow is markedly more sensitive than blood culture (~0.6), matching
the literature (~90-96%). This is what justifies treating bone marrow as the near-reference AND
what makes `Se_BC` identifiable in the paired data. Estimating `Se_BM` (rather than fixing it, as
Arora 2019 did) also relaxes Mogasale's union-complete assumption (the both-negative cell is a free
consequence of `phi_s`, `Se_BC`, `Se_BM`, not assumed all-true-negative).

### 5c. Specificity
`Sp_BC = Sp_BM = 0.99` — a positive culture is essentially always a true case, so confirmed cases
are ~all true typhoid and `pi` drops out of the confirmed-case forward map (matters for absolute
burden, not ratios).

### 5d. PPV of the clinical case definition — hospital vs community
The clinical/suspected definition's PPV is **not one number** — it depends on the population:
- **Hospital (`phi_s`, historic studies): 0.66-0.99, median 0.83.** Hospital-referred suspected
  typhoid patients are heavily enriched for true typhoid.
- **Community surveillance (`pi_o`, outbreaks): Aye 0.26, Kabwama 0.23, Neil 0.65.** Community
  syndromic surveillance casts a wide febrile net, so most "suspected" are not typhoid — PPV is far
  lower. (Neil is higher and is the IP-skewed / ascertainment-changed outbreak.)

**Hospital PPV (~0.8) markedly exceeds community PPV (~0.25)** — a clean, coherent separation
(`figures/fig_final_ppv.png`). This is the central applied finding: **you cannot use a hospital-
derived PPV for community surveillance, or vice versa.**

## 6. Key inferences

1. **Blood culture misses ~40% of true typhoid at typical volumes**, and the miss rate is
   **volume-sensitive** — drawing 10 vs 2 mL raises detection from ~0.52 to ~0.68. Burden/incidence
   corrections must be volume-aware (echoing Antillon), not a flat 0.5.
2. **Bone marrow (~0.90) is the practical near-gold-standard**, but its invasiveness is exactly why
   it is confined to hospital/severe cases — so it cannot be used to anchor community Se directly.
3. **Community-surveillance PPV is low (~0.25) and setting-specific.** For outbreak burden/CEA, the
   suspected-case count must be scaled by a *community* PPV (~0.25 here), not a hospital PPV; using
   the wrong one biases absolute cases by ~3-fold.
4. **The community PPV is not a transferable constant** — it spans 0.23-0.65 across three outbreaks
   and depends on the local true-typhoid : other-febrile ratio; usable only as a wide prior.
5. **Sanity check passes:** community `pi_o` sit below hospital `phi_s`, as theory requires.

## 7. Decisions made along the way (the audit trail)

- **Novelty (Phase 0):** materially additive over Mogasale (pooled Se only), Arora (fixes Se, no
  PPV, no volume), Antillon (volume but no PPV/observation model). See `PHASE0_novelty_memo.md`.
- **Simulation-recovery gate (Phase 1): PASSED** — the model recovers Se_BC, Se_BM, and per-outbreak
  PPV given the anchor; the single-outbreak ridge, the ">1" plug-in resolution, and a conditional-
  dependence stress test are demonstrated. See `PHASE1_RECOVERY.md`.
- **The severe->mild offset `beta` was DROPPED.** A literature search (Zotero + web: Wain
  quantitation, Andrews&Ryan, TSAP/STRATAA/SETA/Pemba/outpatient cohorts) found: no paired
  blood/bone-marrow data in mild cases exists (bone marrow invasive; community surveillance uses
  serology, not culture) -> `beta` is unidentifiable in principle; and severity is **not** a
  supported independent determinant of Se_BC (the drivers are volume, age, illness-week,
  antibiotics — the sign of a severe->mild gap is even ambiguous). So `Se_BC` transfers across
  settings via the **volume model**, not a severity offset. See `SE_BC_MILD_LITSEARCH.md`.
- **Community PPV analysis:** confirmed PPV varies ~3-fold across settings and cannot be a
  transferable constant; the *level* rides on the Se_BC anchor, the *spread* is data-driven. See
  `COMMUNITY_PPV.md`.

## 8. Caveats and limitations

- **Community PPV level depends on the assumed outbreak blood volume** (the one remaining
  assumption, replacing `beta`). Sensitivity: as assumed volume goes 3 -> 10 mL, `pi_o` shifts
  modestly (Aye 0.28->0.24, Neil 0.70->0.60, Kabwama 0.25->0.21) — the qualitative picture (Aye &
  Kabwama ~0.25, Neil ~0.65) is robust. This is a concrete, researchable quantity (unlike a severity
  offset), but is currently assumed.
- **Only 3 clean community points** (Aye, Neil, Kabwama); community PPV is thin and heterogeneous —
  a wide prior, not a precise estimate.
- **Conditional independence** of blood/bone marrow is assumed; unmodelled positive dependence
  (e.g., prior antibiotics failing both) would over-state Se_BC and under-state PPV (Phase-1 stress
  test: Se +0.08, pi -0.08 at cd=0.06).
- **Historic volumes are imperfect** (Hirsowitz imputed; two studies cultured two volumes, modelled
  at one); the volume slope leans on the Antillon prior.
- **Specificity ~1 assumed** (relaxable as a sensitivity).

## 9. How to use these estimates downstream

- **Absolute burden / DALYs / cost:** scale confirmed counts up by 1/Se_BC(volume) to true
  infections; scale suspected counts by the *community* PPV (~0.25) to true cases. Carry the CrIs.
- **Ratio outcomes (R_t, transmission-mode theta, %-reduction):** the constant scaling cancels
  (invariance) — these estimates change little there; they matter for the absolute arm.
- **Transferable:** the Se_BC volume curve and Se_BM. **Not transferable:** any PPV — re-estimate
  `pi_o` per setting from local culture positivity + the volume-appropriate Se_BC.

## 10. Files
`stan/model_final.stan`, `run_final.R`; outputs (gitignored):
`tables/final_parameters.csv`, `tables/final_phi_hospital.csv`, `tables/final_pi_community.csv`,
`figures/fig_final_se_volume.png`, `figures/fig_final_ppv.png`, `outputs/final_draws.rds`.
Supporting memos: `README.md`, `PHASE0_novelty_memo.md`, `PHASE1_RECOVERY.md`, `COMMUNITY_PPV.md`,
`SE_BC_MILD_LITSEARCH.md`.
