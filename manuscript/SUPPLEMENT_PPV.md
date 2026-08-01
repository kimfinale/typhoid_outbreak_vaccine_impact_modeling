# Supplement S: Positive predictive value of clinical typhoid case definitions

> Draft. Every number below is read from a fitted output, not from prose. Sources:
> `latent_class_ppv/tables/{final_parameters,final_pi_community,final_phi_hospital,fit_diagnostics,volume_sensitivity,community_ppv_strata,community_ppv_strata_params}.csv`.
> Models: `latent_class_ppv/stan/{model_final,community_ppv_strata}.stan`.
> Runners: `latent_class_ppv/run_final.R`, `manuscript/fit_ppv_strata.R`.
> **⚑** = open decision. Fit of record: 4 chains, 6000 post-warmup draws, seed 2026.

## S1. Rationale

Outbreak reports count **suspected** cases: patients meeting a syndromic definition
(typically fever ≥3 days plus one or more of headache, abdominal pain, malaise, diarrhoea
or constipation). Vaccine impact, however, accrues only to **true typhoid**. The two differ
by the positive predictive value (PPV) of the case definition, which is not a property of
the definition alone: it depends on how many of the febrile patients presenting in a given
setting actually have typhoid. In a hospital ward receiving referred, severely ill patients
the PPV is high; in community surveillance during an outbreak, where most fevers are
malaria, respiratory or undifferentiated illness, it is low.

PPV is therefore **local** and must be estimated per setting, whereas the diagnostic
accuracy of the tests used to measure it is **transportable**. This supplement estimates
both, and propagates the local quantity into the outbreak-response model.

A note on terminology, because the literature conflates two quantities. Studies commonly
report the fraction of suspected cases that are blood-culture positive and call it PPV.
That is **positivity**, and it equals `π · Se_BC` — the true PPV attenuated by an imperfect
reference test. Since blood culture misses roughly 40% of true typhoid (S7), positivity
systematically understates PPV. We estimate π; we do not treat positivity as π.

## S2. Evidence base and selection

Two targeted syntheses, reported separately because they answer different questions.

### S2.1 Paired blood + bone-marrow studies (test accuracy)

Historic studies culturing **both** blood and bone marrow in the same patients, seeded from
Mogasale et al. (2016). Bone marrow is substantially more sensitive than blood, so the
paired 2×2 identifies blood-culture sensitivity without a gold standard (S5).

**10 studies, 635 patients** (`mogasale2016_paired_bc_bmc_seed.csv`):

| Study | Location | Blood (mL) | BC+BM+ | BC+BM− | BC−BM+ | BC−BM− | N |
|---|---|---|---|---|---|---|---|
| Akoh 1991 | Zaria, Nigeria | 2 | 8 | 3 | 11 | 9 | 31 |
| Hirsowitz 1951 | Transvaal, S. Africa | ⚑ NR | 13 | 0 | 5 | 10 | 28 |
| Farooqui 1991 | Karachi, Pakistan | 5 | 58 | 0 | 30 | 12 | 100 |
| Gasem 1995 | Semarang, Indonesia | ⚑ 3 (also 10) | 39 | 13 | 28 | 0 | 80 |
| Gasem 2002 | Semarang, Indonesia | 9 | 43 | 1 | 10 | 7 | 61 |
| Chaicumpa 1992 | Jakarta, Indonesia | 3 | 17 | 5 | 15 | 15 | 52 |
| Wain 2008 | Dong Thap, Vietnam | ⚑ 5 (also 15) | 53 | 4 | 16 | 30 | 103 |
| Gilman 1975 | Mexico City | 2 | 24 | 1 | 32 | 5 | 62 |
| Guerra-Caceres 1979 | Lima, Peru | 5 | 26 | 0 | 31 | 3 | 60 |
| Vallenas 1985 | Lima, Peru | 3 | 9 | 10 | 24 | ⚑ 15 | 58 |

All counts are transcribed from **Mogasale et al. (2016) Table 3**; the primary papers were
not re-read. Verification against that table (2026-07) resolved one flag and sharpened two.

**Vallenas 1985 — `d = 15` is correct; Mogasale's printed `d = 0` is a typo in the source.**
Provable from their own arithmetic: (i) every other study's cells sum exactly to its
"tested by both" N, while Vallenas gives 9+10+24+0 = 43 against 58 tested; (ii) their totals
row prints Σd = 91, but they report 635 tested and 529 detected (= Σ(a+b+c)), requiring
Σd = 106; (iii) 106 − 91 = 15, exactly the Vallenas shortfall. The error is invisible in
their analysis because `ProCase_BC = (a+b)/(a+b+c)` never uses `d` — but `d` identifies **φ**
here, so it matters to us. Primary source for a belt-and-braces check: *Pediatr Infect Dis*
1985;4:496–8.

**Gasem 1995 (B: "3 and 10 ml") and Wain 2008 (B: "5 and 15 ml") — cells are pooled across
both volumes in the source and cannot be split without the primary papers.** Each currently
enters at the *lower* volume. **This biases α₁ downward with a known sign:** both studies
have among the highest raw sensitivities in the set (0.65, 0.78), so attributing high
sensitivity to low volume flattens the slope. Mogasale note that Wain's 15 mL arm reached
sensitivity "comparable to that of the bone marrow culture". Primaries: *Trop Geogr Med*
1995;47:164–7; *J Infect Dev Ctries* 2008;2:469–74.

**Hirsowitz 1951** — volume recorded as "NP" (not presented) in the source; imputed at the
sample median (3 mL). The primary is a one-page report (*BMJ* 1951;1:862) and may never
state it.

**⚑ Culture media is an unmodelled confounder, plausibly larger than the volume effect.**
Mogasale report Luria-Bertani broth recovering 44/54 (0.81) against Difco Oxgall/Selenite F
19/43 (0.44) — a wider spread than our entire fitted volume range (0.52 → 0.68). Those are
Gasem 2002 (LB, 9 mL; our highest-sensitivity study) and Vallenas (oxgall, 3 mL; among our
lowest). **Volume and media are confounded across these ten studies, and α₁ may be partly
capturing a media effect.** With 10 studies both cannot be separated; this should be stated
as a limitation rather than modelled.

### S2.2 Community and surveillance positivity (local PPV)

Studies reporting **both** a blood-culture-tested denominator and a confirmed numerator
among syndromically-defined suspected cases. This denominator is the binding constraint:
most outbreak reports do not publish it, which is why the community evidence base is small.

**Inclusion** required (i) a syndromic case definition applied at presentation, (ii) an
explicit count blood-cultured, (iii) an explicit count culture-confirmed, and (iv) ≥20
tested.

**Exclusion on principled grounds.** Under the measurement model, positivity is
`θ = π·Se_BC + (1−π)(1−Sp_BC)`, which for π ≤ 1 cannot exceed ≈ Se_BC (≈ 0.62 at 5 mL).
An observed positivity at or above Se_BC is therefore **unattainable for any admissible π**,
and indicates that culture was ordered selectively on patients already likely to be positive
— violating the assumption that the tested are a random subset of the suspected. Such
studies are excluded: they measure a confirmed-case series, not a case-definition PPV.
Excluded on this basis: **Ali 2017** (0.85), **Nimonkar 2022** (1.00), **Srinivasan 2022**
(1.00), **Wang 2022** (1.00).

> ⚑ The rule is currently implemented as `positivity ≥ 0.8`, a rough proxy. It should be
> restated as `positivity ≥ Se_BC`, which is the actual identifying constraint.

⚑ **Misra 2005** (150 tested / 93 confirmed; positivity 0.62) is extracted as usable but
enters no fit. Recommended disposition: **exclude from the community stratum** and, if used
at all, treat as a hospital-ascertainment (φ) point — it is hospital-referred with a highly
restrictive definition (fever ≥8 days), and its positivity ≈ Se_BC implies π ≈ 1, consistent
with that reading. *Decision pending.*

**Stage-1 community set** — restricted to studies with a *community syndromic* definition:

| Study | Setting | Tested | Confirmed | Positivity |
|---|---|---|---|---|
| Aye 2004 | Madaya, Myanmar | 22 | 3 | 0.136 |
| Neil 2012 | Kasese, Uganda | 54 | 19 | 0.352 |
| Kabwama 2017 | Kampala, Uganda | 364 | 56 | 0.154 |

⚑ **Lewis 2005** (Nepal; 104/39) and **Roy 2016** (India; 13/6) are usable but excluded from
stage 1 as not community-syndromic; Roy also fails the ≥20 rule. Lewis is a real excluded
data point and the reason should be stated in the paper, not only in the code.

> **Provenance note.** Neil 2012 is **19/54 blood cultures**. An earlier version of this
> file recorded 63/25, which was the blood cultures (19/54) summed with 6/9 cultures from an
> unlabelled source (blood *or* stool). Because `Se_BC` is specifically blood-culture
> sensitivity, only the blood-culture counts are admissible. Kabwama is 56/364. Both were
> corrected before the fit of record; see S8 for the effect.

### S2.3 Strata set (stage 2)

Stage 2 adds endemic-surveillance studies that culture **all** suspected patients, plus
N'Cho 2019, giving 10 studies across two strata (`community_ppv_strata.csv`):

*Outbreak (selective testing):* Aye 2004 (22/3), Kabwama 2017 (364/56), N'Cho 2019
(286/74), Neil 2012 (54/19).
*Endemic surveillance (all-tested):* Aiemjoy/SEAP 2020 (20 899/2116), Andrews 2018
(4309/87), Meiring 2021 STRATAA (12 082/624), Srikantiah 2006 (1815/90), Thriemer 2012
(2209/46), Voysey 2020 (1615/32).

Two points on this set:

- **SEAP enters as counts, not as its reported PPV.** Aiemjoy et al. report ~10.5% PPV, but
  that figure is blood-culture-referenced — it is positivity (`π·Se_BC`), not π. Entering it
  as π would double-count the reference-test attenuation.
- **N'Cho 2019** contributes the hospital case-management subset (286 tested / 74 confirmed,
  Oct–Dec 2017). The full outbreak (3187 suspected / 191 confirmed) lacks a tested
  denominator and cannot be used.
- ⚑ Strata are assigned in `fit_ppv_strata.R` by **study type**, which overrides the
  `stratum` column in `ppv_positivity_table.csv`. That column is stale and misassigns
  Voysey 2020 to the outbreak stratum. It should be corrected or deleted before anyone
  documents from it.

## S3. Model

### Stage 1 — accuracy and local PPV (`model_final.stan`)

*Paired layer* (conditional-independence latent class; Hui–Walter). For study `s` with
hospital PPV `φ_s`, cells (BC+BM+, BC+BM−, BC−BM+, BC−BM−) ~ Multinomial(N_s, p) with

```
p₁ = φ·Se_BC·Se_BM        + (1−φ)·(1−Sp_BC)·(1−Sp_BM)
p₂ = φ·Se_BC·(1−Se_BM)    + (1−φ)·(1−Sp_BC)·Sp_BM
p₃ = φ·(1−Se_BC)·Se_BM    + (1−φ)·Sp_BC·(1−Sp_BM)
p₄ = φ·(1−Se_BC)·(1−Se_BM)+ (1−φ)·Sp_BC·Sp_BM
```

*Volume sub-model.* `logit Se_BC(v) = α₀ + α₁·log(v) + τ·u`, `u ~ N(0,1)`.

*Community layer.* `k_o ~ Binomial(n_o, π_o·Se_BC(v_o) + (1−π_o)(1−Sp_BC))`.

*Transportable:* α₀, α₁, τ, Se_BM, Sp. *Local:* φ_s (hospital), π_o (community).

### Stage 2 — endemicity gradient (`community_ppv_strata.stan`)

```
logit π_o = μ + β·outbreak_o + σ·z_o ,  z_o ~ N(0,1)
k_o ~ Binomial(n_o, π_o · Se_BC)
```

`μ` = mean logit PPV under endemic surveillance; `μ+β` = under outbreak conditions;
`exp(β)` = the endemicity odds ratio. Generated: `π_surv`, `π_out`, `OR`, and the
posterior-predictive `π_new_out`, `π_new_surv` for an **unobserved** setting.

**Why two stages rather than one joint fit.** Stage 2 takes `Se_BC` as an informative prior
instead of re-estimating it, i.e. cut inference. This is deliberate: it prevents the
outbreak positivity data — which are consistent with a continuum of (π, Se_BC) pairs — from
feeding back and distorting an accuracy parameter that the paired studies identify properly.
The prior is calibrated to stage 1: `Beta(24,15)` has mean 0.615 and sd 0.077 against stage
1's posterior median 0.617 (sd ≈ 0.047), i.e. correctly centred and slightly conservative.
A fully joint model is the natural extension.

## S4. Priors

| Parameter | Prior | Basis |
|---|---|---|
| α₀ | N(0, 1.5) | weakly informative on logit Se_BC |
| **α₁** | **N(0.36, 0.15)** | **Antillón 2018**: 0.51 @2 mL → 0.65 @10 mL ⇒ Δlogit/Δlog v = 0.36 |
| τ | half-N(0, 1) | between-study heterogeneity |
| Se_BM | Beta(12, 2) — mean 0.857 | bone marrow ≫ blood; **breaks label symmetry** (S5) |
| Sp_BC, Sp_BM | Beta(200, 1) ≈ 0.995 | culture is near-perfectly specific |
| φ_s, π_o | Beta(1, 1) | flat; PPV is local and not pooled |
| *Stage 2:* μ, β | N(0, 1.5) | weakly informative |
| *Stage 2:* σ | half-N(0, 1) | between-study heterogeneity |
| *Stage 2:* Se_BC | Beta(24, 15) | calibrated to stage 1 (S3) |

## S5. Identification

Positivity data alone identify only the **product** `π · Se_BC`: any low-PPV/high-sensitivity
combination is observationally equivalent to a high-PPV/low-sensitivity one. This ridge is
the central obstacle, and two features break it:

1. **The paired bone-marrow layer.** With two tests of differing sensitivity applied to the
   same patients, the discordant cells (BC−BM+ vs BC+BM−) carry information about Se_BC
   separately from prevalence — the Hui–Walter argument. The `Se_BM ~ Beta(12,2)` prior
   resolves the label-switching symmetry (the likelihood is invariant to swapping
   "diseased"/"healthy" with complemented accuracies).
2. **The volume prior on α₁.** This anchors the *level* of the Se_BC curve using external
   evidence.

**This is the model's principal dependency and it should be stated plainly.** α₁'s posterior
(0.420, 90% CrI 0.185–0.649; sd ≈ 0.141) is barely narrower than its prior (sd 0.15) — the
10 paired studies reduce the variance by only ≈ 12%. The volume slope is therefore
**substantially imported from Antillón rather than re-estimated here**, and with it a share
of the identification of Se_BC's level, and hence of every π. S8 quantifies the consequence.

## S6. Estimation and diagnostics

Stan/cmdstanr, 4 chains × (1500 warmup + 1500 sampling) = 6000 post-warmup draws,
`adapt_delta = 0.98`, `max_treedepth = 12`, seed 2026.

**Fit of record: 0 divergent transitions, max R̂ = 1.003, min bulk ESS = 1670, min tail
ESS = 1756.** (`fit_diagnostics.csv`)

## S7. Results

### S7.1 Test accuracy (transportable)

| Parameter | Median | 90% CrI |
|---|---|---|
| Se_BC @ 2 mL | 0.52 | 0.44 – 0.61 |
| **Se_BC @ 5 mL** | **0.62** | **0.53 – 0.69** |
| Se_BC @ 10 mL | 0.68 | 0.59 – 0.76 |
| Se_BM | 0.90 | 0.88 – 0.92 |
| Sp_BC | 0.995 | 0.977 – 1.000 |
| α₀ | −0.197 | −0.620 – 0.261 |
| α₁ (volume slope) | 0.420 | 0.185 – 0.649 |
| τ (heterogeneity) | 0.484 | 0.247 – 0.861 |

Blood culture at a routine 5 mL draw detects **≈ 62%** of true typhoid; bone marrow ≈ 90%.
α₁ excludes zero: sensitivity rises with volume, ≈ 0.52 → 0.68 from 2 to 10 mL.

**Coherence with the source review.** Mogasale et al. analyse these same ten studies with a
*composite reference* (positive on blood **or** bone marrow), obtaining `ProCase_BC` = 0.61
(95% CI 0.52–0.70) and `ProCase_BMC` = 0.96 (0.93–0.99). This is **not** external validation
— it is the same data under a different estimator — and the informative contrast is
`ProCase_BMC` **0.96 vs our Se_BM 0.90**. A composite reference misses cases positive on
neither test and therefore *overstates* both sensitivities, as Mogasale concede
("the composite reference standard may have still missed some typhoid fever cases … the
proportion of S. Typhi detected could be an overestimation"). The latent-class model, which
estimates against latent true status rather than a composite, corrects downward — in the
predicted direction and by a plausible magnitude. That agreement of *sign and size* is the
check; numerical similarity of 0.61 and 0.617 is not, and the two are different estimands.

### S7.2 PPV by ascertainment (local)

**Hospital (φ_s), paired studies:** median 0.83, range 0.66 (Hirsowitz 1951) to 0.99
(Gasem 1995).

**Community outbreak (π_o), stage 1:**

| Study | Positivity | π median | 90% CrI |
|---|---|---|---|
| Aye 2004 (Madaya) | 0.136 | 0.26 | 0.08 – 0.61 |
| Neil 2012 (Kasese) | 0.352 | 0.58 | 0.38 – 0.88 |
| Kabwama 2017 (Kampala) | 0.154 | 0.25 | 0.17 – 0.45 |

The gap between hospital (0.83) and community (0.25–0.58) is the quantitative statement of
the opening argument: **the same clinical definition means different things in different
places.** Note that every π exceeds its positivity, as it must — positivity is π attenuated
by Se_BC.

### S7.3 Endemicity gradient (stage 2)

| Quantity | Median | 90% CrI |
|---|---|---|
| π, endemic surveillance | 0.075 | 0.041 – 0.166 |
| π, outbreak | 0.338 | 0.159 – 0.587 |
| **OR, outbreak vs surveillance** | **6.40** | **1.74 – 16.85** |
| σ (between-study) | 0.779 | 0.473 – 1.466 |
| Se_BC (prior-informed) | 0.598 | 0.438 – 0.743 |
| **π_new_out** (predictive, new outbreak) | **0.338** | **0.063 – 0.806** |
| π_new_surv (predictive, new surveillance) | 0.074 | 0.013 – 0.395 |

A syndromic definition applied during an outbreak is **≈ 6 times more likely** (odds) to
identify true typhoid than the same definition applied in routine endemic surveillance —
because an outbreak raises typhoid's share of all fevers. The CrI excludes 1. This is the
result that makes PPV setting-dependent rather than a constant.

**Internal consistency.** Stage 1 and stage 2 estimate π for three shared studies by
different routes (Se_BC estimated vs prior-informed): Aye 0.262 vs 0.278, Neil 0.584 vs
0.537, Kabwama 0.252 vs 0.261. The agreement is a check on the cut-inference design.

## S8. Sensitivity

**Assumed outbreak blood volume** (the single strongest assumption in stage 1; outbreak
reports do not state draw volume). Refit across 3–10 mL (`volume_sensitivity.csv`):

| Assumed volume | Se_BC | π Aye | π Neil | π Kabwama |
|---|---|---|---|---|
| 3 mL | 0.57 | 0.29 | 0.63 | 0.28 |
| **5 mL (base)** | **0.61** | **0.27** | **0.59** | **0.25** |
| 8 mL | 0.66 | 0.25 | 0.56 | 0.24 |
| 10 mL | 0.68 | 0.24 | 0.54 | 0.23 |

π moves inversely to assumed volume (a more sensitive test implies fewer missed cases and
hence a lower PPV to explain the same positivity), but the range is modest — Kabwama
0.23–0.28 across a threefold volume span. **Conclusions are not driven by this assumption.**

**α₁ prior dependence.** Given S5, the honest sensitivity is to α₁'s prior, not only to the
assumed volume. ⚑ *Not yet run: refit under α₁ ~ N(0.36, 0.30) and α₁ ~ N(0, 0.5) to show
what the paired data alone support.* This is the analysis a methods reviewer will ask for.

**Omitted false-positive term in stage 2.** Stage 2 fits `k ~ Binomial(n, π·Se_BC)`, dropping
`(1−π)(1−Sp)`. With Sp = 0.995 this is negligible at high π but not at low π: it inflates
π_surv by ≈ 12% and π_out by ≈ 1.6%. Correcting it would move the OR from 6.40 to ≈ 6.9 —
i.e. **the reported gradient is conservative.**

**Data correction.** Refitting with corrected counts (S2.2) moved π_Neil 0.647 → 0.584 and
π_Kabwama 0.227 → 0.252 — nearly offsetting. The ORI anchor distribution's central π moved
0.366 → 0.356 (−3%), with the spread narrowing (σ 1.01 → 0.81). **Downstream ORI estimates
are materially unchanged.**

## S9. Limitations

1. **Three community anchors.** π for outbreak settings rests on Aye, Neil and Kabwama.
   The 40-PDF extraction added no further community-outbreak π: outbreak reports
   overwhelmingly lack a blood-culture-tested denominator. This is the binding evidence
   constraint, not a modelling choice.
2. **The volume slope is largely imported** (S5). Se_BC's level, and hence π, depends on
   Antillón's external evidence.
3. **Paired studies are old** (1951–2008) and mostly hospital-based. Se_BC is assumed
   transportable to contemporary community outbreaks; only volume modulates it. Blood
   culture technique, prior antibiotic use, and bacterial load may differ.
4. **Prior antimicrobial use is not modelled** and plausibly lowers Se_BC in outbreak
   settings — which would bias π *downward*, making the reported PPVs conservative.
5. **Conditional independence** of blood and bone marrow given true status is assumed, and
   **the source review explicitly warns it may fail** — Mogasale note that bone marrow "might
   have been collected in a later stage … or after antibiotic prescription", and that
   "those who were tested negative for blood culture might have been sampled for bone marrow,
   which would create a blood culture sensitivity underestimation bias." If sampling was
   conditional on the blood result, the paired 2×2 is not a random draw and both Se_BC and
   φ are biased. This is the model's least testable assumption.
6. **Volume is confounded with culture media and era** (S2.1); α₁ is not a clean volume
   effect.
7. **The strata contrast is confounded with geography and era**: outbreak studies are
   African/Asian community responses; surveillance studies are largely Asian programmes.
   The OR should be read as outbreak-vs-surveillance *ascertainment*, not as a clean
   causal effect of outbreak conditions.

## S10. Propagation into the outbreak-response model

Outbreak reports give suspected counts `S(t)`. The ORI model requires true typhoid `I(t)`.
We use the **additive** decomposition

```
S(t) = I(t) + B(t)
```

where `B(t)` is the background of non-typhoidal febrile illness meeting the definition,
taken as constant over the outbreak window and calibrated so that `sum(I) / sum(S) = π`.
The multiplicative alternative (`I(t) = π·S(t)`) is rejected: it assumes non-typhoidal fever
scales with the outbreak, whereas background fever is precisely what does *not* respond to
a typhoid epidemic. Implementation: `renewal/R/ppv.R::ppv_true_incidence`.

At a second additive level, true typhoid is separated into background/common-source and
propagated incidence:

```
I(t) = X(t) + P(t)
P(t) = R_P(t) * sum_s w(s) I(t-s)
```

Vaccination reduces both `X(t)` and `P(t)` directly, and averted propagated cases also
reduce subsequent infectiousness. The structural source fraction θ allocates baseline
incidence inside this recursion. It no longer weights two completed outcomes; the earlier
`θ·static + (1−θ)·renewal` calculation is retained only as a historical comparator.

Consequences, verified in the pipeline:

- Case counts scale with π; **deaths do not** (deaths are observed, not derived), so YLD
  scales with π while YLL is π-invariant.
- Under the primary additive observation model, the other-febrile denominator is fixed, so
  PPV need not cancel from percentage reductions. π-invariance applies only to the
  multiplicative sensitivity model.
- In paired base-case simulations, the additive source-plus-propagated model reduced true
  typhoid by 25.6%, versus 25.9% under the historical θ-weighted calculation; the paired
  difference was −0.5 percentage points (95% simulation interval −2.4 to 0.7).
- `renewal/R/ppv.R` now uses fitted `mu_pi` and `sigma_pi` draws to generate the model's
  posterior predictive PPV for unanchored outbreaks; the earlier ad hoc moment construction
  is no longer used.

## Open items

| # | Item | Where |
|---|---|---|
| 1 | Misra 2005 disposition — recommend exclude from community stratum | S2.2 |
| 2 | Verify Vallenas 1985 cell `d` against source | S2.1 |
| 3 | Restate exclusion rule as `positivity ≥ Se_BC` | S2.2 |
| 4 | Fix/delete stale `stratum` column (misassigns Voysey) | S2.3 |
| 5 | Run α₁ prior-sensitivity refit | S8 |
| 6 | ~~Switch `ppv.R` to fitted posterior-predictive PPV~~ — completed | S10 |
| 7 | State Lewis 2005 exclusion in the paper | S2.2 |
| 8 | Add `(1−π)(1−Sp)` to stage 2, or report the ≈6.9 OR as a sensitivity | S8 |
