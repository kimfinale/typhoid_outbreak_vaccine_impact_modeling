# Case definitions as observation operators

Extracts each outbreak study's **case definition**, classifies it into an
**observation-process type (A–F)**, and designs how each should be fitted — pooled,
harmonized to a common true-case basis, stratified, or excluded — so heterogeneous
definitions are handled coherently. Foundational: it determines which outbreaks enter
which pool in every downstream analysis.

## The observation-operator framing
One latent infection process λ_o(t) per outbreak; each study observes a *different
transform*:

```
observed_{o,t} = rho_o · ∫ λ_o(t − u) · d_o(u) du
```

The case definition sets both the ascertainment scale **rho_o** and the delay
distribution **d_o**. Mishandling it biases (a) absolute/cost outcomes (different rho_o),
(b) the **percent-reduction invariance** — which holds only if **rho_o is time-constant**,
(c) θ / renewal (different d_o distorts curve shape), and (d) DALYs (severity mix).

## Types
A suspected/clinical community · B blood-culture-confirmed (counts are confirmed) ·
C mixed (suspected + confirmed subset) · D perforation/surgical series ·
E Widal-serological · F ambiguous/not stated.

## Key findings (and how they change the pipeline)

- **Muyembe-Tamfum 2009 = Type D (peritonitis/surgical case series)** — 144 perforation
  admissions, the severe ~1–3% tip of the iceberg, long severity-filtered delay. It is
  currently in the **renewal-13**, so its R_t/θ are built on a *perforation* curve, not an
  infection curve. **Recommendation: exclude from the renewal / θ / %-reduction pools (a
  sensitivity case); keep it only for a perforation-mortality (YLL) subgroup** analyzed
  together with Walters-Kasese.
- **Neil 2012 = Type A surveillance, NOT a perforation series** (perforation is one optional
  qualifier). This corrects the prompt's premise. But its **ascertainment changed
  mid-outbreak** (retrospective IP-skewed → prospective enhanced surveillance), so
  **invariance is violated** and the early curve is IP-distorted → include with caution,
  flag the change, consider trimming the IP-skewed early phase.
- **The actual second perforation-biased series is Walters-Kasese** (de facto IP-restricted
  retrospectively, 82% IP). The perforation subgroup is **Muyembe + Walters-Kasese**.
- **Qamar 2018 & Yousafzai 2019 = Type B (culture-confirmed-only counts)** — a *different
  rho scale*; to harmonize their absolute/cost outcomes to suspected-based outbreaks, scale
  *up* by 1/confirmation-fraction (the **opposite direction** to the α suspected→true
  adjustment). Their curves are still usable for %-reduction/θ if culture effort is constant.
- **Ascertainment changed mid-outbreak** for **Neil, Ali, Lutterloh** → percent-reduction
  invariance at risk for those. Imanishi only reduced the *culturing fraction* (not the
  suspected-case definition) → invariance holds.
- **Davis 2018 = Type C** — classified from the real MMWR (the Zotero copy was misfiled):
  780 suspected / 80 confirmed (10.3%), blood/stool culture, city surveillance, no severity
  restriction.

## A and C are fitted identically (fit the suspected series)

Type A (suspected only reported) and Type C (suspected + a confirmed subset) get the **same
fitting treatment: fit the *suspected* series.** The confirmed subset does **not** calibrate
the ascertainment scale rho, because **culture is performed on a *selected* subset of
suspected cases** — so confirmed:suspected conflates culture sensitivity with the
(non-random) decision of whom to test. rho (PPV) instead comes from the α adjustment using
**culture positivity *among the tested* (AI/AG ÷ sensitivity)**, which conditions on testing
and removes the selection factor. So A and C collapse to one "suspected-series" class; only
B (confirmed-only counts) and D (perforation) differ — they have no comparable suspected series.

## Poolability summary (the 13 renewal outbreaks)
- **Suspected-series pool (%-reduction & θ; ρ via α):** Lewis, Polonsky×2, Muti (A) +
  Aye, Kabwama, N'Cho, Davis (C) — fitted identically on the suspected curve.
- **Type B — dynamics OK; harmonize absolute to true-case (opposite of α):** Qamar, Yousafzai.
- **Include with caution (invariance violated):** Neil, Ali.
- **Exclude from renewal/θ (Type D perforation series):** Muyembe-Tamfum 2009.

## Linkages (do not re-run those analyses here)
- **Confirmation adjustment (α):** applies to the suspected-series classes A **and** C
  (suspected→true, using culture positivity AI/AG). Type B needs the *inverse* scaling; Type
  D needs the perforation-rate scale, not α.
- **θ / transmission-mode:** drop Muyembe (perforation curve); caution on Neil (ascertainment).
- **DALY/cost:** Type D carries real mortality (YLL via perforation); Types A/C/B morbidity-dominated.

## Files
`tab_case_definitions.csv` (Steps 1–2), `tab_fitting_design.csv` (Step 3),
`report_case_definitions.qmd`, `run_case_definitions.R`. Run:
`Rscript case_definitions/run_case_definitions.R`. Every definition cited; "not clearly
stated" recorded where vague (Hendriksen, Davis); no definition imputed.
