# Phase 0 — Novelty / worthwhileness memo (for sign-off)

**Author of tracked changes:** Claude. **Status:** GATE — no model code written; awaiting your
approval before Phase 1.

## Question
Is the proposed model materially additive over prior work? Proposed contribution:
*joint estimation of per-outbreak surveillance PPV of the clinical case definition (π_o)
together with a volume-anchored, heterogeneity-aware blood-culture sensitivity (Se_BC),
embedded in a transmission/observation model, with an explicit transportable (Se) vs
non-transportable (PPV) parameter split.*

## What the three papers do / do NOT do (from full-text reads)

| Capability | Mogasale 2016 | Arora 2019 | Antillón 2018 |
|---|---|---|---|
| Estimates Se_BC | pooled point only (66% vs BMC; 61% "proportion detected") | **fixes** it as a prior (~50% FN), not estimated | meta-regresses it on volume (not a latent-class posterior) |
| Se_BC ~ blood volume | no (narrative only) | no (volume extracted, not modelled) | **yes** — the core result |
| Estimates the both-tests-negative-but-truly-positive cell | no (union-complete: both-neg = TN by assumption) | resamples culture-negatives at a **fixed** FN fraction | n/a (single-test outcome) |
| PPV of clinical case definition, P(true typhoid \| suspected) | **no** | **no** (explicitly sidesteps prevalence via BMC priors) | no |
| Per-population / per-outbreak observation operator | no | no (region-level test Se/Sp only) | no |
| Conditional dependence (antibiotic-driven joint failure) | no | narrative only, no covariance term | quantifies antibiotic effect but on Se, not as BC/BMC covariance |
| Embedded in transmission / burden model | no | no | no |

**Arora is the closest prior work, and it is a *test-accuracy* paper** (Se/Sp of
RDTs/serology/Widal/PCR via a latent-class network meta-analysis). It **fixes** blood-culture
sensitivity (~50%) and bone marrow (Se 95%/Sp 99%) as informative priors — it does not estimate
Se_BC as a posterior, does not use volume, and produces **no clinical-definition PPV and no
observation operator**. Its own stated reason for not doing PPV: study-level prevalence "was not
available … since all studies only enrolled patients with suspected typhoid."

## Verdict: **materially additive — YES**, with two honest scoping caveats.

Each proposed piece maps to a gap none of the three fills:

1. **Se_BC as a proper latent-class posterior** that (a) relaxes Mogasale's union-complete
   assumption by *estimating* the both-negative-truly-positive cell, and (b) *estimates* Se_BC
   rather than fixing it (Arora). Gap: real.
2. **Volume-anchored, heterogeneity-aware Se_BC** — reuse Antillón's volume relationship as a
   structural prior so each outbreak's Se_BC is set by its blood volume + between-study τ, not by
   the outbreak's own data. Gap: real (nobody folds volume into a PPV/observation estimator).
3. **Per-outbreak surveillance PPV (π_o)** of the clinical case definition, with the
   **transportable Se / local PPV split**, delivered as an **observation operator** for the
   transmission/burden model. Gap: real and central — this is the genuine novelty.

## Two honest caveats (why the novelty is thinner than a first glance, and how to keep it real)

**(a) Conditional on an anchored Se_BC, π_o is arithmetically shallow.** With Sp_BC≈1,
π_o = positivity / Se_BC — the same `AI/AG ÷ sensitivity` we already use for the α confirmation
adjustment. So the methodological value is *not* the PPV formula; it is: honest **uncertainty
propagation** of Se_BC (volume + heterogeneity) into π_o; the **falsification signal** (plug-in
`positivity/0.6 > 1` becomes a proper posterior that keeps π_o in [0,1] and reallocates mass to
higher Se_BC); and the **transportability logic** made explicit in the prior structure. Pitch the
contribution as the *integrated, uncertainty-honest observation model with a transportable/local
split*, **not** the latent-class machinery itself (which is standard: Hui–Walter, Menten–Lesaffre,
and used by Arora).

**(b) Volume is necessary but not sufficient — antibiotics/duration/age dominate.** Antillón
shows Se_BC is depressed **34% by prior antimicrobials** (up to 73% after 5 days), **31% by
sampling after week 1**, and is **higher in adults (0.74)**. These often differ between the
historic paired studies and outbreak line lists (and are usually **unreported** in outbreaks).
**This is the main threat to Se_BC transportability** — more than volume. The plan should
therefore (i) treat the volume model as *one* covariate on Se_BC with a healthy between-study τ
absorbing the rest, and (ii) headline "unreported antibiotic pretreatment in outbreaks" as the
key sensitivity analysis. The optional conditional-dependence term addresses only part of this.

## Recommended (narrowed, sharper) scope
Proceed, but frame the contribution precisely:
- **Primary novelty:** a volume-and-heterogeneity-anchored Se_BC **posterior** (estimating the
  both-negative cell, unlike Mogasale/Arora) that transfers to outbreaks to yield **per-outbreak
  π_o with propagated uncertainty**, as a plug-in observation operator — the transportable-Se /
  local-PPV split being the organizing idea.
- **Do not oversell** per-outbreak PPV as a deep estimation problem (it is positivity/Se given the
  anchor); sell the uncertainty honesty + the ">1" resolution + burden embedding.
- **Add antibiotic/duration covariates** (or explicit τ + a sensitivity analysis) to Se_BC — a
  volume-only anchor is not defensible given Antillón's own effect sizes.
- Keep the transmission-embedding as the **Phase-2 optional** extension where the anchored Se_BC
  actually pays off (absolute burden/DALYs), since for ratio outcomes (R_t, θ, %-reduction) the
  constant scaling cancels — as already established elsewhere in this repo.

## Seed data extracted (Phase-0 deliverable)
- `data/mogasale2016_paired_bc_bmc_seed.csv` — the 10 paired blood/bone-marrow 2×2 studies
  (a,b,c,d + blood volume). **Flags:** Vallenas 1985 has a+b+c+d = 43 but "tested-both" = 58 →
  d likely 15, not 0 (verify at source). Gasem 1995 and Wain 2008 each cultured **two volumes**
  from the same patients (3 & 10 mL; 5 & 15 mL) — needs a per-arm split or a chosen volume.
- `data/antillon2018_volume_sensitivity.csv` — 25 per-study (volume, Se) rows for the volume
  prior + cross-check against Mogasale's 10.
- **Volume model for the prior:** Antillón's best fit is a log-linear slope-and-intercept
  meta-regression of Se on volume (lowest WAIC); linear check ≈ **+3% sensitivity per mL
  (95% CI 1–6%)** over 1–10 mL; predicted Se 0.51 (2 mL) → 0.56 (5 mL) → 0.65 (10 mL); pooled Se
  0.59 (0.54–0.64); I² = 76%. <!-- MARGIN COMMENT (Claude): exact log-linear intercept a and
  slope b live in Antillón Supplementary Table S2, not in the main-text PDF — PLACEHOLDER, obtain
  before building the volume prior; use the +3%/mL linear value as an interim. -->

## Decision requested
Approve to proceed to **Phase 1** (model build + simulation-recovery gate) under the narrowed
framing above? Specifically confirm: (1) include antibiotic/duration as Se_BC covariates (or τ +
sensitivity)? (2) is the "integrated observation model with transportable/local split" framing the
one you want, rather than claiming novel PPV estimation? (3) resolve the two seed-data flags
(Vallenas d; two-volume studies) now or in Phase 2?
