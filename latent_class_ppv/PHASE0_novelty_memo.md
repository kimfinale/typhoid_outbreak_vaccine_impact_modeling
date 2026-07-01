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
   assumption by *estimating* the both-negative-truly-positive cell, (b) *estimates* Se_BC rather
   than fixing it (Arora), and (c) treats **bone marrow as imperfect** (estimates Se_BM). This last
   point is directly motivated by Antillón, who flags BMC is only ~90% sensitive (<1 mL in 14
   studies) and that their composite "may have missed cases, so BC sensitivity may be
   overestimated." Estimating Se_BM resolves the three-way denominator ambiguity that yields
   66% (vs BMC) / 61% (vs BC∪BMC) / 59% (vs any-site). Gap: real.
2. **Volume-anchored, heterogeneity-aware Se_BC** — reuse Antillón's volume relationship as a
   structural prior so each outbreak's Se_BC is set by its blood volume + between-study τ, not by
   the outbreak's own data. Gap: real (nobody folds volume into a PPV/observation estimator).
3. **Per-outbreak surveillance PPV (π_o)** of the clinical case definition, with the
   **transportable Se / local PPV split**, delivered as an **observation operator** for the
   transmission/burden model. Note Arora *explicitly could not* use prevalence ("all studies only
   enrolled suspected typhoid"); our φ_s (prevalence-among-suspected = hospital PPV) is exactly that
   quantity, estimated freely from the paired 2×2 — we do the thing Arora sidestepped. Gap: real and
   central — this is the genuine novelty.

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
historic paired studies and outbreak line lists (and are usually **unreported** in outbreaks). The
plan should therefore (i) treat the volume model as *one* covariate on Se_BC with a healthy
between-study τ absorbing the rest, and (ii) headline "unreported antibiotic pretreatment in
outbreaks" as a key sensitivity analysis. The optional conditional-dependence term addresses only
part of this.

**(c) The severity/bacteremia gradient — the transfer's biggest threat (from reading Antillón
directly).** **38/40 paired BC/BMC studies are hospitalized** patients (severe, high-bacteremia),
whereas only **1–20% of community-surveillance cases are hospitalized**. So the anchored Se_BC is a
**hospital** sensitivity; the outbreak suspected-case population is milder with lower bacteremia →
its true Se_BC is plausibly **lower**. Since π_o = positivity / Se_BC, transferring the (higher)
hospital Se **over-estimates** outbreak Se and therefore **under-estimates** π_o (a conservative
bias, but a bias). Volume and τ do **not** absorb this — it is a population-setting shift, not
between-study noise, and there is **no community paired BC/BMC data** to estimate it directly. This
is the honest headline limitation and should drive an explicit sensitivity analysis (a
hospital→community Se offset with a wide prior), not be buried.

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
  (a,b,c,d + blood volume), verified against Mogasale Table 3. **Flags:** Vallenas 1985 — the
  published Table 3 prints d (both-negative) = 0, but that is an **erratum**: tested-both = 58 and
  a+b+c = 43, so d must = 15 (every *other* study reconciles as d = tested_both − true_pos). Seed
  uses d=15. Gasem 1995 and Wain 2008 each cultured **two volumes** from the same patients
  (3 & 10 mL; 5 & 15 mL) — needs a per-arm split or a chosen volume. Also note the historic 2×2 d
  cell (both-negative) is the very cell that identifies φ_s once BMC is imperfect — so the Vallenas
  correction is not cosmetic.
- `data/antillon2018_volume_sensitivity.csv` — 25 per-study (volume, Se) rows for the volume
  prior + cross-check against Mogasale's 10.
- **Volume model for the prior:** Antillón's best fit is a log-linear slope-and-intercept
  meta-regression of Se on volume (lowest WAIC); linear check ≈ **+3% sensitivity per mL
  (95% CI 1–6%)** over 1–10 mL; predicted Se 0.51 (2 mL) → 0.56 (5 mL) → 0.65 (10 mL); pooled Se
  0.59 (0.54–0.64); I² = 76%. <!-- MARGIN COMMENT (Claude): exact log-linear intercept a and
  slope b live in Antillón Supplementary Table S2, not in the main-text PDF — PLACEHOLDER, obtain
  before building the volume prior; use the +3%/mL linear value as an interim. -->

## Decision requested (adjusted after reading the three papers in full)
Approve to proceed to **Phase 1** (model build + simulation-recovery gate) under the narrowed
framing above? Four specific choices, in priority order:

1. **The hospital→community Se_BC transfer (caveat c) is the crux — how do we handle it?**
   Options: (a) transfer Se_BC via the volume model + τ only, and report the severity/bacteremia
   gradient as a one-directional bias (π_o under-estimated) — simplest, honest; (b) add an explicit
   hospital→community Se offset with a wide prior and run it as the primary sensitivity analysis
   (recommended); (c) hunt for any community-based paired BC/BMC (or BC-vs-PCR) data to anchor the
   offset (likely none exists — Phase-2 search). My recommendation: **(b)**.
2. **Se_BC covariates:** include antibiotic-pretreatment (± symptom-duration) as covariates, or
   just a wide between-study τ + a sensitivity analysis? (Antillón gives the effect sizes either
   way; outbreaks rarely report these, so the covariate is often latent.)
3. **Framing:** confirm the contribution is pitched as the *integrated, uncertainty-honest
   observation model with a transportable-Se / local-PPV split* (estimating φ_s and Se_BM that
   Arora/Mogasale fixed or assumed away) — **not** as novel PPV arithmetic.
4. **Seed-data:** the Vallenas d=15 correction and the two-volume studies (Gasem 1995, Wain 2008) —
   lock these in now (my proposed values), or defer per-arm splits to Phase 2?
