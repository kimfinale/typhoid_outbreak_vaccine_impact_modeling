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

**(c) The transfer's biggest threat is a BACTERIAL-LOAD composition shift (not a care-setting
label).** Se_BC is mechanistically a function of blood bacterial load — typhoid bacteremia is sparse
(Wain 1998: median **1 CFU/mL**, range <0.3–387), which is *why* volume matters (+3%/mL) and why
Se falls with duration and antibiotics. Load's **measurable** drivers: **volume** (have it),
**symptom duration** (load declines with time, P=0.002 → −31% after wk 1), **age** (children higher,
1.5 vs 0.6 CFU/mL, P=0.008), **antibiotics** (−34%). Crucially, load does **not** map cleanly onto
"severity" or "hospital": Wain ties it to duration/age/shedding/resistance, and children carry
*higher* load yet Antillón found *adults* have higher observed Se (volume confounds). So
"hospital vs community" and even "severe vs mild" are **noisy proxies** for the real axis, which is
bacterial load.

The transfer threat, restated: the 38/40 hospitalized paired studies sample a **load
composition** (illness-stage × age × volume mix) that differs from an outbreak's suspected-case
line list; if the surveillance population sits at lower load, its true Se_BC is lower, so
transferring the paired Se **over-estimates** Se and **under-estimates** π_o. **Plan:** (i) model
Se_BC on the *measurable* load drivers — volume as the in-hand covariate, plus duration/antibiotic
where reported; (ii) represent the residual paired→surveillance **load-composition** difference as
an **explicit offset with a wide prior**, run as the primary sensitivity analysis (per your
steer). There is **no community paired BC/BMC data** to pin this offset, so it stays a stated
structural assumption, not a fitted quantity.

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

## Resolved direction (per user) + concrete Phase-1 plan

**Resolved:** the transfer offset is **bacterial-load-based**, not care-setting-based. First
establish Se_BC ~ load (via volume + severity/duration), then add the explicit load-composition
offset. Mechanism confirmed (Wain 1998 + Antillón), with the honest nuance that load ≠ clean
"severity" (children higher-load; volume confounds) — so the offset is a **load-composition** shift.

**Phase-1 Se_BC sub-model (proposed):**
```
logit(Se_BC,s) = alpha0 + alpha1*g(volume_s) + beta*load_shift_s + u_s,   u_s ~ N(0, tau^2)
  g(volume)       log-linear (Antillón best-fit; +3%/mL linear check as interim)
  load_shift_s    a study-level standardized "load composition" index built from what IS
                  reported (mean/median symptom duration; age mix; inpatient fraction), = 0 for
                  the hospital paired-study baseline
  beta            effect of a lower-load population on logit(Se); INFORMATIVE-but-wide prior
                  (sign constrained: lower load -> lower Se), the primary sensitivity knob
For an outbreak o:  Se_BC,o = inv_logit(alpha0 + alpha1*g(volume_o) + beta*load_shift_o + u_o)
  with load_shift_o set from the outbreak's reported severity/duration mix (or a prior if
  unreported), u_o ~ N(0, tau^2) from the same population.
```
The severe/mild separation enters as `load_shift` (a composition covariate), and the
hospital→surveillance transfer is `beta * (load_shift_o - 0)` with a wide prior — exactly the
explicit offset you asked for, but on the load axis.

**Still-open choices for your quick confirm before I write Stan:**
1. **`load_shift` construction** — build it from reported duration + age + inpatient-fraction as a
   crude standardized index (my proposal), or keep it a single binary hospital-paired(0)/
   surveillance(1) indicator with a wide `beta`? (The index is better if the paired studies report
   enough; several report age/duration, few report inpatient fraction — I'll use what's there.)
2. **`beta` prior** — how strong? Proposal: weakly-informative, sign-constrained so a lower-load
   surveillance population cannot raise Se, centered on a modest reduction with wide tails, and
   swept in the sensitivity analysis.
3. **Seed data** — lock Vallenas d=15 and pick volumes for the two-volume studies now (Gasem 1995
   → model both 3 & 10 mL arms; Wain 2008 → both 5 & 15 mL arms, since the paper reports per-volume)?

Answer these three and I build Phase 1 (Stan model + simulation-recovery: parameter recovery, the
single-outbreak ridge demo, the ">1" scenario, and the conditional-dependence stress test), then
**stop at the recovery gate** before any real data.
