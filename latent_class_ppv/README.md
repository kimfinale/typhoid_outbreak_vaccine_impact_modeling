# Joint Se_BC + clinical-definition PPV for typhoid (latent-class, volume-anchored)

A Bayesian latent-class model that estimates, from real data: blood-culture sensitivity `Se_BC` as
a function of **blood volume**, bone-marrow-culture sensitivity `Se_BM`, and the positive predictive
value (PPV) of the clinical/suspected case definition under **hospital** (`phi_s`) and
**community-surveillance** (`pi_o`) ascertainment — with an explicit **transportable (test) vs
non-transportable (population) parameter split**, built to plug into the outbreak
observation/transmission model as an observation operator.

Status: **Phase 0 (novelty) done; Phase 1 (simulation-recovery) PASSED; Phase 2 (real-data fit)
DONE.** The final results and the full write-up are in **`FINAL_REPORT.md`**. This README is the
front door / map of the module.

---

## 1. Why this exists (the problem)

For the *tested* subset of suspected patients, with blood-culture specificity ~ 1:

```
positivity  =  k / n  =  pi * Se_BC          (Sp_BC ~ 1)

  pi     = P(true typhoid | suspected)  = PPV of the clinical case definition   [WANT, local]
  Se_BC  = blood-culture sensitivity     = test property                         [WANT, transportable]
  k, n   = culture-positive / cultured   (observed from line lists)
```

Two facts drive the design:

1. **Non-identifiable from one test.** The likelihood sees only the *product* `pi * Se_BC`. A
   single outbreak identifies the product, not the parts. Naive `positivity / 0.6` can exceed 1
   (a falsification signal, not a valid PPV).
2. **Transportability asymmetry.** `Se_BC` is a *test* property (transportable; depends mainly on
   blood volume, plus age/illness-week/antibiotics) and can be shared. `pi` is a *population*
   property (depends on local prevalence and case-definition specificity; ~2-4% in low-endemicity
   surveillance, much higher in outbreaks) and must stay local.

The identifying leverage: **historic studies that cultured the same patients by both blood AND bone
marrow**. Two conditionally-independent, highly-specific tests are jointly identified (Hui-Walter),
so these pin `Se_BC` and `Se_BM`; that anchored `Se_BC` then lets each outbreak's positivity
identify its own `pi`.

## 2. What is novel (vs the three anchor papers)

- **Mogasale 2016** gives a *pooled* Se_BC only, under a union-complete assumption (both-negative =
  true-negative). We *estimate* the both-negative-truly-positive cell (bone marrow imperfect).
- **Arora 2019** (closest) is a test-accuracy network meta-analysis; it *fixes* Se_BC and Se_BM as
  priors, ignores volume, and could not use prevalence ("all studies enrolled suspected"). We
  estimate that prevalence: `phi_s` (hospital PPV) and `pi_o` (community PPV).
- **Antillon 2018** established Se_BC vs blood volume but did not fold it into a PPV/observation
  estimator. We reuse its volume relationship as a structural prior.

Novel contribution: *a volume-anchored `Se_BC` posterior that transfers to outbreaks to yield
per-outbreak surveillance PPV, as a plug-in observation operator, organized by a transportable-Se /
local-PPV split* — and, as the applied headline, the **hospital-vs-community PPV gap**. Honestly:
conditional on the anchored Se_BC, `pi = positivity/Se_BC` is arithmetically shallow — the value is
uncertainty propagation, the ">1" resolution, and the transportability framing.

## 3. Data

```
Historic paired studies s = 1..S      (blood + bone marrow, same patients):
  N_s        enrolled & tested by both
  a,b,c,d    = (BC+BM+, BC+BM-, BC-BM+, BC-BM-)      N_s = a+b+c+d
  volume_s   blood volume (mL)

Outbreak single-test studies o = 1..O (blood culture only):
  n_o        cultured;  k_o  of those BC+;   volume_o  blood volume (mL, assumed for outbreaks)
```

Real data: the 10 Mogasale paired 2x2 + volumes (Vallenas both-negative cell corrected to 15;
Hirsowitz volume imputed) in `data/mogasale2016_paired_bc_bmc_seed.csv`; the 25 Antillon volume/Se
points (`data/antillon2018_volume_sensitivity.csv`, volume-slope prior); and the community outbreak
positivity (`data/community_surveillance_ppv.csv`, tight syndromic set: Aye, Neil, Kabwama).

## 4. The model (final — `stan/model_final.stan`)

```
Se sub-model:   logit(Se_BC) = alpha0 + alpha1*log(volume) + tau*u,   u ~ Normal(0,1)
                  alpha0 = level, alpha1 = volume slope (Antillon-informed), tau = heterogeneity

Historic LCM (conditional-independence, marginalised over status):
  p_a = phi*Se_BC*Se_BM        + (1-phi)*(1-Sp_BC)*(1-Sp_BM)
  p_b = phi*Se_BC*(1-Se_BM)    + (1-phi)*(1-Sp_BC)*Sp_BM
  p_c = phi*(1-Se_BC)*Se_BM    + (1-phi)*Sp_BC*(1-Sp_BM)
  p_d = phi*(1-Se_BC)*(1-Se_BM)+ (1-phi)*Sp_BC*Sp_BM
  (a,b,c,d) ~ Multinomial(N, (p_a,p_b,p_c,p_d))

Outbreak (single-test):
  k_o ~ Binomial(n_o, pi_o*Se_BC(vol_o) + (1-pi_o)*(1-Sp_BC))
```

Priors: `alpha1 ~ N(0.36,0.15)` (Antillon volume slope), `Se_BM ~ Beta(12,2)` (~0.86, favours high
to break LCM label-switch), `Sp ~ Beta(200,1)` (~0.995), `phi_s, pi_o ~ Beta(1,1)` **independent
(NOT pooled)**. There is **no severity/`mild`/`beta` term** — see Section 7.

## 5. Assumptions (read before using)

1. **Blood-culture specificity ~ 1** (a positive culture is a true case; confirmed cases ~ all true
   typhoid). Relaxable as a sensitivity.
2. **Conditional independence** of blood and bone marrow given true status (baseline). Unmodelled
   positive dependence (prior antibiotics failing both) over-states Se_BC / under-states PPV; the
   optional covariance term (Phase-1 model) stress-tests it.
3. **Bone marrow is imperfect** (Se_BM estimated ~0.90) — not a gold standard. This is what lets us
   estimate the both-negative-truly-positive cell and identify Se_BC.
4. **Se_BC transfers via blood VOLUME** (+ residual heterogeneity `tau`). Severity is NOT modelled
   (Section 7).
5. **Community-outbreak blood volume is assumed** (default 5 mL; swept 3-10 mL) — the one remaining
   assumption, replacing the dropped severity offset. Concrete and researchable per outbreak.
6. **Prevalences are local**: `phi_s` (hospital PPV) and `pi_o` (community PPV) are independent per
   study/outbreak — never pooled, never transferred.

## 6. What transfers and what does not

- **Transferred (shared):** the Se sub-model `alpha0, alpha1, tau` and `Se_BM`.
- **Local (never transferred):** every PPV — `phi_s`, `pi_o`. Hospital PPVs are a validation range
  only (community PPVs should sit below them), never a prior for outbreaks.

## 7. Identifiability, and why the severity offset was DROPPED

- **The ridge:** one outbreak alone identifies only `pi*Se` (Phase-1: pi sd 0.26, Se sd 0.26, but
  product sd 0.02). The historic paired anchor is what separates them.
- **No severity offset.** An earlier design had a severe->mild offset `beta` on logit(Se_BC). A
  literature search (Zotero + web; `SE_BC_MILD_LITSEARCH.md`) found **no paired blood/bone-marrow
  data in mild/community cases can exist** (bone marrow is invasive; community typhoid surveillance
  uses serology, not culture), so `beta` is unidentifiable *in principle*; and severity is **not** a
  supported independent determinant of Se_BC (the drivers are volume, age, illness-week, antibiotics
  — the sign of a severe->mild gap is even ambiguous). So Se_BC transfers across settings via the
  **volume model**, and the outbreak Se level rides on the assumed outbreak volume (swept), not a
  severity parameter.

## 8. Results (real-data fit — see `FINAL_REPORT.md` for detail)

Converged: 0 divergences, max R-hat 1.002, min ESS 1826.

| Quantity | Estimate (median [90% CrI]) |
|---|---|
| Se_BC at 2 / 5 / 10 mL | **0.52 / 0.62 / 0.68** (slope 0.42 per log-mL) |
| Se_BM (bone marrow) | **0.90 [0.87, 0.93]** |
| Sp_BC | 0.99 |
| Hospital PPV `phi_s` (10 studies) | **0.66-0.99** (median 0.83) |
| Community PPV `pi_o` | Aye **0.26**, Kabwama **0.23**, Neil **0.65** |

**Headline:** the clinical definition's PPV is population-specific — **hospital ~0.8 vs community
~0.25** (a ~3-fold gap; `figures/fig_final_ppv.png`); and Se_BC is volume-dependent
(`figures/fig_final_se_volume.png`). These reproduce Mogasale/Antillon (Se_BC ~0.6, volume gradient)
and literature Se_BM (~0.9) from an independent latent-class fit; sanity check (community pi <
hospital phi) passes.

Phase-1 simulation-recovery (the gate, `PHASE1_RECOVERY.md`): PASSED — recovery of Se_BC/Se_BM/PPV,
the ridge demo, the ">1" plug-in resolution, and a conditional-dependence stress test.

## 9. Honest caveats

- **Community PPV level depends on the assumed outbreak blood volume** (swept 3-10 mL -> pi_o robust:
  Aye 0.28->0.24, Neil 0.70->0.60, Kabwama 0.25->0.21). This is the one remaining assumption.
- **Only 3 clean community points** — community PPV is thin and heterogeneous (a wide prior, not a
  precise estimate; ~3-fold spread across settings, `COMMUNITY_PPV.md`).
- **Conditional independence assumed** (dependence over-states Se, under-states PPV).
- **Historic volumes imperfect** (one imputed; two studies cultured two volumes at one modelled
  value); the volume slope leans on the Antillon prior.
- **Where it pays off:** for ratio outcomes (R_t, theta, %-reduction) the constant scaling cancels;
  Se_BC and pi_o matter for **absolute burden / DALYs / cost**.

## 10. Files & how to run

```
latent_class_ppv/
  README.md                  <- this file (map)
  FINAL_REPORT.md            <- the comprehensive results write-up  [START HERE]
  PHASE0_novelty_memo.md     <- novelty verdict vs the 3 anchor papers
  PHASE1_RECOVERY.md         <- simulation-recovery gate report
  COMMUNITY_PPV.md           <- community-PPV analysis (not a constant; spread vs Se_BC)
  SE_BC_MILD_LITSEARCH.md    <- literature search -> the beta-drop verdict
  stan/model_final.stan      <- FINAL model (volume-anchored, beta-free)
  stan/bc_bmc_ppv.stan       <- Phase-1 gate model (has the optional/legacy beta term)
  stan/community_ppv.stan    <- standalone community-PPV model
  run_final.R                <- consolidate real data + fit + inferences + figures  [MAIN]
  run_phase1_recovery.R      <- simulation-recovery gate
  analyze_community_ppv.R    <- community-PPV analysis
  config.yml                 <- Phase-1 sim/gate config
  data/                      <- Mogasale paired, Antillon volume/Se, community positivity
  figures/ tables/ outputs/  <- generated (gitignored)
```

Needs R 4.5.3 + cmdstanr/cmdstan 2.39. Final fit (~2-3 min): `Rscript latent_class_ppv/run_final.R`.
Recovery gate (~5-7 min): `Rscript latent_class_ppv/run_phase1_recovery.R`.

## 11. Roadmap

- **Phase 0** — novelty memo. DONE.
- **Phase 1** — model + simulation-recovery gate. DONE, PASS.
- **Phase 2** — real-data fit (beta-free, volume-anchored). DONE — see `FINAL_REPORT.md`.
- **Next (optional):** fill in real per-outbreak blood volumes (removes the last assumption); add
  more community-surveillance points; embed Se_BC/PPV in the DALY/cost arm; potential standalone
  paper on the hospital-vs-community PPV gap + volume-dependent Se_BC.
