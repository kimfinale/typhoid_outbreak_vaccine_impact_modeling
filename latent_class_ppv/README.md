# Joint Se_BC + clinical-definition PPV for typhoid (latent-class, load-anchored)

A Bayesian latent-class model that estimates, **per outbreak**, the positive predictive value
(PPV) of the clinical/suspected case definition — i.e. `pi_o = P(true typhoid | suspected)` —
jointly with a **blood-volume- and severity-anchored blood-culture sensitivity `Se_BC`**, with an
explicit **transportable (test) vs non-transportable (population) parameter split**, built to plug
into the outbreak observation/transmission model as an observation operator.

Status: **Phase 0 (novelty) done; Phase 1 (simulation-recovery) PASSED; Phase 2 (real data) awaiting
sign-off.** This README is the front door; see `PHASE0_novelty_memo.md` and `PHASE1_RECOVERY.md`
for the detailed memos.

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
2. **Transportability asymmetry.** `Se_BC` is a *test* property (roughly transportable; depends on
   blood volume, bacterial load, antibiotics) and can be pooled. `pi` is a *population* property
   (depends on local prevalence and case-definition specificity; ~2-4% in low-endemicity
   Pemba/Zanzibar, much higher in outbreaks) and must stay local.

The identifying leverage comes from **historic studies that cultured the same patients by both
blood AND bone marrow**. Two conditionally-independent, highly-specific tests are jointly
identified (Hui-Walter), so these studies pin `Se_BC` (and `Se_BM`); that anchored `Se_BC` then lets
each outbreak's positivity identify its own `pi`.

## 2. What is novel (vs the three anchor papers)

- **Mogasale 2016** gives a *pooled* Se_BC only, under a union-complete assumption (both-negative =
  true-negative). We *estimate* the both-negative-truly-positive cell (bone marrow imperfect).
- **Arora 2019** (closest) is a test-accuracy network meta-analysis; it *fixes* Se_BC (~50%) and
  Se_BM as priors, ignores volume, and — its own words — could not use prevalence ("all studies
  enrolled suspected"). We estimate that prevalence: `phi_s` (hospital PPV) and `pi_o` (surveillance
  PPV).
- **Antillon 2018** established Se_BC vs blood volume but did not fold it into a PPV/observation
  estimator. We reuse its volume relationship as a structural prior.

Novel contribution: *a volume-and-load-anchored `Se_BC` **posterior** that transfers to outbreaks
to yield **per-outbreak surveillance PPV with propagated uncertainty**, as a plug-in observation
operator, organized by a transportable-Se / local-PPV split.* Honestly: conditional on an anchored
Se_BC, `pi = positivity/Se_BC` is arithmetically shallow — the value is uncertainty propagation, the
">1" resolution, and the transportability logic, **not** the latent-class machinery itself.

## 3. Data schema

```
Historic paired studies s = 1..S      (blood + bone marrow, same patients):
  N_s        enrolled & tested by both
  a,b,c,d    = (BC+BM+, BC+BM-, BC-BM+, BC-BM-)      N_s = a+b+c+d
  volume_s   blood volume (mL);  mild_s = 0 (severe/hospital baseline)

Outbreak single-test studies o = 1..O (blood culture only):
  n_o        cultured;  k_o  of those BC+
  volume_o   blood volume (mL);  mild_o = 1 (surveillance = mild, by default)
```

Seed historic data (`data/`): the 10 Mogasale paired 2x2 + volumes (Vallenas d corrected to 15;
Gasem 1995 & Wain 2008 modelled as both volume arms) and the 25 Antillon volume/Se points.

## 4. The model

```
Se sub-model (bacterial-LOAD axis):
  logit(Se_BC) = alpha0 + alpha1*log(volume) + beta*mild + tau*u,   u ~ Normal(0,1)
     alpha0,alpha1  volume model (Antillon log-linear form)
     beta <= 0      severe->mild offset (lower load -> lower Se)
     tau            between-study heterogeneity

Historic likelihood (conditional-independence latent class, marginalised over status):
  p_a = phi*Se_BC*Se_BM       + (1-phi)*(1-Sp_BC)*(1-Sp_BM)
  p_b = phi*Se_BC*(1-Se_BM)   + (1-phi)*(1-Sp_BC)*Sp_BM
  p_c = phi*(1-Se_BC)*Se_BM   + (1-phi)*Sp_BC*(1-Sp_BM)
  p_d = phi*(1-Se_BC)*(1-Se_BM)+ (1-phi)*Sp_BC*Sp_BM
  (a,b,c,d) ~ Multinomial(N, (p_a,p_b,p_c,p_d))

Outbreak likelihood (single-test):
  theta_o = pi_o*Se_BC_o + (1-pi_o)*(1-Sp_BC)
  k_o ~ Binomial(n_o, theta_o)

Optional conditional-dependence (stress test): a fixed covariance cd among the truly infected
shifts mass onto the agreement cells (models antibiotic-driven joint BC/BMC failure).
```

Priors: `alpha0 ~ N(0,1.5)`, `alpha1 ~ N(0,1)`, `beta ~ half-N(0,1)` (<=0), `tau ~ half-N(0,1)`,
`Se_BM ~ Beta(90,8)` (~0.92, imperfect), `Sp_BC,Sp_BM ~ Beta(200,1)` (~0.995), `phi_s, pi_o ~
Beta(1,1)` **independent per study/outbreak (NOT pooled)**.

## 5. Assumptions (read before using)

1. **Blood-culture specificity ~ 1** (Sp_BC, Sp_BM ~ 0.995): a positive culture is a true case, so
   confirmed cases are essentially all true typhoid. Relaxable (priors provided) as a sensitivity.
2. **Conditional independence** of blood and bone marrow given true status (baseline). Violated by
   prior antibiotics / low load failing both together; the optional covariance term stress-tests it
   (see caveat below).
3. **Bone marrow is imperfect** (Se_BM ~ 0.92, estimated) — not a gold standard. This is what lets
   us estimate the both-negative-truly-positive cell.
4. **Se_BC transfers via the load axis** (volume + severity + heterogeneity tau); everything not
   captured by volume/severity is absorbed by tau.
5. **The severe->mild offset `beta` is a FIXED assumption, not a fitted parameter** (see 6/7). All
   paired studies are severe and outbreaks are single-test, so no data identifies it.
6. **Surveillance = mild (lower load) by default.** Care setting ("hospital vs community") is *not*
   the axis — hospitals see mild outpatients too; the axis is bacterial load, proxied by severity +
   volume (Wain 1998: blood load median ~1 CFU/mL, falls with duration, higher in children).
7. **Prevalences are local**: `phi_s` (hospital PPV) and `pi_o` (surveillance PPV) are independent
   per study/outbreak — never pooled, never transferred.

## 6. What transfers and what does not

- **Transferred (shared across datasets):** the Se sub-model `alpha0, alpha1, beta, tau` and
  `Se_BM`. Only these connect historic studies to outbreaks.
- **Local (never transferred):** every prevalence/PPV — `phi_s` and `pi_o`. Historic hospital PPVs
  are a *validation range* only (surveillance PPVs should sit below them), never a prior for
  outbreaks.

## 7. Identifiability (proven in Phase 1)

- **The ridge:** one outbreak alone identifies only `pi*Se`. In simulation: `pi` sd 0.26, `Se` sd
  0.26, but product `theta` sd **0.02**. The historic anchor is what separates them.
- **The offset is not identifiable:** a *free* `beta` posterior lands at -1.29 (90% CrI excludes the
  truth -0.50) — it tracks the prior/boundary, not data. So `beta` is **fixed and swept**, and
  `pi_o` is reported *conditional on the offset* with a sensitivity band (pi_o moves ~0.09 per 0.5
  logit of assumed offset).

## 8. Results — Phase 1 simulation-recovery: PASS

| Check | Result |
|---|---|
| Recovery, offset LOCKED at truth | coverage(alpha0,alpha1,tau,Se_BM + pi_o) = **1.00**; median\|bias\| pi = 0.04; phi_s coverage 1.00 |
| Convergence | 0 divergences, max R-hat 1.008, min ESS 554 |
| Free `beta` (non-identified) | median -1.29 (CrI -1.58..-0.86), excludes truth -0.50 |
| Offset sweep (beta = 0 / -0.5 / -1.0) | mean pi_o = 0.33 / 0.42 / 0.53 |
| ">1" plug-in (positivity 0.62) | plug-in 1.03 > 1; Bayesian pi median **0.89**, P(pi>1)=0, Se ~0.69 |
| Conditional-dependence stress (cd=0.06) | Se +0.078, Se_BM 0.92->0.98, pi -0.076 |

Figures in `figures/`: `fig_recovery_{scalars,pi}.png`, `fig_offset_sweep.png`, `fig_ridge.png`,
`fig_gt1.png`.

## 9. Honest caveats

- **`pi_o` is only as good as the assumed offset** — the dominant uncertainty is the fixed
  severe->mild `beta`, reported as a sensitivity band, because no mild/community paired BC/BMC data
  exists to estimate it.
- **Load ≠ a clean "severity" label** — children carry *higher* load yet adults show higher
  observed Se (volume confounds); severity/volume are proxies for the true load axis.
- **Conditional dependence biases the direction we care about** — unmodelled positive BC/BMC
  correlation (prior antibiotics) over-states Se and under-states PPV; keep the covariance term as a
  sensitivity check, and headline "unreported antibiotic pretreatment in outbreaks."
- **Where it pays off:** for ratio outcomes (R_t, theta, %-reduction) the constant scaling cancels;
  the anchored Se_BC and pi_o matter for **absolute burden / DALYs / cost**.

## 10. Files & how to run

```
latent_class_ppv/
  README.md                     <- this file
  PHASE0_novelty_memo.md        <- novelty verdict vs the 3 anchor papers
  PHASE1_RECOVERY.md            <- recovery-gate report
  config.yml                    <- truth, sim design, priors, gate thresholds
  stan/bc_bmc_ppv.stan          <- the model
  run_phase1_recovery.R         <- simulate + fit + recovery + ridge + >1 + dependence + GATE
  data/mogasale2016_paired_bc_bmc_seed.csv    <- 10 paired 2x2 + volumes (seed)
  data/antillon2018_volume_sensitivity.csv    <- 25 volume/Se points (volume prior)
  figures/ tables/ outputs/     <- generated (gitignored)
```

Run the recovery gate (needs R 4.5.3 + cmdstanr/cmdstan 2.39; ~5-7 min):

```
Rscript latent_class_ppv/run_phase1_recovery.R
```

## 11. Roadmap

- **Phase 0** — novelty memo. DONE (materially additive; scope narrowed to the load-anchored
  observation model).
- **Phase 1** — model + simulation-recovery gate. DONE, **PASS**.
- **Phase 2** — real data (awaiting sign-off): fit the 10 Mogasale paired studies + outbreak line
  lists with the offset fixed & swept; report per-outbreak `pi_o` + offset band; validate
  surveillance `pi_o` < hospital `phi_s`. Optional: embed in the renewal/burden model for absolute
  DALYs/cost.
