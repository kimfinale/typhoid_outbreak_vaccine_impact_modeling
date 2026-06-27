# Rigorous R_t analysis — propagated uncertainty + full Bayesian model

Extends the renewal (R_t) model to account, with literature-grounded ranges, for
the three uncertainties that actually move the estimand:

1. **Generation interval** — μ_g, CV, and a chronic-carrier long tail.
2. **Asymptomatic / unreported-but-transmitting infection** — enters through the
   transmission-relevant efficacy ψ_T = ψ·(1 − φ_a), φ_a = f_a·c_a/(f_a·c_a+(1−f_a)).
3. **PPV of clinical "suspected" cases** — time-varying PPV reshapes R_t; a
   *constant* PPV is provably invariant (tested).

Parameter ranges and citations: `params_uncertainty.csv` + `EVIDENCE.md`.

## How to run

```sh
# from the repo root (after renewal/ is built)
Rscript uncertainty/run_uncertainty.R
```

Outputs to `uncertainty/{tables,figures,outputs}` (all gitignored, regenerated):
- `tab_mc_impact.csv`, `tab_mc_pooled.csv` — per-outbreak & pooled % reduction with 95% CrI.
- `tab_prcc.csv` + `fig_prcc_tornado` — global sensitivity (partial rank correlation).
- `rt_ribbon.csv` + `fig_Rt_uncertainty` — R_t with the full uncertainty band over the data.
- `fig_forest_uncertainty` — per-outbreak % reduction, full CrI vs point estimate.
- `tab_bayes_impact.csv` + `fig_bayes_vs_mc` — full Bayesian model vs exact+MC.

## Two inference engines

- **Primary — exact reconstruction + Monte Carlo:** `R_t = I_t/Λ_t` (self-consistent),
  with N draws over {GI mean/CV, carrier tail, asymptomatic f_a/c_a, PPV level+trend,
  culture sensitivity, ψ}. Fast, reuses the validated engine.
- **Extension — full Bayesian renewal model** (`R/bayes.R`, self-contained, no Stan):
  latent log R_t random walk + Poisson renewal likelihood, with the GI mean sampled
  so GI uncertainty is integrated natively; componentwise adaptive Metropolis. A
  deliberately different reconstruction (joint smoothing + full posterior) to test
  robustness of conclusions to the inference approach.

## Headline findings (base scenario τ=8 wk, 80% coverage)

- Pooled % reduction: **point 36.0% → full-uncertainty MC 29.0% (95% CrI 18.8–36.3)**;
  full Bayesian 35.4% (matches the point estimate, isolating the reconstruction layer).
- **Dominant uncertainty is the asymptomatic-transmission fraction (f_a, c_a) and ψ —
  not the generation interval.** GI shape is secondary; carrier tail and PPV level are minor.
- **Level-invariance confirmed:** PPV level and culture sensitivity have ≈0 PRCC on the
  percent reduction (they rescale absolute true-case numbers, which are non-identifiable
  from one outbreak, but not the relative estimand).

## Honesty notes

- No direct empirical typhoid generation interval exists; it is derived from natural
  history (incubation + shedding), anchored above by the 4-week endemic-model GI.
- Constant reporting / PPV / asymptomatic fraction are invariant for % reduction (tested);
  rigor is concentrated on the components that genuinely change R_t or impact.
- Bayesian uses base ψ_T (no asymptomatic discount) so the MC-vs-Bayesian gap is exactly
  the asymptomatic structural uncertainty, not a reconstruction artifact.
