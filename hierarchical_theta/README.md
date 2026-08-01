# Hierarchical Bayesian renewal-with-source model (source fraction θ)

Attempts to estimate each outbreak's **source fraction θ_o** (1 = source-dominated;
0 = propagated-dominated) by fitting a hierarchical
renewal-with-source model to all outbreaks jointly (partial pooling), in
**R + Stan (cmdstanr)**. Tests honestly whether epidemic-curve *shape* recovers the
transmission mode that field investigators establish — without baking the article
classification into the prior.

## Model (`stan/renewal_source.stan`)

```
I[o,t] ~ NegBinomial2( mu[o,t], phi )
mu[o,t] = s[o,t]  +  R_tr[o] · Λ[o,t]
Λ[o,t]  = Σ_{k≥1} w[k] · I[o,t−k]            (conditions on OBSERVED lagged incidence)
```
- **s** = exogenous common-source term: a single smooth gamma pulse × amplitude. The
  single-pulse constraint is the identifiability lever — it stops the source absorbing
  arbitrary curve structure; for a continuous source the pulse *width* grows.
- **R_tr[o] ∈ (0,1)** = subcritical per-outbreak transmission (propagated channel).
- **θ_o = 1 − R_tr[o]** (verified on simulation). `logit(θ_o)` is partially pooled,
  non-centred: neutral model → pooled to a grand mean; informed model → pooled to a
  **mode-specific** mean (A/B/C covariate).
- **GI mean** is a parameter (lognormal prior ~2 weeks) so its uncertainty propagates
  into θ; `gi_cv` fixed at 0.6.

## Gated workflow (honesty)

1. **Phase 1 — synthetic identifiability gate** (`run_phase1_gate.R`). Simulate ~18
   outbreaks across true θ ∈ {0.1,…,0.9}, fit the **neutral** model, check recovery
   (90% CrI coverage ≥ 0.8, |median bias| ≤ 0.12, rank cor ≥ 0.7, converged).
   - **PASS** → `GATE_PASSED.txt`, proceed to Phase 2.
   - **FAIL** → `IDENTIFIABILITY_FAILED.md`, **STOP** — do not report real-data θ point
     estimates; the article classification is primary, the model a loose check.
2. **Phase 2 — real-data fits** (only if gate passes): M1 **neutral** + M2 **article-
   informed**, 4 chains, full diagnostics (R̂ < 1.01, ESS, divergences), PPCs.
3. **Phase 3 — comparison, if the gate passes:** M1 posterior θ vs the article θ band
   (convergent validity), rank concordance, and τ_θ (pooling). Any impact analysis must
   allocate source and propagated incidence inside the additive renewal recursion.
   The earlier post hoc θ-weighted static-plus-renewal bracket is historical only.

## Why this differs from `uncertainty/`'s Bayesian model
That module's self-contained MCMC reconstructs R_t (no source term). This is a
*generative* renewal+source model fit jointly across outbreaks to estimate θ — it needs
Stan (the renewal convolution can't be expressed in brms).

## Result: gate FAILED (see `IDENTIFIABILITY_FAILED.md`)

The synthetic gate failed across two specifications (neutral and source-leaning priors):
θ's **rank** is recoverable (Spearman ≈ 0.85) but its **level** is not (90% CrI coverage
0.50, median bias ≈ −0.18 to −0.20), even on data simulated from the model itself, with clean
convergence (R̂≈1.00). Mechanism: conditioning Λ on observed incidence lets the transmission
term mimic a common-source hump. **Phase 2/3 was therefore not run** — real-data θ point
estimates are not reported. The article classification (`source_decomposition/`) is primary;
this model is an ordinal cross-check only.

## Reproduce
```sh
Rscript hierarchical_theta/run_phase1_gate.R     # gate first
Rscript hierarchical_theta/run_phase2_3.R        # only if GATE_PASSED.txt exists
```
Requires cmdstanr + a compiled cmdstan (RTools45 on Windows). Seeds (R + Stan) are set
and printed; Stan files, fits, and `sessionInfo()` are saved. Outputs gitignored.
