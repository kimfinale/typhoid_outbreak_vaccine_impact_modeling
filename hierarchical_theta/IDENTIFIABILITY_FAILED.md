# Identifiability gate FAILED — θ level is not recoverable from curve shape

The synthetic identifiability gate (Phase 1) **failed**, robustly, across two model
specifications. Per the gated workflow we **do not fit real-data θ point estimates**;
the article classification (`source_decomposition/`, Stream A) remains primary, and this
Bayesian model is at most a loose **ordinal** consistency check.

## What was tested
Simulate 18 outbreaks across true θ ∈ {0.1, 0.3, 0.5, 0.7, 0.9} **from the generative
model itself** (best case for recovery), fit the neutral hierarchical model, and check
recovery of θ. PASS required 90% CrI coverage ≥ 0.80, |median bias| ≤ 0.12, rank cor ≥ 0.70,
and convergence.

## Result — both specifications fail the same way

| Specification | coverage | median bias | rank cor | divergences | max R-hat |
|---|---|---|---|---|---|
| (1) neutral wide prior (theta ~ logit-N(0, 1.5)) | 0.50 | -0.176 | 0.84 | 0 | 1.006 |
| (2) source-leaning prior (theta ~ 0.70) + broad source pulse | 0.50 | -0.202 | 0.85 | 6 | 1.004 |

The MCMC **converged** in both (R-hat ~ 1.00) — this is a statistical identifiability limit,
not a sampling or technical failure. The signature is consistent:

- **Rank is recoverable** (rho ~ 0.84-0.85): the model can tell a more-propagated outbreak
  from a more-source-driven one in *order*.
- **Level is not** (coverage 0.50, bias ~ -0.18 to -0.20): theta is systematically
  **under-estimated** — the model attributes too much of the curve to transmission, even
  for truly source-dominated outbreaks, and even under a source-leaning prior.

## Why (mechanism)
The renewal term Lambda[o,t] = sum_k w[k]*I[o,t-k] conditions on **observed lagged
incidence**. A common-source hump (a source pulse convolved with the incubation period) is
itself autocorrelated, so recent cases mechanically predict present cases — and the
transmission term R_tr*Lambda can mimic a source-driven curve. The source-vs-transmission
split is therefore weakly identified: the likelihood resolves the degeneracy toward
transmission (theta low). A single smooth source pulse and a domain-justified source-leaning
prior reduce but do not remove this; only the *ordinal* contrast survives.

## Recommendation (the honest fallback)
1. **Article classification is primary** — the field-epidemiology mode (water source tested,
   spatial clustering, secondary attack, remediation) in `source_decomposition/` Stream A
   identifies the mode that curve shape cannot.
2. Use this model, if at all, as an **ordinal** cross-check (the rank of theta), not for
   theta point estimates or the impact bracket. The static<->renewal bracket should use the
   article theta bands (as in `source_decomposition/`), with theta treated as a range.
3. A fundamentally different design — a fully latent-infection state-space model that does
   not condition Lambda on observed incidence — *might* improve level identification but is a
   much larger build with its own priors-do-the-work risks; not pursued here.

## Reproduce
`Rscript hierarchical_theta/run_phase1_gate.R` -> `tables/synthetic_recovery.csv`,
`figures/fig_synthetic_recovery.png`, `outputs/phase1_gate.rds`. Phase 2/3
(`run_phase2_3.R`) refuses to run without `outputs/GATE_PASSED.txt`.
