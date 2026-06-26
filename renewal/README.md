# Renewal-equation (R_t) model for typhoid ORI

A dynamic complement to the manuscript's static proportional-reduction model.
The renewal model reconstructs the instantaneous reproduction number R_t from each
observed outbreak, reduces transmission under a vaccination counterfactual, and
propagates the renewal equation forward so that **indirect (herd) protection emerges
endogenously** — no fixed indirect-effectiveness (η) term is supplied. The only thing
that differs from the static model is how cases-averted is computed; downstream
DALY/cost machinery is shared (R/2-functions.R).

## How to run

```sh
# from the repo root
Rscript renewal/run_analysis.R
```

Requires R ≥ 4.2 with `yaml` (and `EpiEstim` for the Phase-2 R_t cross-check).
All parameters live in `renewal/config.yml` (sourced once; seed printed at startup).

## Status

**Complete.** Core engine, Sobol uncertainty pipeline, DALY/cost integration,
validation tests, all figures/tables, and the rendered Quarto report. All 14
validation checks pass; self-consistency holds to ~1e-13. `Rscript renewal/run_analysis.R`
reproduces every output (and renders `report.html` if Quarto is installed) in ~1 min.

## Outputs

`tables/`: `tab_resolution.csv`, `tab_selfconsistency.csv`, `tab_eta_eff.csv`,
`tab_impact_summary.csv`, `summary_delay_coverage.csv`, `pooled_estimates.csv`,
`sens_gi.csv`, `sens_psiT.csv`.
`figures/` (PNG+PDF ≥300 dpi): `fig_Rt_panels`, `fig_forest_pctreduction`,
`fig_amplification_vs_delay`, `fig_eta_eff`, `fig_policy_grid`, `fig_gi_sensitivity`,
`fig_psiT_sensitivity`, `fig_curves_examples`.
`outputs/`: `sessionInfo.txt`, `results.rds`. `report.html`: assembled report.

## Layout

```
renewal/
├── config.yml            # all parameters (manuscript Table 1) + paths + seed
├── run_analysis.R        # top-level reproducible driver (Phase 1 + Phase 2 + report)
├── report.qmd            # Quarto report assembling all outputs
├── R/
│   ├── gi.R              # discretize_gi(), GI mass diagnostics
│   ├── renewal_core.R    # reconstruct_Rt, static/renewal counterfactual, impact, eta_eff
│   ├── data_prep.R       # load, exclusion, resolution rule, daily→weekly aggregation
│   ├── scenario.R        # Sobol draws + static-vs-renewal scenario runner
│   ├── cost_daly.R       # cost/DALY loaders + REUSE of R/2-functions.R add_cea_results
│   ├── summarise.R       # median+UI, pooled/pop-weighted, eta_eff table
│   ├── epiestim_rt.R     # EpiEstim R_t cross-check (plots only)
│   └── figures.R         # all ggplot figure builders
├── tests/test_renewal.R  # 5 validation tests (+ extras), base-R harness
├── tables/  figures/  outputs/   # generated (gitignored)
```

## Required external inputs (in `data/`)

Cost/DALY needs `country_costs.xlsx`, an IMF CPI export, `wpp_life_expectancy_*.csv`,
and `GDP_WorldBank.xls`. If any are missing the cost arm is skipped automatically and
the transmission outputs are unaffected. **Note:** the provided IMF file
(`imf-dm-export-20241028.xls`) has no `World` row, so the Sub-Saharan Africa CPI series
is used as the inflation deflator (most outbreaks are in SSA; non-SSA outbreaks use it
as an approximation). Swap in a World-series export to change this.

## Key assumptions

- **Incidence series:** `total_cases` relabeled `suspected_cases` (matches R/1-data.R,
  where raw `suspected_cases` is frequently NA). Static and renewal use the identical series.
- **Generation interval:** gamma, base-case mean 14 d (CV ≈ 0.6); the first weekly bin
  (0–7 d) is zeroed then renormalized (incubation ≥ 1 wk, so no within-week transmission).
  Sensitivity over μ_g ∈ {7,14,21,28} d.
- **Resolution rule:** renewal is applied only where μ_g/Δ ≥ 2 — the 13 daily/weekly
  outbreaks. Daily series are aggregated to weekly; fortnightly (μ_g/Δ = 1) and monthly
  (μ_g/Δ ≈ 0.47, ~95% of GI mass in one bin) are excluded and remain in the static analysis.
- **Vaccine on transmission:** R_t^v = [1 − c(t)·π·ψ_T]·R̂_t from t* onward; direct effect
  only (π·ψ_T). Base ψ = 0.83 (matches the static Sobol pipeline, not parameters.csv's 0.87).
  ψ_T (transmission-relevant VE) defaults to ψ with a sensitivity floor of 0.6·ψ.

## Validation tests

1. Self-consistency: max|I^v(ψ=0) − I^obs| < 1e-8 for every outbreak.
2. GI normalization: Σ w_s = 1 for every μ_g.
3. Static identity: renewal with transmission feedback OFF and η = 0 equals the static
   model with η = 0 to machine precision (documents static = renewal − herd feedback + fixed η).
4. Resolution assert: μ_g/Δ ≥ 2 for all included outbreaks.
5. Monotonicity: earlier intervention ⇒ weakly greater cases averted.
