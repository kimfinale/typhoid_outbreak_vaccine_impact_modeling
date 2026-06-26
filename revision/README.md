# Typhoid ORI manuscript-revision outputs

Generates the numerical results, tables, and figures that fill the 10 margin
comments (c0–c9) in `papers/manuscript_tracked_revisions.docx`, plus a single
**results manifest** mapping every pending value to its computed result.

Reuses the renewal engine (`renewal/R/*`) and the manuscript CEA functions
(`R/2-functions.R`). Static and renewal share identical inputs and downstream
cost/DALY machinery — only cases-averted differs.

## How to run

```sh
# from the repo root (after renewal/ is built)
Rscript revision/run_revision_outputs.R
```

Produces `revision/outputs/revision_values.csv` (the manifest), the supporting
tables/figures, runs the guard tests, and renders `revision_manifest.html`.

## Comment → deliverable map

| Comment | Issue | Output |
|---|---|---|
| c0 | Exclusion breakdown (20/1/1/2=24 vs submitted 22) | `tab_reconciliations.csv` |
| c1 | Renewal-model results | `renewal/` outputs (companion module) |
| c2 | Indirect-VE coverage dependence (negative >67%) | `tab_ive_coverage.csv`, `fig_ive_coverage.png` |
| c3 | Confirmation adjustment (α); % reduction invariant | `tab_confirmation_adjusted.csv` |
| c4 | 1-week cost/DALY (not the duplicated 8-week value) | `tab_ce_decomposition.csv` (manifest row) |
| c5 | Outbreak span / per-year denominator (26,298/18=1,461) | `tab_reconciliations.csv` |
| c6 | Attack-rate units (per 100) | `tab_reconciliations.csv` |
| c7 | Impossible 96% 1-week extrapolation, recomputed | `tab_global_extrapolation.csv` |
| c8 | Cost-effectiveness honesty; deaths/YLL/YLD; spatial | `tab_ce_decomposition.csv`, `tab_spatial_targeting.csv`, figs |
| c9 | "95% CI" → percentile/IQR; pop-weighted; forest | `tab_intervals.csv`, `tab_weighted_impact.csv`, `fig_forest_pctreduction.png` |

## Guard tests (run automatically)

1. **Invariance** — % reduction identical before/after α (Part D).
2. **Impossibility** — averted ≤ post-intervention cases; static median ≤ P_max.
3. **Self-consistency** — renewal(ψ=0) reproduces observed (<1e-8).
4. **Static identity** — renewal(feedback off, η=0) == static(η=0).
5. **Resolution** — renewal only where μ_g/Δ ≥ 2.

## Key assumptions & sources

- **s_culture = 0.61** (95% CI 0.52–0.70) — Mogasale et al. (2016) *Ann Clin Microbiol Antimicrob* 15:32.
- **Opportunity-cost threshold = 0.5× GDP/capita** — Ochalek, Lomas & Claxton (2018) *BMJ Glob Health*; plus 1× and 3× GDP (WHO-CHOICE).
- **IVE algebra** — `IVE(f)=(OVE−f·TVE)/(1−f)`, Khanam et al. (2023) Vi-TT trial OVE=0.57, TVE=0.85 → zero at f≈67%.
- **Spatial φ** — φ_cases=0.80 captured in φ_pop∈{0.25–0.50}; **illustrative only**, anchored loosely to Harare (Poncin 2022). Not a fitted estimate.
- **Span** — data 2000–2018; manuscript "18 years"; denominator 26,298/18 = 1,461/yr (daily/weekly analysis set).
- **GI / ψ / inflation** — as in `renewal/` (mean 14 d, drop-first bin; ψ=0.83; SSA CPI deflator since the IMF file lacks a World row).

**Honesty note:** Part E is a parameterized illustrative sensitivity analysis, not a fitted estimate. Part D uses per-outbreak confirmed/suspected where available (15/19) with a pooled fallback.
