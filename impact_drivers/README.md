# Drivers of ORI impact

Which outbreak/environment characteristics — **knowable at the time of
vaccination** — predict the magnitude of outbreak-response immunization impact?
Combines the 13 real outbreaks (empirical, exploratory) with a renewal-model
**simulation experiment** (confounding-free, large n), then overlays the real
outbreaks on the simulated surface. Reuses the renewal engine, scenario runner,
and cost/DALY machinery.

## How to run

```sh
Rscript impact_drivers/run_impact_drivers.R   # ~50s; renders report.html if quarto present
```

## Design principles
- **Prospective predictors only:** computed from the first τ weeks + static
  context (early growth rate, R_t at τ, cumulative cases/attack-rate by τ, log
  population). Retrospective quantities (final size, peak, duration) are never
  predictors.
- **Three outcomes, ranked:** cases averted per 1,000 doses (primary), % case
  reduction, cases averted.
- **Simulation-anchored:** functional forms and the characteristic×delay
  interaction come from the simulation; the 13 real outbreaks validate them.

## Phases / outputs
1. **Features** → `tab_feature_outcomes.csv` (13 outbreaks × 7 delays).
2. **Empirical association** → `tab_associations.csv` (Spearman + bootstrap CI,
   ranked), `tab_mediation.csv`, `fig_ranking_*`, `fig_assoc_by_delay`.
3. **Simulation** (1,500 synthetic outbreaks) → `sim_outbreaks.csv`,
   `tab_sim_importance.csv`.
4. **Interaction** → `tab_delay_sensitivity.csv`, `fig_interaction`.
5. **Decision support** → `fig_impact_surface` (the impact map), `fig_validation`,
   `decision_rule.txt` + `fig_decision_tree`, LOO-CV.

## Findings
- **R_t at the time of vaccination (epidemic phase) is the dominant,
  mechanistically-grounded driver.** Early growth rate alone does not predict
  impact at a fixed late delay because fast outbreaks burn out sooner; what
  matters is the state (R_t / remaining burden) at the decision time.
- **Delay matters most for still-growing outbreaks:** % reduction falls ~2.7
  points per week of delay when R_t>1.2 at vaccination, vs ~0 when declining.
- Cumulative attack rate by τ (how much has already burned) is the strongest
  empirical correlate of % reduction (ρ ≈ −0.90); R_t-at-τ correlates with the
  efficiency metric (cases/1,000 doses, ρ ≈ +0.70).

## Limitations
n = 13 empirically (hypothesis-generating; simulation carries the functional
forms; LOO-CV skill is weak at this n); prospective-only predictors; reconstructed
early R_t is uncertain (winsorized for display); detection/selection bias; urban/
rural dropped; population missing for 4 outbreaks (cases/1,000 uses the rest).
