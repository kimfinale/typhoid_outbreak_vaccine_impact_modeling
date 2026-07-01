# Community-surveillance PPV from the outbreak Summary sheet

Can we estimate the clinical-definition PPV under community surveillance from the suspected-case
attack-rate + culture columns, and is it transferable? Data: Summary sheet cols AF (suspected),
AG (culture-tested), AI (culture-confirmed), AN (attack rate %). Build + fit:
`Rscript latent_class_ppv/analyze_community_ppv.R` -> `data/community_surveillance_ppv.csv`,
`tables/community_ppv_estimates.csv`, `figures/fig_community_ppv.png`.

## IMPORTANT: what is and isn't estimated
Each point gives ONE number, positivity = k/n = **pi * Se_BC** (a product; one equation, two
unknowns). So from these single-test points **Se_BC is NOT estimated** - it is fixed by an
informative prior (community/mild ~0.55). Given that assumed Se_BC, each positivity implies a
pi = positivity/Se_BC; partial pooling then summarises those as a mean + spread. **Se_BC is
estimated only in the separate historic paired blood+bone-marrow layer** (`bc_bmc_ppv.stan`),
where two tests per patient make it identifiable (Hui-Walter).

## Dataset & tightening
9 points have positivity + attack rate. Excluded:
- **4 confirmed-only / selective** (positivity >= 0.8 => pi > 1, or tested >= suspected):
  Ali 2017, Nimonkar 2022, Srinivasan 2022, Wang 2022.
- **2 not genuinely community-syndromic** (case defs are hospital/loose, not comparable):
  **Lewis 2005** (hospital admission; "admitting physician suspects typhoid"; no criteria) and
  **Roy 2016** (no case definition; presented with fever to a medical college / nursing homes).

The case definitions are heterogeneous - stringency ranges from "fever" (Roy) / "physician
suspects" (Lewis) to strict multi-symptom + malaria-negative (Kabwama), and catchment from a
village (Aye) to hospital admissions (Lewis/Roy). Tellingly the data run OPPOSITE to stringency
(strictest Kabwama -> lowest PPV; loosest Roy -> highest), so positivity is dominated by true
prevalence and testing-selection, not a shared definitional PPV. Hence the tightening.

**Tight community-syndromic set (explicit criteria):** Aye 2004, Neil 2012, Kabwama 2017.
positivity: Aye 0.136, Kabwama 0.140, Neil 0.397. The two clean community points (Aye, Kabwama)
agree at ~0.14 -> pi ~ 0.25; Neil (0.40) is the outlier and is IP-skewed / ascertainment-changed.

## Estimate (tight set; community Se_BC ~ Beta(12,10), mean 0.545)
| Quantity | Estimate (90% CrI) |
|---|---|
| typical community PPV `pi_mean` | **0.43 (0.20-0.74)** |
| between-setting spread `sigma` (logit) | **1.05 (0.48-2.05)** (data-driven) |
| predictive PPV, new outbreak `pi_new` | **0.43 (0.08-0.90)** = the transferable prior |

## Level vs spread across the Se_BC anchor (answers "does the spread move with Se_BC?")
The LEVEL tracks Se_BC; the SPREAD is nearly invariant (drifts only near the c_max=0.40 floor):

```
Se_BC~    pi_mean(LEVEL)   sigma(SPREAD)
0.45          0.47             1.15
0.50          0.45             1.09
0.55          0.42             1.01
0.60          0.39             1.03
0.66          0.35             1.00
0.75          0.33             0.92
```

Mechanism: pi = positivity/Se_BC rescales all pi by 1/Se_BC. The ORDER of outbreaks is exactly
invariant to Se_BC; on the logit scale the spread is ~constant, stretching slightly as Se_BC
approaches the largest positivity (top pi -> 1).

## Findings (honest)
1. **PPV is NOT a transferable constant** - even after tightening, the predictive PPV for a new
   outbreak is 0.43 with a 0.08-0.90 interval. Consistent-looking case definitions do not make PPV
   consistent; local true-typhoid : other-febrile ratio and testing-selection dominate.
2. **The attack rate is not a usable covariate** (no signal; denominator inconsistent across
   studies - e.g. Neil AN 8.09% vs 0.086% implied by suspected/population).
3. **Level needs an external Se_BC; spread is data-driven and Se_BC-robust.**

## Bearing on the severe->mild offset (beta)
These points ARE the community/mild layer measured directly, so Phase 2 can estimate community
pi from real positivity instead of a non-identifiable hospital->community `beta`. What remains
irreducible is the absolute Se_BC level (the ridge), pinned by the historic paired layer + a prior.
