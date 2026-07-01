# Community-surveillance PPV from the outbreak Summary sheet

Answers: can we estimate the clinical-definition PPV under community surveillance from the
suspected-case attack-rate + culture columns, and is it transferable?

Data: `data/Typhoid_Outbreak_Time_Series_2000_2022_Summary.csv` cols AF (suspected), AG
(culture-tested), AI (culture-confirmed), AN (attack rate %). Build + fit:
`Rscript latent_class_ppv/analyze_community_ppv.R` -> `data/community_surveillance_ppv.csv`,
`tables/community_ppv_estimates.csv`, `figures/fig_community_ppv.png`.

## Dataset
9 points have both culture positivity (AI/AG) and an attack rate. **4 excluded** as
confirmed-only / selectively-tested (positivity >= 0.8 => pi > 1 impossible; or tested >=
suspected = a confirmed-case series): Ali 2017, Nimonkar 2022, Srinivasan 2022, Wang 2022.
**5 clean syndromic community points:** Aye 2004, Lewis 2005, Neil 2012, Roy 2016, Kabwama 2017.

## Estimate (partial-pooled; community Se_BC ~ Beta(12,10), mean 0.545)
| Quantity | Estimate (90% CrI) |
|---|---|
| typical community PPV `pi_mean` | **0.56 (0.32-0.82)** |
| between-setting spread `sigma` (logit) | **1.11 (0.58-2.11)** (data-driven) |
| predictive PPV, new outbreak `pi_new` | **0.56 (0.12-0.94)** = the transferable prior |

Level sensitivity to the Se_BC anchor: `pi_mean` = 0.58 / 0.51 / 0.47 at Se_BC = 0.50 / 0.60 / 0.66.
Per-outbreak pi (Se_BC~0.55): Aye ~0.25, Kabwama ~0.26, Lewis ~0.69, Neil ~0.73, Roy ~0.85.

## Findings (honest)
1. **PPV is NOT a transferable constant.** Positivity spans 0.14-0.46 (**3.4-fold**); the
   predictive PPV for a new outbreak is 0.56 but with a **0.12-0.94** interval. Consistent case
   definitions do not make PPV consistent - the local true-typhoid : other-febrile ratio varies.
   It is transferable only as a **wide prior distribution**, not a usable point.
2. **The attack rate cannot rescue it.** No signal (Spearman 0.20, n=5) and its denominator is
   inconsistent across studies (Neil AN 8.09% vs 0.086% implied by suspected/population) -> not a
   clean covariate.
3. **The LEVEL still needs an external Se_BC** (positivity = pi*Se identifies only the product);
   the SPREAD is data-driven.

## Bearing on the severe->mild offset (beta)
These 5 points ARE the community/mild layer, measured directly. So Phase 2 can estimate community
pi_o **directly from these positivity data** (partially pooled), anchoring Se_BC from the historic
paired studies - instead of positing a non-identifiable hospital->community `beta`. The community
data supply what `beta` was guessing; only the absolute Se_BC level remains an external anchor.
