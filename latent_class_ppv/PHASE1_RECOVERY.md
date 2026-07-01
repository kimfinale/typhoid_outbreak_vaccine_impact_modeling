# Phase 1 — simulation-recovery gate: PASS (recovery GIVEN the offset)

**Status:** GATE PASSED. No real data has been touched. Awaiting sign-off before Phase 2.
Reproduce: `Rscript latent_class_ppv/run_phase1_recovery.R` (writes `outputs/GATE_PASSED.txt`,
`tables/`, `figures/`).

## The model (as approved)
Historic paired blood/bone-marrow studies -> a conditional-independence latent-class 2x2 that
estimates Se_BC **and** Se_BM (bone marrow imperfect) and each study's hospital PPV phi_s.
Outbreak single-test data -> each surveillance PPV pi_o, anchored by Se_BC. Se sub-model on the
**bacterial-load axis**: `logit(Se_BC) = alpha0 + alpha1*log(volume) + beta*mild + tau*u`, with
`mild` = severe(0)/mild(1) and `beta<=0` the severe->mild offset. Transportable
(alpha0, alpha1, beta, tau, Se_BM) vs local (phi_s, pi_o).

## Headline result: the offset is NOT identifiable — so it is fixed, not fitted
All historic paired studies are severe (mild=0) and outbreaks are single-test (Se/pi
confounded), so **no data identifies `beta`**. The gate therefore certifies recovery *given*
the offset, and separately quantifies the offset's non-identifiability:

| Check | Result |
|---|---|
| **Recovery, offset LOCKED at truth** | coverage(alpha0,alpha1,tau,Se_BM + pi_o) = **1.00**; max\|bias\| shared (logit) = 0.21; median\|bias\| pi = **0.04**; phi_s coverage = 1.00 |
| Convergence | 0 divergences, max R-hat 1.008, min ESS 554 |
| **Free `beta` (try to fit it)** | posterior median **-1.29** (90% CrI -1.58..-0.86) — *excludes* the truth -0.50; posterior ~ prior/geometry => **UNIDENTIFIED**. Resulting pi bias +0.22 (vs +0.03 when locked) |
| **Offset sensitivity sweep** | mean pi_o = 0.33 / 0.42 / 0.53 at beta = 0 / -0.5(truth) / -1.0 -> **pi_o moves ~0.09 per 0.5 logit of assumed offset** |

Interpretation: **the model recovers every identifiable quantity, and pi_o is pinned once the
offset is fixed. The severe->mild offset is the one structural assumption that data cannot
supply** — exactly caveat (c). Report pi_o *conditional on the offset*, with the sweep as the
primary sensitivity band.

## The three required demonstrations
- **Identifiability ridge** (`fig_ridge.png`): one outbreak, no historic anchor -> pi sd 0.26,
  Se sd 0.26, but **product theta sd only 0.02** — only pi*Se is identified; the anchor is what
  breaks the ridge.
- **">1" plug-in scenario** (`fig_gt1.png`): a high-yield outbreak with positivity 0.62 ->
  plug-in 0.62/0.60 = **1.03 > 1** (invalid). The Bayesian pi stays in [0,1]: **median 0.89,
  P(pi>1)=0**, with Se reallocated up to ~0.69. The pathology is dissolved.
- **Conditional-dependence stress** (simulate correlated BC/BMC failure cd=0.06, fit
  independence): Se(severe,5mL) biased **+0.078** (0.66 -> 0.74), Se_BM biased up (0.92 -> 0.98),
  **pi biased -0.076**. So unmodelled positive dependence (e.g., prior antibiotics failing both
  cultures together) **over-states Se and under-states PPV** — a caveat, and a reason to keep the
  optional covariance term as a sensitivity check.

## Files
`stan/bc_bmc_ppv.stan`, `config.yml`, `run_phase1_recovery.R`;
outputs (gitignored): `tables/recovery.csv`, `tables/offset_sweep.csv`,
`figures/fig_recovery_{scalars,pi}.png`, `fig_offset_sweep.png`, `fig_ridge.png`, `fig_gt1.png`,
`outputs/phase1.rds`, `outputs/GATE_PASSED.txt`.

## Decision requested
Recovery is acceptable (coverage 1.00 given the offset; identifiable structure works; the one
non-identified quantity is transparently the fixed offset). **Approve Phase 2 (real data):**
load the 10 Mogasale paired 2x2 (+ the two-volume arms, Vallenas d=15), attach the outbreak line
lists (n_o, k_o, volume), fit with the offset **fixed and swept**, and report per-outbreak
surveillance PPV pi_o with the offset-sensitivity band + a sanity check that surveillance PPVs
sit below hospital phi_s. Or adjust the offset prior/sweep grid first.
