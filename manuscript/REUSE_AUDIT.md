# REUSE AUDIT — infrastructure to reuse vs evidence base to rebuild

The `latent_class_ppv/` module contains validated analytical **machinery**. This run reuses the
machinery and **rebuilds the evidence base** via a systematic search, then **re-derives every number**.
No prior point estimate is carried forward. (There is no `manuscript/` draft yet — the referenced
"typhoid_diagnostic_accuracy / 2026-07-02 manuscript" does not exist in this repo; the reusable assets
are the `latent_class_ppv/` files below.)

## REUSE AS-IS — method, not answer (extend, do not rewrite)

| File | What it provides | Reuse plan |
|---|---|---|
| `latent_class_ppv/stan/model_final.stan` | Hui–Walter paired-culture LCM (Se_BC, Se_BM, Sp pinned) + volume-anchored Se_BC sub-model + outbreak binomial → hospital φ_s / community π_o; transportable-Se / local-PPV split | **Extend** the *structure*; refit on the systematically-assembled data. Add per-definition layer only if fallback K cleared. |
| `latent_class_ppv/stan/bc_bmc_ppv.stan` | Recovery-gate model; carries the **optional conditional-dependence** term (`cd_infected`) and a **legacy `beta`** (severity offset, DROPPED — keep off) | Reuse the **conditional-dependence** machinery for the sensitivity analysis. `beta` stays disabled (see `SE_BC_MILD_LITSEARCH.md`). |
| `latent_class_ppv/stan/community_ppv.stan` | Standalone partial-pooled community-PPV | Reuse if a pooled community summary is needed. |
| `latent_class_ppv/run_phase1_recovery.R` | **The simulation-recovery GATE** (coverage/bias/convergence; ridge, ">1", dependence-stress demos) | Reuse the **recovery-gate pattern**. *(Prompt calls this `run_recovery_final.R` — it does not exist by that name; this is the file.)* |
| `latent_class_ppv/run_final.R` | Consolidated real-data fit + inferences + figures + **volume-sensitivity sweep** | Reuse. *(Prompt's `run_sensitivity.R` does not exist; the sensitivity logic lives HERE.)* |
| `latent_class_ppv/ppv_spectrum.R` | PPV = f(definition Se/Sp, prevalence) surface + empirical anchors | Reuse; **re-run with the systematically-assembled Se/Sp inputs.** |
| `latent_class_ppv/data/*.csv` | prior seed data (Mogasale paired, Antillón volume/Se, community positivity, clinical Se/Sp) | **Seed / cross-check only** — the study set of record must come from the systematic search. |

Sound analytical concepts to keep: Hui–Walter identification from paired imperfect tests; volume-anchored
Se_BC; **Sp ≈ 1 pinned** (a positive S. Typhi culture is definitionally a true case); the single-outbreak
**product ridge** (π·Se_BC) broken by the paired-culture anchor; transportable-Se / local-PPV split; the
mandatory recovery gate.

## REBUILD SYSTEMATICALLY — the DATA / evidence base

- **The study set and its cells.** The current data (Mogasale's 10 paired 2×2; the outbreak-positivity
  points Aye/Neil/Kabwama; the clinical Se/Sp from Thriemer/Storey/Dong/Hosoglu) came from **targeted,
  hand-picked extraction**. For a PRISMA-DTA systematic review these must be re-assembled from a
  **registered protocol + database searches + dual screening** (the human runs the searches; see the
  HARD GATE). The prior CSVs become a *candidate pool* to seed extraction and cross-check the search —
  **not the search of record.**
- **The entire PRISMA-DTA apparatus** (protocol, per-layer strings, screening log, flow diagram,
  QUADAS-2): none exists — build it.

## RE-DERIVE — every number

The recovery gate re-runs; the model refits on the expanded data; estimates may move. **Do NOT carry
forward** the prior targeted-analysis numbers (Se_BC 0.52/0.62/0.68, Se_BM 0.90, hospital φ ~0.83,
community π ~0.25, the specific PPV-surface values). They are superseded pending the systematic evidence
base, and any movement is reported plainly.

## Rebuild targets (checklist)
1. Registered protocol + per-layer search strings (Phase 2). 
2. Candidate pool CSV (Phase 2, seed/cross-check only).
3. Screening log + included PDFs (**human**, at the gate).
4. Systematically-extracted 2×2 / positivity / Se-Sp dataset + QUADAS-2 (Phase 3–4).
5. Refit + recovery gate + PPV surface on the new data (Phase 5).
6. PRISMA-DTA flow from the human's counts; manuscript (Phase 6).
