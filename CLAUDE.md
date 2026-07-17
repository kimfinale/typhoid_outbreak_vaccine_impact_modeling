# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> `AGENTS.md` is the Codex-facing equivalent and covers overlapping ground. It is untracked, so
> **this file must stand alone** — don't turn its content into a pointer. If you change a rule that
> lives in both, update both. Where the two disagree on a Zotero misfiling, this file is the more
> current record.

This repo backs a scientific manuscript on typhoid outbreak-response immunization (ORI).
Counts in an older README, a manuscript PDF, or a rendered report are **not** presumed current —
verify against the active inputs and fitted outputs before carrying a number into prose.

## AI collaboration, attribution, and Git workflow

Claude Code and Codex may both be used on this repository. Attribution is reliable only when each
agent works in its own branch and worktree. Do not allow two AI agents to edit the same branch or
worktree concurrently.

### Before editing

1. Read the repository instructions that apply to the requested files.
2. Report the current branch and worktree path.
3. List the files expected to change and give a brief implementation plan.
4. Inspect the working tree and preserve all pre-existing changes. If another agent appears to be
   active in the same worktree, or the requested files contain unattributed changes, stop and ask
   the user before editing.

For new tasks, use an agent-specific task branch and preferably a separate worktree:

- Codex: `codex/<task-name>`
- Claude Code: `claude/<task-name>`

Never work directly on `main`. Do not create or switch branches or worktrees when doing so would
disrupt uncommitted work; report the condition and ask the user how to proceed.

### During editing

- Work on one task only and make minimal, focused changes.
- Modify only the files named in the pre-edit scope unless a newly discovered dependency is required;
  explain any scope expansion before making it.
- Preserve existing style and avoid formatting-only rewrites.
- Do not change APIs, equations, scientific assumptions, parameter definitions, units, or default
  values unless the user explicitly requests the change.
- Distinguish bug fixes, numerical improvements, statistical changes, and scientific/modeling
  changes in progress updates and the final report.

### Commits and attribution

A request to edit or run analyses does not authorize a commit. Commit only when the user explicitly
asks to commit or to finish and commit. Do not push, open a pull request, amend, squash, rebase, or
rewrite another agent's commits unless explicitly requested.

Every AI-created commit must be focused and include the appropriate trailer without changing the
repository or global Git identity:

```text
AI-Agent: Codex
```

or

```text
AI-Agent: Claude Code
```

Do not remove or alter the other agent's attribution trailer. When no commit was authorized, state
that the changes remain uncommitted; exclusive worktree ownership is then the attribution record.

### Validation and handoff

After editing:

1. Run tests or validation proportional to the change.
2. Show `git diff --stat` for the task files and inspect the substantive diff.
3. Summarize every modified file, including regenerated artifacts.
4. Report tests and results, consequential warnings, limitations, and remaining concerns.
5. If a commit was authorized, report its hash and confirm whether task-related changes remain
   uncommitted. If no commit was authorized, explicitly say so.

For review-only requests, do not edit unless the user subsequently asks for implementation. Separate
definite problems from optional suggestions, and focus on correctness, numerical stability,
statistical validity, scientific assumptions, missing tests, and unintended side effects.

## Commands

R 4.5.3 and Quarto are on the System PATH. Run everything from the repo root.

Each module is reproduced by its own top-level driver — read the module README before running
scripts out of order:

```powershell
Rscript renewal/run_analysis.R              # renewal R_t model + Sobol
Rscript latent_class_ppv/run_final.R        # Se_BC + PPV latent-class fit (stage 1)
Rscript uncertainty/run_uncertainty.R
Rscript impact_drivers/run_impact_drivers.R
Rscript source_decomposition/run_source_decomposition.R
Rscript case_definitions/run_case_definitions.R
Rscript rtr_anchors/run_rtr_anchors.R
Rscript revision/run_revision_outputs.R
Rscript hierarchical_theta/run_phase1_gate.R   # gate; run_phase2_3.R only if the gate passes
```

`R/` has no driver — source `1-data.R` → `2-functions.R` → `3-model.R` in order, then
`figures.R` / `tables.R`.

### PPV rebuild chain (ordered — later steps consume earlier outputs)

After **any** PPV input or model change, run the whole chain; do not render from stale draws:

```powershell
Rscript latent_class_ppv/run_final.R        # regenerates final_draws.rds + tables/figures
Rscript renewal/run_analysis.R              # consumes new draws -> renewal/tables/tab_ppv_effect.csv
Rscript revision/run_selection_sensitivity.R # propagates six selection scenarios + denominator-frame audit
quarto render revision/ppv_supplement.qmd --to html
quarto render revision/ppv_supplement.qmd --to docx --output ppv_supplement_selection_sensitivity.docx
Move-Item ppv_supplement_selection_sensitivity.docx revision/ppv_supplement_selection_sensitivity.docx -Force
```

If stage-2 PPV strata results change, also run `Rscript manuscript/fit_ppv_strata.R`.

Validate after a rebuild: fit log reads Neil `54/19` and Kabwama `364/56`; `final_draws.rds` is
newer than the corrected input; `tab_ppv_effect.csv` is newer than the posterior; check
`latent_class_ppv/tables/fit_diagnostics.csv` for divergences / R-hat / ESS. A correct input CSV
does **not** imply a current posterior — compare timestamps or refit.
The canonical PPV DOCX is `revision/ppv_supplement_selection_sensitivity.docx`; the older
`revision/ppv_supplement.docx` is superseded.

## Repo orientation

Analysis modules are self-contained, each with a README and gitignored `figures/`/`tables/`/
`outputs/`. Those directories existing does not prove they match current inputs.

| Module | What it is |
|---|---|
| `R/` | Original manuscript pipeline: static model, DALY/cost, figures. `R/2-functions.R` holds the shared CEA functions (`compute_yld/yll`, `add_cea_results`). |
| `renewal/` | Renewal-equation (R_t) ORI model + Sobol. Indirect effects **emerge** from transmission feedback — not an imposed fixed parameter. |
| `latent_class_ppv/` | Bayesian latent-class fit for blood-culture sensitivity `Se_BC` (transportable, volume-anchored) and clinical case-definition PPV (local). `FINAL_REPORT.md` is the write-up. |
| `hierarchical_theta/` | Hierarchical renewal-with-source model for per-outbreak source fraction θ (R + Stan). **See `IDENTIFIABILITY_FAILED.md` before reusing** — the phase-1 gate is the entry point. |
| `source_decomposition/` | Common-source / propagated / mixed classification + the static↔renewal bracket. |
| `case_definitions/` | Case-definition harmonization and study-level ascertainment cautions. |
| `rtr_anchors/` | Literature priors for transmission R_tr. |
| `uncertainty/` | Rigorous R_t analysis: GI / asymptomatic / PPV propagation + Bayesian. |
| `impact_drivers/` | Which outbreak characteristics predict ORI impact. |
| `manuscript/` | Manuscript-facing methods, supplements, evidence extraction, screening rubrics, POC analyses. Many loose scripts — no single driver. |
| `revision/` | **Destination for new deliverables.** Current rendered manuscript/supplement HTML, DOCX, PDF + run logs. |
| `papers/` | Manuscript PDFs + Supplementary Table 1 (outbreak → primary-report mapping). |

Key identifiability point behind `latent_class_ppv/`: one test sees only the *product*
`positivity = pi * Se_BC`, so a single outbreak cannot separate PPV from sensitivity
(naive `positivity / 0.6` can exceed 1 — that is a falsification signal, not a PPV). Studies that
cultured the same patients by **both blood and bone marrow** are the Hui-Walter leverage that pins
`Se_BC`.

## Reference PDFs live in the user's Zotero library (on disk)

Full-text sources (outbreak reports, vaccine trials, methods refs) are stored by Zotero, **not** in
this repo:

```
C:\Users\jonghoon.kim\Zotero\storage\<8-CHAR-KEY>\<Author> - <Year> - <Title…>.pdf
```

- ~2,000 PDFs, one folder per item (random 8-char key, e.g. `2ARGT7WZ`). `Read` opens PDFs directly.
- **Find a paper by filename** (author–year is in the filename):
  `Get-ChildItem "C:\Users\jonghoon.kim\Zotero\storage" -Recurse -Filter *.pdf | Where-Object { $_.Name -match "Kabwama" }`
- Metadata export: `references/Typhoid modeling paper.rdf` (titles/authors/some abstracts;
  **no attached PDFs** — `references/files/` is empty).
- Outbreak → primary-report mapping: manuscript **Supplementary Table 1** in `papers/`.
  19 outbreaks ≈ 16 papers (Walters / Polonsky / Imanishi each cover two).

**Caveats when using the Zotero PDFs:**
- **Verify the content matches the filename** — always, on open. Status as of 2026-07-16:
  Davis 2018 **fixed**; **Poncin is two separate, correctly-linked papers** (cholera + typhoid
  — never a misfiling); **Valenciano 2000 fixed — use `ED36QT97`** (verified: AJE 152(10):934–9,
  floating restaurant, Seine). **Do NOT use `JEMF2WZ7` or `AW7IPIVI`** — both still contain
  Borgdorff 2001 on *M. tuberculosis* and should be deleted. The cause was publisher-side:
  Valenciano is AJE **152**(10):934–9, Borgdorff is **154**(10):934–43 (same journal, issue and
  start page, different volume), and OUP's URL `academic.oup.com/aje/article/152/10/934/55560`
  serves the Borgdorff PDF — so two independent downloads a decade apart were identically wrong.
  A working copy exists, so this class of error is recoverable, but **never trust a same-journal
  same-page PDF without opening it**.
- Filenames also state location only loosely — confirm against the paper. Known fix:
  **Neil 2012 is Kasese, Uganda** (intestinal-perforation epidemic), not DRC.
- **Counts need the specimen checked, not just the number.** Neil 2012 is **19/54 blood
  cultures**; the 63/25 once in the data was blood summed with 6/9 cultures from an
  unlabelled (blood *or* stool) source. `Se_BC` is blood-culture-specific, so pooled
  blood/stool/urine counts are inadmissible. This is the single most common defect in the
  outbreak literature — see `latent_class_ppv/data/appiah2020_unique_outbreak_ppv_audit.csv`.
- **Davis 2018 vs N'Cho 2019 (verified 2026-07-17).** Correct MMWR copies: `EDH27L72` = Davis 2018
  (Harare **2016–2017**, MMWR 67:342–3) and `RAYYTM4G` = N'Cho 2019 (Harare **2017–2018**, MMWR
  68[2]:44–5). Folder `CTL54YK3`, filename-labelled "Davis 2018", actually **contains the N'Cho
  2019 paper** — relabel/delete it; never cite `CTL54YK3` as Davis. Same-city, adjacent-season pair
  that swap easily.
- **N'Cho 2019's 74/286 is not blood-disaggregated.** The primary says "cultures were processed for
  286 (49%) inpatients; 74 (26%) yielded Typhi" under a blood/stool/rectal-swab confirmation
  definition — a mixed-specimen positivity, not a verified blood-culture pair. It is **kept in the
  main model** (real tested/positive pair in the 583-patient frame) but flagged for a
  leave-N'Cho-out sensitivity check. The other Harare reports (Polonsky 2014, Muti 2014, Imanishi
  2014, Davis 2018, Hechaichi 2023) all pool blood with stool/urine/bone marrow and are excluded on
  that basis — see the 2026-07-17 reconstruction verdicts in `revision/build_merged_ppv_study_audit.R`.
- For a multi-PDF extraction pass, parallelize with subagents (one paper each, Read +
  a structured rubric) — see the transmission-mode classification in `source_decomposition/`.

## Count conventions and interpretation rules

These are the errors that actually recur in this project:

- **Store `tested` / `confirmed`; write proportions in prose as positive/tested.** Corrected PPV
  anchors: Neil 2012 = **19/54** (not 25/63); Kabwama 2017 = **56/364** (not 51/364).
- **Blood-culture positivity ≠ PPV of a clinical case definition.** Positivity ≈ `PPV × Se_BC` when
  specificity ≈ 1. Keep suspected cases, culture-confirmed cases, and latent true typhoid distinct
  in code, tables, and prose.
- **State the target population.** PPV among *selectively cultured* suspected cases is not PPV among
  *all* suspected cases — assess verification/selection bias.
- **`Se_BC` is only conditionally transportable** — blood volume, media, antibiotics, illness timing,
  and lab practice matter. `pi` (PPV) is local and does not transport.
- **The static model is a conservative bracket** when indirect effects are minimal; the renewal model
  *generates* indirect effects through transmission feedback. Never describe renewal indirect effects
  as an imposed fixed parameter.
- **Multiplicative vs additive PPV observation models are structural alternatives** unless one is
  empirically identified. Under an additive background, PPV need not cancel from percentage reductions.
- **Outbreak-vs-surveillance contrasts are not causal endemicity effects** — case definitions, testing
  selection, geography, era, and lab practice all differ.
- **Match interval labels to computation.** Stage 1 commonly reports 90% CrIs; the current stage-2
  runner uses 2.5th/97.5th percentiles, i.e. 95% CrIs.
- Active stage-1 PPV model input: `latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv`
  (community layer); `community_surveillance_ppv.csv` now supplies blood volumes to the audit build.
- For load-bearing claims, check the primary paper even when an extracted table exists. If a cell is
  **reconstructed** rather than reported, label it so and document the arithmetic and uncertainty.

When reporting work: lead with the scientific or reproducibility outcome; name inputs, commands run,
validation results, and deliverable paths; distinguish files you edited from regenerated artifacts.
If an analysis is still limited by selection bias, unidentified structure, or missing source
verification, say so rather than presenting the output as definitive. Report consequential warnings
rather than suppressing them.

## Environment notes

- R 4.5.3 at `C:\Program Files\R\R-4.5.3\bin\Rscript.exe` (on System PATH).
- `data.table::fread` **segfaults under Rscript here** — use base `read.csv` / `readRDS` / `readxl`
  in headless scripts. Run R via PowerShell; for multi-line code use a `.R` file, not `-e`.
- Canonical data: `data/Typhoid_Outbreak_Time_Series_2000_2022_{Summary,Timeseries}.csv`
  (human-readable headers; summary uses `(1)/(2)` study-id naming → normalized by
  `renewal/R/data_prep.R::read_summary` / `read_timeseries`). Don't add a second normalization
  scheme without a clear need.
- Stan modules (`latent_class_ppv/`, `hierarchical_theta/`) use **cmdstanr**.
- The working tree is **often intentionally dirty** — preserve unrelated changes, and don't commit,
  push, or discard without being asked.
- Edit sources (`.R`, `.qmd`, `.md`, `.csv`); regenerate HTML/DOCX/PDF rather than hand-editing
  rendered output. UTF-8 — check math symbols, en dashes, and accented author names after rendering.
