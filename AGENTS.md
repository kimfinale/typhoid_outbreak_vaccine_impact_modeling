# AGENTS.md

## Scope and purpose

These instructions apply to the entire repository. Work here supports a scientific manuscript on
typhoid outbreak-response immunization (ORI), so prioritize reproducibility, source provenance,
and consistency between data, fitted artifacts, manuscript text, and rendered outputs.

Do not assume that counts or terminology in an older README, manuscript PDF, or generated report
are current. The repository is under active revision; verify claims against the active input files
and fitted outputs before carrying them into prose.

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

## Repository orientation

- `R/`: original static ORI model, DALY/cost calculations, figures, and tables.
  `R/2-functions.R` contains shared cost-effectiveness functions.
- `renewal/`: semi-mechanistic renewal-equation model in which indirect effects emerge through
  transmission feedback. The reproducible driver is `renewal/run_analysis.R`.
- `latent_class_ppv/`: latent-class model for blood-culture sensitivity and the PPV of suspected
  clinical case definitions. The main driver is `latent_class_ppv/run_final.R`.
- `manuscript/`: manuscript-facing methods, supplements, evidence extractions, and development
  analyses.
- `revision/`: current rendered manuscript and supplement deliverables. Put newly requested HTML,
  DOCX, PDF, and associated run logs here unless the user specifies another destination.
- `source_decomposition/`: common-source, propagated, and mixed outbreak classification and the
  static/renewal interpretation.
- `case_definitions/`: case-definition harmonization and study-level ascertainment cautions.
- `uncertainty/` and `impact_drivers/`: uncertainty analyses and predictors of ORI impact.
- `papers/`: manuscript PDFs and supplementary tables used to map outbreaks to primary reports.

## Source papers and provenance

Full-text papers are stored in the user's Zotero library, not in this repository:

```text
C:\Users\jonghoon.kim\Zotero\storage\<8-CHAR-KEY>\<Author> - <Year> - <Title>.pdf
```

Find a paper recursively by author or title, for example:

```powershell
Get-ChildItem "C:\Users\jonghoon.kim\Zotero\storage" -Recurse -Filter *.pdf |
  Where-Object { $_.Name -match "Kabwama" }
```

The Zotero metadata export is `references/Typhoid modeling paper.rdf`; it contains bibliographic
metadata but not the attached PDFs. The outbreak-to-primary-report mapping is in the manuscript
supplement under `papers/`.

Always confirm the title, population, location, specimen type, numerator, and denominator inside
the paper. Do not trust filenames alone. Known traps:

- Davis 2018 vs N'Cho 2019 (verified 2026-07-17): the correct MMWR copies are `EDH27L72` (Davis
  2018, Harare 2016-2017) and `RAYYTM4G` (N'Cho 2019, Harare 2017-2018). `CTL54YK3`, labelled
  "Davis 2018", actually contains the N'Cho 2019 paper, and `U68ZE8ZZ` is a yellow-fever paper.
  Never cite `CTL54YK3` or `U68ZE8ZZ` as Davis; do not trust a Davis/N'Cho filename without opening
  the PDF.
- Neil 2012 concerns Kasese, Uganda, not the Democratic Republic of the Congo.
- Distinguish blood-culture counts from stool, bone-marrow, serology, TUBEX, or unlabeled cultures.
- Record whether counts refer to specimens, tested patients, unique confirmed patients, or all
  suspected cases.

For load-bearing claims, check the primary paper even when an extracted table is available. If a
cell is reconstructed rather than directly reported, label it as reconstructed and document the
arithmetic and uncertainty. When explicitly asked to parallelize a large PDF extraction, assign one
paper per agent and use a shared structured extraction rubric.

## Canonical data and count conventions

The canonical outbreak files are:

```text
data/Typhoid_Outbreak_Time_Series_2000_2022_Summary.csv
data/Typhoid_Outbreak_Time_Series_2000_2022_Timeseries.csv
```

Study identifiers are normalized by `renewal/R/data_prep.R::read_summary()` and
`read_timeseries()`. Do not create a second normalization scheme without a clear need.

For diagnostic-yield data, store columns as `tested` and `confirmed` but write proportions in prose
as positive/tested. The corrected PPV anchors are:

- Neil 2012: `tested = 54`, `confirmed = 19`; positivity is **19/54**, not 25/63.
- Kabwama 2017: `tested = 364`, `confirmed = 56`; positivity is **56/364**, not 51/364.

The active stage-1 PPV model input is
`latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv` (community layer);
`latent_class_ppv/data/community_surveillance_ppv.csv` now supplies blood volumes to the audit
build. Do not infer that a saved posterior is current merely because the CSV is correct; compare
timestamps or regenerate the fit.

## Scientific interpretation rules

- Keep **blood-culture positivity** distinct from the **PPV of a clinical case definition**.
  Positivity is approximately `PPV * blood-culture sensitivity` when specificity is near one.
- State the target population explicitly. PPV among selectively cultured suspected cases is not
  automatically PPV among all suspected cases; assess verification/selection bias.
- Treat blood-culture sensitivity as only conditionally transportable. Blood volume, culture media,
  antibiotics, illness timing, laboratory practice, and study heterogeneity may matter.
- Distinguish suspected cases, culture-confirmed cases, and latent true typhoid throughout code,
  tables, and prose.
- Use the outbreak vocabulary `common-source`, `propagated`, and `mixed` consistently. Explain any
  finer point-source or continuous-source distinction when it affects classification.
- The static model provides a conservative bracket when indirect transmission effects are minimal.
  The renewal model generates indirect effects through transmission feedback; do not describe them
  as an imposed fixed indirect-effect parameter.
- Treat multiplicative and additive PPV observation models as structural alternatives unless one is
  empirically identified. Under an additive background, PPV need not cancel from percentage
  reductions.
- Do not interpret outbreak-versus-surveillance contrasts as causal endemicity effects when case
  definitions, testing selection, geography, era, and laboratory practice also differ.
- Match interval labels to computation. Stage 1 commonly reports 90% credible intervals; the
  current stage-2 runner uses 2.5th and 97.5th percentiles, i.e. 95% credible intervals.

## Reproducible workflows

Run commands from the repository root. R 4.5.3 and Quarto are available on the system PATH.
`data.table::fread()` has segfaulted under headless `Rscript` in this environment; prefer base
`read.csv()`, `readRDS()`, or `readxl`. Use an `.R` script for substantial multi-line R code.

### Rebuild the PPV supplement

After any PPV input or model change, run the full dependency chain:

```powershell
Rscript latent_class_ppv/run_final.R
Rscript renewal/run_analysis.R
Rscript revision/run_selection_sensitivity.R
quarto render revision/ppv_supplement.qmd --to html
quarto render revision/ppv_supplement.qmd --to docx --output ppv_supplement_selection_sensitivity.docx
Move-Item ppv_supplement_selection_sensitivity.docx revision/ppv_supplement_selection_sensitivity.docx -Force
```

The first command must regenerate `final_draws.rds`, parameter tables, diagnostics, sensitivity
tables, and PPV figures. The renewal run must then consume the new draws and regenerate
`renewal/tables/tab_ppv_effect.csv` before rendering the supplement.
The sensitivity driver must regenerate the six scenario-specific renewal results,
`renewal/tables/tab_ppv_selection_sensitivity.csv`, and the PPV-to-renewal denominator-frame audit.
Treat `revision/ppv_supplement_selection_sensitivity.docx` as the canonical DOCX; the older
`revision/ppv_supplement.docx` is superseded.

If stage-2 PPV strata results change, also run:

```powershell
Rscript manuscript/fit_ppv_strata.R
```

### Minimum validation after a PPV rebuild

- Confirm the fit log reads Neil `54/19` and Kabwama `364/56` in tested/confirmed columns.
- Confirm `latent_class_ppv/outputs/final_draws.rds` is newer than the corrected input.
- Inspect `latent_class_ppv/tables/fit_diagnostics.csv` for divergences, R-hat, and effective sample
  sizes.
- Confirm the renewal validation tests pass and `tab_ppv_effect.csv` is newer than the posterior.
- Confirm both rendered supplements contain the corrected counts and current posterior summaries.
- Treat warnings separately from failures; report consequential warnings rather than suppressing
  them silently.

For other modules, use their top-level `run_*.R` driver and read the module README before invoking
individual scripts out of order.

## Editing and generated artifacts

- Preserve unrelated working-tree changes. This repository is often intentionally dirty.
- Do not overwrite source manuscripts or user documents unless the requested workflow requires it.
- Edit source `.R`, `.qmd`, `.md`, `.csv`, or configuration files; regenerate HTML/DOCX/PDF outputs
  from those sources rather than hand-editing rendered files.
- Keep manuscript-facing deliverables and useful run/render logs in `revision/`.
- Generated module directories such as `figures/`, `tables/`, and `outputs/` may be gitignored. Their
  presence does not prove they correspond to current inputs.
- Use UTF-8 and check mathematical symbols, en dashes, multiplication signs, and accented author
  names after rendering.
- Do not commit, push, or discard user changes unless explicitly asked.

## Reporting completed work

Lead with the scientific or reproducibility outcome. Name the inputs used, commands run, validation
results, and deliverable paths. Clearly distinguish files edited by the agent from regenerated
artifacts. If an analysis remains limited by selection bias, unidentified structure, missing source
verification, or an unresolved manuscript decision, state that directly rather than presenting the
output as definitive.
