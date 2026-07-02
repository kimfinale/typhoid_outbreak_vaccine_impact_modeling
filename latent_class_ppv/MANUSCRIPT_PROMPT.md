# Claude Code task — draft a manuscript: accuracy of typhoid case-ascertainment tools (clinical case definitions, blood culture, bone-marrow culture) via a Bayesian latent-class synthesis

**One-line intent.** Produce a submission-ready **draft manuscript + reproducible code + a references.bib**
that consolidates, in one Bayesian latent-class framework, the **sensitivity, specificity, and PPV**
of the tools used to ascertain typhoid fever cases — the **clinical/suspected case definition**,
**blood culture**, and **bone-marrow culture** — across **hospital** and **community** settings, and
quantifies how a case definition's **PPV varies across the definition spectrum and prevalence**.
Antibody-based rapid tests are OUT of scope (cite Arora 2019).

This runs autonomously for a while. Work on a branch, commit at every checkpoint, and **stop-and-note**
rather than fabricate when data or convergence fail.

---

## 0. Environment & orientation (read first)
- R 4.5.3 at `C:\Program Files\R\R-4.5.3\bin\Rscript.exe`; cmdstanr + cmdstan 2.39 installed. Run R via
  PowerShell with script files (not `-e`); `data.table::fread` segfaults under Rscript — use `read.csv`.
- Reference PDFs are in the user's **Zotero** library on disk (see `CLAUDE.md` for the path pattern and
  the **misfiling caveat** — always verify a PDF's content matches its filename).
- **Build on existing assets — do NOT restart.** The `latent_class_ppv/` module already contains most of
  this analysis. Read these before doing anything:
  - `FINAL_REPORT.md` (fitted results), `README.md` (model + assumptions + estimand split),
    `PHASE0_novelty_memo.md` (novelty vs Mogasale/Arora/Antillon), `PHASE1_RECOVERY.md` (the recovery
    gate), `COMMUNITY_PPV.md`, `SE_BC_MILD_LITSEARCH.md` (why the severity offset was dropped),
    `PPV_SPECTRUM.md` (PPV vs definition x prevalence).
  - `stan/model_final.stan` (the fitted model), `run_final.R`, `ppv_spectrum.R`.
  - `data/`: `mogasale2016_paired_bc_bmc_seed.csv`, `antillon2018_volume_sensitivity.csv`,
    `community_surveillance_ppv.csv`, `clinical_definition_accuracy.csv`.
  - Broader repo: `case_definitions/` (case-definition observation-operator typology), the outbreak
    Summary/Timeseries data under `data/`.

## 1. Scope decisions (FIXED — do not relitigate)
- **IN:** (a) clinical/syndromic case-definition accuracy (Se_clin, Sp_clin) and PPV; (b) blood-culture
  sensitivity as a function of **blood volume** + specificity; (c) bone-marrow-culture sensitivity; (d)
  the **hospital vs community PPV gap**; (e) the **PPV-vs-definition-spectrum x prevalence** surface.
- **OUT:** antibody/serologic rapid tests (Widal / TUBEX / Typhidot) as objects of estimation. **Cite
  Arora 2019** (PLoS NTD 13:e0007303) for their accuracy and use those values only as *inputs* where a
  serologic "probable" tier illustrates PPV. Do not re-estimate them.
- **Estimand framing (keep explicit):** transportable *test* properties (Se_BC(volume), Se_BM, Sp) vs
  *local* population PPVs (phi_s hospital, pi_o community). PPV is never a single transferable number.

## 2. Data to consolidate (with a data-inventory gate)
Use parallel subagents for extraction (one paper each, structured rubric), as the module already did.
1. **Paired blood + bone-marrow (culture-accuracy layer).** Start from Mogasale's 10 studies (already
   extracted). Systematically look for **additional** paired blood+bone-marrow studies (Zotero + web),
   extract the 2x2 (BC+BM+, BC+BM-, BC-BM+, BC-BM-) + blood volume. Correct known issues (Vallenas
   both-negative = 15; impute/flag missing volumes).
2. **Clinical-definition accuracy layer (the expansion).** Find studies that culture a **defined febrile
   denominator** and report **clinical-definition status** (so a definition-vs-culture 2x2 exists).
   Anchors: Thriemer 2012 (Pemba, WHO suspected 0.83/0.36 and probable 0.41/0.997), Dong 2010 (Hechi),
   Storey 2015 (serology Se/Sp), Hosoglu 2006, Ross & Abraham 1987. Extract each 2x2 or reported Se/Sp +
   reference standard + setting + definition wording.
3. **Outbreak/community positivity (PPV/observation layer).** The tight community-syndromic set (Aye,
   Neil, Kabwama) + any additional clean syndromic community points from the Summary sheet.
4. **Antibody Se/Sp** — from Arora 2019 only (do not re-extract primaries).
Write a **data-inventory table** (source, layer, n, what was extracted, Zotero-vs-web, confidence) and
sanity-check counts before modelling.

## 3. Model
- **Culture + PPV layers:** use `stan/model_final.stan` as-is (volume-anchored, beta-free): historic
  paired 2x2 -> Se_BC(volume), Se_BM, Sp; outbreak positivity -> phi_s, pi_o.
- **Optional clinical-definition layer:** if >=3-4 definition-vs-culture 2x2 studies are found, extend
  the LCM so the **clinical definition is a third imperfect test** applied to the febrile denominator
  (blood culture is imperfect, so this is a latent-class, not a fixed-reference, problem) -> estimate
  Se_clin, Sp_clin within the framework. If data are too thin, fall back to a **bivariate (Reitsma)
  meta-analysis** of Se_clin/Sp_clin citing the primaries, and say so plainly.
- **MANDATORY gate:** a **simulation-recovery** check before any real-data claim (reuse the
  `run_phase1_recovery.R` pattern for whatever structure you fit). Respect PASS/FAIL; report divergences,
  max R-hat, min ESS; **do not present a non-converged fit**.

## 4. Analyses to finalize
- Refit the consolidated model; report Se_BC at 2/5/10 mL, Se_BM, Sp, phi_s, pi_o (and Se_clin/Sp_clin
  if fitted). Reproduce Mogasale/Antillon/Thriemer where overlapping (validation) and note divergences.
- Finalize the **PPV-vs-definition-spectrum x prevalence** surface (`ppv_spectrum.R`) with consolidated
  inputs; overlay the empirical anchors (Pemba WHO defs; our outbreak PPVs on the syndromic curve).
- Sensitivity analyses: conditional dependence (antibiotic-driven joint culture failure); specificity
  relaxation; the outbreak blood-volume assumption; the **Widal-specificity-degrades-in-endemic-areas**
  caveat for the serologic-tier PPV.

## 5. Manuscript
- **Format:** Quarto (`manuscript/typhoid_diagnostic_accuracy.qmd`) -> PDF + docx. Target a
  methods-applied journal (PLoS NTD / J Infect Dis / BMC Infect Dis) — pick one and format loosely to it.
- **Structure:**
  - *Abstract* (structured: background/methods/results/conclusions with the headline numbers).
  - *Introduction*: the layered ascertainment problem (`observed = rho * (lambda * d)`); what Mogasale
    (pooled Se), Arora (test accuracy, fixes culture Se, no PPV), and Antillon (volume) did, and the gap
    we fill (joint volume-anchored culture Se + Se_BM + clinical-definition PPV, transportable/local
    split, hospital-vs-community).
  - *Methods*: data sources + inventory; the latent-class model(s) + priors + estimands; the identifiability
    structure (the product ridge; culture pins the parts; PPV level rides on the Se anchor); the recovery
    gate; why a severity offset was dropped (cite the lit search).
  - *Results*: Se_BC(volume); Se_BM; specificity; clinical Se/Sp (fitted or meta-analytic); hospital phi_s
    vs community pi_o; the PPV surface. Tables + figures below.
  - *Discussion*: implications for burden estimation / CEA (the confirmation adjustment is not one number),
    for surveillance (a high-PPV definition needs a serologic tier, not more clinical criteria), the
    **mild-case blind spot** (no paired culture data can exist in the community — serology-based
    surveillance), and limitations.
  - *Data & code availability* (point to the repo).
- **Tables:** (T1) parameter estimates with 90% CrI; (T2) definition-accuracy spectrum (Se/Sp + source);
  (T3) PPV by definition x prevalence.
- **Figures:** (F1) Se_BC vs volume with historic points; (F2) hospital vs community PPV forest; (F3) PPV
  spectrum surface; (F4) study/data flow (PRISMA-ish) or a study-level forest. Regenerate all from code.
- **Every number must be traceable to a code output; every citation must be real** (assemble
  `manuscript/references.bib`).

## 6. Honesty & quality rules (hard constraints)
- **Never fabricate** a number, result, or citation. Mark "not reported" and data gaps explicitly.
- **Distinguish** quantities *estimated from our data* from those *taken from the literature* (and which).
- Report the **identifiability caveats** honestly (product ridge; PPV level = Se-prior-dependent; the
  spread/ordering are data-driven). Present the severity-offset drop and its evidence.
- **Verify Zotero PDFs** match filenames before quoting (known misfiling; see CLAUDE.md). Prefer the
  user's Zotero copies; for web sources, record the DOI/PMC ID and flag them as web-retrieved.
- If a planned analysis isn't supported by data, **downgrade the claim** (meta-analysis / literature
  value / "not estimable") rather than overreach.
- Reproduce known results where they overlap as **validation**; state divergences.

## 7. Deliverables
- `manuscript/typhoid_diagnostic_accuracy.qmd` + rendered `.pdf`/`.docx`.
- `manuscript/references.bib` (all cited sources; Zotero-vs-web flagged in a comment).
- Consolidated data CSVs + a single reproducible driver (`manuscript/run_manuscript.R`) that regenerates
  every table/figure and the fit.
- `manuscript/README.md` (how to reproduce; the data-inventory; open gaps/TODOs for the human author).

## 8. Checkpoints (commit at each; clear messages)
1. **Data inventory** complete -> commit the inventory table + consolidated CSVs; sanity-check counts.
2. **Recovery gate** run -> commit diagnostics; if not converged, STOP and write a NEEDS_REVIEW note.
3. **Real-data fit** done -> commit estimates + figures.
4. **Manuscript draft** done -> self-review against Section 6; commit; write an open-issues/TODO list for
   the human (missing data, weak claims, decisions needed).
Use parallel subagents / a workflow for the literature extraction and verification passes where it saves
wall-clock. Keep outputs gitignored except curated data, code, and the manuscript.

## 9. What "done" looks like
A coherent draft the human can read end-to-end: headline estimates (Se_BC ~0.52-0.68 by volume, Se_BM
~0.90, Sp ~0.99, hospital PPV ~0.8 vs community ~0.25, the PPV-spectrum surface), each backed by a table
or figure regenerated from committed code, each citation in the .bib, the honesty caveats stated, and a
TODO list of what still needs a human decision before submission.
