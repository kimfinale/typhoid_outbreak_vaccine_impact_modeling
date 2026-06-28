# CLAUDE.md — typhoid outbreak ORI modeling

## Reference PDFs live in the user's Zotero library (on disk)

The full-text source papers (outbreak reports, vaccine trials, methods refs) are
stored locally by Zotero, **not** in this repo:

```
C:\Users\jonghoon.kim\Zotero\storage\<8-CHAR-KEY>\<Author> - <Year> - <Title…>.pdf
```

- ~2,000 PDFs across one folder per item (random 8-char key, e.g. `2ARGT7WZ`).
- The `Read` tool reads PDFs directly — pass the full path.
- **Find a paper by filename** (author–year is in the filename), e.g. search recursively:
  `Get-ChildItem "C:\Users\jonghoon.kim\Zotero\storage" -Recurse -Filter *.pdf | Where-Object { $_.Name -match "Kabwama" }`
- The Zotero metadata export (BibTeX/RDF) is at `references/Typhoid modeling paper.rdf`
  (titles/authors/some abstracts; **no attached PDFs** — `references/files/` is empty).
- The outbreak → primary-report mapping is in the manuscript **Supplementary Table 1**
  (`papers/Typhoid Fever Modeling Manuscript-supp.*`). 19 outbreaks ≈ 16 papers
  (Walters / Polonsky / Imanishi each cover two outbreaks).

**Caveats when using the Zotero PDFs:**
- **Verify the content matches the filename** — some items are misfiled. Known bad:
  `storage\U68ZE8ZZ\Davis et al. - 2018 …Typhoid… Harare.pdf` actually contains a
  *Yellow Fever in Travelers to Brazil* paper. Always confirm title/topic after opening.
- Filenames also state location only loosely — confirm against the paper. Known fix:
  **Neil 2012 is Kasese, Uganda** (intestinal-perforation epidemic), not DRC.
- For a multi-PDF extraction pass, parallelize with subagents (one paper each, Read +
  a structured rubric) — see the transmission-mode classification in `source_decomposition/`.

## Repo orientation (analysis modules; each self-contained, gitignored outputs)

- `R/` — original manuscript pipeline (static model, DALY/cost, figures). `R/2-functions.R`
  holds the reused CEA functions (`compute_yld/yll`, `add_cea_results`).
- `renewal/` — renewal-equation (R_t) ORI model + Sobol pipeline (`run_analysis.R`).
- `revision/` — manuscript-revision outputs filling margin comments c0–c9.
- `uncertainty/` — rigorous R_t analysis (GI/asymptomatic/PPV propagation + Bayesian).
- `impact_drivers/` — which outbreak characteristics predict ORI impact.
- `source_decomposition/` — transmission-mode classification (common-source vs propagated)
  + static↔renewal bracket.

## Environment notes

- R 4.5.3 at `C:\Program Files\R\R-4.5.3\bin\Rscript.exe` (on System PATH).
- `data.table::fread` segfaults under Rscript here — use base `read.csv` / readxl in
  headless scripts. Run R via PowerShell; for multi-line code use a `.R` file, not `-e`.
- Canonical data: `data/Typhoid_Outbreak_Time_Series_2000_2022_{Summary,Timeseries}.csv`
  (human-readable headers; summary uses `(1)/(2)` study-id naming → normalized by
  `renewal/R/data_prep.R::read_summary` / `read_timeseries`).
- Each module's `run_*.R` reproduces its outputs; figures/tables/outputs are gitignored.
