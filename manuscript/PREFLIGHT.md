# PREFLIGHT — environment

**Note on machine.** This run is executing on the **local development machine** (Windows 11), not a
fresh remote. The toolchain is already provisioned and was exercised throughout the prior
`latent_class_ppv/` analysis, so most preflight items are already green. (The remote-machine preflight
in `MANUSCRIPT_PROMPT_v2.md` still applies if/when this is re-run elsewhere.)

| Item | Status |
|---|---|
| R | **4.5.3** at `C:\Program Files\R\R-4.5.3\bin\Rscript.exe` (on PATH). OK (>= 4.4). |
| C++ toolchain / RTools | RTools45 present; used to compile all prior Stan models. OK. |
| cmdstanr + cmdstan | **cmdstan 2.39.0** installed; used successfully for `model_final.stan`, recovery gate, etc. OK. |
| R packages (core) | cmdstanr, posterior, ggplot2, dplyr, tidyr, yaml, scales, knitr, rmarkdown — present (used in prior runs). OK. |
| R packages (DTA meta-analysis) | **meta / metafor / mada — NOT yet verified/installed.** Only needed IF the clinical-layer bivariate/HSROC is fit (Phase 5, post-gate, and only if the fallback threshold K is cleared). **Deferred** — install at Phase 5 if triggered. |
| Quarto | Available (prior `.qmd` reports rendered, e.g. `case_definitions/report_case_definitions.qmd`). OK. |
| Zotero PDFs | **Local library present** (~2,000 PDFs under `C:\Users\jonghoon.kim\Zotero\storage\`). So **Zotero-local mode**, not web-only — better than the fresh-remote case. Per CLAUDE.md, verify each PDF's content matches its filename (known misfiling). |

**Env notes (Windows):** invoke R via script files (not `-e`); `data.table::fread` segfaults under
Rscript — use `read.csv`. Run R via PowerShell.

**SETUP_NEEDED:** none for Phases 0–4. The only deferred item is `meta`/`metafor`/`mada`, installed on
demand at Phase 5 if the clinical bivariate/HSROC layer is triggered.

**Path decision:** review-specific data (candidate pool, QUADAS-2, extraction flags, consolidated CSVs,
included PDFs) live under **`manuscript/data/`** to avoid polluting the canonical outbreak `data/` at
repo root. `MANUSCRIPT_PROMPT_v2.md` writes `data/…`; read that as `manuscript/data/…` here.
