# QUADAS-3 first-pass appraisal

This directory records the 2026-08-21 AI-prefilled, estimate-level QUADAS-3
assessment for the current PPV evidence:

- 26 JHK-approved clinical-definition PPV estimates; and
- nine JHK-approved paired blood/bone-marrow culture estimates.

The appraisal uses the current QUADAS-3 version 1.2 structure. It is not final:
all rows require JHK confirmation and an independent human review before the
evidence freeze. A high judgment does not automatically exclude an estimate
and is not converted into a numerical study weight.

The authoritative working review records are in the external staging bundle:

```text
C:\Users\jonghoon.kim\Workspace\typhoid_ms_20260804
```

The clinical rows were generated from the current human adjudication log,
observation extraction, and model-facing audit in that bundle. The paired rows
were generated from `latent_class_ppv/data/mogasale2016_paired_bc_bmc_seed.csv`.
The staging-bundle amendment log records PPV-009: JHK confirmed targeted efforts
to identify additional paired studies and decided that additional Scopus or Web
of Science searching is not required. The evidence should be described as
prior-review seeding plus targeted supplementary searching, not as a newly
completed exhaustive four-database search.

## Files

- `quadas3_first_pass.xlsx`: user-facing workbook.
- `clinical_ppv_quadas3_first_pass.csv`: 26 clinical estimate-level rows.
- `paired_bc_bm_quadas3_first_pass.csv`: nine paired estimate-level rows.
- `quadas3_first_pass_summary.csv`: domain-level judgment counts.
- `quadas3_signaling_questions.csv`: signaling-question reference.
- `quadas3_tailoring_guidance.md`: synthesis questions, ideal trials, and
  interpretation rules.
- `build_quadas3_first_pass_workbook.ps1`: recreates the workbook from the CSV
  and Markdown files in this directory.

To recreate the workbook before human edits:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File revision/quadas3/build_quadas3_first_pass_workbook.ps1 -Overwrite
```

Do not run that command after reviewers enter decisions unless their editable
fields have first been synchronized to the source CSVs.

## Modeling consequence

The first pass does not change the configured PPV likelihood inputs. QUADAS-3
is used to define sensitivity analyses and interpretation, not to create an
ad hoc quality score or automatic numerical weight. Any input exclusion or
alternative fit should follow human-confirmed consensus judgments and be
recorded as an explicit analysis decision.

