# QUADAS-3 first-pass appraisal: tailoring and review instructions

Version: first pass, 2026-08-21  
Status: AI-prefilled; JHK confirmation and an independent human review remain required.

## Why QUADAS-3 is used

QUADAS-3 is the current recommended QUADAS version. It assesses risk of bias
and applicability at the level of each accuracy estimate rather than assigning
one judgment to an entire report. The official tool has four domains:
Participants, Index test, Target condition, and Analysis. Applicability is
assessed for the first three domains, followed by an overall judgment.

Official sources:

- University of Bristol QUADAS-3 page:
  https://www.bristol.ac.uk/population-health-sciences/projects/quadas/quadas-3/
- QUADAS-3 tool structure and version 1.2:
  https://www.bristol.ac.uk/population-health-sciences/projects/quadas/quadas-3/quadas-3-tool/
- Whiting et al. 2026, Annals of Internal Medicine:
  https://doi.org/10.7326/ANNALS-25-02104

This project does not contain conventional diagnostic-accuracy studies alone.
The tool is therefore tailored transparently to two related evidence layers.
The signaling questions retain the QUADAS-3 structure, while the rationales
state how each question was interpreted for this model.

## Synthesis question 1: clinical-definition PPV

Among people meeting a suspected typhoid clinical definition during a bounded
community outbreak in a World Bank FY2027 low- or lower-middle-income setting,
what proportion truly has *Salmonella Typhi* infection?

Ideal test-accuracy trial:

- Participants: a prospective, consecutive or population-complete cohort of
  all people meeting a prespecified suspected-case definition in the intended
  outbreak population.
- Index test: the prespecified suspected clinical case definition, applied and
  recorded before culture results are known.
- Target condition: latent true *S. Typhi* infection, assessed uniformly in
  every participant with a sufficiently sensitive reference strategy.
- Analysis: unique participants with aligned suspected, tested, and positive
  counts; missing outcomes and verification selection are handled explicitly;
  overlapping cohorts are not counted as independent.

Important adaptation: blood culture is highly specific but imperfectly
sensitive. The first pass therefore rates the Target condition domain high risk
for all raw clinical-yield estimates. The downstream latent model adjusts for
blood-culture sensitivity, but that model-based adjustment does not erase the
primary-study reference-standard and verification limitations.

## Synthesis question 2: paired blood and bone-marrow culture

Among clinically suspected typhoid patients receiving paired blood and
bone-marrow cultures, what is blood-culture sensitivity as a function of blood
volume and between-study heterogeneity?

Ideal test-accuracy trial:

- Participants: a prospective consecutive cohort representative of the full
  suspected-typhoid spectrum.
- Index test: blood culture collected before treatment, with blood volume,
  timing, media, incubation, and interpretation reported.
- Target condition: latent true *S. Typhi* infection informed by both cultures
  and any suitable additional reference information; bone-marrow culture is not
  assumed to be an infallible gold standard.
- Analysis: one paired 2x2 table per independent cohort, one row per unique
  participant, explicit missingness, and an analysis that acknowledges the
  absence of a perfect reference standard and possible conditional dependence.

Important adaptation: the paired latent-class model intentionally does not
treat bone-marrow culture as perfect. A high Target condition risk judgment is
therefore an expected structural warning, not an instruction to discard all
paired studies. It should motivate sensitivity analyses and cautious
interpretation.

## How the first pass was produced

- The 26 clinical rows come from the current JHK-approved observation log and
  its model-facing audit.
- The nine paired rows come from the JHK-approved paired 2x2 seed table.
- Directly verified counts and current selection/denominator notes were used.
- `NI` means that the current audit did not contain enough information to infer
  an answer. It was used instead of guessing.
- No PPV, renewal, or ORI model was run.

The paired first pass distinguishes direct primary-table checks already
recorded for Gasem 2002 and Wain 2008 from rows inherited through the
Mogasale-derived seed. The latter are not rejected; their primary-method detail
is marked insufficient where it has not been re-extracted in the current audit.

## How JHK should review the workbook

For each row:

1. Review the signaling answers and domain rationales against the primary
   paper or the source audit.
2. Edit any answer or judgment that is not supported.
3. Enter `agree_with_first_pass`, `revise`, `needs_full_text`, or
   `not_applicable` in `jhk_decision`.
4. Enter an ISO date (`YYYY-MM-DD`) and explain any revision in `jhk_notes`.
5. Leave the independent-review fields for a second human reviewer.

After the independent reviewer completes the same exercise, resolve
disagreements and set `consensus_status` to `agreed` or
`resolved_after_discussion`. Until then, these files are not the final
systematic-review risk-of-bias assessment.

## Interpretation guardrails

- A high risk-of-bias judgment does not automatically exclude an estimate.
- Applicability and risk of bias are separate judgments.
- Blood-culture positivity is not identical to PPV of the clinical definition.
- Selected verification can bias positivity even when the positive/tested pair
  itself is arithmetically correct.
- Expanded-setting studies remain outside the target synthesis even if their
  internal risk of bias is low.
- The paired search should be described as prior-review seeding plus JHK's
  targeted supplementary attempts, not as a newly completed exhaustive search
  across PubMed, Embase, Scopus, and Web of Science.

