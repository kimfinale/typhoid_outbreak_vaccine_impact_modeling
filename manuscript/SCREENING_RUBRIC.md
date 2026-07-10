# Title/abstract screening rubric (pre-specified, before screening)

Operationalizes PROTOCOL §2–3 into title/abstract-stage decisions. **Conservative for recall:** the TA
stage only removes *clear* non-matches; anything plausibly eligible advances to full-text. CC's screen is
**advisory** (per the division-of-labor); the human adjudicates includes + maybes and spot-checks excludes.
Decision ∈ {include, maybe, exclude}; every record also gets a **layer tag** and a one-line reason.

## INCLUDE / MAYBE (advance to full-text) if the record plausibly involves typhoid/enteric fever AND any of:
- **Layer 1 (culture accuracy):** paired **blood + bone-marrow** culture, or **blood-culture sensitivity /
  yield / detection / proportion detected**, or blood-culture vs another reference.
- **Layer 2 (clinical-definition accuracy):** a **clinical/suspected case definition** evaluated vs culture;
  **clinical prediction rule**; **sensitivity/specificity of clinical features**; **predictive value** of a
  definition.
- **Layer 3 (surveillance/outbreak positivity):** an **outbreak** or **surveillance/population-based** study
  reporting **suspected AND confirmed (culture) counts**, **positivity / proportion confirmed**, or an
  **attack rate** with a case definition.

Tag the likely layer(s). If it fits none cleanly but mentions typhoid + culture/definition/positivity →
**maybe**.

## EXCLUDE at TA stage ONLY if clearly ineligible (state which rule):
1. **Not typhoid/enteric fever / not *S.* Typhi or Paratyphi** (other pathogen; keep only if typhoid is a
   reported comparator with extractable counts).
2. **Not human** (animal, environmental, in-vitro only).
3. **No diagnostic / positivity / case-definition / culture content** — e.g., pure **genomics/AMR-mechanism**,
   **vaccine efficacy/immunogenicity without positivity or a suspected denominator**, **cost/economic-only**,
   **risk-factor/WASH-only**, **molecular/pathogenesis**, **treatment trial without diagnostic counts**.
4. **Case report / small case series with no denominator** (< 10 suspected; no k/n).
5. **Review / editorial / commentary / news / protocol / erratum** → **exclude as data** but flag
   `citation_chase = yes` if it is a systematic review/meta-analysis or an outbreak catalogue (Khaliq-type)
   whose reference list should be harvested.
6. **Special-population-only** cohort (HIV/cancer/returning travellers/chronic carriers only).
7. **Wrong language** (not English/French/Spanish/Chinese) — mark `language_exclude`.

## Ambiguity rule
When genuinely uncertain → **maybe** (advance). Do **not** exclude on suspicion. Selective-testing /
sampling-representativeness and exact definition class are assessed at **full-text/extraction**, not here —
do not exclude a positivity study at TA for a subset-testing concern.

## Output (per record)
`pmid, title, decision {include|maybe|exclude}, layer {L1|L2|L3|multiple|none}, reason, citation_chase, language_exclude`
→ `manuscript/data/screening_ta_pubmed.csv`. Includes + maybes → full-text retrieval (human).
