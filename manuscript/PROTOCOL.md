# Protocol — systematic review + latent-class synthesis of the PPV of typhoid clinical case definitions (with blood/bone-marrow culture sensitivity as machinery)

Registerable protocol (target: **PROSPERO**; note PROSPERO accepts DTA reviews). Reporting: **PRISMA-DTA**
(McInnes 2018). Risk of bias: **QUADAS-2** (Whiting 2011). This protocol is **pre-committed**; deviations
will be reported. It follows the Wiens 2023 (PLoS Med) template adapted to typhoid.

## 1. Objective / estimand structure (PICO-style; FIXED)

**Primary estimand.** The positive predictive value of the clinical/suspected typhoid case definition —
`π = P(true typhoid | meets suspected-case definition)` — **by setting** (hospital-referred vs
community/surveillance) and **as a function of prevalence** (typhoid among suspected/febrile presenters).
"True typhoid" = latent infection status inferred from imperfect culture(s), not a naive culture-positive.

**Supporting estimands (machinery, not headline).** Blood-culture sensitivity `Se_BC` as a function of
blood volume; bone-marrow-culture sensitivity `Se_BM`; culture specificity `Sp ≈ 1` (**pinned**, not
estimated — a positive *S.* Typhi/Paratyphi culture is definitionally a true case). These make `π`
identifiable from positivity (`positivity = π · Se_BC`).

- **Population:** patients meeting a stated clinical/suspected typhoid definition (febrile, care-seeking);
  and, for the culture layer, patients cultured by both blood and bone marrow.
- **Index:** the clinical/suspected case definition (primary); blood culture, bone-marrow culture (machinery).
- **Reference:** culture (imperfect) → latent true-infection status via the latent-class model.
- **Out of scope as objects of estimation:** antibody RDTs (Widal/TUBEX/Typhidot) — cite Arora 2019;
  used only as illustrative serologic-tier inputs on the PPV surface.

## 2. Evidence streams & inclusion criteria

Three streams; a study may contribute to more than one.

1. **Culture-accuracy (paired BC/BM).** Same patients cultured by **both** blood and bone marrow, with a
   recoverable 2×2 (BC±/BM±) or per-test detection counts, ideally + blood volume. *No date restriction*
   (paired-culture methodology is stable; Mogasale/Antillón include pre-2000 studies).
2. **Clinical-definition accuracy.** Culture applied to a **defined febrile denominator**, allowing a
   **definition-vs-culture 2×2** or reported Se/Sp of the definition; OR **positivity** (k culture-confirmed
   / n meeting the definition) with a **stated case definition** + **setting**. *Samples ≥ 1 Jan 2000.*
3. **Surveillance/outbreak positivity (the spine).** Outbreak reports and surveillance cohorts reporting
   **suspected AND confirmed (culture) counts**, with a stated case definition and setting label. *≥ 2000.*

## 3. Exclusion criteria (Wiens-style, adapted to typhoid)

1. **No stated case definition**, or a definition **not specific for suspected typhoid** (accept
   fever [± duration] with typhoid-compatible features, incl. the WHO suspected definition; **exclude**
   "any acute febrile illness / any fever / any diarrhoea" with no typhoid-specific criteria).
2. **Epidemiological-link case selection** (cases chosen by link to other cases or an environmental source)
   — for the positivity/PPV estimand only (biases PPV upward); such studies may still inform transmission
   context but are excluded from π/φ.
3. **Special populations only** (HIV/cancer, chronic carriers, returning travellers, pregnancy cohorts).
4. **< 10 suspected cases cultured** (unstable positivity).
5. **Language** other than English, French, Spanish, Chinese (state any not screenable).
6. **Reviews/editorials/modelling-only** without extractable primary counts (harvested for citation
   chasing, not included as data).

## 4. Setting-labeling rule (pre-committed)

Every primary-evidence study is labeled: **hospital-referred** (inpatient / tertiary referral) /
**community-or-surveillance** (population-based, outpatient, community/outbreak field surveillance) /
**mixed-facility** (facility-based spanning in- and out-patients, or explicitly both).
- If a study reports **stratified** positivity by setting/age/time, **each stratum is a separate
  observation** (Wiens approach). Priority order for multiply-stratified studies: setting → age →
  antibiotic use → year → geography.
- **mixed-facility** is a **third stratum**, never forced into hospital or community. If genuinely
  unstratifiable, it is **kept in the overall summary + prevalence surface but excluded from the
  hospital-vs-community contrast**, with the reason recorded. No ad hoc adjudication later.

### 4.1 Sampling-representativeness of the cultured subset (pre-committed) — the key positivity modifier

Positivity is computed among the **cultured subset**, and `positivity / Se_BC` recovers the PPV of the
*full* suspected definition **only if culturing was representative among suspected cases**. `Se_BC`
corrects for the test *missing* true cases; it does **not** correct for *who was selected to be tested*.
So every positivity observation is additionally labeled by how its cultured subset was chosen
(extraction field `sampling_representativeness`, with the verbatim description recorded):

- **reference tier — ALL suspected cultured:** positivity is the direct definition-level PPV × Se_BC.
- **high — random / systematic (e.g., "every k-th case") / consecutive** sampling of suspected cases:
  treated as representative.
- **low — convenience / clinician-selected / sicker-first / failing-antimalarials / unreported:**
  positivity is expected **biased upward** (tested subset enriched for true typhoid) and **cannot be
  de-biased analytically**.

Rules (pre-committed): (i) `sampling_representativeness` is an extracted variable and a **covariate/
stratifier** in the positivity model; (ii) the **primary analysis restricts to reference + high tiers**,
with low-tier studies added only in a sensitivity analysis (and reported separately) to show the
direction/size of the selection bias; (iii) studies culturing **all** suspected are the reference-quality
tier for the outbreak/community PPV. This mirrors Wiens 2023's sampling-quality (high/low) handling.

### 4.2 Case-definition classification (pre-committed) — the PPV-spectrum axis

Suspected-case definitions are heterogeneous and non-standardized (Khaliq 2021), and their **breadth
drives specificity and hence PPV** (broad C1/C2 vs serology-tier C5 differ by ~0.6 in specificity — see
`latent_class_ppv/PPV_SPECTRUM.md`). Each study's clinical/suspected definition is classified into a
pre-committed class **and** coded on discriminating features (verbatim wording recorded).

**Primary class (syndromic-breadth spectrum, broad → strict):**
- **C1 Bare fever** — fever alone (± duration), no required additional feature.
- **C2 Fever + ≥1 typhoid symptom** — WHO-suspected style (fever ± duration + ≥1 GI/abdominal/systemic).
- **C3 Fever + multiple required symptoms** — stricter syndromic (e.g., fever ≥3 d + headache +
  abdominal pain + ≥2 GI).
- **C4 Syndromic + malaria-exclusion** — C2/C3 plus negative malaria test and/or antimalarial
  non-response (endemic-Africa tightening).
- **C5 Probable (syndromic + serology)** — fever + a positive rapid serology (Widal/TUBEX/Typhidot).

**Non-syndromic types (classified but handled SEPARATELY — not on the suspected-PPV spectrum):**
- **D Perforation / surgical series** — cases defined by intestinal perforation/peritonitis (severity-
  filtered, e.g., Muyembe); excluded from the suspected-definition PPV, reported separately.
- **X Confirmed-only** — "cases" are culture-confirmed by construction (no suspected denominator for a
  PPV, e.g., Qamar/Yousafzai XDR counts); contributes to the culture layer, not the suspected-PPV estimand.

**Coded discriminating features (reproducible stringency proxy):** fever-duration threshold
(none / ≥3 d / ≥7 d); number of required additional symptoms (0 / 1 / ≥2); malaria-exclusion (Y/N);
serologic tier (Y/N); entry point (self-referral vs screened). These map each definition onto the
PPV-surface specificity axis.

**Use (pre-committed):** `definition_class` is a **stratifier/covariate** in the positivity/PPV model and
the **primary axis of the PPV × prevalence surface**; positivity is **not pooled across classes** without
accounting for class. Ambiguous/unclassifiable definitions are flagged for human adjudication
(`data/EXTRACTION_FLAGS.md`), not silently binned.

## 5. "Overall PPV" definition (pre-commit ONE — Wiens-style)

**Overall = the pooled transportable ingredients (definition Se/Sp, positivity) with PPV always shown as a
curve over prevalence; hospital and community reported within-stratum.** Any single marginal PPV figure is
explicitly labeled a **prevalence-weighted summary of the study mix** (with I² and the full stratified
estimates alongside) — **never** presented as a single context-free "the PPV of typhoid clinical
definitions." We report **overall summary + stratified + prevalence surface together** (Wiens' "headline
overall, immediately qualified"). We do **not** compute a burden-correction multiplier (see framing memo).

## 6. Search strategy (per layer, per database; run by the HUMAN)

No DTA methodological/search filter (poor sensitivity per Cochrane DTA guidance) — accept a large pull.
Databases: **PubMed, Embase, Scopus, Web of Science** (distinct syntaxes below). Date limits per stream
(§2). Deduplicate across databases (e.g., in Covidence/Rayyan).

### Layer 1 — culture accuracy (paired blood + bone marrow)
- **PubMed:** `(typhoid[tiab] OR "enteric fever"[tiab] OR "Salmonella Typhi"[tiab] OR "Salmonella Paratyphi"[tiab] OR "Salmonella typhi"[MeSH]) AND ("bone marrow"[tiab] AND ("blood culture"[tiab] OR "blood cultures"[tiab])) AND (sensitiv*[tiab] OR yield[tiab] OR detect*[tiab] OR positiv*[tiab] OR recover*[tiab])`
- **Embase:** `('typhoid fever'/exp OR 'Salmonella typhi'/exp OR typhoid:ti,ab OR 'enteric fever':ti,ab) AND ('bone marrow culture':ti,ab OR ('bone marrow':ti,ab AND 'blood culture':ti,ab)) AND (sensitivity:ti,ab OR yield:ti,ab OR detect*:ti,ab OR positiv*:ti,ab)`
- **Scopus:** `TITLE-ABS-KEY((typhoid OR "enteric fever" OR "Salmonella Typhi" OR "Salmonella Paratyphi") AND "bone marrow" AND ("blood culture" OR "blood cultures") AND (sensitiv* OR yield OR detect* OR positiv* OR recover*))`
- **Web of Science:** `TS=((typhoid OR "enteric fever" OR "Salmonella Typhi") AND "bone marrow" AND "blood culture*" AND (sensitiv* OR yield OR detect* OR positiv*))`

### Layer 2 — clinical case-definition accuracy
- **PubMed:** `(typhoid[tiab] OR "enteric fever"[tiab] OR "Salmonella Typhi"[tiab]) AND ("case definition"[tiab] OR "clinical diagnosis"[tiab] OR "clinical features"[tiab] OR "clinical prediction"[tiab] OR "prediction rule"[tiab] OR "suspected typhoid"[tiab] OR "predictive value*"[tiab] OR (sensitiv*[tiab] AND specific*[tiab])) AND ("blood culture"[tiab] OR "culture-confirmed"[tiab] OR "culture confirmed"[tiab] OR "bone marrow"[tiab])`
- **Embase:** `(typhoid:ti,ab OR 'enteric fever':ti,ab) AND ('case definition':ti,ab OR 'clinical diagnosis':ti,ab OR 'clinical feature':ti,ab OR 'prediction rule':ti,ab OR 'suspected typhoid':ti,ab OR 'predictive value':ti,ab OR (sensitivity:ti,ab AND specificity:ti,ab)) AND ('blood culture':ti,ab OR 'culture confirmed':ti,ab OR 'bone marrow':ti,ab)`
- **Scopus:** `TITLE-ABS-KEY((typhoid OR "enteric fever") AND ("case definition" OR "clinical diagnosis" OR "clinical features" OR "prediction rule" OR "suspected typhoid" OR "predictive value" OR (sensitiv* AND specific*)) AND ("blood culture" OR "culture confirmed" OR "bone marrow"))`
- **Web of Science:** `TS=((typhoid OR "enteric fever") AND ("case definition" OR "clinical diagnosis" OR "clinical feature*" OR "prediction rule" OR "suspected typhoid" OR "predictive value*") AND ("blood culture*" OR "culture confirmed" OR "bone marrow"))`

### Layer 3 — outbreak / surveillance positivity (PRIORITY — the spine)
- **PubMed:** `(typhoid[tiab] OR "enteric fever"[tiab] OR "Salmonella Typhi"[tiab]) AND (outbreak*[tiab] OR surveillance[tiab] OR "population-based"[tiab] OR incidence[tiab] OR "case-based"[tiab]) AND (positiv*[tiab] OR "proportion confirmed"[tiab] OR "confirmed cases"[tiab] OR "suspected cases"[tiab] OR "attack rate"[tiab] OR "blood culture"[tiab] OR "culture-confirmed"[tiab])`
- **Embase:** `(typhoid:ti,ab OR 'enteric fever':ti,ab) AND (outbreak*:ti,ab OR surveillance:ti,ab OR 'population based':ti,ab OR incidence:ti,ab) AND (positiv*:ti,ab OR 'proportion confirmed':ti,ab OR 'confirmed case*':ti,ab OR 'suspected case*':ti,ab OR 'attack rate':ti,ab OR 'blood culture':ti,ab)`
- **Scopus:** `TITLE-ABS-KEY((typhoid OR "enteric fever" OR "Salmonella Typhi") AND (outbreak* OR surveillance OR "population-based" OR incidence) AND (positiv* OR "proportion confirmed" OR "confirmed cases" OR "suspected cases" OR "attack rate" OR "blood culture" OR "culture-confirmed"))`
- **Web of Science:** `TS=((typhoid OR "enteric fever" OR "Salmonella Typhi") AND (outbreak* OR surveillance OR "population-based" OR incidence) AND (positiv* OR "confirmed case*" OR "suspected case*" OR "attack rate" OR "blood culture*"))`

### Citation chasing (legitimate SR method)
Forward/backward citation-chase: **Thriemer 2012** (Pemba WHO-definition Se/Sp), **Aiemjoy 2020**
(model-based multi-country), **Mogasale 2016** (its 10 paired BC/BM studies), **Antillón 2018** (its 25
volume studies), **Khaliq 2021** (its 68 outbreak case-definition studies), **SEAP/Carey 2020**, and the
**Wiens 2023** method chain (Wang–Lin–Nelson 2020).

## 7. Overlap / double-counting rule (pre-committed)
Where samples overlap (e.g., **Wain 2008 / Wain 2001 / Gasem** Dong-Thap–Vietnam overlaps; multi-copy
Zotero entries), **prefer the larger / primary / more-representative sample; exclude the overlapping
one(s); record the decision** in `manuscript/data/EXTRACTION_FLAGS.md`. Successive distinct outbreaks in
one location (e.g., Davis 2016–17 vs N'Cho 2017–18, Harare) are **separate**, not overlaps.

## 8. Risk of bias — QUADAS-2 plan
Four domains with typhoid-specific signalling questions:
- **Patient selection:** consecutive/random vs convenience; case-control avoided; epi-link exclusion applied.
- **Index test (case definition / culture):** definition pre-specified & applied uniformly; blood volume
  reported; antibiotic pre-treatment.
- **Reference standard (culture — imperfect):** blood vs bone marrow; volume; timing vs illness week.
- **Flow & timing:** all suspected cultured vs a selected subset (selection bias — a key concern);
  same reference for all.
Concerns map to **sensitivity analyses** (exclude high-RoB; exclude convenience-sampled; exclude
selected-testing subsets).

## 9. Pre-committed clinical-layer fallback (registered decision)
**Set K = 4** (minimum for a stable bivariate DTA meta-analysis). *If fewer than 4 studies with a clean
definition-vs-culture 2×2 on a defined febrile denominator survive screening, do NOT fit a clinical
latent-class or bivariate/HSROC layer.* Instead report the clean 2×2(s) individually and characterize PPV
with the **prevalence-varying Bayes surface** from external Se/Sp inputs (Thriemer, Storey, etc.), and
state this downgrade as **pre-specified** (in `CLINICAL_LAYER_DECISION.md`).

## 10. Analysis plan
- **Reuse** `latent_class_ppv/stan/model_final.stan` structure: paired-culture latent class (Hui–Walter) →
  `Se_BC(volume)`, `Se_BM`, `Sp` pinned; outbreak positivity → hospital `φ_s` / community `π_o`;
  transportable-Se / local-PPV split.
- **Recovery gate FIRST** (reuse `run_phase1_recovery.R` pattern): coverage/bias/convergence over
  replicates; respect PASS/FAIL.
- **Positivity synthesis (spine):** per-setting `π`/`φ` with a study random effect, adjusting for
  Se_BC and including **`definition_class`** (§4.2) and **`sampling_representativeness`** (§4.1) as
  covariates; positivity is **not pooled across definition classes** without accounting for class;
  **primary fit restricted to reference + high representativeness tiers.** Optional
  **conditional-dependence** sensitivity (Wang–Lin–Nelson 2020; Dendukuri & Joseph 2001).
- **Clinical Se/Sp:** bivariate (Reitsma 2005) / HSROC (Rutter & Gatsonis 2001) **iff K cleared**; else
  the fallback surface.
- **PPV × prevalence × definition surface** (reuse `ppv_spectrum.R`) with empirical anchors overlaid.
- **Sensitivity:** add low-representativeness studies (quantify the upward selection bias vs the primary
  reference+high fit); conditional dependence; Sp-prior relaxation; blood-volume sweep (3–10 mL); the
  endemic-Widal-specificity caveat for any serologic-tier PPV; overlap-swap.
- Report divergences, max R-hat, min ESS; regenerate all figures from code.

## 11. Deviations
Any change to this protocol after registration is logged with date + rationale in
`manuscript/PROTOCOL_DEVIATIONS.md`.
