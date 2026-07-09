# Methods scoping memo

Purpose: justify (and cite) the synthesis method per estimand-family, confirm the review frameworks,
position against the two nearest precedents (Arora, Wiens), and record the existing-review landscape.
Two Phase-1 reads inform this: a full-text read of Wiens 2023 and a web landscape/gap scan. Honesty
note: one premise in `MANUSCRIPT_PROMPT_v2.md` needed correcting — see §5.

## 1. Review frameworks (confirmed, with a hybrid caveat)

- **Reporting: PRISMA-DTA** (McInnes et al. 2018, JAMA 319:388) — the DTA extension of PRISMA. Fixed by
  the prompt and appropriate for the two **accuracy** layers (culture Se/Sp; clinical-definition Se/Sp).
- **Risk of bias: QUADAS-2** (Whiting et al. 2011, Ann Intern Med 155:529) — for the accuracy layers.
- **Honest caveat (hybrid genre).** The **spine** of this paper — setting-stratified *positivity* /
  PPV of the suspected-case definition — is not a classical Se/Sp estimand, it is a
  proportion/prevalence synthesis. The direct precedent (Wiens 2023) reported under **generic PRISMA
  2020, not PRISMA-DTA, and did not use QUADAS-2.** We nonetheless adopt PRISMA-DTA (superset of PRISMA)
  + QUADAS-2 because our culture and clinical-definition layers are genuine DTA; for the positivity
  layer we report the PRISMA flow and apply QUADAS-2 domains where they map (patient selection, flow/
  timing), noting where a domain is not applicable. This hybrid is stated in Methods, not hidden.

## 2. Existing-review landscape & the gap (Phase-1 scan)

- **Category (a) — diagnostic-TEST accuracy: densely reviewed** (no novelty available here): Storey 2015
  (PLoS ONE, composite reference), Wijedoru 2017 (Cochrane DTA, RDTs), Arora 2019 (network LCM),
  Mogasale 2016 (BC detection proportion 61%), Antillón 2018 (BC Se vs blood volume), + immunodiagnostic
  reviews.
- **Category (b) — PPV of the clinical/suspected case definition, or setting-stratified
  suspected→confirmed positivity: EMPTY.** No published systematic review and no indexed PROSPERO
  protocol targets P(true typhoid | suspected) or the hospital-vs-community positivity relationship.
  - **Khaliq 2021** (J Glob Health 11:04031) is the nearest: it catalogued outbreak case definitions
    (only 15/68 studies even stated one), found **none reported PPV/positivity**, and explicitly did not
    synthesize it — i.e., it *documents the void* and recommends standardization.
  - **SEAP 2020** (Clin Infect Dis 71(S3):S257) is the nearest *primary* PPV: clinical case definition
    PPV ≈ **10.5% (7.2–15.0)** across Bangladesh/Nepal/Pakistan — a single project, not a review, not
    hospital-vs-community stratified. A useful external anchor/prior, not prior art that closes the gap.
  - The hospital (~30–73%) vs community/surveillance (~2–5%) positivity spread exists only as scattered
    figures in incidence/CFR papers, **never pooled**.
- **Residual uncertainty (FLAG FOR HUMAN):** the live PROSPERO registry is a JS single-page app that
  could not be queried programmatically. Recommend one manual PROSPERO search ("typhoid case
  definition", "typhoid positive predictive value", "enteric fever blood culture positivity") to fully
  close the category-(b) novelty claim before submission.

**Positioning sentence:** *diagnostic-test accuracy in typhoid is exhaustively reviewed, but no
systematic review synthesizes the PPV of the clinical/suspected case definition or the setting-
stratified positivity relationship; the one review to touch outbreak case definitions (Khaliq 2021)
found it un-synthesized and called for exactly this.*

## 3. Estimand-family → synthesis method (with citations)

| Estimand | Target | Reference problem | Synthesis (method + cite) |
|---|---|---|---|
| Blood culture | **Se only** (Sp≈1 **pinned**) | no gold standard; bone marrow is a better-but-imperfect reference | conditional-independence **latent class on paired BC/BM** (Hui & Walter 1980, Biometrics 36:167) + **volume model** (Antillón 2018 prior); optional conditional-dependence (Wang, Lin & Nelson 2020; Dendukuri & Joseph 2001) as sensitivity |
| Bone marrow | **Se only** (Sp≈1 pinned) | it is the anchor | same paired studies |
| Clinical case definition | **Se, Sp → PPV(prevalence, setting)** | needs culture on a defined febrile denominator; definitions are not one exchangeable test | **bivariate** (Reitsma 2005, J Clin Epidemiol 58:982) or **HSROC** (Rutter & Gatsonis 2001, Stat Med 20:2865) **iff ≥K studies**; else the **prevalence-varying Bayes PPV surface** from external Se/Sp (pre-specified fallback, §PROTOCOL) |
| Setting-stratified **positivity** (the spine) | community π vs hospital φ | positivity = π·Se_BC; Se_BC from the paired layer | **positivity synthesis** in the same framework: per-setting π backed out via the anchored Se_BC (as in `model_final.stan`); a Wiens-style hierarchical model with a study random effect + test-accuracy adjustment is the template |
| Antibody RDTs | — | — | **out of scope**; cite Arora 2019 (values as illustrative inputs only) |

Rationale: Se/Sp of imperfect tests without a gold standard is the Hui–Walter/latent-class problem; the
paired BC/BM design is what makes Se_BC identifiable (breaks the single-test π·Se ridge). Standard DTA
pooling (bivariate/HSROC) applies to the clinical-definition Se/Sp *if enough clean definition-vs-culture
2×2s survive screening*; otherwise PPV is characterized as a prevalence surface (the fallback). The
positivity spine follows Wiens' template but adapted to typhoid's invasive, insensitive, volume-dependent
confirmation.

## 4. Positioning vs the two nearest precedents

- **Arora 2019** (PLoS NTD 13:e0007303): Bayesian latent-class **network** meta-analysis of typhoid
  *test* accuracy (RDTs/serology/Widal/PCR vs a modeled bone-marrow reference). It **fixes** culture Se
  as a prior, ignores **blood volume**, evaluates **no clinical case definition**, and produces **no
  PPV/prevalence** quantity. We differ by estimating culture Se as a volume-anchored posterior and by
  targeting the **clinical-definition PPV by setting** — which Arora does not touch.
- **Wiens 2023** (PLoS Med 20:e1004286): the **methodological template** — "proportion of suspected
  [cholera] cases that are true infections" via a hierarchical conditional-dependence LCM (Wang, Lin &
  Nelson 2020) + a positivity meta-GLM, stratified by setting (outbreak vs non-outbreak). We adopt the
  estimand and model family and extend to typhoid's culture-specific machinery (bone-marrow anchor,
  blood-volume Se model). **Structural distinction:** cholera confirmation is a stool RDT; typhoid
  confirmation is **invasive, insensitive, and volume-dependent** (blood/bone-marrow culture), which
  makes the suspected count more load-bearing and the Se-identification problem harder — the reason a
  typhoid-specific synthesis is warranted.

## 5. Framing correction (honesty note — diverges from the prompt's premise)

`MANUSCRIPT_PROMPT_v2.md` §"READ THIS FIRST" states Wiens "reframed **surveillance interpretation, not
burden**." The full-text read shows this is **only half right**: Wiens makes **burden a co-primary
motivation** ("suspected cases may overestimate incidence ~2-fold") and does **not** disclaim it; its
caveats are narrower (it corrects only *over*counting, not missed cases; warns against applying the
global number locally; scope = medically-attended). Per the prompt's own honesty rule ("do not
attribute a burden-disclaimer to it"), we will **not** misattribute a disclaimer to Wiens.

**We keep the typhoid paper's surveillance-not-burden framing, but justify it on typhoid-specific
grounds (which is a cleaner differentiator):** cholera surveillance counts **suspected** cases, so
positivity legitimately bears on cholera burden (hence Wiens' burden framing); the dominant **typhoid**
burden pipeline (SEAP/GBD) instead starts from a **culture-confirmed numerator** and adjusts *upward*
for BC sensitivity, blood volume, care-seeking and BC coverage — it never multiplies suspected cases by
a positivity fraction, so the clinical-definition PPV **is not a term in typhoid burden estimation.**
Our contribution therefore lands on **surveillance interpretation, outbreak case-counting and response
sizing (incl. TCV campaign triggers), and testing strategy** — and we explicitly do **not** claim to
revise confirmation-based typhoid burden. This positioning is stronger than the prompt's original
(mis)attribution and is stated as our own argument, not Wiens'.

## 6. Flags for the human
1. **Manual PROSPERO search** to fully close the category-(b) gap claim (see §2).
2. The **PRISMA-DTA + QUADAS-2 vs generic-PRISMA** choice (§1) is defensible but is a reviewer-facing
   decision — confirm at write-up.
3. `meta`/`metafor`/`mada` install deferred to Phase 5 (only if the clinical bivariate/HSROC layer is
   triggered).
