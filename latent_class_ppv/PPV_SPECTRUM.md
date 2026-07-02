# PPV across the spectrum of typhoid case definitions

Answers: how does the PPV of a typhoid case definition change as the definition goes from broad
(fever) to strict (fever + serology), and across settings from endemic to outbreak? Because a
cross-outbreak comparison confounds the definition with local prevalence, we instead compute the
Bayes surface with **external** definition accuracy and vary prevalence explicitly.

```
PPV(definition, prevalence) = prev*Se_clin / [prev*Se_clin + (1-prev)*(1-Sp_clin)]
  prev    = typhoid prevalence among febrile presenters (the population the definition screens)
  Se_clin, Sp_clin = sensitivity/specificity of the case definition vs culture (literature)
```

Build: `Rscript latent_class_ppv/ppv_spectrum.R` -> `figures/fig_ppv_spectrum.png`,
`tables/ppv_spectrum_table.csv`. Inputs in `data/clinical_definition_accuracy.csv`.

## The definition spectrum (Se/Sp from the literature)

| definition (broad -> strict) | Se_clin | Sp_clin | tier | source |
|---|---|---|---|---|
| Bare fever >=3d | 0.80 | 0.20 | syndromic (assumed) | Storey 2015 (illustrative) |
| WHO suspected (fever + symptom) | 0.83 | **0.36** | syndromic | Thriemer 2012 (Pemba) / Storey 2015 |
| Prediction rule (Hosoglu; +Widal+WBC) | 0.84 | 0.82 | syndromic + lab | Hosoglu 2006 |
| Fever + single serology (TUBEX) | 0.75 | 0.88 | + serology | Storey 2015 (TUBEX pooled) |
| WHO probable (fever + dual Widal) | 0.41 | **0.997** | + serology | Thriemer 2012 (Pemba) |
| Fever + triple serology | 0.39 | 1.00 | + serology | Dong 2010 (Hechi) |

**The structural fact:** clinical criteria alone cap specificity at ~**0.20-0.36**; only a
**serologic tier** pushes Sp to ~0.85-1.0. (Storey serology meta-analysis: Widal 0.83, TUBEX 0.88,
Typhidot 0.80 specificity.)

## PPV by definition x prevalence

| definition | Sp | endemic 0.03 | moderate 0.10 | outbreak 0.30 | high 0.50 |
|---|---|---|---|---|---|
| Bare fever | 0.20 | 0.03 | 0.10 | 0.30 | 0.50 |
| WHO suspected | 0.36 | 0.04 | 0.13 | 0.36 | 0.56 |
| Hosoglu rule | 0.82 | 0.13 | 0.34 | 0.67 | 0.82 |
| Fever + TUBEX | 0.88 | 0.16 | 0.41 | 0.73 | 0.86 |
| WHO probable | 0.997 | 0.81 | 0.94 | 0.98 | 0.99 |
| Fever + triple serology | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

## Key findings

1. **A syndromic-only definition cannot give high PPV — its low specificity is the ceiling.**
   With Sp ~0.2-0.36, PPV ~= the prevalence: **3-4% in endemic surveillance, and even in a large
   outbreak (prev 0.30) only ~0.36.** No amount of clinical wording fixes this; the false-positive
   febrile pool dominates.
2. **The PPV lift requires a serologic tier.** WHO probable (fever + dual Widal, Sp 0.997) gives
   PPV **0.81 at 3% prevalence** vs 0.04 for WHO suspected at the same prevalence - a ~20-fold lift
   from adding serology, not from clinical criteria. Prediction rules / single serology (Sp
   0.82-0.88) are intermediate (~2-3x lift).
3. **Prevalence dominates for syndromic definitions.** The same WHO-suspected definition yields PPV
   0.04 (endemic) -> 0.56 (high outbreak) - which is exactly why outbreak PPV >> endemic PPV, and
   why our outbreak `pi_o` (0.23-0.65) sit far above endemic surveillance positivity (STRATAA
   3-8%).
4. **Internal validation.** Two anchors land on the curves: the Pemba WHO-suspected (2.7%) and
   WHO-probable (73.1%) PPVs (Thriemer 2012) sit on their respective lines at ~2% prevalence; and
   our independently-estimated outbreak PPVs (Aye 0.26, Kabwama 0.23, Neil 0.65) sit on the
   syndromic (WHO-suspected) curve at implied typhoid-among-febrile prevalences of **0.19-0.59** -
   plausible outbreak values, and consistent with the latent-class fit.

## Practical implications

- **For surveillance/confirmation:** to get a usable-PPV case definition without culture, you must
  add a rapid serologic test (a "probable" tier). Tightening clinical criteria alone (more
  symptoms, malaria-negativity) moves Sp only from ~0.20 to ~0.36 - not enough.
- **For burden/CEA:** the suspected-case -> true-case scaling (the confirmation adjustment / PPV) is
  not a fixed number - it must be read off this surface using the *specific definition* and the
  *setting prevalence*. Using a hospital/endemic value in an outbreak (or vice versa) biases
  absolute cases severalfold (cf. `FINAL_REPORT.md`: hospital PPV ~0.8 vs community ~0.25).
- **Relation to definition sensitivity:** the same Se_clin/Sp_clin inputs answer the companion
  question - a syndromic definition is sensitive (Se ~0.83) but non-specific; the serologic tiers
  trade sensitivity (Se falls to ~0.4) for specificity (Sp -> ~1). High PPV comes at the cost of
  missing cases.

## Caveats

- **Serologic-tier PPVs assume the serology's specificity holds**, but Widal/TUBEX specificity is
  setting-dependent and **degrades in endemic areas** (background antibody) - so the "probable"
  curves are optimistic where typhoid is highly endemic.
- Se/Sp are from **heterogeneous settings and reference standards** (mostly vs imperfect blood
  culture, which modestly understates Se); the "bare fever" point is an illustrative assumption.
- **Prediction rules (Hosoglu, Ross-Abraham) embed a Widal term**, so they are syndromic+serology in
  operation, not purely clinical.
- **"Prevalence among febrile presenters"** is the operative x-axis; it is itself setting-specific
  and not always directly measured (here inferred for our outbreaks from the syndromic curve).

## Sources
Thriemer et al. 2012 PLoS One 7:e51823 (Pemba WHO defs); Storey et al. 2015 PLoS One 10:e0142364
(diagnostic-accuracy meta-analysis; serology Se/Sp); Hosoglu et al. 2006 Trans R Soc Trop Med Hyg
100:1068; Ross & Abraham 1987 Trans R Soc Trop Med Hyg 81:374; Dong et al. 2010 PLoS One
(Hechi combined definitions); Khaliq et al. 2021 (case-definition spectrum, qualitative).
