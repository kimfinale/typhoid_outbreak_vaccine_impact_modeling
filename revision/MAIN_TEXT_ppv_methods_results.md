# Main-text Methods and Results — PPV adjustment for the ORI impact analysis

*Draft for the main manuscript. Numbers are from the current fit (nine paired
studies, Vallenas excluded; eight community outbreaks). Detail and equations are
in the Supplement (S1–S5). Figures 1–3 are in `revision/main_text_figures/`.*

---

## Methods

### Positive predictive value of the suspected case definition

Outbreak time series used in the outbreak-response immunization (ORI) analysis
count **suspected** typhoid cases, defined for sensitive detection and therefore
including febrile illnesses of other cause. We defined the positive predictive
value of the suspected case definition in setting *o* as
$\pi_o = P(\text{true typhoid} \mid \text{suspected})$ — the prevalence of true
*Salmonella* Typhi infection among suspected cases — and propagated its
uncertainty into estimates of true cases averted, morbidity, treatment cost, and
cases averted per dose. This is a setting-specific observation model for the ORI
analysis, not a systematic review of case definitions.

### Latent-class estimation of blood-culture sensitivity

Because blood culture is highly specific but incompletely sensitive, culture
positivity understates true typhoid, while a syndromic suspected count overstates
it. We reconciled the two with a Bayesian latent-class model. Nine historical
studies that cultured the same patients by both blood and bone marrow (a
Hui–Walter design) identified blood-culture sensitivity ($Se_{BC}$), bone-marrow
sensitivity ($Se_{BM}$), and culture specificity without treating either culture
as a perfect gold standard. Blood-culture sensitivity was modelled as a function
of cultured blood **volume**, $\operatorname{logit}(Se_{BC,s}) = \alpha_0 +
\alpha_1 \log V_s + \tau u_s$, with an Antillón (2018)-informed prior on the
volume slope; an earlier severity term was dropped because it was not identifiable
from the available hospital-based paired evidence (Supplement S4).

### Community PPV and verification-bias selection

For each community outbreak we observed the number of suspected cases *N*, the
number blood-cultured *n*, and the number culture-positive *k*. PPV among **all**
suspected cases was estimated hierarchically,
$\operatorname{logit}(\pi_o) = \mu_\pi + \sigma_\pi z_o$, jointly with blood-culture
sensitivity from the paired layer. Because outbreaks preferentially culture
sicker or more typhoid-like patients, a selection sub-model allowed the
probability of culture to differ between true-typhoid and non-typhoid suspected
cases, with the enrichment scaled by a prespecified per-study audit grade (low,
moderate, major). The three observed counts (culture-positive, culture-negative,
untested) followed a multinomial likelihood; because population PPV and selection
are not separately identified by these counts alone, the selection layer is
governed by structured priors and is treated as a sensitivity dimension rather
than a data-driven estimate.

### Evidence selection

Community observations were drawn from 107 candidate outbreak records assembled
from two systematic reviews (Koh 2025; Appiah 2020) and audited at full text
(Supplement S2.3). Admissibility required a syndromic suspected-case denominator
defined independently of culture, a **blood-culture-specific** tested count (not
pooled across specimens), a positive count from the same subset, and outbreak
independence. Eight outbreaks met all criteria; the remainder were excluded
(most commonly for confirmed-only case series or pooled-specimen denominators),
retained only for sensitivity analysis, or held pending source retrieval. Blood
volume was unreported in all eight community reports and was set to 5 mL, varied
from 3 to 10 mL as the principal structural sensitivity.

### Propagation and sampling

Per-outbreak, per-draw PPV values allocated suspected incidence additively into
true typhoid and other febrile illness. True typhoid was then allocated, within a
shared renewal recursion, into background/common-source and propagated components.
Vaccination acted directly on both typhoid components, while propagated cases also
generated transmission feedback; other febrile illness was held fixed. The former
post hoc θ-weighted static-plus-renewal calculation was retained only as a paired
historical comparator. Deaths and years of life lost, already defined on true
infections, were left unscaled. The PPV model was fitted
in Stan (four chains, 1,500 warm-up and 1,500 retained iterations each). We report
posterior medians and 90% credible intervals (CrIs), and we assessed robustness
with the blood-volume sweep, six selection-prior/audit-grade scenarios, and a
leave-one-outbreak-out check for the sole specimen-ambiguous anchor.

## Results

### Blood culture detects about two-thirds of true typhoid (Figure 1)

Estimated blood-culture sensitivity was **0.65 (90% CrI 0.56–0.73) at 5 mL**,
rising from 0.56 (0.46–0.66) at 2 mL to 0.72 (0.62–0.80) at 10 mL, and was
consistently above the conventional 0.5 correction factor. Bone-marrow sensitivity
was 0.92 (0.90–0.94) and culture specificity was near unity (0.996). Because a
single blood culture misses roughly a third of true cases, raw culture positivity
materially understates the proportion of suspected patients with true typhoid.

### Community PPV is low and highly heterogeneous (Figure 2)

Across the eight outbreaks, the fitted PPV of the suspected case definition ranged
from **0.10 to 0.86**, with a population median of **0.31 (0.16–0.52)** and a
large between-outbreak logit-scale standard deviation of 1.35 (0.85–2.03) that was
clearly updated from its prior. Culture positivity was not interchangeable with
PPV: imperfect sensitivity inflates the implied PPV, while preferential culture of
likely-typhoid patients deflates it, and the two adjustments partly offset in the
most verification-biased outbreaks. Most outbreaks — including the large Kabwama
(0.15) and Neil (0.38) series — had a minority of suspected cases representing
true typhoid.

### Additive propagation reveals a larger true than surveillance-visible ORI effect (Figure 3)

In the primary two-level additive model, ORI averted a pooled median **1,846 true
typhoid cases** (95% simulation interval 644–3,742), reducing true typhoid by
**25.6%** but suspected-case surveillance by only **7.0%** because other febrile
illness was unchanged. The historical θ-weighted calculation gave **25.9%**; the
paired additive-minus-weighted difference was **−0.5 percentage points**
(−2.4 to 0.7). A separate PPV-selection sensitivity conditional on the
pure-renewal workflow bracketed cases averted at **1,511–2,277**; it is not an
uncertainty interval for the primary source-plus-propagated model. The main PPV fit sampled well
(0 divergent transitions, maximum $\widehat R = 1.004$, minimum bulk effective
sample size ≈ 1,700).

---

## Figure captions

**Figure 1. Blood-culture sensitivity increases with cultured blood volume and
exceeds the conventional 0.5 assumption.** Posterior median (line) and 90% CrI
(band) of $Se_{BC}$ from the latent-class model; grey points are raw
blood-culture detection proportions in the nine historical paired studies. The
orange point marks the 5 mL estimate used for community outbreaks; the dashed line
is the 0.5 correction factor commonly assumed in burden studies.

**Figure 2. True-typhoid PPV of the suspected case definition is low and varies
widely across outbreaks.** Filled points and bars are the fitted PPV among all
suspected cases (posterior median, 90% CrI), coloured by prespecified
verification-adjustment grade; hollow points are raw culture positivity *k/n*. The
shaded band and dotted line show the population PPV (median and 90% CrI).

**Figure 3. Two-level additive ORI model versus the historical θ-weighted
calculation.** (A) Paired true-typhoid cases averted with 95% simulation intervals.
(B) Percent reduction in true typhoid versus the reduction visible in suspected-case
surveillance for each structural formulation.
