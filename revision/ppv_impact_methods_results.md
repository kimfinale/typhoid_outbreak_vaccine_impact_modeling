# Positive predictive value of outbreak case definitions and its effect on measured vaccine impact

*Draft Methods and Results sections. Numbers are drawn from the current fitted outputs
(`latent_class_ppv/tables/final_parameters.csv`, `final_pi_community.csv`,
`renewal/tables/tab_ppv_effect.csv`, `revision/tables/outbreak_allocation_brackets_production.csv`,
`revision/tables/clause_model/`). Stage-1 latent-class intervals are 90% credible intervals (CrI);
counts are written positive/tested.*

---

## Methods

### 2.x Blood-culture sensitivity and community positive predictive value

During a typhoid outbreak the quantity that surveillance reports — the count of *suspected* cases
meeting a clinical case definition — is separated from the epidemiologically relevant quantity — the
count of true *Salmonella* Typhi infections — by two distinct filters. The first is the **positive
predictive value (PPV, π)** of the clinical case definition: the fraction of suspected cases that are
truly typhoid. The second is the **sensitivity of blood culture (Se_BC)**: the fraction of true
typhoid cases that a single blood culture detects. A single stream of outbreak data observes only
the *product* of these two quantities. Writing the blood-culture positivity among cultured suspected
cases as *positivity* and assuming near-perfect culture specificity,

$$\text{positivity} \;\approx\; \pi \times \text{Se}_{BC},$$

so that positivity, PPV, and sensitivity lie on a confounded ridge: a single outbreak cannot
separate a low-PPV case definition from an insensitive culture. Dividing a naïve positivity by an
assumed sensitivity (e.g. `positivity / 0.6`) can return a value exceeding one, which is a
falsification signal rather than an estimate of PPV.

We break this identifiability problem with a Bayesian latent-class (Hui–Walter) model that couples
two evidence streams. **Historic paired-specimen studies**, in which the *same* patients were
cultured by both blood and bone marrow, provide the leverage that pins Se_BC: bone marrow is a
near-gold-standard specimen (Se_BM), so the discordant cells of the blood × bone-marrow table
identify blood-culture sensitivity directly. Because sensitivity depends on inoculated blood volume,
media, prior antibiotics and illness timing, we model Se_BC as **volume-dependent and only
conditionally transportable**,

$$\text{logit}\,\text{Se}_{BC,s} \;=\; \alpha_0 + \alpha_1 \log v_s + \tau\,u_s, \qquad u_s\sim\mathcal N(0,1),$$

with $v_s$ the blood volume (mL) and $u_s$ a study random effect. The paired studies contribute a
4-cell multinomial per study built from Se_BC, Se_BM, and specificities (Sp_BC, Sp_BM). PPV, by
contrast, is a **local** property of each outbreak's case definition and setting and is not assumed
to transport.

**Community outbreak layer.** Each outbreak *o* contributes suspected cases $N_o$, of which $n_o$
were blood-cultured and $k_o$ were culture-positive. Selective culturing is common — clinicians
preferentially culture patients who look more like typhoid — so the tested subset is not a random
sample of the suspected. We therefore model the culture-selection process explicitly as a
verification-bias multinomial over three outcomes (culture-positive, culture-negative, untested),
with culture-selection probabilities $q_0$ and $q_1$ for non-typhoid and true-typhoid suspects and a
selectivity contrast $\delta = \text{logit}\,q_1 - \text{logit}\,q_0 = \mu_\delta + \beta_{\text{bias}} B_o + \sigma_\delta z_o$,
where $B_o\in\{0,1,2\}$ grades each study's documented verification bias. Community PPV is given a
hierarchical prior,

$$\text{logit}\,\pi_o \;=\; \mu_\pi + \sigma_\pi\, z_o, \qquad z_o\sim\mathcal N(0,1),$$

so that the eight anchor outbreaks that supply a genuine blood-culture-specific tested/positive pair
inform both their own $\pi_o$ and the population distribution $(\mu_\pi,\sigma_\pi)$. For outbreaks
without an admissible blood-culture pair we assign the **posterior-predictive population PPV**
$\pi_{\text{pop}} = \text{logit}^{-1}(\mu_\pi + \sigma_\pi z)$, which propagates the full
between-outbreak heterogeneity rather than a point value.

**Admissibility of the anchor counts.** Blood-culture sensitivity is specimen-specific, so only
blood-culture tested/positive pairs are admissible; counts that pool blood with stool, urine, rectal
swab, or bone marrow are excluded. Two recurring corrections in the primary literature are applied:
Neil 2012 (Kasese, Uganda) contributes **19/54** blood cultures — not the 25/63 that arises from
adding an unlabelled blood-or-stool source — and Kabwama 2017 contributes **56/364**. N'Cho 2019
(74/286) is retained as a real tested/positive pair in the 583-patient frame but flagged as a
mixed-specimen confirmation definition for a leave-one-out check; the remaining Harare reports, which
pool specimens, are excluded on that basis.

The model is fit in Stan (cmdstanr, four chains). We report posterior medians and 90% CrIs and
monitor divergences, $\hat R$, and effective sample size.

### 2.y Does a case-definition covariate improve PPV assignment?

Community PPV plausibly rises with the stringency of the clinical case definition — in particular
when the definition requires *excluding alternative diagnoses* — which raises the question of whether
unanchored outbreaks should receive a covariate-adjusted PPV rather than the pooled population value.
With only eight anchors this must be tested, not assumed. We coded one pre-specified, objectively
extractable covariate — the presence of an **alternative-diagnosis exclusion clause** in the case
definition (present for 3 of 8 anchors) — and entered it with a **monotone, non-negative** increment
that encodes the one certainty (an exclusion clause cannot lower PPV) and shrinks to pooling when the
data are silent:

$$\text{logit}\,\pi_o \;=\; \mu_\pi + \delta\,\text{clause}_o + \sigma_\pi z_o, \qquad \delta\sim\text{half-Normal}(0,s),\ \delta\ge 0.$$

A flag recovers the exact pooled model ($\delta$ removed) for a like-for-like comparison. The two
models are compared by leave-one-anchor-out predictive accuracy (WAIC-based expected log predictive
density, ELPD) on the eight community observations, and $\delta$ is examined for prior sensitivity
across $s\in\{0.25,0.5,1.0\}$.

### 2.z Renewal-equation model of outbreak-response vaccine impact

Vaccine impact is estimated with a semi-mechanistic **renewal-equation** model in which indirect
(transmission) effects **emerge** from feedback rather than being imposed as a fixed parameter. For
each outbreak we reconstruct the instantaneous reproduction number from the incidence $I_t$ (under the
primary additive observation model this is the PPV-de-backgrounded true-typhoid series; Section 2.w)
using the Cori estimator, $R_t = I_t / \Lambda_t$, where $\Lambda_t = \sum_{u\ge 1} w_u I_{t-u}$
is the total infectiousness and $w$ is the discretised generation interval (Gamma, mean 14 days,
CV 0.6; weekly time step). Outbreak-response immunisation reduces transmission multiplicatively from
the protection week $t_{\text{eff}}$ onward, $R_t \to R_t\,(1 - c\,\psi_T)$, with coverage
$c = 0.80$, transmission efficacy $\psi_T \in [0.6\psi,\ \psi]$, and direct efficacy $\psi = 0.83$;
the counterfactual epidemic is then propagated forward through the same renewal recursion, so that
averted infections include both directly and indirectly protected cases. Base-case protection is
reached at week $t_{\text{eff}} = 13$ (an 8-week campaign-request delay plus a 14-day campaign and a
21-day immunity-onset lag), reflecting realistic operational timing. Of the 19 reconstructed
outbreaks, the 13 with sufficient temporal resolution ($\mu_g/\Delta \ge 2$) enter the renewal
analysis.

### 2.w From true to observed impact: the PPV bridge and its weekly-allocation robustness

The vaccine acts on true typhoid infections, but surveillance counts suspected cases, so the bridge
between them is an explicit observation model — and that choice is precisely the choice of *which
curve $R_t$ is reconstructed from*. Two structurally distinct models bracket the possibilities. Under
a **multiplicative** model the true series is a constant fraction of the suspected, $I_t = \pi S_t$;
here $R_t$ and every *percentage* reduction are PPV-invariant (the constant $\pi$ cancels in
$R_t = I_t/\Lambda_t$ and in any ratio), so reconstructing on the suspected curve is legitimate and
the observed reduction equals the true reduction. Under an **additive** model — **our primary
specification** — the suspected series is true typhoid on a non-typhoid febrile **background**,
$S_t = I_t + B_t$; we de-background the curve ($B$ taken constant, $I_t = \max(S_t - B, 0)$
renormalised so $\sum_t I_t = \pi\sum_t S_t$), **reconstruct $R_t$ from the de-backgrounded
true-typhoid incidence $I_t$**, and let the vaccine act on $I_t$ only. Because the vaccine cannot
reduce non-typhoid fever, $B$ is unchanged, and the same absolute number of true infections averted is
diluted against it, so the reduction surveillance *sees* is smaller than the true reduction. These are
genuine structural alternatives: reconstructing $R_t$ from the suspected curve rather than the
de-backgrounded one is not a neutral preprocessing choice — it *is* the multiplicative assumption.

A subtlety is that the cumulative PPV pins the *total* true burden, $\sum_t I_t = \pi \sum_t S_t$, but
**not** the weekly allocation $I_t$. To show that our conclusions do not rest on an arbitrary choice
of $I_t$, we bracket each outbreak over five epidemiologically motivated weekly allocations —
multiplicative ($I_t=\pi S_t$), additive excess over a flat background, a declining/front-loaded
time-varying PPV, a renewal-consistent smoothing, and a reporting-shock spike-then-background — and
restrict each outbreak to the subset of allocations consistent with its transmission mode
(common-source, propagated, or mixed), case definition, and documented reporting artefacts.

Throughout the bridge we use a **background-fixed counterfactual**: because outbreak-response
immunisation cannot avert non-typhoid fever, the non-typhoid background $B_t = S_t - I_t$ is held
fixed between factual and counterfactual worlds. Under this counterfactual the reduction that
surveillance observes and the true reduction are related exactly by the PPV,

$$\text{suspected reduction} \;=\; \pi \times \text{true reduction},$$

**independently of how the true cases are allocated across weeks.** The allocation choice therefore
moves the *magnitude* of averted true infections but never the PPV dilution factor that separates
observed from true impact — the central qualitative result is allocation-invariant by construction.

---

## Results

### 3.x Blood-culture sensitivity and community PPV

The paired-specimen anchors identify a blood-culture sensitivity that rises with inoculated volume:
Se_BC is **0.56 (90% CrI 0.46–0.66)** at 2 mL, **0.65 (0.56–0.73)** at 5 mL, and **0.72 (0.62–0.80)**
at 10 mL, with a positive log-volume slope $\alpha_1 = 0.42\ (0.18{-}0.66)$. Bone-marrow sensitivity
is high (Se_BM 0.92, 0.90–0.94) and blood-culture specificity near-perfect (Sp_BC 0.996), confirming
that observed positivity is essentially $\pi\times\text{Se}_{BC}$ and that the naïve divide-by-0.6
shortcut is not admissible.

Community PPV is **low and highly heterogeneous** across the eight anchor outbreaks
(Table 1), ranging from 0.10 (Dhadwal 2008, 5/92) to 0.86 (Hu 2022, 24/31). The fitted
between-outbreak SD is large ($\sigma_\pi = 1.35$ on the logit scale, 0.85–2.03), and the
posterior-predictive **population PPV is 0.31 (0.16–0.52)** — i.e. under a typical outbreak case
definition roughly two of every three suspected cases are not typhoid. The verification-bias term is
positive ($\beta_{\text{bias}} = 0.39$, 0.05–0.94), confirming that clinicians preferentially culture
patients more likely to be truly infected, so that raw positivity overstates PPV where culturing was
selective.

**Table 1. Community-PPV anchors (culture-positive/tested; posterior PPV).**

| Outbreak | Positive/tested | Positivity | PPV π (90% CrI) |
|---|---|---|---|
| Dhadwal 2008 | 5/92 | 0.05 | 0.10 (0.04–0.24) |
| Kabwama 2017 | 56/364 | 0.15 | 0.15 (0.04–0.42) |
| Makungo 2020 | 7/65 | 0.11 | 0.16 (0.07–0.39) |
| Aye 2004 | 3/22 | 0.14 | 0.19 (0.06–0.47) |
| N'Cho 2019 | 74/286 | 0.26 | 0.32 (0.19–0.58) |
| Neil 2012 | 19/54 | 0.35 | 0.38 (0.14–0.75) |
| Swaddiwudhipong 2001 | 5/10 | 0.50 | 0.60 (0.30–0.89) |
| Hu 2022 | 24/31 | 0.77 | 0.86 (0.69–0.97) |
| **Population (posterior-predictive)** | — | — | **0.31 (0.16–0.52)** |

### 3.y A case-definition covariate does not improve PPV assignment

The exclusion-clause covariate is directionally correct but does not earn its place. Under the
reference prior the clause raises PPV (odds ratio 1.66, 90% CrI 1.06–3.19; group PPV 0.28 [0.14–0.49]
without a clause vs 0.40 [0.20–0.64] with one), yet it does **not** improve out-of-sample prediction:
the leave-one-anchor-out ELPD difference versus pooling is $+0.26 \pm 0.23$ — a statistical tie —
and the residual heterogeneity barely moves ($\sigma_\pi$ 1.35 → 1.28). The *predictive* PPV for a
new outbreak is essentially indistinguishable between groups (0.28 [0.03–0.80] without vs
0.40 [0.06–0.88] with a clause). The clause effect is moreover prior-driven: the posterior $\delta$
tracks the prior scale almost linearly (median 0.21, 0.51, 1.14 for $s = 0.25, 0.50, 1.0$), the
signature of a covariate that eight anchors — three of them carrying the clause, one of which
(Kabwama, 0.15) is a within-group counterexample — cannot pin. We therefore retain the **pooled
population PPV as the primary assignment** for unanchored outbreaks and report the clause model only
as a pre-specified sensitivity analysis; reassigning unanchored outbreaks by clause status shifts
their PPV modestly (0.28 or 0.40 vs the pooled 0.31) and does not change the aggregate impact.

### 3.z Vaccine impact: large true reduction, small surveillance-visible reduction

Across the renewal-eligible outbreaks, outbreak-response immunisation at realistic timing averts on
the order of **1,570 true typhoid infections**, a **≈30% reduction in true typhoid** burden
(pooled median 29.8%). The same intervention reduces the **surveillance-visible (suspected)** count
by only **≈6%** (6.3%). The gap is the PPV bridge in action: because the vaccine cannot avert the
non-typhoid febrile background, and because PPV is low (population ≈0.31, and lower still in the
large outbreaks that dominate the pooled denominator), the true reduction is diluted roughly
**five-fold** by the time it reaches routine surveillance. A programme that genuinely averts a third
of typhoid cases will register as a single-digit percentage dip in a suspected-case time series.

### 3.w Robustness to weekly allocation and to operational timing

The qualitative result is invariant to how true cases are distributed within an outbreak, exactly as
the background-fixed identity predicts. Bracketing each outbreak over its epidemiologically plausible
weekly allocations at production timing, aggregate true typhoid averted lies in the envelope
**1,085–1,541** across the 13 renewal outbreaks (26,298 suspected cases; 6,678 true typhoid; mean
PPV 0.25), a **16.2–23.1% true reduction** that surveillance would see as only **4.1–5.9%**. Because
the counterfactual holds the non-typhoid background fixed, suspected reduction equals PPV times true
reduction in every allocation, so the allocation choice moves the magnitude of the envelope but never
the observed-to-true ratio.

The envelope sits below the 29.8% headline true reduction for a concrete, reportable reason:
**production timing** — an 8-week request delay plus campaign and immunity-onset lags — places
protection at week 13, after three short outbreaks (including two of the larger ones) have already
ended, so those contribute essentially no averted cases. Faster deployment is the single largest
lever on realised impact. Whichever timing and allocation are assumed, however, both the primary
renewal engine (29.8% true → 6.3% observed) and the allocation-bracket analysis (16–23% true →
4–6% observed) agree on the load-bearing message: **the reduction in true typhoid is several times
larger than the reduction a suspected-case surveillance system will record, and PPV is the factor
that separates them.** This has a direct programmatic implication — evaluations that judge
outbreak-response immunisation on suspected-case counts will systematically understate its true
impact, and culture-confirmation (or an explicit PPV correction) is required to recover it.
