# Literature anchors for the transmission reproduction number (R_tr)

Mines the outbreak reports + the broader typhoid household-transmission literature for
field measurements of **person-to-person transmission**, which load on R_tr *independently
of the epidemic-curve shape* — the only information that can break the source/transmission
confound that left the source fraction θ unidentified (`hierarchical_theta/`).

## Why this works (and its limit)
Curve-derived R is confounded with the source pulse (same ridge). Household secondary
attack rates, contact-tracing secondary counts, and case-control attributable fractions are
observed *separately* from the incidence trajectory, so they identify R_tr along the
direction the curve cannot. **But:** household measurements capture only **short-cycle
(within-household)** transmission. The model's R_tr lumps in **long-cycle (community,
water-mediated)** transmission — which is the dominant propagated channel in these typhoid
outbreaks *and* the part actually confounded with the source. No field measurement isolates
it. So the anchors constrain the **minor** channel; they are unlikely to fully resolve the
identifiability problem.

## Independent vs confounded (the key distinction)
- **Independent anchors:** household SAR, secondary-case counts, attributable fraction,
  proportion-secondary, contact-derived R. These observe transmission directly.
- **Confounded estimates:** growth-rate / curve-derived R — **discarded as anchors**.

## Quantity → R_tr mapping (transparent, kept separate from the raw values)
- Household secondary cases per index ≈ SAR × (mean household size − 1) → household R_tr.
- Proportion-secondary p_sec → one-generation household R_tr ≈ p_sec/(1−p_sec).
- Attributable fraction or proportion-secondary → prior on (1−θ) = R_tr (household channel).
- Contact OR with **no denominator** → ordinal only (transmission present vs negligible);
  widen the prior, do not pin the level.
- No outbreak-specific data → pooled prior from `tab_typhoid_transmission_priors.csv`
  (flagged as shared, not outbreak-specific).

## Findings — coverage is thin (as expected)
Of the 13 renewal outbreaks: **5/13** report any independent p2p quantity; only **2/13**
(Polonsky's two suburbs, prop-secondary 5.3% = 135/2570) give a level-anchoring quantity
*with a denominator*; **3/13** give ordinal contact ORs (Aye 22, Muti 8.34, Qamar 3.81 — no
denominator); **8/13** fall back to the pooled Vellore prior. Michel 2005 (excluded from the
pipeline) gives a clean *lower* anchor: explicit **0% secondary**, R_tr≈0. No outbreak
reports a true household SAR with a denominator, so the proposed binomial likelihood term has
no data to attach to.

Pooled prior (Srinivasan/Marudaiappan et al. 2021, Vellore): **7/56 (12.5%)** households had a
blood-culture-confirmed secondary case; 12/172 (7%) contacts shed S. Typhi; serial interval
~12 d. → one-generation household R_tr ≈ 0.14.

## Honest bottom line
The independent evidence that exists points to **low household transmission** (~5–12.5%),
consistent with source/community-dominated typhoid. It constrains a household sub-component
of R_tr but **not** the community-level R_tr that drives the confound. Adding these as
likelihood terms / priors is worth trying in `hierarchical_theta/` (then re-run the Phase-1
gate), but should be expected to tighten only the household channel — it is unlikely on its
own to make the total θ identifiable.

## Files
`tab_Rtr_anchors.csv` (per-outbreak extraction), `tab_typhoid_transmission_priors.csv`
(pooled literature), `tab_Rtr_priors_for_model.csv` (per-outbreak constraint for the model),
`run_rtr_anchors.R` (coverage + figure + render), `report_Rtr_anchors.qmd`.
Run: `Rscript rtr_anchors/run_rtr_anchors.R`.
