# Assemble the self-contained catch-up HTML (embeds figures as base64).
#   Rscript manuscript/build_summary_html.R <out.html>
setwd(Sys.getenv("RENEWAL_ROOT", "."))
out <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out)) out <- "manuscript/summary.html"
fig <- "manuscript/figures"
b64 <- function(f) jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size))
img <- function(f, alt) sprintf('<img src="data:image/png;base64,%s" alt="%s">', b64(file.path(fig, f)), alt)

par <- read.csv("latent_class_ppv/tables/final_parameters.csv", check.names = FALSE)
gp <- function(nm) { r <- par[par$param == nm, ]; sprintf("%.2f <span class=ci>(%.2f–%.2f)</span>", r$med, r$`lo.5.`, r$`hi.95.`) }

css <- '
:root{--paper:#FBFAF7;--card:#FFFFFF;--ink:#1A1D1E;--muted:#5B6566;--line:#E6E2D9;
--accent:#0E6E78;--accent-soft:#E1EEEF;--amber:#B85C1A;
--serif:"Charter","Iowan Old Style","Palatino Linotype","Book Antiqua",Georgia,serif;
--sans:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
--mono:ui-monospace,"SF Mono","Cascadia Code",Menlo,Consolas,monospace;}
@media (prefers-color-scheme:dark){:root{--paper:#14171A;--card:#1B1F22;--ink:#ECEEEA;--muted:#9AA6A7;
--line:#2C3234;--accent:#54B7C2;--accent-soft:#16302F;--amber:#E0904E;}}
:root[data-theme="light"]{--paper:#FBFAF7;--card:#FFFFFF;--ink:#1A1D1E;--muted:#5B6566;--line:#E6E2D9;
--accent:#0E6E78;--accent-soft:#E1EEEF;--amber:#B85C1A;}
:root[data-theme="dark"]{--paper:#14171A;--card:#1B1F22;--ink:#ECEEEA;--muted:#9AA6A7;
--line:#2C3234;--accent:#54B7C2;--accent-soft:#16302F;--amber:#E0904E;}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);line-height:1.6;
-webkit-font-smoothing:antialiased}
.wrap{max-width:920px;margin:0 auto;padding:56px 24px 80px}
.eyebrow{font-family:var(--mono);font-size:12px;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);margin:0 0 10px}
h1{font-family:var(--serif);font-weight:600;font-size:2.5rem;line-height:1.1;margin:0 0 12px;text-wrap:balance}
h2{font-family:var(--serif);font-weight:600;font-size:1.6rem;margin:0 0 6px;text-wrap:balance}
.lede{font-size:1.12rem;color:var(--muted);max-width:64ch;margin:0 0 8px}
.rail{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px;margin:32px 0 8px}
.stat{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 16px}
.stat .n{font-family:var(--serif);font-size:1.7rem;font-weight:600;color:var(--accent);font-variant-numeric:tabular-nums}
.stat .l{font-family:var(--mono);font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);margin-top:2px}
section{margin-top:56px;border-top:1px solid var(--line);padding-top:36px}
.kicker{font-family:var(--mono);font-size:12px;letter-spacing:.12em;text-transform:uppercase;color:var(--amber)}
figure{margin:26px 0;background:var(--card);border:1px solid var(--line);border-radius:14px;padding:14px;overflow:hidden}
figure img{width:100%;height:auto;display:block;border-radius:8px}
figcaption{font-size:.94rem;color:var(--muted);margin-top:12px;padding:0 4px}
figcaption b{color:var(--ink);font-weight:600}
.tbl-wrap{overflow-x:auto;margin:22px 0}
table{border-collapse:collapse;width:100%;font-size:.95rem}
caption{text-align:left;font-family:var(--mono);font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);padding-bottom:8px}
th,td{text-align:left;padding:9px 14px;border-bottom:1px solid var(--line)}
thead th{font-family:var(--mono);font-size:11px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);font-weight:600}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums}
.ci{color:var(--muted);font-size:.86em}
.note{background:var(--accent-soft);border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin-top:20px;font-size:.94rem}
.note b{color:var(--accent)}
footer{margin-top:56px;border-top:1px solid var(--line);padding-top:22px;color:var(--muted);font-size:.85rem}
code{font-family:var(--mono);font-size:.88em;background:var(--accent-soft);padding:1px 6px;border-radius:5px}
'

html <- paste0('<div class="wrap">',
'<p class="eyebrow">Typhoid · diagnostics + vaccine impact</p>',
'<h1>PPV of the clinical case definition &amp; outbreak-response immunization</h1>',
'<p class="lede">A results catch-up: Bayesian latent-class estimates of typhoid diagnostic accuracy and the positive predictive value (PPV) of clinical case definitions, propagated into a renewal-equation model of reactive vaccination impact.</p>',
'<div class="rail">',
'<div class="stat"><div class="n">0.52–0.68</div><div class="l">Blood-culture Se (2–10 mL)</div></div>',
'<div class="stat"><div class="n">0.90</div><div class="l">Bone-marrow Se</div></div>',
'<div class="stat"><div class="n">0.23–0.65</div><div class="l">Community PPV π</div></div>',
'<div class="stat"><div class="n">~43%</div><div class="l">ORI true reduction (base)</div></div>',
'<div class="stat"><div class="n">~4×</div><div class="l">Suspected overcount vs true</div></div>',
'</div>',

'<section><p class="kicker">Part 1</p><h2>The latent-class model</h2>',
'<p class="lede">There is no gold standard for typhoid. A Bayesian latent-class model (<code>model_final.stan</code>) estimates test accuracy and case-definition PPV jointly, from two data layers that share a common set of <em>transportable</em> test parameters.</p>',
'<figure>', img("fig0_latent_class_model.png", "latent-class model schematic"),
'<figcaption><b>Model structure.</b> <b>Left</b> &mdash; 10 historic studies with paired blood + bone-marrow culture give 2&times;2 counts modelled as a Hui&ndash;Walter conditional-independence latent-class multinomial; this identifies blood-culture sensitivity, bone-marrow sensitivity, and the hospital PPV φ. <b>Right</b> &mdash; outbreak suspected&rarr;confirmed counts are binomial in <code>π·Se_BC + (1&minus;π)(1&minus;Sp)</code>, identifying the community PPV π. <b>Bridge</b> &mdash; <code>logit&nbsp;Se_BC = α0 + α1·log(volume) + τu</code> (Antillón-informed slope) makes sensitivity transportable by blood volume, so the paired layer informs the outbreak layer.</figcaption></figure>',
'<div class="note"><b>Why the split matters.</b> <b>Transportable</b> test parameters (Se_BC(volume), Se_BM, Sp) are reused across settings. <b>Local</b> population parameters &mdash; hospital φ and community π &mdash; are the PPV of the clinical case definition, and they are <em>not</em> one transferable constant: the same definition means P(true typhoid | suspected) ≈ 0.9 in a referral hospital but ≈ 0.2 in community surveillance. Identifiability comes from the paired bone-marrow layer plus the volume prior breaking the π&middot;Se ridge.</div>',
'<h2 style="margin-top:40px">PPV &amp; culture accuracy &mdash; estimates</h2>',
'<figure>', img("fig1a_culture_accuracy.png", "culture accuracy"),
'<figcaption><b>Blood-culture sensitivity rises with cultured volume</b> (0.52 &rarr; 0.68 from 2 to 10 mL), well below bone-marrow sensitivity (~0.90); specificity is pinned near 1.0. Orange points are the model estimates at 2/5/10 mL with 90% CrIs; grey points are Antillón 2018 study-level observations.</figcaption></figure>',
'<figure>', img("fig1b_ppv_by_setting.png", "ppv by setting"),
'<figcaption><b>PPV is not a constant &mdash; it depends on setting.</b> Hospital-referred patients (φ, green) yield PPVs of 0.66–0.99; community/outbreak surveillance (π, purple) yields 0.23–0.65. The same clinical definition means very different things in different places.</figcaption></figure>',
'<figure>', img("fig1c_ppv_spectrum.png", "ppv spectrum"),
'<figcaption><b>PPV surface across case definitions and prevalence.</b> Syndromic definitions (bare fever, WHO suspected) cap specificity low, so their PPV stays modest except at high prevalence; adding serology or a prediction rule lifts PPV sharply.</figcaption></figure>',
'<div class="tbl-wrap"><table><caption>Table 1 &mdash; Latent-class diagnostic accuracy (median, 90% CrI)</caption>',
'<thead><tr><th>Parameter</th><th>Estimate</th></tr></thead><tbody>',
'<tr><td>Blood-culture sensitivity @ 2 mL</td><td>', gp("Se_BC_2mL"), '</td></tr>',
'<tr><td>Blood-culture sensitivity @ 5 mL</td><td>', gp("Se_BC_5mL"), '</td></tr>',
'<tr><td>Blood-culture sensitivity @ 10 mL</td><td>', gp("Se_BC_10mL"), '</td></tr>',
'<tr><td>Bone-marrow sensitivity</td><td>', gp("Se_BM"), '</td></tr>',
'<tr><td>Blood-culture specificity</td><td>', gp("Sp_BC"), '</td></tr>',
'<tr><td>Volume slope α<sub>1</sub> (logit per log-mL)</td><td>', gp("alpha1"), '</td></tr>',
'</tbody></table></div></section>',

'<section><p class="kicker">Data provenance</p><h2>Positivity extraction &mdash; 40 local PDFs</h2>',
'<p class="lede">To test whether the community-PPV dataset could be expanded, suspected&rarr;tested&rarr;confirmed counts were extracted from 40 outbreak/surveillance PDFs. Implied PPV = positivity / Se_BC.</p>',
'<figure>', img("fig1d_ppv_extraction.png", "ppv extraction"),
'<figcaption><b>The extraction validates the model and separates two strata.</b> <b>Outbreak</b> studies with selective testing (triangles = the three already in the model) sit at π ≈ 0.2&ndash;0.6; the newly-quantified <b>endemic-surveillance</b> stratum (all-tested cohorts: STRATAA, Srikantiah, Thriemer, Voysey) sits ~5&ndash;10&times; lower at π ≈ 0.03&ndash;0.08.</figcaption></figure>',
'<figure>', img("fig1e_ppv_strata.png", "ppv endemicity strata fit"),
'<figcaption><b>Formal two-stratum fit of the endemicity gradient</b> (10 studies incl. the SEAP outpatient cohort, 20,899 tested). Outbreak π = <b>0.34</b> [0.16&ndash;0.59] versus endemic-surveillance π = <b>0.075</b> [0.04&ndash;0.17] &mdash; an odds ratio of <b>6.4</b> [1.7&ndash;16.9]. The clinical case definition is ~6&times; more predictive of true typhoid during an outbreak than in routine surveillance (CrI excludes 1); SEAP sits at the top of the surveillance stratum as a high-burden setting.</figcaption></figure>',
'<div class="note"><b>PPV vs reported &ldquo;PPV&rdquo; &mdash; the key correction.</b> SEAP reports a clinical-definition PPV of 10.5%, but that is <em>blood-culture-referenced</em>: P(culture+ | suspected) = positivity = π&middot;Se_BC, not our latent-true π. Entering SEAP&rsquo;s counts (20,899 tested / 2,116 <em>S.</em> Typhi) the model recovers π &asymp; <b>0.16</b> &mdash; ~1.6&times; higher than the reported 10.5%, because blood culture misses ~40% of true cases. That gap is exactly what the latent-class model corrects. <b>ORI unchanged:</b> no clean new <em>outbreak</em> points were recoverable and the fitted outbreak π (0.34) sits inside the existing 0.23&ndash;0.65 range. Model <code>community_ppv_strata.stan</code>; tables <code>community_ppv_strata.csv</code>, <code>ppv_positivity_table.csv</code>.</div></section>',
'<section><p class="kicker">Part 2</p><h2>Outbreak-response immunization (ORI)</h2>',
'<p class="lede">Outbreak surveillance counts <em>suspected</em> cases. Under an additive model <code>S(t) = I(t) + B(t)</code>, only the true typhoid signal I(t) responds to vaccination &mdash; the non-typhoid background B(t) does not.</p>',
'<figure>', img("fig2a_suspected_vs_true.png", "suspected vs true"),
'<figcaption><b>Observed suspected cases vs model-predicted true typhoid.</b> Grey is the suspected surveillance curve; red is the de-backgrounded true typhoid I(t). Where PPV is low (Kabwama, Aye: π≈ 0.23–0.26) most of the curve is non-typhoid background; where it is high (Neil: π = 0.65) most is true typhoid.</figcaption></figure>',
'<figure>', img("fig2b_ori_impact_grid.png", "ori impact grid"),
'<figcaption><b>Reactive vaccination impact is dominated by timing.</b> Pooled true-typhoid reduction falls from ~80% when protection arrives in week 1–2 to ~35% by week 12; higher coverage helps, but a late campaign at 90% still underperforms an early one at 65%.</figcaption></figure>',
'<figure>', img("fig2c_cases_averted_grid.png", "cases averted"),
'<figcaption><b>True typhoid cases averted</b> across the same timing × coverage grid, summed over the developing-country outbreak cohort.</figcaption></figure>',
'<div class="note"><b>The surveillance-dilution result.</b> Because the vaccine acts only on true typhoid, the reduction a surveillance system <em>observes</em> in suspected counts is diluted by roughly π. A campaign that truly cuts typhoid 34% in Kampala (Kabwama, π=0.23) shows only an ~8% drop in suspected cases; in Kasese (Neil, π=0.65) a true 97% reads as 63%. Ignoring PPV overcounts averted cases nearly 4× (cohort: 8,498 suspected-naive vs 2,330 true typhoid).</div>',
'</section>',

'<section><p class="kicker">Sources</p><h2>Definitive model &amp; data sources</h2>',
'<div class="tbl-wrap"><table><caption>PPV modeling</caption><thead><tr><th>Role</th><th>Source</th></tr></thead><tbody>',
'<tr><td>Model</td><td><code>latent_class_ppv/stan/model_final.stan</code> (companion: <code>community_ppv.stan</code>)</td></tr>',
'<tr><td>Fit / driver</td><td><code>latent_class_ppv/run_final.R</code> (cmdstanr, 4 chains, seed 2026)</td></tr>',
'<tr><td>Paired BC/BM data</td><td><code>data/mogasale2016_paired_bc_bmc_seed.csv</code> &middot; Mogasale 2016</td></tr>',
'<tr><td>Outbreak positivity</td><td><code>data/community_surveillance_ppv.csv</code> &middot; Aye/Neil/Kabwama</td></tr>',
'<tr><td>Volume-Se prior</td><td><code>data/antillon2018_volume_sensitivity.csv</code> &middot; Antillón 2018</td></tr>',
'<tr><td>Estimates out</td><td><code>final_parameters.csv</code>, <code>final_pi_community.csv</code>, <code>final_phi_hospital.csv</code>, <code>final_draws.rds</code></td></tr>',
'</tbody></table></div>',
'<div class="tbl-wrap"><table><caption>ORI modeling</caption><thead><tr><th>Role</th><th>Source</th></tr></thead><tbody>',
'<tr><td>Engine</td><td><code>renewal/R/renewal_core.R</code> (R_t + counterfactual), <code>scenario.R</code>, <code>ppv.R</code> (π propagation)</td></tr>',
'<tr><td>Driver / config</td><td><code>renewal/run_analysis.R</code> &middot; <code>renewal/config.yml</code></td></tr>',
'<tr><td>CEA (reused)</td><td><code>R/2-functions.R</code> (DALY/cost)</td></tr>',
'<tr><td>Outbreak time series</td><td><code>data/Typhoid_Outbreak_Time_Series_2000_2022_{Timeseries,Summary}.csv</code></td></tr>',
'<tr><td>PPV posterior input</td><td><code>latent_class_ppv/outputs/final_draws.rds</code></td></tr>',
'</tbody></table></div>',
'<p class="lede" style="font-size:.95rem">Methods lineage: Hui &amp; Walter 1980 (latent class); Mogasale 2016, Antillón 2018 (culture accuracy); Cori 2013 / Wallinga&ndash;Teunis (renewal R_t). Full write-up: <code>manuscript/MODELS_AND_SOURCES.md</code>.</p></section>',
'<footer>Latent-class PPV model (paired blood/bone-marrow culture; volume-anchored Se) &middot; renewal-equation ORI model with PPV posterior propagation (additive S=I+B regime). Estimates are advisory / proof-of-concept; the ORI cohort is the 13 resolution-eligible outbreaks from the 2000–2022 time-series dataset. Background B modeled constant &mdash; a data-derived B(t) is a one-function swap.</footer>',
'</div>')

writeLines(paste0("<style>", css, "</style>", html), out)
cat("wrote", out, "(", round(file.info(out)$size/1024), "KB )\n")
