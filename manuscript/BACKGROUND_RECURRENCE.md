# Background incidence & outbreak recurrence in the tf outbreak studies

For modelling **TCV's long-term impact beyond the ongoing outbreak** (reducing endemic
background incidence + averting future recurrences for the protection window). Extracted
from 20 source PDFs of the tf outbreak dataset. Table: `data/tf_background_recurrence.csv`.

## Headline
- **Pre-outbreak background *incidence rate* is rarely reported** at the outbreak locality
  (9/20 give any incidence figure; most are *regional/national*, not a local pre-outbreak baseline).
- **Recurrence is well documented** (17/20): several ORI sites are *repeat-outbreak* regions.
- One study empirically shows the **beyond-outbreak vaccine effect** the user hypothesises.

## 1. Background (pre-outbreak / endemic) incidence — what the papers give
| Study (region) | Reported incidence | Type |
|---|---|---|
| **Limpitikul (Songkhla, Thailand)** | **1.8/100k/yr pre-outbreak** → 25.2 during → back to baseline after | **local pre-outbreak baseline (cleanest)** |
| **Vellore, India (SEFI cohort)** | **1,977/100k person-yrs** (95% UI 1,740–2,236), children 6mo–15y | cohort endemic incidence (SEFI, recovered from web) |
| Hechaichi (Tunisia) | 0.32/100k national; ~35 cases/yr | national (low-endemic) |
| Scobie (Fiji) | 33–52/100k/yr national; up to 136–1,052 Northern div | national + subnational |
| Qamar (Pakistan) | 413/100k (2–4y), 573/100k (5–15y) — Ochiai 2008 | regional (children) |
| Yousafzai (Pakistan/Karachi) | up to 1,000/100k child-yrs | regional (children) |
| Lutterloh (Malawi/Mozambique) | 10–100/100k/yr | national (WHO band) |
| Neil (Uganda) | Africa-wide 13–845/100k/yr | continental range |
| Imanishi (Zimbabwe) | sub-Saharan Africa 724.6/100k | regional |
| Makungo (South Africa) | national counts 6,000 (1985) → 200 (2002) | national trend |

Srinivasan (Vellore) defines the outbreak against a **background baseline** (mean + 2·SD of
prior-year cases under continuous SEFI surveillance) but does not print the rate here.

## 2. Outbreak recurrence / frequency — the strong signal
| Region | Recurrence pattern |
|---|---|
| **Harare, Zimbabwe** | **Annual since 2010** (per Poncin 2022, the TCV-campaign paper): typhoid 2010, 2011–12, 2012, 2016–17 (860 cases), 2017–18 (>4,300), 2018–19 (1,967) + cholera 2008/09/10; rainy-season (Oct–Apr), water-shortage-driven |
| **Delmas, South Africa** | 1993 (>1,000) → 2005 (>600), same strain (12-yr gap) + isolated 2007/2009 |
| Kasese/Bundibugyo, Uganda | historically recurrent water/food-borne (cholera+typhoid), 2008–2011 |
| **Fiji** | ≥13 outbreaks since 2005, cyclone-linked |
| Pakistan (Sindh) | endemic MDR since 1980s; XDR clonal expansion Nov 2016 → 486 cases by Dec 2017 |
| Myanmar | ~4 major water-borne outbreaks/decade since 1989 |
| Zambia/Zimbabwe/DRC | concurrent regional outbreaks from Nov 2010 (H58), persisting ~33 months |
| Songkhla, Thailand | 2 waves (2009–10, 2010–11), 5-mo remission, then back to endemic norm |

## 3. Direct beyond-outbreak vaccine evidence (the payoff)
**Scobie 2014 (Fiji, Vi-PS post-cyclone campaign):** high-vaccination subdivision incidence
**IRR = 0.23 (~77 % reduction)** while low/no-vaccination subdivisions rose **2–8×**. Empirical
support that a campaign suppresses endemic/background incidence for the protection window
(not only the index outbreak).

**Poncin 2022 (Harare 2019, TCV outbreak-response campaign):** the actual TCV campaign in the
hyper-recurrent Harare setting — 373,027 targeted across 9 suburbs, **85.4 % coverage** (72.4 %
under-5, 97.1 % ages 5–15). It reports **no post-campaign follow-up incidence** — i.e. the
long-term / beyond-outbreak impact was never measured. **That is exactly the gap this modelling
would fill:** given Harare's annual recurrence, projecting averted future outbreaks over the
TCV protection window is the missing quantification.

## Implications for the long-term-impact model
1. **The assumption is supported** — recurrence is real (Harare near-annual; Delmas; Fiji;
   Pakistan) and Scobie shows a sustained ~77 % post-campaign incidence drop.
2. **Background-incidence input** (rarely in the outbreak papers): use, in order,
   (a) the paper's own pre-outbreak baseline where given (Limpitikul 1.8; Hechaichi 0.32;
   Scobie 33–52/100k); (b) **the endemic-surveillance stratum already extracted** for the PPV
   work (STRATAA, SEAP, Srikantiah, Andrews, Thriemer, Hancuh — these are endemic incidence/
   positivity for the same regions); (c) GBD 2017 country estimates.
3. **Recurrence frequency** for a "future outbreaks averted" term is highly region-specific:
   ~annual (Harare) to decadal (Delmas 12 yr). Model per-site, not one global rate.

## Corrected sources (mislabelled Zotero items — now resolved via web)
- **Davis 2018** — real paper = MMWR "Typhoid Fever Outbreak, Harare, Oct 2016–Mar 2017"
  (PMID 29565843, PMC5868204); 860 cases (780 suspected / 80 confirmed / 4 deaths). The Zotero
  file was a duplicate of N'Cho 2019.
- **Poncin 2022** — real paper = "Implementation of an outbreak response vaccination campaign
  with typhoid conjugate vaccine – Harare, Zimbabwe, 2019" (PMC9379662), **not** the cholera
  paper in Zotero and **not** Zambia. Data extracted above.
- **SEFI/Vellore incidence** — recovered: 1,977/100k person-yrs (Srinivasan ref 4 = SEFI, John
  et al. 2021, NEJM/JID). Added to the background table.
- Still open: **"Lutterloh 2012"** is the **Malawi–Mozambique** border outbreak (Neno/Tsangano),
  mislabelled as Uganda in the tf dataset (correct paper, wrong tf label).

Corrected PDFs to add to Zotero: Davis → https://pmc.ncbi.nlm.nih.gov/articles/PMC5868204/ ;
Poncin → https://pmc.ncbi.nlm.nih.gov/articles/PMC9379662/ (Europe PMC render endpoint refused
automated download; use the PMC page).
