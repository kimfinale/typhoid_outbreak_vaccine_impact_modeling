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
| **Harare, Zimbabwe** | **Near-annual**: cholera 2008/09/10, typhoid 2010, 2011–12, 2012, 2016–17, 2017–18 ("4th enteric epidemic since 2008"); water-shortage-driven |
| **Delmas, South Africa** | 1993 (>1,000) → 2005 (>600), same strain (12-yr gap) + isolated 2007/2009 |
| Kasese/Bundibugyo, Uganda | historically recurrent water/food-borne (cholera+typhoid), 2008–2011 |
| **Fiji** | ≥13 outbreaks since 2005, cyclone-linked |
| Pakistan (Sindh) | endemic MDR since 1980s; XDR clonal expansion Nov 2016 → 486 cases by Dec 2017 |
| Myanmar | ~4 major water-borne outbreaks/decade since 1989 |
| Zambia/Zimbabwe/DRC | concurrent regional outbreaks from Nov 2010 (H58), persisting ~33 months |
| Songkhla, Thailand | 2 waves (2009–10, 2010–11), 5-mo remission, then back to endemic norm |

## 3. Direct beyond-outbreak vaccine evidence (the payoff)
**Scobie 2014 (Fiji, Vi-PS post-cyclone campaign):** high-vaccination subdivision incidence
**IRR = 0.23 (~77 % reduction)** while low/no-vaccination subdivisions rose **2–8×**. This is
empirical support for the modelling assumption that a campaign suppresses endemic/background
incidence for the protection window (not only the index outbreak).

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

## Data caveats (mislabelled Zotero items found during extraction)
- "Davis 2018" Zotero PDF is actually the **N'Cho 2019** report (duplicate/misfiled).
- "Lutterloh 2012" is the **Malawi–Mozambique** border outbreak (Neno/Tsangano), not Uganda.
- "Poncin 2018" PDF is a **cholera** paper (wrong file for the typhoid entry).
