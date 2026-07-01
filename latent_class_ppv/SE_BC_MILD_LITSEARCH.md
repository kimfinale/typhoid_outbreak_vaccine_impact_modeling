# Does blood-culture sensitivity differ in MILD cases? (literature search -> `beta` verdict)

Question driving this: the model had a severe->mild offset `beta` on logit(Se_BC). Is it
identifiable / supported by data, or should it be dropped? Searched the Zotero library
(surveillance cohorts, diagnostics review, quantitation) + web (paired-in-mild; load-vs-severity).

## What we looked for, and what exists

| Target | Would it identify `beta`? | Found? |
|---|---|---|
| Paired BC + bone-marrow culture in mild/outpatient/community cases | **yes** | **NO — none exists** |
| Blood-culture Se stratified by severity (outpatient vs inpatient) | partly | **NO direct comparison found** |
| Bacterial load (CFU/mL) vs disease severity | mechanistic | **weak / mostly absent** |
| Community-surveillance culture positivity (setting PPV) | grounds the level | yes (but see below) |

## Findings (sources retrieved)

**1. No paired BC/BMC data in mild cases — structurally, and confirmed.**
Antillon 2018: **38/40 paired studies were hospitalized inpatients**; bone-marrow aspiration's
invasiveness is the reason. Even Wain 2001's "81 nonsevere" cases were *hospitalized uncomplicated*
inpatients (severe = perforation/encephalopathy), not ambulatory, and it did **not** split Se by
severity. **No study anywhere compares Se_BC between outpatient and inpatient typhoid.** And
community typhoid surveillance uses **serology** (anti-Vi seroconversion, STRATAA), not blood
culture — you cannot blood-culture a random community sample — so the mild-population Se_BC is
**structurally unmeasurable**. `beta` cannot be identified from data that could exist.

**2. Severity is not a supported independent determinant of Se_BC.** The definitive quantitation
studies (Wain 1998, blood; Wain 2001, marrow) **excluded severe cases** and found blood load
median ~1 CFU/mL driven by **age (children higher), illness week (wk1 1.7 -> wk4 0.3), MDR, and
fecal excretion — not complications**. Wain 1998 explicitly: typhoid pathology "cannot be explained
by differences in the numbers of organisms circulating in the blood." Only *marrow* load weakly
tracked fever-clearance time; *blood* load showed no clinical-severity correlation. Andrews & Ryan
2015 (review) ties Se to load/week/age/antibiotics, **not severity**, and notes hospitalized
patients are *more* likely pre-treated with antibiotics (which *lowers* yield) — so even the SIGN
of a severe->mild gap is ambiguous.

**3. The one severity-direction signal is really a load/fever proxy.** Andrews et al. JID 2021
(children): blood-culture positivity rose with **peak fever OR 3.77**, age OR 1.96, blood volume OR
2.82; prior antibiotics non-significant. Higher fever ~ higher bacteremia ~ higher yield — i.e. the
same load mechanism as volume/duration, not a distinct "severity" axis, and it's positivity (pre-
test prob) not a paired-sensitivity estimate.

**4. Endemic-surveillance culture PPV is low and similar across settings** — and far below outbreak
PPV. Typhoidal positivity among blood-cultured febrile: SETA (severe-inclusive) **1.3%**, Mahende
(outpatient children) **2.1%**, Pemba (low-endemicity) **2.1%**, STRATAA **3.2-8.2%**, SEAP
outpatients **~10.1%** — versus **14-46%** in our outbreaks. The mild/outpatient settings do NOT
show a lower per-culture typhoid yield than severe-inclusive ones. PPV is **context-driven
(outbreak >> endemic)**, not a severity property — and not transferable as a constant.

**5. What IS solidly quantified (and modelable):** Se_BC ~ **60%** (Mogasale 61-66%, SETA 50-70%,
Pemba 50%); driven by **blood volume** (age-differentiated: children 1-3 mL, adults 8-10 mL),
illness **week**, and **prior antibiotics**. These are identifiable covariates (volume already in
the model via Antillon).

## Verdict: drop `beta`

- **Unidentifiable in principle** (target #1 data cannot exist).
- **Unsupported as an axis** — severity is not a data-backed independent determinant of Se_BC;
  what drives Se_BC is volume, age, illness-week, and antibiotics, all of which act through
  bacterial load and are (partly) modelable without a severity term.
- **Ambiguous even in sign** (milder -> lower load -> lower Se, but hospitalized -> more prior
  antibiotics -> lower Se).

**Recommendation:** remove the severe->mild `beta` offset. Transfer Se_BC across settings via the
**volume model** (identified, Antillon-grounded), optionally add **prior-antibiotic fraction** and
**illness-week** as covariates where reported. Community/outbreak Se_BC = the volume-appropriate
value; carry the residual population uncertainty in the between-study `tau` + a **Se_BC prior
sensitivity sweep** (as the recovery gate already does). Do NOT model a distinct severity offset.

## Sources
Antillon 2018 JID 218(S4):S255 (PMC6226661); Wain 1998 JCM 36:1683 (PMC104900); Wain 2001 JCM
39:1571 (PMC87972); Andrews & Ryan 2015 Vaccine 33:C8; Andrews et al. JID 2021 224(S5):S484;
Mogasale 2016 Ann Clin Microbiol Antimicrob 15:32; von Kalckreuth 2016 TSAP CID; Meiring 2021
STRATAA; Marks 2022 SETA; Mahende 2015; Thriemer 2012 Pemba; Carey/Barkume 2020-22 SEAP.
