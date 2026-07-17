#!/usr/bin/env Rscript

# Merge three discovery routes for outbreak PPV evidence:
# 1. outbreak studies represented in Koh et al. (2025),
# 2. full texts excluded from Koh because no usable time series was available,
# 3. publications represented in Appiah et al. (2020) but not as distinct Koh
#    time series.
#
# One row is retained per audited study/publication label. Source membership is
# many-to-many and is represented by flags rather than duplicated rows.

koh_csv <- "latent_class_ppv/data/community_outbreak_ppv_eligibility_audit.csv"
excluded_csv <- "latent_class_ppv/data/covidence_no_timeseries_fulltext_ppv_audit.csv"
appiah_csv <- "latent_class_ppv/data/appiah2020_not_in_koh2025_fulltext_ppv_audit.csv"
koh_summary_csv <- "data/Typhoid_Outbreak_Time_Series_2000_2022_Summary.csv"
live_ppv_csv <- "latent_class_ppv/data/community_surveillance_ppv.csv"
output_csv <- Sys.getenv("MERGED_PPV_OUTPUT",
                         "latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv")
se_draws_rds <- "latent_class_ppv/outputs/final_draws.rds"

koh <- read.csv(koh_csv, check.names = FALSE, stringsAsFactors = FALSE)
excluded <- read.csv(excluded_csv, check.names = FALSE, stringsAsFactors = FALSE)
appiah <- read.csv(appiah_csv, check.names = FALSE, stringsAsFactors = FALSE)
koh_summary <- read.csv(koh_summary_csv, check.names = FALSE, stringsAsFactors = FALSE)
live_ppv <- read.csv(live_ppv_csv, check.names = FALSE, stringsAsFactors = FALSE)

audit_cols <- names(koh)

clean <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

key <- function(x) {
  x <- tolower(trimws(clean(x)))
  gsub("[[:space:]]+", " ", x)
}

for (nm in audit_cols) {
  koh[[nm]] <- clean(koh[[nm]])
  excluded[[nm]] <- clean(excluded[[nm]])
  appiah[[nm]] <- clean(appiah[[nm]])
}

# The pre-existing Koh PPV audit contained only 23 articles selected for an
# earlier positivity exercise, not Koh's complete article inventory. Add the
# 14 omitted Koh articles explicitly. Polonsky and Imanishi each contribute
# multiple outbreak rows in the canonical summary; their total/article row is
# used here so source inventory is counted at publication level.
missing_koh_articles <- c(
  "Al-Sanouri 2008", "Muehlen 2007", "Bayram 2011", "Hendriksen 2015",
  "Polonsky 2014 (total)", "Muti 2014", "Imanishi 2014 (total)",
  "Cherian 2015", "Burnsed 2018", "Davis 2018", "Hechaichi 2023",
  "N'Cho 2019", "Bano-Zaidi 2018", "Poncin 2022"
)

make_blank_rows <- function(ids) {
  z <- as.data.frame(setNames(replicate(length(audit_cols), rep("", length(ids)),
                                        simplify = FALSE), audit_cols),
                     stringsAsFactors = FALSE, check.names = FALSE)
  z$study <- ids
  z$summary_row_id <- paste0("Koh inventory: ", ids)
  z$legacy_attack_rate_pass <- "FALSE"
  z$current_stage1 <- "FALSE"
  z$denominator_unit <- "patients"
  z$verification_bias_risk <- "unassessable"
  z$source_status <- "canonical Koh summary reviewed; PPV fields require primary-paper audit"
  z$simple_model_status <- "not eligible pending compatible blood-culture denominator"
  z$verification_model_status <- "source reconstruction required"
  z$recommended_use <- "not yet: source reconstruction required"
  z$independence_group <- ids
  z
}

koh_added <- make_blank_rows(missing_koh_articles)
summary_ids <- koh_summary[["Study ID"]]
for (i in seq_len(nrow(koh_added))) {
  j <- match(koh_added$study[i], summary_ids)
  if (is.na(j)) next
  koh_added$location[i] <- koh_summary$Location[j]
  koh_added$outbreak_year[i] <- paste(unique(na.omit(c(koh_summary$`Outbreak start year`[j],
                                                        koh_summary$`Outbreak end year`[j]))),
                                      collapse = "-")
  koh_added$case_definition_summary[i] <- koh_summary$`Case definitions`[j]
}

set_added <- function(study, ...) {
  i <- match(study, koh_added$study)
  values <- list(...)
  for (nm in names(values)) koh_added[[nm]][i] <<- as.character(values[[nm]])
}

set_added("Muehlen 2007", audited_suspected_testing_frame="Six outbreak-linked suspects",
          audited_blood_tested="6", audited_blood_positive="4", audited_positivity="0.667",
          tested_fraction="1.000", specimen_scope="blood/stool confirmation; blood-flow detail requires recheck",
          denominator_compatibility="Small fully tested outbreak frame; confirm that the denominator and numerator are both blood-specific rather than pooled blood/stool",
          verification_bias_risk="high: exposure-defined outbreak frame; internal PPV may be valid but transportability is limited",
          simple_model_status="not primary until blood-specific flow is verified",
          verification_model_status="sensitivity candidate after specimen-scope verification",
          recommended_use="sensitivity only; small n is not an exclusion criterion",
          audit_notes="Observed 4/6 is not arithmetically inadmissible. Remaining concern is specimen/denominator compatibility and target-population transport, not sample size or positivity alone.")
set_added("Bayram 2011", denominator_compatibility="Duplicate report of the Van/Anatolia outbreak represented by Aypak 2010",
          simple_model_status="not eligible as an independent observation",
          verification_model_status="duplicate outbreak/report",
          recommended_use="exclude duplicate; retain Aypak 2010 counts",
          independence_group="Van/Anatolia 2008 outbreak",
          audit_notes="Explicit active decision: Aypak supplies the recoverable 34/154 blood-culture pair; Bayram is not an independent outbreak.")
set_added("N'Cho 2019", summary_suspected="583", audited_suspected_testing_frame="583 hospitalized patients in the Oct-Dec 2017 case-management review subset",
          audited_blood_tested="286", audited_blood_positive="74", audited_positivity="0.259",
          tested_fraction="0.489", specimen_scope="cultures processed for the audited hospital subset; specimen not disaggregated in the primary (blood/stool/rectal swab), likely predominantly blood in an inpatient setting",
          testing_selection_mechanism="Hospital-referred and selectively cultured subset; full outbreak had 3,187 suspected and 191 confirmed cases",
          denominator_compatibility="Compatible within the 583-patient hospital subset; not representative of the full community outbreak",
          verification_bias_risk="high",
          source_status="primary extraction in manuscript/data/ppv_positivity_table.csv; current stage-2 strata member",
          simple_model_status="not eligible without selection adjustment",
          verification_model_status="candidate with major verification adjustment",
          recommended_use="verification-adjusted candidate",
          independence_group="Harare 2017-2018 N'Cho-Poncin overlap",
          audit_notes="Current stage-2 anchor: 74/286. This hospital subset lies within the Oct 2017-Jun 2018 outbreak summarized retrospectively by Poncin 2022 (4,330 cases; 250 blood-or-stool confirmed). Retain N'Cho as the sole PPV observation; Poncin has no tested denominator. Specimen caveat (primary verified 2026-07-17, Zotero RAYYTM4G, MMWR 68[2]:44-45): the paper states 'cultures were processed for 286 (49%) inpatients; 74 (26%) yielded Typhi' under a blood/stool/rectal-swab confirmation definition, so 74/286 is not disaggregated to blood culture. Retained in the main model as a real tested/positive pair within the 583-patient frame; a leave-N'Cho-out sensitivity check is recommended given the specimen ambiguity.")
set_added("Poncin 2022",
          audited_suspected_testing_frame="Primary report covers the subsequent Aug 2018-Feb 2019 season (1,967 cases) and cites the prior Oct 2017-Jun 2018 season (4,330 cases) that contains N'Cho's cohort",
          specimen_scope="suspected and confirmed surveillance counts; historical confirmation pooled blood and stool",
          temporal_alignment="Historical 2017-2018 counts overlap N'Cho; primary campaign-response series is the subsequent 2018-2019 season",
          denominator_compatibility="No blood-culture-tested denominator for either season; cannot supply a PPV observation",
          verification_bias_risk="unassessable",
          source_status="primary Poncin 2022 full text reviewed (Zotero P3VSMY2P) and reconciled against N'Cho 2019",
          simple_model_status="not eligible: no blood-specific tested denominator",
          verification_model_status="not a PPV observation; historical overlap with N'Cho explicitly linked",
          recommended_use="exclude from PPV likelihood; contextual/time-series source only",
          independence_group="Harare 2017-2018 N'Cho-Poncin overlap",
          audit_notes="Do not combine Poncin's 250 blood-or-stool confirmations among 4,330 cases with N'Cho 74/286. N'Cho is the only compatible tested/positive pair. Poncin's 1,967 cases are a later season and still lack a tested denominator.")

koh <- rbind(koh, koh_added[, audit_cols])

# Koh is preferred for its canonical 23-row outbreak audit; the no-time-series
# audit is next because it contains the latest primary-paper corrections;
# Appiah contributes studies not already present in either set.
combined <- rbind(
  cbind(discovery_source = "Koh 2025 review", koh[, audit_cols]),
  cbind(discovery_source = "Koh 2025 exclusion: no time-series data", excluded[, audit_cols]),
  cbind(discovery_source = "Appiah 2020 review", appiah[, audit_cols])
)
combined$study_key <- key(combined$study)

all_keys <- unique(combined$study_key)
rows <- vector("list", length(all_keys))

for (i in seq_along(all_keys)) {
  candidates <- combined[combined$study_key == all_keys[i], , drop = FALSE]
  chosen <- candidates[1, , drop = FALSE]

  # Fill any blanks from lower-priority versions of the same audited study.
  for (nm in audit_cols) {
    if (!nzchar(chosen[[nm]][1])) {
      value <- candidates[[nm]][nzchar(candidates[[nm]])]
      if (length(value)) chosen[[nm]][1] <- value[1]
    }
  }

  sources <- unique(candidates$discovery_source)
  chosen$source_koh2025 <- "Koh 2025 review" %in% sources
  chosen$source_koh_excluded_no_timeseries <- "Koh 2025 exclusion: no time-series data" %in% sources
  chosen$source_appiah2020 <- "Appiah 2020 review" %in% sources
  chosen$source_membership <- paste(sources, collapse = "; ")
  rows[[i]] <- chosen
}

out <- do.call(rbind, rows)

# ---------------------------------------------------------------------------
# Full-text reconstruction pass (2026-07-17). Each record below was audited
# against its primary paper to test for a blood-culture-specific
# suspected/tested/positive triple, replacing the earlier "pending source
# reconstruction" placeholder with a resolved status and reason. None yielded
# an admissible community-PPV pair. Qamar 2018 already carries a valid
# active-surveillance pair (170/70) and is routed to sensitivity-only with the
# resistance-only-numerator caveat; Al-Sanouri 2008 has a partial blood pair
# but only an inferred tested denominator, so it is documented as a flagged
# sensitivity candidate without a fabricated pair. Wording in the three
# role-driving fields deliberately avoids "reconstruction/required/pending" so
# each record resolves to a documented exclusion (or sensitivity) rather than
# the pending placeholder.
# ---------------------------------------------------------------------------
reconstruction_resolved <- list(
  "Holt 2011" = c(
    recommended_use = "exclude: trial-enriched estimand, not community PPV",
    verification_model_status = "treatment-trial confirmation rate; not a suspected-case PPV",
    simple_model_status = "not eligible: enrolment admitted culture-confirmed patients",
    audit_notes = "Source trial Dolecek 2008 (PLoS ONE 3:e2188) full text audited 2026-07-17: 358 randomised on clinical suspicion, 288 blood-or-marrow confirmed (287 per-protocol, incl. 5 S. Paratyphi A); eligibility explicitly admitted already-culture-confirmed patients, so the ~80 percent confirmation fraction is trial-enriched, not community PPV. Blood volume 5-8 mL adults / 3-5 mL children. The prior 287/358 summary figure is inadmissible."),
  "Scobie 2014" = c(
    recommended_use = "exclude: confirmed-only surveillance, no suspected-case tested denominator",
    verification_model_status = "laboratory-confirmed case series; not a PPV",
    simple_model_status = "not eligible: 1449 counts blood-positive cases, not people tested",
    audit_notes = "Full text audited 2026-07-17: 1582 confirmed (1560 with specimen data): 1449 blood, 108 stool, 3 other. 1449 is a positive count, not a tested denominator. Suspected-case data sit only in unavailable Supplemental Tables 2-3; the Table 3 blood-culture positivity uses a lab-throughput denominator, not a suspected frame."),
  "Qamar 2018" = c(
    recommended_use = "sensitivity only (flagged): active-surveillance pair with resistance-only numerator",
    verification_model_status = "sensitivity candidate; all-Typhi positive count not separable",
    simple_model_status = "eligible for sensitivity only: k counts ceftriaxone-resistant Typhi",
    audit_notes = "Full text audited 2026-07-17: nested active-surveillance arm gives 170 suspected / 170 blood-cultured / 70 positive (positivity 41 percent), a clean culture-independent frame. But the 70 counts ceftriaxone-resistant Typhi only (all-Typhi at least 70, unreported); 89 percent prior antibiotics, median 14-day illness, paediatric low blood volume (1.5-2.5 mL). Do not divide 486 or 70 by 586."),
  "Nimonkar 2022" = c(
    recommended_use = "exclude: blood-positive count unreported, serology-inclusive case set",
    verification_model_status = "probable plus confirmed case set; not a blood-culture pair",
    simple_model_status = "not eligible: 75 cultures taken, positive yield never stated",
    audit_notes = "Full text audited 2026-07-17: 75 cases = probable (Widal at least 1:80) plus confirmed; 75 blood cultures taken but the positive count is never reported; at least 6 cases diagnosed by TyphiDot serology. The 75/75 summary is a non-informative artifact, not positivity 1.0."),
  "Srinivasan 2022" = c(
    recommended_use = "exclude: confirmed-only, number blood-cultured not reported",
    verification_model_status = "cohort confirmed-case series; denominator absent",
    simple_model_status = "not eligible: k=51 firm but n tested is unreported",
    audit_notes = "Full text audited 2026-07-17: 51 blood-culture-confirmed pediatric cases (clean, blood-specific), but the number of fever at least 3-day children actually blood-cultured is never reported. The SEFI cohort (5705-6760) and area populations are attack-rate denominators, not the cultured-suspected frame."),
  "Al-Sanouri 2008" = c(
    recommended_use = "sensitivity only (flagged): partial blood pair, tested denominator inferred not reported",
    verification_model_status = "sensitivity candidate; n tested not directly reported",
    simple_model_status = "not eligible as a clean pair: hospital severity-selected, denominator inferred",
    audit_notes = "Full text audited 2026-07-17: outbreak N=83 clinically suspected (fever at least 39C plus abdominal pain plus diarrhoea, admitted), k=20 blood-culture-positive (blood-specific; the pooled 24 also includes 5 faecal). n tested not stated, only inferable near 83 from the authors 24/83=29 percent. Hospital severity-selected. Do not use the pooled 24/83."),
  "Hendriksen 2015" = c(
    recommended_use = "exclude: genomics/isolate-characterisation paper, no suspected denominator",
    verification_model_status = "stored-isolate WGS study; no testing flow",
    simple_model_status = "not eligible: no clinical suspected-case denominator",
    audit_notes = "Full text audited 2026-07-17: WGST of 94 stored isolates (86 blood, 8 stool) from a 2040-case surveillance total. Isolates are a storage-survival convenience sample, not a cultured-suspected denominator. No PPV recoverable."),
  "Polonsky 2014 (total)" = c(
    recommended_use = "exclude: pooled-specimen confirmation, no blood tested denominator",
    verification_model_status = "descriptive-epi report; specimen pooled",
    simple_model_status = "not eligible: 62 confirmed pool blood/bone marrow/bowel fluid/stool",
    audit_notes = "Full text audited 2026-07-17: 3795 suspected, 62 confirmed by isolation from blood/bone marrow/bowel fluid/stool (pooled); number cultured never reported. Same 2011-2012 Harare epidemic as Muti 2014 and Imanishi 2014 (non-independent)."),
  "Muti 2014" = c(
    recommended_use = "exclude: pooled-specimen confirmation (stool/urine/blood), no blood tested denominator",
    verification_model_status = "outbreak investigation; specimen pooled",
    simple_model_status = "not eligible: 24 positives pooled across stool/urine/blood",
    audit_notes = "Full text audited 2026-07-17: 1078 suspected, 24 positive for S. Typhi in stool, urine or blood (pooled, stool-dominant); number blood-cultured never reported. Same 2011-2012 Harare epidemic as Polonsky and Imanishi (non-independent)."),
  "Imanishi 2014 (total)" = c(
    recommended_use = "exclude: WASH-survey paper, pooled specimen, time-varying sampling",
    verification_model_status = "household WASH survey; not the case-count source",
    simple_model_status = "not eligible: pooled confirmation, no fixed tested denominator",
    audit_notes = "Full text audited 2026-07-17: 4181 suspected / 52 confirmed by blood, stool or urine culture (pooled); culturing dropped from all high-risk suspects to 1-in-5 then 1-in-10, so no fixed denominator. Same 2011-2012 Harare epidemic as Polonsky and Muti (non-independent)."),
  "Burnsed 2018" = c(
    recommended_use = "exclude: WGS/descriptive, contact-defined probables, no cultured-suspected denominator",
    verification_model_status = "molecular plus descriptive; no testing flow",
    simple_model_status = "not eligible: probable cases are epi-linked contacts, not a cultured frame",
    audit_notes = "Full text audited 2026-07-17: 38 cases (25 culture-confirmed: 12 blood, 9 stool, 4 both; 13 probable contacts). US point-source Marshallese-community outbreak; confirmation required a matching PFGE pattern and there is no cultured-suspected denominator."),
  "Davis 2018" = c(
    recommended_use = "exclude: pooled blood/stool confirmation, no blood-culture tested denominator",
    verification_model_status = "short-form surveillance report; specimen pooled",
    simple_model_status = "not eligible: 80 confirmed pool blood/stool; no tested count",
    audit_notes = "Correct primary audited 2026-07-17 (Zotero EDH27L72, MMWR 67:342-3): separate 2016-2017 Harare outbreak; 780 suspected, 80 confirmed by blood or stool; the only culture number is 73 isolates AST-tested from blood and stool. No blood-specific tested denominator. NB the folder CTL54YK3 mislabels an N'Cho 2019 copy as Davis 2018."),
  "Hechaichi 2023" = c(
    recommended_use = "exclude: pooled-specimen confirmation, no blood tested denominator",
    verification_model_status = "outbreak investigation; specimen pooled",
    simple_model_status = "not eligible: 74 confirmed pool blood/bone marrow/bowel fluid/stool",
    audit_notes = "Full text audited 2026-07-17: 187 met suspected criteria (74 confirmed plus 28 probable plus 85 suspected); confirmation pooled blood/bone marrow/bowel fluid/stool, not disaggregated; culture-negative suspects reclassified as healthy. No blood tested denominator."),
  "Bano-Zaidi 2018" = c(
    recommended_use = "exclude: severity-selected referral cohort, blood-tested denominator unreported",
    verification_model_status = "referral-hospital inpatient series; no tested count",
    simple_model_status = "not eligible: k=46 blood-positive but n tested unreported",
    audit_notes = "Full text audited 2026-07-17: 110 hospitalised (34 percent with severe complications, 81 percent ill at least 2 weeks); confirmed by blood (46), bone marrow (6) or stool (4), plus 10 clinical-perforation and 44 Widal over 1:320 probables. k=46 is blood-specific but the blood-cultured denominator is not reported; severity/hospital-selected frame.")
)
for (nm in names(reconstruction_resolved)) {
  i <- which(key(out$study) == key(nm))
  if (!length(i)) next
  vals <- reconstruction_resolved[[nm]]
  for (f in names(vals)) out[[f]][i] <- vals[[f]]
}

out$blood_volume_ml <- NA_real_
out$blood_volume_source <- "not reported/retrieved"
if ("blood_volume_ml" %in% names(live_ppv)) {
  j <- match(key(out$study), key(live_ppv$study))
  live_volume <- suppressWarnings(as.numeric(live_ppv$blood_volume_ml[j]))
  have_volume <- !is.na(j) & !is.na(live_volume) & live_volume > 0
  out$blood_volume_ml[have_volume] <- live_volume[have_volume]
  out$blood_volume_source[have_volume] <- "active community PPV input; primary-source value"
}

# Main-model inclusion requires a usable blood-specific observation and either
# current stage-1 status, complete verification, or an explicit verification-
# adjusted candidacy assessment. Sensitivity-only observations remain excluded
# from the main model but retain their role below.
recommendation <- paste(out$recommended_use, out$verification_model_status,
                        out$simple_model_status, sep = " | ")
main_candidate <- out$current_stage1 == "TRUE" |
  grepl("base candidate|primary new candidate|verification-adjusted candidate",
        recommendation, ignore.case = TRUE)

has_pair <- nzchar(out$audited_blood_tested) & nzchar(out$audited_blood_positive)
tested <- suppressWarnings(as.numeric(out$audited_blood_tested))
positive <- suppressWarnings(as.numeric(out$audited_blood_positive))
valid_pair <- has_pair & !is.na(tested) & !is.na(positive) & tested > 0 &
  positive >= 0 & positive <= tested

out$observed_positivity <- ifelse(valid_pair, positive / tested, NA_real_)
if (!requireNamespace("posterior", quietly = TRUE)) stop("Package 'posterior' is required")
se_draws <- posterior::as_draws_matrix(readRDS(se_draws_rds))
se5 <- se_draws[, "Se_BC_5mL"]
se10 <- se_draws[, "Se_BC_10mL"]
se_bc_median_5ml <- unname(stats::median(se5))
se_bc_lower90_5ml <- unname(stats::quantile(se5, 0.05))
se_bc_upper90_5ml <- unname(stats::quantile(se5, 0.95))
out$se_bc_reference_median_5ml <- se_bc_median_5ml
out$se_bc_reference_90cri <- sprintf("%.3f-%.3f", se_bc_lower90_5ml, se_bc_upper90_5ml)
posterior_predictive_upper_tail <- function(y, n, se) {
  mean(stats::pbinom(y - 1, size = n, prob = se, lower.tail = FALSE))
}
out$pp_max_pi_upper_tail_5ml <- NA_real_
out$pp_max_pi_upper_tail_10ml <- NA_real_
out$pp_max_pi_upper_tail_study_volume <- NA_real_
for (i in which(valid_pair)) {
  out$pp_max_pi_upper_tail_5ml[i] <- posterior_predictive_upper_tail(positive[i], tested[i], se5)
  out$pp_max_pi_upper_tail_10ml[i] <- posterior_predictive_upper_tail(positive[i], tested[i], se10)
  if (!is.na(out$blood_volume_ml[i]) && out$blood_volume_ml[i] > 0) {
    se_study <- plogis(se_draws[, "alpha0"] + se_draws[, "alpha1"] * log(out$blood_volume_ml[i]))
    out$pp_max_pi_upper_tail_study_volume[i] <-
      posterior_predictive_upper_tail(positive[i], tested[i], se_study)
  }
}
out$ppv_likelihood_admissibility <- ifelse(
  !valid_pair, "not assessable: no valid blood-culture pair",
  ifelse(out$pp_max_pi_upper_tail_5ml < 0.05,
         "posterior-predictive conflict at 5 mL reference (upper-tail p<0.05)",
         ifelse(out$pp_max_pi_upper_tail_5ml < 0.10,
                "mild posterior-predictive tension at 5 mL reference (0.05<=p<0.10)",
                "no posterior-predictive conflict at maximum pi"))
)

# Correct the n>=20 inversion: small n controls precision, not validity.
swad <- key(out$study) == "swaddiwudhipong 2001"
main_candidate[swad & valid_pair] <- TRUE
out$recommended_use[swad] <- "main candidate: low-bias near-complete verification; small n handled by likelihood"
out$verification_model_status[swad] <- "candidate with negligible verification adjustment"
out$simple_model_status[swad] <- "eligible; small n implies uncertainty, not exclusion"

# Remove the obsolete hard sample-size rationale wherever it survives from an
# earlier audit. This does not automatically promote a study: setting,
# specimen compatibility, and selection can still restrict it to sensitivity.
small_n_legacy <- grepl("fails n>=20|n<20|too few tested", out$recommended_use,
                        ignore.case = TRUE)
out$recommended_use[small_n_legacy & grepl("Kato", out$study, ignore.case = TRUE)] <-
  "sensitivity only: returning-traveler point-source estimand; small n handled by likelihood"

out$include_for_main_ppv_model <- main_candidate & valid_pair
recommendation <- paste(out$recommended_use, out$verification_model_status,
                        out$simple_model_status, sep = " | ")
out$ppv_model_role <- ifelse(
  out$include_for_main_ppv_model,
  ifelse(grepl("complete verification|base candidate", recommendation, ignore.case = TRUE),
         "main model: directly observed PPV candidate",
         "main model: verification-adjusted PPV candidate"),
  ifelse(grepl("sensitivity", recommendation, ignore.case = TRUE) & valid_pair,
         "sensitivity analysis only",
         ifelse(grepl("reconstruction|required|pending", recommendation, ignore.case = TRUE),
                "not yet: source reconstruction or retrieval required", "exclude from PPV modeling"))
)

# Separate the magnitude of verification correction from inclusion status.
out$verification_adjustment_grade <- ifelse(
  !out$include_for_main_ppv_model, "not applicable",
  ifelse(grepl("moderate", out$verification_bias_risk, ignore.case = TRUE), "moderate",
         ifelse((!is.na(suppressWarnings(as.numeric(out$tested_fraction))) &
                   abs(suppressWarnings(as.numeric(out$tested_fraction)) - 1) < 1e-8) |
                  grepl("(^|[^a-z])low([^a-z]|$)", out$verification_bias_risk,
                        ignore.case = TRUE), "negligible/low", "major"))
)
out$ppv_model_role[out$include_for_main_ppv_model &
                     out$verification_adjustment_grade == "negligible/low"] <-
  "main model: negligible/low verification adjustment"
out$ppv_model_role[out$include_for_main_ppv_model &
                     out$verification_adjustment_grade == "moderate"] <-
  "main model: moderate verification adjustment"
out$ppv_model_role[out$include_for_main_ppv_model &
                     out$verification_adjustment_grade == "major"] <-
  "main model: major verification adjustment"

reason_part <- function(label, value) ifelse(nzchar(value), paste0(label, value), "")
justification <- mapply(
  function(denom, bias, simple, verify, use, notes) {
    parts <- c(reason_part("Denominator: ", denom),
               reason_part("Verification bias: ", bias),
               reason_part("Simple-model assessment: ", simple),
               reason_part("Verification-model assessment: ", verify),
               reason_part("Recommended use: ", use),
               reason_part("Notes: ", notes))
    paste(parts[nzchar(parts)], collapse = " | ")
  },
  out$denominator_compatibility, out$verification_bias_risk,
  out$simple_model_status, out$verification_model_status,
  out$recommended_use, out$audit_notes,
  USE.NAMES = FALSE
)
out$modeling_justification <- justification

count_discrepancy <- function(summary_value, audited_value) {
  s <- suppressWarnings(as.numeric(summary_value))
  a <- suppressWarnings(as.numeric(audited_value))
  !is.na(s) & !is.na(a) & s != a
}
out$superseded_summary_count_discrepancy <-
  count_discrepancy(out$summary_blood_tested, out$audited_blood_tested) |
  count_discrepancy(out$summary_blood_positive, out$audited_blood_positive)
out$modeling_count_fields <- "audited_blood_tested; audited_blood_positive"
names(out)[names(out) == "summary_suspected"] <- "superseded_summary_suspected"
names(out)[names(out) == "summary_blood_tested"] <- "superseded_summary_blood_tested"
names(out)[names(out) == "summary_blood_positive"] <- "superseded_summary_blood_positive"
names(out)[names(out) == "summary_positivity"] <- "superseded_summary_positivity"

# Remove invalid legacy byte sequences inherited from the canonical summary so
# the consolidated CSV is valid UTF-8 on Windows.
for (nm in names(out)) {
  if (is.character(out[[nm]])) out[[nm]] <- iconv(out[[nm]], from = "", to = "UTF-8", sub = "")
}

front <- c("study", "location", "outbreak_year", "source_koh2025",
           "source_koh_excluded_no_timeseries", "source_appiah2020",
           "source_membership", "include_for_main_ppv_model", "ppv_model_role",
           "verification_adjustment_grade", "observed_positivity",
           "blood_volume_ml", "blood_volume_source",
           "ppv_likelihood_admissibility", "se_bc_reference_median_5ml",
           "se_bc_reference_90cri", "pp_max_pi_upper_tail_5ml",
           "pp_max_pi_upper_tail_10ml", "pp_max_pi_upper_tail_study_volume",
           "modeling_justification")
drop <- c("discovery_source", "study_key")
out <- out[, c(front, setdiff(names(out), c(front, drop))), drop = FALSE]

write.csv(out, output_csv, row.names = FALSE, na = "", fileEncoding = "UTF-8")

cat(sprintf("Wrote %d unique audited studies to %s\n", nrow(out), output_csv))
cat("Source membership counts:\n")
print(c(Koh_2025 = sum(out$source_koh2025),
        Koh_excluded = sum(out$source_koh_excluded_no_timeseries),
        Appiah_2020 = sum(out$source_appiah2020)))
cat("Model-role counts:\n")
print(sort(table(out$ppv_model_role), decreasing = TRUE))
