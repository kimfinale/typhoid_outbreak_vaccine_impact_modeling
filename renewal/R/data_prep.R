# =============================================================================
# Data preparation for the renewal analysis.
#
# Reuses the existing pipeline's conventions so static and renewal share inputs:
#   * incidence = `total_cases` RELABELED `suspected_cases` (as in R/1-data.R,
#     where raw `suspected_cases` is frequently NA);
#   * the same 24-study eligibility exclusion list as R/1-data.R (-> 19 outbreaks).
#
# Renewal runs only on outbreaks whose GI spans several reporting intervals
# (resolution rule mu_g/Delta >= 2): daily (aggregated to weekly) + weekly = 13.
#
# Base R loaders only (data.table::fread segfaults under Rscript in this env).
# Pure-ish: reads files named in cfg, returns data; no hidden global mutation.
# =============================================================================

# 24-study eligibility exclusion (verbatim from R/1-data.R:114-121).
RENEWAL_EXCLUSIONS <- c(
  "Cherian 2015", "Wang 2022 (2)", "Imanishi 2014 (total)", "Polonsky 2014 (total)",
  "Limpitikul 2014", "Limpitikul 2014 (1)", "Limpitikul 2014 (2)", "Poncin 2022",
  "Scobie 2014", "Muehlen 2007", "Bayram 2011", "Keddy 2011", "Michel 2005",
  "Burnsed 2018", "Hu 2022", "Bano-Zaidi 2018", "Wang 2022 (1)", "Makungo 2020",
  "Al-Sanouri 2008", "Hechaichi 2023", "Holt 2011", "Roy 2016", "Nimonkar 2022",
  "Srinivasan 2022"
)

# Analysis time step (weeks) by native reporting resolution. Daily series are
# aggregated up to weekly, so their analysis step is 1 week; fortnightly/monthly
# cannot be refined below their native bin.
STEP_WEEKS_BY_UNIT <- c(daily = 1, weekly = 1, fortnightly = 2, monthly = 30 / 7)

# --- Summary normalizer ------------------------------------------------------
# The canonical summary (Typhoid_Outbreak_Time_Series_2000_2022_Summary.csv) uses
# human-readable headers and its own study-id naming. read_summary() maps both to
# the snake_case columns and the ts_v2 study-id names the pipeline expects, so all
# downstream code is unchanged. (The legacy typhoid_outbreak_summary_v1.csv passes
# through untouched.)
SUMMARY_COLMAP <- c(
  "Study ID" = "study_id", "Outbreak start year" = "outbreak_start_year",
  "Outbreak end year" = "outbreak_end_year", "Location" = "location",
  "Country 3-letter ISO" = "country_iso", "WHO region" = "WHO_region",
  "AMR status" = "amr_status", "Time unit" = "time_unit", "Start date" = "start_date",
  "End date" = "end_date", "Peak date" = "peak_date", "Intervention date" = "intervention_date",
  "Total suspected cases" = "tot_suspected_cases",
  "Lab tested cases by blood/bone marrow culture" = "lab_tested_cases",
  "Total confirmed cases" = "tot_confirmed_cases",
  "Confirmed cases by blood/ bone marrow culture" = "lab_confirmed_cases",
  "Total probable cases" = "tot_probable_cases", "Number of hospitalized" = "hospitalized",
  "Number of complications" = "complicated", "Total deaths" = "tot_deaths",
  "Attack rate (%)" = "attack_rate", "CFR (%)" = "cfr", "Population" = "population")

# new-file study id -> ts_v2 canonical name (matched by cumulative cases / dates).
STUDY_ID_CROSSWALK <- c(
  "Polonsky 2014 (1)" = "Polonsky 2014, Dzivaresekwa",
  "Polonsky 2014 (2)" = "Polonsky 2014_Kuwadzana",
  "Imanishi 2014 (1)" = "Imanishi 2014, Dzivaresekwa",
  "Imanishi 2014 (2)" = "Imanishi 2014, Kuwadzana",
  "Walters 2014 (1)" = "Walters 2014, Kasese",
  "Walters 2014 (2)" = "Walters 2014, Bundibugyo",
  "Wang 2022" = "Wang 2022 (1)",          # excluded; map to an excluded canonical name
  "Limpitikul 2014 (total)" = "Limpitikul 2014")  # aggregate; excluded base name

read_summary <- function(path) {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("Study ID" %in% names(raw)) {                     # new canonical format
    raw <- raw[, names(raw) %in% names(SUMMARY_COLMAP), drop = FALSE]
    names(raw) <- SUMMARY_COLMAP[names(raw)]
    raw <- raw[!is.na(raw$study_id) & trimws(raw$study_id) != "", ]
    raw$study_id <- trimws(raw$study_id)
    hit <- raw$study_id %in% names(STUDY_ID_CROSSWALK)
    raw$study_id[hit] <- STUDY_ID_CROSSWALK[raw$study_id[hit]]
  }
  raw
}

# Time-series normalizer: the canonical Timeseries file uses human-readable
# headers but the same (ts_v2) study-id naming, so only column renaming is needed.
# "No. of Patients" is the count the static model uses (old total_cases).
TS_COLMAP <- c(
  "Study ID" = "study_id", "start date" = "start_date", "end date" = "end_date",
  "No. of Patients" = "total_cases", "Suspected" = "suspected_cases",
  "Confirmed" = "confirmed_cases", "Probable" = "probable_cases",
  "Date of onset (T/F)" = "date_onset", "Date of confirmation (T/F)" = "date_confirmation")

read_timeseries <- function(path) {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if ("Study ID" %in% names(raw)) {                     # new canonical format
    raw <- raw[, names(raw) %in% names(TS_COLMAP), drop = FALSE]
    names(raw) <- TS_COLMAP[names(raw)]
    raw <- raw[!is.na(raw$study_id) & trimws(raw$study_id) != "", ]
    raw$study_id <- trimws(raw$study_id)
  }
  raw
}

load_inputs <- function(cfg) {
  ts  <- read_timeseries(cfg$paths$ts)
  dat <- read_summary(cfg$paths$summary)
  # incidence series the static model uses: total_cases relabeled suspected_cases
  ts$incidence <- suppressWarnings(as.numeric(ts$total_cases))
  list(ts = ts, summary = dat)
}

# Aggregate a daily outbreak's rows to ISO weeks (Monday week start), summing
# incidence; returns an ordered numeric vector of weekly counts.
aggregate_daily_to_weekly <- function(df) {
  d <- as.Date(df$start_date)
  monday <- d - (as.integer(format(d, "%u")) - 1L)   # back up to Monday
  agg <- tapply(df$incidence, monday, sum, na.rm = TRUE)
  as.numeric(agg[order(as.Date(names(agg)))])
}

# Build the weekly incidence vector for one outbreak given its native unit.
outbreak_weekly_series <- function(df, time_unit) {
  df <- df[order(as.Date(df$start_date)), ]
  if (identical(time_unit, "daily")) aggregate_daily_to_weekly(df)
  else as.numeric(df$incidence)                       # weekly: already weekly
}

# Resolution diagnostic table for ALL eligible outbreaks (19), plus the included
# weekly incidence series for the 13 that pass the mu_g/Delta rule.
prep_outbreaks <- function(cfg, inputs = NULL) {
  if (is.null(inputs)) inputs <- load_inputs(cfg)
  ts  <- inputs$ts
  dat <- inputs$summary

  dat <- dat[!dat$study_id %in% RENEWAL_EXCLUSIONS, ]   # -> 19 outbreaks
  mu_g <- cfg$gi$mean_days
  sd_g <- mu_g * cfg$gi$cv
  min_ratio <- cfg$gi$resolution_min_ratio

  res <- data.frame(study_id = character(), time_unit = character(),
                    delta_weeks = numeric(), mu_g_over_delta = numeric(),
                    pct_gi_first_bin = numeric(), included = character(),
                    n_weeks = integer(), stringsAsFactors = FALSE)
  series <- list(); meta <- list()

  for (sid in dat$study_id) {
    tu <- dat$time_unit[dat$study_id == sid][1]
    step_w <- unname(STEP_WEEKS_BY_UNIT[tu])
    if (is.na(step_w)) step_w <- NA_real_
    ratio <- mu_g / (step_w * cfg$step_days)               # mu_g(days) / Delta(days)
    first_bin <- gi_first_bin_mass(mu_g, sd_g, step_days = step_w * cfg$step_days,
                                   max_steps = cfg$gi$max_steps)
    inc_rows <- ts[ts$study_id == sid, ]
    n_weeks <- NA_integer_
    incl <- !is.na(ratio) && ratio >= min_ratio
    if (incl && nrow(inc_rows) > 0) {
      vec <- outbreak_weekly_series(inc_rows, tu)
      n_weeks <- length(vec)
      series[[sid]] <- vec
      drow <- dat[dat$study_id == sid, ]
      meta[[sid]] <- list(
        study_id = sid, time_unit = tu,
        country_iso = trimws(drow$country_iso[1]),   # summary CSV has trailing spaces
        population = suppressWarnings(as.numeric(drow$population[1])),
        who_region = drow$WHO_region[1],
        amr_status = drow$amr_status[1],
        year = suppressWarnings(as.integer(substr(drow$start_date[1], 1, 4))),
        tot_deaths = suppressWarnings(as.numeric(drow$tot_deaths[1])),
        tot_cases = sum(vec, na.rm = TRUE),
        intervention_week = .intervention_week(drow, cfg)
      )
    }
    res <- rbind(res, data.frame(
      study_id = sid, time_unit = tu, delta_weeks = step_w,
      mu_g_over_delta = round(ratio, 3), pct_gi_first_bin = round(100 * first_bin, 1),
      included = ifelse(incl, "Y", "N"), n_weeks = n_weeks,
      stringsAsFactors = FALSE))
  }
  res <- res[order(res$included == "N", res$time_unit, res$study_id), ]
  list(resolution = res, series = series, meta = meta)
}

# Observed community-intervention timing (weeks from outbreak start), if present.
.intervention_week <- function(row, cfg) {
  idate <- row$intervention_date[1]; sdate <- row$start_date[1]
  if (is.null(idate) || is.na(idate) || idate %in% c("", "-")) return(NA_real_)
  as.numeric(as.Date(idate) - as.Date(sdate)) / cfg$step_days
}
