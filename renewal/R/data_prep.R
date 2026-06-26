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

load_inputs <- function(cfg) {
  ts  <- read.csv(cfg$paths$ts,      stringsAsFactors = FALSE, check.names = FALSE)
  dat <- read.csv(cfg$paths$summary, stringsAsFactors = FALSE, check.names = FALSE)
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
      meta[[sid]] <- list(
        study_id = sid, time_unit = tu,
        country_iso = dat$country_iso[dat$study_id == sid][1],
        population = suppressWarnings(as.numeric(dat$population[dat$study_id == sid][1])),
        who_region = dat$WHO_region[dat$study_id == sid][1],
        amr_status = dat$amr_status[dat$study_id == sid][1],
        intervention_week = .intervention_week(dat[dat$study_id == sid, ], cfg)
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
