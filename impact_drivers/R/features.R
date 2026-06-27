# =============================================================================
# Prospective predictors of ORI impact, computed AT the time of vaccination.
#
# Every feature here is knowable from data available by decision week tau (the
# first tau weeks of the curve + static context) -- NO retrospective quantities
# (final size, peak, duration) are used as predictors. R_t comes from the exact
# reconstruction (renewal_core); the early growth rate is a fixed per-outbreak
# descriptor (log-linear slope over the first k weeks).
# =============================================================================

# Early exponential growth rate: slope of log(incidence) over the first k weeks.
early_growth_rate <- function(inc, k) {
  kk <- min(k, length(inc))
  if (kk < 3) return(NA_real_)
  wk <- seq_len(kk); y <- log(inc[seq_len(kk)] + 0.5)
  tryCatch(unname(coef(stats::lm(y ~ wk))[2]), error = function(e) NA_real_)
}

# Build the per-(outbreak, delay) prospective feature table.
compute_features <- function(prep, rcfg, cfg) {
  w <- gi_from_config(rcfg)
  k <- cfg$growth_window_weeks
  rows <- list()
  for (sid in names(prep$series)) {
    inc <- prep$series[[sid]]; Tn <- length(inc); m <- prep$meta[[sid]]
    Rt <- reconstruct_Rt(inc, w)$Rt
    r_early <- early_growth_rate(inc, k)
    for (tau in cfg$delays_weeks) {
      tt <- min(tau, Tn)                              # decision-time index (prospective)
      rt_pre <- Rt[seq_len(tt)]
      rows[[length(rows) + 1]] <- data.frame(
        study_id = sid, tau = tau,
        # --- dynamic predictors (known by week tau) ---
        r_early = r_early,
        doubling_time = if (is.finite(r_early) && r_early > 0) log(2) / r_early else NA_real_,
        Rt_at_tau = Rt[tt],
        Rt_mean_pre = mean(rt_pre, na.rm = TRUE),
        Rt_slope_pre = if (sum(is.finite(rt_pre)) >= 3)
          tryCatch(unname(coef(stats::lm(rt_pre ~ seq_along(rt_pre)))[2]), error = function(e) NA_real_) else NA_real_,
        phase_growing = as.integer(isTRUE(Rt[tt] > 1)),
        cum_cases_by_tau = sum(inc[seq_len(tt)], na.rm = TRUE),
        cum_AR_by_tau = 100 * sum(inc[seq_len(tt)], na.rm = TRUE) / m$population,
        incidence_at_tau = inc[tt],
        # --- static context ---
        population = m$population, log_pop = log10(m$population),
        who_region = m$who_region, amr_status = m$amr_status,
        time_unit = m$time_unit, year = m$year,
        # --- retrospective descriptors (NOT predictors; for validation only) ---
        tot_cases = m$tot_cases, duration_weeks = Tn,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
