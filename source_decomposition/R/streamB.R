# =============================================================================
# Stream B — data-driven source fraction theta_data from the incidence series.
#
# IDENTIFIABILITY CAVEAT: a single incidence curve cannot definitively separate
# "exogenous common-source pulse convolved with the incubation period" from
# "propagated transmission" -- they can generate identical curves (the same
# non-identifiability the manuscript notes). Stream B therefore yields SUGGESTIVE
# diagnostics, not a definitive theta; it is cross-validated against the paper
# evidence (Stream A).
#
# Common-source fit: incidence ~ (short exogenous source) convolved with the
# incubation-period distribution. A good low-width fit + no multi-generational
# growth / resurgence -> common-source (theta -> 1). Sustained growth over many
# incubation widths and/or multiple waves -> propagated component (theta down).
# =============================================================================

# Weekly incubation-period pmf over lags 0,1,2,... (gamma, includes first week).
incubation_pmf <- function(mean_days, sd_days, step_days = 7, max_lag = 12) {
  shape <- (mean_days / sd_days)^2; rate <- mean_days / sd_days^2
  edges <- (0:(max_lag + 1)) * step_days
  m <- diff(pgamma(edges, shape = shape, rate = rate))
  m / sum(m)
}

# Best common-source fit: a rectangular source of width 1..max_width starting at
# any week, convolved with the incubation pmf, scaled by least squares. Returns
# the best R^2 (how well a short exogenous source explains the whole curve).
common_source_fit <- function(inc, f_inc, max_width = 6) {
  Tn <- length(inc); L <- length(f_inc)
  conv <- function(src) {
    out <- numeric(Tn)
    for (tau in which(src > 0)) for (d in 0:(L - 1)) {
      t <- tau + d; if (t <= Tn) out[t] <- out[t] + src[tau] * f_inc[d + 1]
    }
    out
  }
  sstot <- sum((inc - mean(inc))^2); best <- -Inf; best_par <- c(NA, NA)
  for (t0 in seq_len(Tn)) for (wdt in seq_len(max_width)) {
    src <- numeric(Tn); idx <- t0:min(t0 + wdt - 1, Tn); src[idx] <- 1
    shape <- conv(src); den <- sum(shape^2); if (den <= 0) next
    A <- sum(inc * shape) / den; pred <- A * shape
    r2 <- 1 - sum((inc - pred)^2) / sstot
    if (is.finite(r2) && r2 > best) { best <- r2; best_par <- c(t0, wdt) }
  }
  list(r2 = max(best, 0), t0 = best_par[1], width = best_par[2])
}

# Propagation diagnostics + theta_data for one outbreak.
theta_data_one <- function(inc, w_gi, cfg, mean_days) {
  Tn <- length(inc)
  f_inc <- incubation_pmf(mean_days, cfg$incubation$sd_days)
  inc_sd_weeks <- cfg$incubation$sd_days / 7
  Rt <- reconstruct_Rt(inc, w_gi)$Rt

  peak <- which.max(inc)
  growth_weeks <- max(peak - 1, 0)
  growth_ratio <- growth_weeks / max(inc_sd_weeks, 0.5)        # rise length vs incubation SD
  # epidemic waves: local maxima above 30% of the peak, separated by >= GI (~2 wk)
  thr <- 0.3 * max(inc)
  loc_max <- which(c(FALSE, diff(sign(diff(inc))) < 0, FALSE) & inc > thr)
  n_waves <- if (length(loc_max) == 0) 1L else {
    keep <- loc_max[c(TRUE, diff(loc_max) >= 2)]; max(length(keep), 1L) }
  # resurgence: upward crossings of R_t = 1 after week 2 beyond the first
  rtv <- Rt[is.finite(Rt)]; up <- sum(diff(rtv > 1) == 1)
  n_resurg <- max(up - 1L, 0L)
  cs <- common_source_fit(inc, f_inc)

  # Propagation score from the two most principled signals: R_t resurgence
  # (upward re-crossings of 1 beyond the first) and poor common-source fit. The
  # growth-duration and raw wave count are reported as diagnostics only -- a
  # multi-week rise occurs in BOTH modes, so it is too weak to drive theta (the
  # identifiability limit). theta_data is a SUGGESTIVE cross-check, not definitive.
  prop <- min(0.7, 0.12 * n_resurg) + 0.5 * (1 - cs$r2)
  theta <- max(0.05, min(0.95, 1 - prop))
  list(theta = theta, cs_r2 = cs$r2, growth_weeks = growth_weeks,
       growth_ratio = growth_ratio, n_waves = n_waves, n_resurg = n_resurg,
       peak_week = peak)
}

# theta_data for all outbreaks, with an interval from the incubation-mean range.
stream_B <- function(prep, rcfg, cfg) {
  w_gi <- gi_from_config(rcfg)
  means <- c(cfg$incubation$mean_days, cfg$incubation$mean_days_range)
  do.call(rbind, lapply(names(prep$series), function(sid) {
    inc <- prep$series[[sid]]
    base <- theta_data_one(inc, w_gi, cfg, cfg$incubation$mean_days)
    ths <- vapply(means, function(md) theta_data_one(inc, w_gi, cfg, md)$theta, numeric(1))
    data.frame(study_id = sid,
               theta_data = round(base$theta, 2),
               theta_data_lo = round(min(ths), 2), theta_data_hi = round(max(ths), 2),
               cs_r2 = round(base$cs_r2, 2), growth_weeks = base$growth_weeks,
               growth_ratio = round(base$growth_ratio, 2), n_waves = base$n_waves,
               n_resurg = base$n_resurg,
               mode_data = ifelse(base$theta >= 0.7, "A_common_source",
                           ifelse(base$theta >= 0.3, "C_mixed", "B_propagated")),
               stringsAsFactors = FALSE)
  }))
}
