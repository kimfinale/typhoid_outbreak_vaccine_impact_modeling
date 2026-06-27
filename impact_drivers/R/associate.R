# =============================================================================
# Phase 2 — empirical association (the real outbreaks; exploratory, small n).
# Single-predictor Spearman screening with bootstrap CIs, ranked per outcome,
# plus a mechanistic mediation check (does R_t-at-vaccination explain impact,
# and do contextual variables add anything beyond it?).
# =============================================================================

suppressMessages({ library(dplyr) })

# Spearman rho + percentile bootstrap CI for one predictor/outcome pair.
.boot_spearman <- function(x, y, B = 2000) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 5 || sd(x) == 0 || sd(y) == 0) return(c(rho = NA, lo = NA, hi = NA, n = n))
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  bs <- replicate(B, { i <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman")) })
  c(rho = rho, lo = unname(quantile(bs, 0.025, na.rm = TRUE)),
    hi = unname(quantile(bs, 0.975, na.rm = TRUE)), n = n)
}

# Screen all predictors against all outcomes at a fixed delay (default tau = 8).
associate_screen <- function(dat, predictors, outcomes, tau = 8) {
  d <- dat[dat$tau == tau, ]
  out <- list()
  for (o in outcomes) for (p in predictors) {
    s <- .boot_spearman(d[[p]], d[[o]])
    out[[length(out) + 1]] <- data.frame(
      outcome = o, predictor = p, tau = tau,
      rho = unname(s["rho"]), ci_lo = unname(s["lo"]), ci_hi = unname(s["hi"]),
      n = unname(s["n"]), stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  res %>% group_by(outcome) %>% arrange(desc(abs(rho)), .by_group = TRUE) %>% ungroup()
}

# How the leading predictors' associations change across delay (R_t-at-tau is
# time-varying, so its predictive value should shift with tau).
associate_by_delay <- function(dat, predictors, outcome, taus) {
  do.call(rbind, lapply(taus, function(tau) {
    d <- dat[dat$tau == tau, ]
    do.call(rbind, lapply(predictors, function(p) {
      s <- .boot_spearman(d[[p]], d[[outcome]], B = 500)
      data.frame(outcome = outcome, predictor = p, tau = tau,
                 rho = unname(s["rho"]), stringsAsFactors = FALSE)
    }))
  }))
}

# Mechanistic mediation: incremental variance explained beyond R_t-at-vaccination.
# (1) outcome ~ R_t_at_tau ; (2) + log_pop ; (3) + cumulative AR. Reports R^2 and
# the increment, i.e. whether context adds signal once R_t is known.
mediation_check <- function(dat, outcome, tau = 8) {
  d <- dat[dat$tau == tau, ]
  d <- d[is.finite(d[[outcome]]) & is.finite(d$Rt_at_tau), ]
  r2 <- function(form) tryCatch(summary(lm(form, data = d))$r.squared, error = function(e) NA)
  m1 <- r2(reformulate("Rt_at_tau", outcome))
  m2 <- r2(reformulate(c("Rt_at_tau", "log_pop"), outcome))
  m3 <- r2(reformulate(c("Rt_at_tau", "log_pop", "cum_AR_by_tau"), outcome))
  data.frame(outcome = outcome, tau = tau, n = nrow(d),
             R2_Rt = round(m1, 3), R2_Rt_pop = round(m2, 3), R2_Rt_pop_AR = round(m3, 3),
             incremental_pop = round(m2 - m1, 3), incremental_AR = round(m3 - m2, 3),
             stringsAsFactors = FALSE)
}
