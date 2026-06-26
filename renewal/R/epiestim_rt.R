# =============================================================================
# EpiEstim R_t cross-check — for the R_t PLOTS and credible intervals only.
# The counterfactual point estimates use the exact self-consistent
# reconstruction (renewal_core.R); EpiEstim's smoothed R_t is NOT used for
# propagation (smoothing would break self-consistency).
# =============================================================================

# Cori (instantaneous) R_t with a sliding weekly window. GI mean/sd in WEEKS.
# Returns data.frame(week, mean_r, lo, hi) aligned to the window end, or NULL.
epiestim_rt <- function(incidence, cfg, window = 3L) {
  if (!requireNamespace("EpiEstim", quietly = TRUE)) return(NULL)
  inc <- pmax(round(as.numeric(incidence)), 0)
  Tn <- length(inc)
  if (Tn < window + 2L || sum(inc) < 5) return(NULL)
  mean_si <- cfg$gi$mean_days / cfg$step_days          # weeks
  std_si  <- (cfg$gi$mean_days * cfg$gi$cv) / cfg$step_days
  t_start <- seq(2L, Tn - window + 1L)
  t_end   <- t_start + window - 1L
  res <- tryCatch(
    EpiEstim::estimate_R(
      incid = inc,
      method = "parametric_si",
      config = EpiEstim::make_config(list(mean_si = mean_si, std_si = std_si,
                                          t_start = t_start, t_end = t_end))),
    error = function(e) NULL)
  if (is.null(res)) return(NULL)
  r <- res$R
  data.frame(week = r$t_end,
             mean_r = r$`Mean(R)`,
             lo = r$`Quantile.0.025(R)`,
             hi = r$`Quantile.0.975(R)`,
             stringsAsFactors = FALSE)
}
