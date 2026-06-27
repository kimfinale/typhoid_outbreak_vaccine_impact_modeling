# =============================================================================
# Monte Carlo propagation of GI / asymptomatic / PPV uncertainty through the
# exact R_t reconstruction and the vaccine counterfactual. Reuses renewal_core.
#
# Per draw and outbreak: build the mixture GI kernel, reshape incidence by the
# (possibly time-varying) PPV, reconstruct R_t exactly, derive psi_T from the
# asymptomatic split, propagate the renewal counterfactual, record impact.
# =============================================================================

# t_eff (first protected week) for the base scenario, deterministic.
base_t_eff <- function(rcfg, cfg) {
  round(cfg$scenario$tau_weeks + rcfg$timing$immunity_onset_days / rcfg$step_days +
        round(rcfg$timing$campaign_duration_days / rcfg$step_days / 2))
}

run_mc <- function(prep, rcfg, cfg, params) {
  te <- base_t_eff(rcfg, cfg); cov <- cfg$scenario$coverage
  sids <- names(prep$series); nd <- nrow(params)
  rt_mats <- setNames(lapply(sids, function(s) {
    matrix(NA_real_, nrow = length(prep$series[[s]]), ncol = nd) }), sids)
  imp <- vector("list", nd)

  for (i in seq_len(nd)) {
    p <- params[i, ]
    w <- discretize_gi_mix(p$gi_mean, p$gi_cv, p$pi_c, p$c_c, p$mu_carrier,
                           step_days = rcfg$step_days, drop_first = rcfg$gi$drop_first_week)
    pT <- psi_T_from_asymptomatic(p$psi, p$f_a, p$c_a)
    rows <- lapply(sids, function(s) {
      inc <- prep$series[[s]]; Tn <- length(inc)
      icorr <- ppv_reshape(inc, p$ppv_trend)
      Rt <- reconstruct_Rt(icorr, w)$Rt
      rt_mats[[s]][, i] <<- pmin(Rt, 8)
      occurred <- te <= (Tn - rcfg$tail_guard_weeks)
      if (!occurred) { rho <- 0; av <- 0 } else {
        cf <- renewal_counterfactual(icorr, w, tau = cfg$scenario$tau_weeks, t_eff = te,
                                     pi = cov, psi_T = pT$psi_T,
                                     c_shape = rcfg$timing$c_shape, feedback = TRUE)
        im <- impact_measures(icorr, cf$incidence_v)
        rho <- im$pct_reduction
        av <- (rho / 100) * sum(inc)            # cases averted on observed (suspected) scale
      }
      data.frame(draw = i, study_id = s, rho = rho, averted = av,
                 psi_T = pT$psi_T, phi_a = pT$phi_a, stringsAsFactors = FALSE)
    })
    imp[[i]] <- do.call(rbind, rows)
  }
  list(impact = do.call(rbind, imp), rt = rt_mats, t_eff = te,
       params = params, sids = sids)
}

# Per-outbreak R_t quantile ribbon from the MC draws.
rt_ribbon <- function(mc) {
  do.call(rbind, lapply(mc$sids, function(s) {
    M <- mc$rt[[s]]
    data.frame(study_id = s, week = seq_len(nrow(M)),
               rt_median = apply(M, 1, median, na.rm = TRUE),
               rt_lo = apply(M, 1, quantile, 0.025, na.rm = TRUE),
               rt_hi = apply(M, 1, quantile, 0.975, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
}

# Impact summary: per outbreak and pooled (median + 95% credible interval).
summarise_mc <- function(mc) {
  q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
  per <- mc$impact %>% dplyr::group_by(study_id) %>%
    dplyr::summarise(rho_median = median(rho), rho_lo = q(rho, 0.025), rho_hi = q(rho, 0.975),
                     averted_median = median(averted), averted_lo = q(averted, 0.025),
                     averted_hi = q(averted, 0.975), .groups = "drop")
  pooled_by_draw <- mc$impact %>% dplyr::group_by(draw) %>%
    dplyr::summarise(rho = median(rho), averted = sum(averted), .groups = "drop")
  pooled <- data.frame(
    rho_median = median(pooled_by_draw$rho),
    rho_lo = q(pooled_by_draw$rho, 0.025), rho_hi = q(pooled_by_draw$rho, 0.975),
    averted_median = median(pooled_by_draw$averted),
    averted_lo = q(pooled_by_draw$averted, 0.025), averted_hi = q(pooled_by_draw$averted, 0.975))
  list(per_outbreak = per, pooled = pooled, pooled_by_draw = pooled_by_draw)
}

# Partial rank correlation coefficients: influence of each input parameter on the
# pooled percent reduction (global sensitivity ranking).
prcc <- function(params, y, inputs = c("gi_mean","gi_cv","pi_c","c_c","mu_carrier",
                                       "f_a","c_a","ppv_trend","psi")) {
  R <- as.data.frame(lapply(params[inputs], rank)); R$y <- rank(y)
  sapply(inputs, function(v) {
    others <- setdiff(inputs, v)
    rx <- residuals(lm(reformulate(others, response = v), data = R))
    ry <- residuals(lm(reformulate(others, response = "y"), data = R))
    cor(rx, ry)
  })
}
