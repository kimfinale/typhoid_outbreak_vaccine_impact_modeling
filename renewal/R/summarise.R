# =============================================================================
# Summarisation: collapse Sobol draws to median + 95% uncertainty intervals,
# compute pooled and population-weighted estimates, and the eta_eff diagnostic.
# =============================================================================

suppressMessages({ library(dplyr) })

# Per-row percent reduction, then median + 95% UI by scenario cell.
summarise_draws <- function(df_cea) {
  df_cea$pct_reduction <- ifelse(df_cea$s_ch_tot > 0,
                                 100 * df_cea$s_ch_averted_tot / df_cea$s_ch_tot, 0)
  q <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE))
  df_cea %>%
    group_by(study_id, model, tau, vacc_cov) %>%
    summarise(
      country = first(country), population = first(population),
      year = first(year), s_ch_tot = first(s_ch_tot),
      ori_occurred = any(ori_occurred),
      pct_reduction_median = median(pct_reduction, na.rm = TRUE),
      pct_reduction_low = q(pct_reduction, 0.025), pct_reduction_high = q(pct_reduction, 0.975),
      s_ch_averted_median = median(s_ch_averted_tot, na.rm = TRUE),
      s_ch_averted_low = q(s_ch_averted_tot, 0.025), s_ch_averted_high = q(s_ch_averted_tot, 0.975),
      death_averted_median = median(death_averted_tot, na.rm = TRUE),
      s_ch_averted_amr_median = median(s_ch_averted_amr_tot, na.rm = TRUE),
      s_ch_averted_mdr_median = median(s_ch_averted_mdr_tot, na.rm = TRUE),
      s_ch_averted_fqns_median = median(s_ch_averted_fqns_tot, na.rm = TRUE),
      s_ch_averted_xdr_median = median(s_ch_averted_xdr_tot, na.rm = TRUE),
      daly_averted_median = if ("daly_averted" %in% names(df_cea)) median(daly_averted, na.rm = TRUE) else NA_real_,
      cost_per_daly_averted_median = if ("cost_per_daly_averted" %in% names(df_cea)) median(cost_per_daly_averted, na.rm = TRUE) else NA_real_,
      cost_per_daly_averted_low = if ("cost_per_daly_averted" %in% names(df_cea)) q(cost_per_daly_averted, 0.025) else NA_real_,
      cost_per_daly_averted_high = if ("cost_per_daly_averted" %in% names(df_cea)) q(cost_per_daly_averted, 0.975) else NA_real_,
      case_averted_per_1000_median = if ("case_averted_per_1000_OCV" %in% names(df_cea)) median(case_averted_per_1000_OCV, na.rm = TRUE) else NA_real_,
      .groups = "drop")
}

# Pooled (median of per-outbreak medians) and population-weighted % reduction,
# per (model, tau, vacc_cov). Used for the forest plot reference lines.
pooled_estimates <- function(summ) {
  summ %>% group_by(model, tau, vacc_cov) %>%
    summarise(
      pooled_pct_median = median(pct_reduction_median, na.rm = TRUE),
      popwt_pct = {
        wv <- population; pv <- pct_reduction_median
        ok <- is.finite(wv) & is.finite(pv)
        if (any(ok)) sum(wv[ok] * pv[ok]) / sum(wv[ok]) else NA_real_
      },
      total_averted_median = sum(s_ch_averted_median, na.rm = TRUE),
      .groups = "drop")
}

# eta_eff per outbreak at base timing/coverage: the fixed indirect VE the STATIC
# formula would need to reproduce the renewal TOTAL reduction. Computed per draw
# (using that draw's t_eff and renewal pct), then summarised to median + UI.
eta_eff_table <- function(prep, cfg, amr_props, tau = cfg$timing$delay_base_weeks,
                          cov = cfg$vaccine$coverage_base, psi_T_frac = 1) {
  w <- gi_from_config(cfg)
  params <- build_param_sets(cfg)
  out <- lapply(names(prep$series), function(sid) {
    inc <- prep$series[[sid]]; Tn <- length(inc); m <- prep$meta[[sid]]
    tail_guard <- Tn - cfg$tail_guard_weeks
    vals <- vapply(seq_len(nrow(params)), function(i) {
      p <- params[i, ]; psi_T <- psi_T_frac * p$dve
      te <- t_eff_weeks(tau, p$immuno_delay, p$campaign_duration, cfg$step_days)
      if (te > tail_guard) return(NA_real_)
      te_int <- max(round(te), 1L)
      ren <- renewal_counterfactual(inc, w, tau, te_int, cov, psi_T,
                                    c_shape = cfg$timing$c_shape, feedback = TRUE)
      rho <- impact_measures(inc, ren$incidence_v)$pct_reduction
      eta_eff(inc, te_int, rho, cov, psi_T)
    }, numeric(1))
    q <- function(p) as.numeric(stats::quantile(vals, p, na.rm = TRUE))
    data.frame(study_id = sid, country = m$country_iso, population = m$population,
               eta_eff_median = median(vals, na.rm = TRUE),
               eta_eff_low = q(0.025), eta_eff_high = q(0.975),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
