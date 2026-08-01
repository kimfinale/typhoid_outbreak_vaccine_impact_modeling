# =============================================================================
# Scenario runner: static vs renewal across Sobol draws x delay x coverage.
#
# Sobol design matches R/3-model.R exactly (dim 3: dve, campaign_duration,
# immuno_delay; n=200; seed in config) so static and renewal run on IDENTICAL
# parameter draws. Coverage and ORI delay (tau) are scenario dimensions, as in
# the existing pipeline. The ONLY difference between the two models is how
# cases-averted is computed (renewal_counterfactual vs static_counterfactual).
#
# Per-outbreak deaths and AMR-subtype cases averted are obtained by scaling the
# total cases averted by the outbreak's death-per-case and resistant-proportion
# ratios (constant across weeks), matching the static pipeline's construction.
# =============================================================================

# AMR resistant-fraction lookup keyed by study_id (proportions of tested cases).
load_amr_props <- function(cfg) {
  amr <- read.csv(cfg$paths$amr, stringsAsFactors = FALSE, check.names = FALSE)
  num <- function(x) suppressWarnings(as.numeric(x))
  tested <- num(amr$tot_tested_amr)
  data.frame(
    study_id = amr$study_id,
    prop_amr       = num(amr$tot_amr_cases) / tested,
    prop_resistant = num(amr$resistant) / tested,
    prop_mdr       = num(amr$mdr) / tested,
    prop_fqns      = num(amr$fqns) / tested,
    prop_xdr       = num(amr$xdr) / tested,
    stringsAsFactors = FALSE
  )
}

# Sobol parameter sets — identical construction to R/3-model.R.
build_param_sets <- function(cfg) {
  if (!requireNamespace("randtoolbox", quietly = TRUE) ||
      !requireNamespace("TruncExpFam", quietly = TRUE) ||
      !requireNamespace("truncnorm", quietly = TRUE))
    stop("Need randtoolbox, TruncExpFam, truncnorm (as in R/3-model.R).")
  sob <- randtoolbox::sobol(n = cfg$n_sobol, dim = 3, scrambling = 0, seed = cfg$seed)
  data.frame(
    dve = truncnorm::qtruncnorm(sob[, 1], a = 0, b = 1,
                                mean = cfg$vaccine$psi_mean, sd = cfg$vaccine$psi_sd),
    campaign_duration = TruncExpFam::qtruncgamma(
      sob[, 2], a = cfg$timing$campaign_duration_sobol[1],
      b = cfg$timing$campaign_duration_sobol[2], shape = 16, scale = 0.5, rate = 2),
    immuno_delay = stats::qunif(sob[, 3], min = cfg$timing$immunity_onset_sobol[1],
                                max = cfg$timing$immunity_onset_sobol[2])
  )
}

# t_eff (first protected week) for a draw — identical formula to R/3-model.R:156.
t_eff_weeks <- function(tau, immuno_delay_days, campaign_duration_days, step_days = 7) {
  tau + immuno_delay_days / step_days + round(campaign_duration_days / step_days / 2)
}

# Run all scenarios. Returns a tidy long data frame at the
# (study_id, tau, coverage, draw, model) level with averted quantities and the
# metadata columns add_cea_results() needs.
#
# psi_T_frac scales transmission VE (psi_T = frac * dve); gi_mean_days overrides
# the GI for sensitivity. eta is applied to the STATIC model only.
run_scenarios <- function(prep, cfg, amr_props,
                          coverage = cfg$vaccine$coverage_base,
                          delays   = cfg$timing$delay_base_weeks,
                          psi_T_frac = 1, gi_mean_days = cfg$gi$mean_days,
                          eta = cfg$vaccine$eta_static, pi_post = NULL,
                          ppv_regime = if (is.null(cfg$ppv$regime)) "multiplicative" else cfg$ppv$regime) {
  w <- gi_from_config(cfg, mean_days = gi_mean_days)
  params <- build_param_sets(cfg)
  ndraw <- nrow(params)
  # PPV (pi) draws per (study, draw); the additive OBSERVATION regime allocates
  # suspected incidence into true typhoid plus a fixed other-febrile residual.
  pim <- if (!is.null(pi_post))
    build_pi_matrix(pi_post, names(prep$series), ndraw, seed = cfg$seed + 777L) else NULL
  additive <- !is.null(pim) && identical(ppv_regime, "additive")
  N <- length(prep$series) * length(delays) * length(coverage) * ndraw * 2L

  # Preallocate columns (filled by index; one data.frame built at the end).
  chr <- function() character(N); dbl <- function() numeric(N)
  C <- list(study_id = chr(), country = chr(), year = integer(N), who_region = chr(),
            population = dbl(), tau = dbl(), vacc_cov = dbl(), draw = integer(N),
            model = chr(), psi = dbl(), psi_T = dbl(), eta = dbl(), gi_mean = dbl(),
            t_eff = dbl(), ori_occurred = logical(N), s_ch_tot = dbl(),
            post_int_cases = dbl(),
            s_ch_averted_tot = dbl(), death_averted_tot = dbl(),
            s_ch_averted_amr_tot = dbl(), s_ch_averted_resistant_tot = dbl(),
            s_ch_averted_mdr_tot = dbl(), s_ch_averted_fqns_tot = dbl(),
            s_ch_averted_xdr_tot = dbl(),
            pi_ppv = dbl(), suspected_tot = dbl(), obs_pct_reduction = dbl())
  k <- 0L

  for (sid in names(prep$series)) {
    inc <- prep$series[[sid]]; m <- prep$meta[[sid]]; Tn <- length(inc)
    ap <- amr_props[amr_props$study_id == sid, ]
    pr <- function(col) if (nrow(ap) == 1 && is.finite(ap[[col]])) ap[[col]] else 0
    props <- c(amr = pr("prop_amr"), resistant = pr("prop_resistant"), mdr = pr("prop_mdr"),
               fqns = pr("prop_fqns"), xdr = pr("prop_xdr"))
    death_per_case <- if (isTRUE(m$tot_cases > 0) && is.finite(m$tot_deaths))
      m$tot_deaths / m$tot_cases else 0
    tail_guard <- Tn - cfg$tail_guard_weeks   # outbreak "missed" if t_eff beyond this
    sum_inc <- sum(inc)                        # observed suspected total (dilution denom)

    for (tau in delays) for (cov in coverage) for (i in seq_len(ndraw)) {
      psi   <- params$dve[i]
      psi_T <- psi_T_frac * psi
      te <- t_eff_weeks(tau, params$immuno_delay[i], params$campaign_duration[i], cfg$step_days)
      te_int <- max(round(te), 1L); occurred <- te <= tail_guard

      # ADDITIVE observation regime models TRUE typhoid T(t) from S(t)=T(t)+F(t).
      pi_i    <- if (!is.null(pim)) pim[i, sid] else NA_real_
      inc_i   <- if (additive) ppv_true_incidence(inc, pi_i) else inc
      s_tot_i <- if (additive) sum(inc_i) else m$tot_cases
      dpc_i   <- if (additive) { s <- sum(inc_i)
                   if (is.finite(m$tot_deaths) && s > 0) m$tot_deaths / s else 0 } else death_per_case

      post_int <- if (occurred && te_int <= Tn) sum(inc_i[te_int:Tn]) else 0
      if (!occurred) { a_ren <- 0; a_sta <- 0 } else {
        ren <- renewal_counterfactual(inc_i, w, tau = tau, t_eff = te_int, pi = cov,
                                      psi_T = psi_T, c_shape = cfg$timing$c_shape, feedback = TRUE)
        a_ren <- impact_measures(inc_i, ren$incidence_v)$averted
        a_sta <- impact_measures(inc_i, static_counterfactual(inc_i, te_int, P_halloran(cov, psi, eta)))$averted
      }

      for (j in 1:2) {
        mdl <- if (j == 1L) "renewal" else "static"
        a   <- if (j == 1L) a_ren else a_sta
        k <- k + 1L
        C$study_id[k] <- sid; C$country[k] <- m$country_iso; C$year[k] <- m$year
        C$who_region[k] <- m$who_region; C$population[k] <- m$population
        C$tau[k] <- tau; C$vacc_cov[k] <- cov; C$draw[k] <- i; C$model[k] <- mdl
        C$psi[k] <- psi; C$psi_T[k] <- psi_T; C$eta[k] <- if (j == 2L) eta else NA_real_
        C$gi_mean[k] <- gi_mean_days; C$t_eff[k] <- te_int; C$ori_occurred[k] <- occurred
        C$s_ch_tot[k] <- s_tot_i; C$post_int_cases[k] <- post_int
        C$s_ch_averted_tot[k] <- a
        C$death_averted_tot[k] <- a * dpc_i
        C$s_ch_averted_amr_tot[k] <- a * props["amr"]
        C$s_ch_averted_resistant_tot[k] <- a * props["resistant"]
        C$s_ch_averted_mdr_tot[k] <- a * props["mdr"]
        C$s_ch_averted_fqns_tot[k] <- a * props["fqns"]
        C$s_ch_averted_xdr_tot[k] <- a * props["xdr"]
        C$pi_ppv[k] <- pi_i; C$suspected_tot[k] <- sum_inc
        C$obs_pct_reduction[k] <- if (additive && sum_inc > 0) 100 * a / sum_inc else NA_real_
      }
    }
  }
  df <- as.data.frame(C, stringsAsFactors = FALSE)
  # MULTIPLICATIVE regime: scale case columns SUSPECTED -> TRUE typhoid (pi-invariant %reduction).
  if (!is.null(pim) && !additive) {
    cc <- intersect(.ppv_case_cols, names(df))
    df[cc] <- df[cc] * df$pi_ppv
    df$obs_pct_reduction <- ifelse(df$s_ch_tot > 0, 100 * df$s_ch_averted_tot / df$s_ch_tot, NA_real_)
  }
  df
}
