# =============================================================================
# Phase 3 — simulation experiment (the workhorse that beats small n).
#
# The renewal model is used as a GENERATIVE engine: synthetic outbreaks with
# CONTROLLED characteristics (initial R0 / early growth, GI, decline timing,
# population) are forward-simulated, vaccinated at each delay, and their impact
# measured. Because we set the drivers, attribution is confounding-free and n is
# unlimited. The real 13 outbreaks are later overlaid to check they lie on the
# simulated surface.
# =============================================================================

suppressMessages({ library(dplyr) })

# Forward-simulate the renewal equation with a prescribed R_t path (deterministic).
.sim_incidence <- function(Rt_path, w, Tn, seed = c(2, 3, 5)) {
  inc <- numeric(Tn); inc[seq_along(seed)] <- seed; K <- length(w)
  for (t in (length(seed) + 1L):Tn) {
    s_max <- min(K, t - 1L)
    inc[t] <- Rt_path[t] * sum(w[1:s_max] * inc[t - (1:s_max)])
  }
  inc
}

# Logistic R_t path: starts at R0 (sets early growth), crosses 1 near t_cross,
# settles at R_end (response/source remediation).
.rt_path <- function(Tn, R0, R_end, t_cross, steep) {
  t <- seq_len(Tn)
  R_end + (R0 - R_end) / (1 + exp(steep * (t - t_cross)))
}

# Run the LHS simulation sweep. Returns one row per (sim, delay) with the SAME
# prospective features as the real data and the three impact outcomes.
run_simulation <- function(rcfg, cfg, n_sim = 1500, Tn = 45, sim_delta = 4) {
  cov <- cfg$coverage; psi_T <- rcfg$vaccine$psi_mean; k <- cfg$growth_window_weeks
  sob <- randtoolbox::sobol(n = n_sim, dim = 6, scrambling = 1, seed = cfg$seed)
  P <- data.frame(
    R0      = 1.3 + sob[, 1] * (3.5 - 1.3),
    R_end   = 0.30 + sob[, 2] * (0.85 - 0.30),
    t_cross = 4 + sob[, 3] * (18 - 4),
    steep   = 0.4 + sob[, 4] * (1.2 - 0.4),
    gi_mean = 10 + sob[, 5] * (28 - 10),
    log_pop = 4 + sob[, 6] * (6.5 - 4))
  delays <- cfg$delays_weeks
  rows <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    p <- P[i, ]
    w <- discretize_gi(mean_days = p$gi_mean, sd_days = p$gi_mean * rcfg$gi$cv,
                       step_days = rcfg$step_days, max_steps = rcfg$gi$max_steps,
                       drop_first = rcfg$gi$drop_first_week)
    Rt_path <- .rt_path(Tn, p$R0, p$R_end, p$t_cross, p$steep)
    inc <- .sim_incidence(Rt_path, w, Tn)
    if (!all(is.finite(inc)) || sum(inc) <= 0) next
    Rt_hat <- reconstruct_Rt(inc, w)$Rt
    pop <- 10^p$log_pop
    r_early <- { kk <- min(k, Tn); wk <- seq_len(kk)
      unname(coef(stats::lm(log(inc[wk] + 0.5) ~ wk))[2]) }
    di <- lapply(delays, function(tau) {
      tt <- min(tau, Tn); te <- tau + sim_delta
      if (te > Tn - rcfg$tail_guard_weeks) { rho <- 0; av <- 0 } else {
        cf <- renewal_counterfactual(inc, w, tau = tau, t_eff = te, pi = cov,
                                     psi_T = psi_T, c_shape = "step", feedback = TRUE)
        im <- impact_measures(inc, cf$incidence_v); rho <- im$pct_reduction; av <- im$averted
      }
      data.frame(sim_id = i, tau = tau,
                 r_early = r_early, Rt_at_tau = Rt_hat[tt],
                 cum_AR_by_tau = 100 * sum(inc[seq_len(tt)]) / pop,
                 log_pop = p$log_pop, gi_mean = p$gi_mean,
                 pct_reduction = rho, cases_averted = av,
                 case_averted_per_1000 = ifelse(pop * cov > 0, 1000 * av / (pop * cov), NA),
                 stringsAsFactors = FALSE)
    })
    rows[[i]] <- do.call(rbind, di)
  }
  do.call(rbind, rows)
}

# Univariate GAM deviance-explained per predictor = global importance ranking
# (large-n, confounding-free). Plus a joint GAM for the response surface.
surface_importance <- function(sim, outcome,
                               preds = c("Rt_at_tau","r_early","cum_AR_by_tau","log_pop","tau","gi_mean")) {
  d <- sim[is.finite(sim[[outcome]]), ]
  imp <- sapply(preds, function(p) {
    f <- as.formula(sprintf("%s ~ s(%s, k=5)", outcome, p))
    tryCatch(summary(mgcv::gam(f, data = d))$dev.expl, error = function(e) NA)
  })
  data.frame(outcome = outcome, predictor = names(imp),
             dev_explained = round(as.numeric(imp), 3),
             stringsAsFactors = FALSE) %>% arrange(desc(dev_explained))
}

# Joint surface model for prediction (impact map / interaction).
fit_surface_model <- function(sim, outcome) {
  d <- sim[is.finite(sim[[outcome]]), ]
  mgcv::gam(as.formula(sprintf(
    "%s ~ ti(Rt_at_tau, k=5) + ti(tau, k=5) + ti(Rt_at_tau, tau, k=4) + s(log_pop, k=5)",
    outcome)), data = d)
}
