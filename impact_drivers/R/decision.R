# =============================================================================
# Phase 4-5 — characteristic x delay interaction, impact map, decision rule.
# =============================================================================

suppressMessages({ library(dplyr) })

# Delay-sensitivity by early growth: on the simulation, the slope d(impact)/d(tau)
# within strata of R_t-at-vaccination. Tests "every week counts more for fast
# (still-growing) outbreaks".
delay_sensitivity <- function(sim, outcome = "pct_reduction") {
  d <- sim[is.finite(sim[[outcome]]), ]
  d$rt_stratum <- cut(d$Rt_at_tau, breaks = c(-Inf, 0.8, 1.2, 2, Inf),
                      labels = c("<0.8 (declining)", "0.8-1.2", "1.2-2", ">2 (explosive)"))
  d %>% filter(!is.na(rt_stratum)) %>% group_by(rt_stratum) %>%
    summarise(slope_per_week = tryCatch(unname(coef(lm(.data[[outcome]] ~ tau))[2]),
                                        error = function(e) NA),
              mean_impact = mean(.data[[outcome]], na.rm = TRUE), n = n(), .groups = "drop")
}

# Predict the real outbreaks' impact from the simulation-trained surface model
# (external validation: do the real outbreaks lie on the simulated surface?).
predict_real_from_sim <- function(surface_model, real_dat, outcome) {
  d <- real_dat[is.finite(real_dat[[outcome]]) & is.finite(real_dat$Rt_at_tau), ]
  d$log_pop[!is.finite(d$log_pop)] <- median(d$log_pop, na.rm = TRUE)
  d$predicted <- as.numeric(predict(surface_model, newdata = d))
  data.frame(study_id = d$study_id, tau = d$tau, actual = d[[outcome]],
             predicted = d$predicted, stringsAsFactors = FALSE)
}

# Honest small-n skill: leave-one-OUTBREAK-out CV using the single mechanistic
# predictor (R_t-at-vaccination) on the real data, at a fixed delay.
loo_cv_real <- function(real_dat, outcome, predictor = "Rt_at_tau", tau = 8) {
  d <- real_dat[real_dat$tau == tau, ]
  d <- d[is.finite(d[[outcome]]) & is.finite(d[[predictor]]), c("study_id", predictor, outcome)]
  names(d) <- c("study_id", "x", "y")
  pred <- sapply(seq_len(nrow(d)), function(i) {
    fit <- lm(y ~ x, data = d[-i, ]); predict(fit, newdata = d[i, ])
  })
  ss_res <- sum((d$y - pred)^2); ss_tot <- sum((d$y - mean(d$y))^2)
  list(table = data.frame(study_id = d$study_id, actual = d$y, loo_pred = as.numeric(pred)),
       loo_r2 = 1 - ss_res / ss_tot,
       loo_cor = suppressWarnings(cor(d$y, pred)), n = nrow(d))
}

# Interpretable decision rule: regression tree on the SIMULATION (large n) using
# only prospective features, predicting the primary efficiency outcome.
fit_decision_tree <- function(sim, outcome = "case_averted_per_1000",
                              preds = c("Rt_at_tau", "r_early", "cum_AR_by_tau", "tau", "log_pop")) {
  d <- sim[is.finite(sim[[outcome]]), ]
  rpart::rpart(reformulate(preds, outcome), data = d,
               control = rpart::rpart.control(maxdepth = 3, cp = 0.01, minbucket = 50))
}
