# =============================================================================
# Full Bayesian renewal model (self-contained; no Stan/compilation required).
#
# Latent log R_t random walk (joint temporal smoothing) + Poisson renewal
# likelihood I_t ~ Poisson(R_t * Lambda_t), with the generation-interval mean
# treated as a sampled parameter so GI uncertainty is integrated natively.
# Componentwise adaptive Metropolis: each log R_t enters only its own Poisson
# term plus random-walk neighbours, so updates are local and cheap; gi_mean is
# updated jointly (it changes every Lambda_t).
#
# This is a deliberately different reconstruction from the exact estimator
# (it imposes a smoothing prior and returns a full joint posterior), so the
# comparison probes whether conclusions are robust to the inference approach.
# =============================================================================

# Build the (acute) weekly GI kernel from a mean, reusing the renewal discretizer.
.bayes_w <- function(gi_mean, cfg, rcfg)
  discretize_gi(mean_days = gi_mean, sd_days = gi_mean * rcfg$gi$cv,
                step_days = rcfg$step_days, max_steps = rcfg$gi$max_steps,
                drop_first = rcfg$gi$drop_first_week)

# Poisson log-likelihood (dropping the constant lgamma term).
.pois_ll <- function(y, mu) ifelse(mu > 0, y * log(mu) - mu, ifelse(y == 0, 0, -Inf))

# Fit one outbreak. Returns posterior draws of R_t (matrix kept x T) and gi_mean.
bayes_fit_outbreak <- function(inc, cfg, rcfg, gi_lo = 10, gi_hi = 28,
                               n_iter = 1500, n_warmup = 500, n_chains = 2,
                               rw_sd = 0.4, thin = 2) {
  Tn <- length(inc)
  sdlog <- (log(gi_hi) - log(gi_lo)) / (2 * 1.96)
  lprior_gi <- function(gm) dlnorm(gm, log(rcfg$gi$mean_days), sdlog, log = TRUE)
  lprior_rw <- function(lr) sum(dnorm(diff(lr), 0, rw_sd, log = TRUE)) +
    dnorm(lr[1], 0, 1, log = TRUE)

  chains <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    gi <- rcfg$gi$mean_days
    w <- .bayes_w(gi, cfg, rcfg); lam <- total_infectiousness(inc, w)
    Rt0 <- ifelse(lam > 0, inc / lam, 1); Rt0[!is.finite(Rt0) | Rt0 <= 0] <- 0.5
    lr <- log(pmax(Rt0, 0.05))
    step_lr <- rep(0.3, Tn); step_gi <- 0.08
    acc_lr <- rep(0, Tn); acc_gi <- 0; n_since <- 0
    keep_lr <- list(); keep_gi <- numeric()

    for (it in seq_len(n_iter)) {
      # --- update each log R_t (local: own Poisson term + RW neighbours) ------
      for (t in seq_len(Tn)) {
        prop <- lr[t] + rnorm(1, 0, step_lr[t])
        ll_cur <- if (lam[t] > 0) .pois_ll(inc[t], exp(lr[t]) * lam[t]) else 0
        ll_pro <- if (lam[t] > 0) .pois_ll(inc[t], exp(prop) * lam[t]) else 0
        pr_cur <- (if (t > 1) dnorm(lr[t] - lr[t - 1], 0, rw_sd, log = TRUE) else dnorm(lr[t], 0, 1, log = TRUE)) +
                  (if (t < Tn) dnorm(lr[t + 1] - lr[t], 0, rw_sd, log = TRUE) else 0)
        pr_pro <- (if (t > 1) dnorm(prop - lr[t - 1], 0, rw_sd, log = TRUE) else dnorm(prop, 0, 1, log = TRUE)) +
                  (if (t < Tn) dnorm(lr[t + 1] - prop, 0, rw_sd, log = TRUE) else 0)
        if (log(runif(1)) < (ll_pro + pr_pro) - (ll_cur + pr_cur)) {
          lr[t] <- prop; acc_lr[t] <- acc_lr[t] + 1
        }
      }
      # --- update gi_mean (changes all Lambda) --------------------------------
      gi_p <- exp(log(gi) + rnorm(1, 0, step_gi))
      if (gi_p >= gi_lo && gi_p <= gi_hi) {
        w_p <- .bayes_w(gi_p, cfg, rcfg); lam_p <- total_infectiousness(inc, w_p)
        ll_c <- sum(.pois_ll(inc[lam > 0], exp(lr[lam > 0]) * lam[lam > 0]))
        ll_p <- sum(.pois_ll(inc[lam_p > 0], exp(lr[lam_p > 0]) * lam_p[lam_p > 0]))
        if (log(runif(1)) < (ll_p + lprior_gi(gi_p)) - (ll_c + lprior_gi(gi))) {
          gi <- gi_p; w <- w_p; lam <- lam_p; acc_gi <- acc_gi + 1
        }
      }
      n_since <- n_since + 1
      # --- adapt step sizes during warmup (target ~0.3 acceptance) -----------
      if (it <= n_warmup && n_since == 100) {
        step_lr <- step_lr * exp((acc_lr / 100 - 0.3))
        step_gi <- step_gi * exp((acc_gi / 100 - 0.3))
        acc_lr[] <- 0; acc_gi <- 0; n_since <- 0
      }
      if (it > n_warmup && (it %% thin == 0)) {
        keep_lr[[length(keep_lr) + 1]] <- exp(lr); keep_gi <- c(keep_gi, gi)
      }
    }
    chains[[ch]] <- list(Rt = do.call(rbind, keep_lr), gi = keep_gi)
  }
  list(Rt = do.call(rbind, lapply(chains, `[[`, "Rt")),
       gi = unlist(lapply(chains, `[[`, "gi")))
}

# Propagate the posterior R_t draws through the vaccine counterfactual.
# Returns posterior % reduction (one per kept draw) at the base scenario.
bayes_impact <- function(inc, fit, cfg, rcfg, te) {
  cov <- cfg$scenario$coverage; psi_T <- rcfg$vaccine$psi_mean
  Tn <- length(inc)
  if (te > (Tn - rcfg$tail_guard_weeks)) return(rep(0, nrow(fit$Rt)))
  vapply(seq_len(nrow(fit$Rt)), function(d) {
    w <- .bayes_w(fit$gi[d], cfg, rcfg)
    Rt <- fit$Rt[d, ]
    c_t <- as.numeric(seq_len(Tn) >= te)
    Rt_v <- Rt * (1 - c_t * cov * psi_T)
    vax <- inc
    for (t in seq_len(Tn)) {
      if (t < te) next
      s_max <- min(length(w), t - 1L)
      lam_t <- if (s_max >= 1L) sum(w[1:s_max] * vax[t - (1:s_max)]) else 0
      vax[t] <- Rt_v[t] * lam_t
    }
    100 * sum(inc - vax) / sum(inc)
  }, numeric(1))
}
