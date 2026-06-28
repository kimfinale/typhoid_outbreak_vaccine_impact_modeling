# =============================================================================
# Helpers for the hierarchical renewal-with-source model: Stan data assembly,
# synthetic data generation (matching the generative model), recovery metrics.
# =============================================================================

# Weekly generation-interval kernel (gamma CDF differences, renormalized).
gi_w <- function(gi_mean, gi_cv, K) {
  shape <- 1 / gi_cv^2; rate <- shape / gi_mean
  w <- pgamma(1:K, shape, rate) - pgamma(0:(K - 1), shape, rate)
  w / sum(w)
}

# Build the Stan data list from a named list of incidence vectors.
# modes: optional integer vector (1=A,2=B,3=C) aligned to series; informed flag.
build_stan_data <- function(series, cfg, modes = NULL, informed = 0L) {
  O <- length(series); Tn <- vapply(series, length, integer(1))
  I <- as.integer(round(unlist(series, use.names = FALSE)))
  start <- as.integer(c(1, head(cumsum(Tn), -1) + 1))
  logA_center <- vapply(series, function(x) log(mean(x) + 1), numeric(1))
  if (is.null(modes)) modes <- rep(1L, O)
  list(O = O, N = length(I), I = I, Tn = as.integer(Tn), start = start,
       logA_center = as.numeric(logA_center),
       K = as.integer(cfg$gi$max_lag_weeks), gi_cv = cfg$gi$cv,
       gi_mean_prior_mu = cfg$gi$mean_prior_log_mu,
       gi_mean_prior_sd = cfg$gi$mean_prior_log_sd,
       informed = as.integer(informed), n_modes = 3L,
       mode = as.integer(modes))
}

# Simulate synthetic outbreaks under the generative model with KNOWN theta.
# Returns series (integer vectors), true_theta, true_Rtr, modes, gi_mean.
simulate_synthetic <- function(cfg, gi_mean = 2.0) {
  set.seed(cfg$seed)
  K <- cfg$gi$max_lag_weeks; w <- gi_w(gi_mean, cfg$gi$cv, K)
  grid <- cfg$gate$theta_grid; n <- cfg$gate$n_synth
  thetas <- rep_len(grid, n)
  series <- vector("list", n); modes <- integer(n)
  for (o in seq_len(n)) {
    th <- thetas[o]; Rtr <- 1 - th
    Tn <- sample(14:40, 1)
    loc <- runif(1, 3, max(4, Tn / 3)); wid <- runif(1, 2, 5)
    A <- exp(runif(1, log(20), log(400)))          # amplitude (total source scale)
    phi <- runif(1, 4, 15)
    psh <- (loc / wid)^2; prt <- loc / wid^2
    p <- dgamma(1:Tn, psh, prt); p <- p / sum(p)
    inc <- numeric(Tn)
    for (t in seq_len(Tn)) {
      lam <- 0; kmax <- min(K, t - 1)
      if (kmax >= 1) lam <- sum(w[1:kmax] * inc[t - (1:kmax)])
      mu <- A * p[t] + Rtr * lam
      inc[t] <- rnbinom(1, mu = max(mu, 1e-6), size = phi)
    }
    series[[o]] <- inc
    modes[o] <- if (th >= 0.66) 1L else if (th <= 0.34) 2L else 3L   # A / B / C
  }
  names(series) <- paste0("synth_", seq_len(n))
  list(series = series, true_theta = thetas, true_Rtr = 1 - thetas,
       modes = modes, gi_mean = gi_mean)
}

# Recovery metrics: posterior theta draws (matrix draws x O) vs true theta.
recovery_metrics <- function(theta_draws, true_theta) {
  med <- apply(theta_draws, 2, median)
  lo <- apply(theta_draws, 2, quantile, 0.05)
  hi <- apply(theta_draws, 2, quantile, 0.95)
  covered <- (true_theta >= lo) & (true_theta <= hi)
  list(
    per = data.frame(true_theta = true_theta, post_median = med,
                     cri90_lo = lo, cri90_hi = hi, covered = covered),
    coverage = mean(covered),
    median_bias = median(med - true_theta),
    rank_cor = suppressWarnings(cor(med, true_theta, method = "spearman")))
}
