# =============================================================================
# Sample the uncertainty parameters (ranges from uncertainty/params_uncertainty.csv,
# evidence in uncertainty/EVIDENCE.md). Returns n draws as a data frame.
# psi (direct VE) is drawn from the renewal config's truncated-normal, matching
# the static/renewal Sobol pipeline.
# =============================================================================

sample_uncertainty_params <- function(cfg, rcfg, n = cfg$n_mc) {
  p <- read.csv(cfg$params_csv, stringsAsFactors = FALSE)
  rng <- function(nm) { r <- p[p$param == nm, ]; c(low = r$low, high = r$high, base = r$base) }
  u <- function(nm) { r <- rng(nm); stats::runif(n, r["low"], r["high"]) }

  gm <- rng("gi_mean")
  # lognormal for the GI mean, clamped to its plausible bounds
  sdlog <- (log(gm["high"]) - log(gm["low"])) / (2 * 1.96)
  gi_mean <- pmin(pmax(stats::rlnorm(n, meanlog = log(gm["base"]), sdlog = sdlog),
                       gm["low"]), gm["high"])
  sc <- rng("culture_sensitivity")
  s_cult <- pmin(pmax(stats::rnorm(n, sc["base"], (sc["high"] - sc["low"]) / (2 * 1.96)),
                      sc["low"]), sc["high"])

  data.frame(
    draw = seq_len(n),
    gi_mean = as.numeric(gi_mean),
    gi_cv = u("gi_cv"),
    pi_c = u("carrier_tail_weight"),
    c_c = u("carrier_rel_infectiousness"),
    mu_carrier = u("carrier_gi_mean"),
    f_a = u("asymptomatic_fraction"),
    c_a = u("asymptomatic_rel_infectiousness"),
    ppv_level = u("ppv_level"),
    ppv_trend = u("ppv_trend"),
    s_cult = as.numeric(s_cult),
    psi = truncnorm::rtruncnorm(n, a = 0, b = 1,
                                mean = rcfg$vaccine$psi_mean, sd = rcfg$vaccine$psi_sd)
  )
}
