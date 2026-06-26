# =============================================================================
# Generation-interval (GI) kernel.
#
# The discretized GI {w_s} (s >= 1, sum w_s = 1) is the ONLY transmission
# assumption the renewal model makes, so it is parameterized explicitly and
# subjected to sensitivity analysis (mu_g in 7-28 days).
#
# Pure functions, no global state.
# =============================================================================

# Discretize a continuous gamma GI onto a time step (default weekly).
#   w_s proportional to G(s*step) - G((s-1)*step), s >= 1, w_0 = 0, then renormalize.
# `drop_first` additionally zeroes the first whole bin (0-1 step): for a weekly
# step this removes within-week onward transmission, defensible because the
# typhoid incubation period alone is >= 1 week (validated prototype choice).
#
# Args:
#   mean_days, sd_days : GI mean and SD in DAYS.
#   step_days          : reporting/analysis step in days (7 = weekly).
#   max_steps          : GI support in number of steps.
#   drop_first         : zero the first step bin then renormalize.
# Returns: numeric vector w of length max_steps, sum(w) == 1, w indexed by s=1..max_steps.
discretize_gi <- function(mean_days = 14, sd_days = 8.4, step_days = 7,
                          max_steps = 20, drop_first = TRUE) {
  stopifnot(mean_days > 0, sd_days > 0, step_days > 0, max_steps >= 2)
  shape <- (mean_days / sd_days)^2
  rate  <- mean_days / sd_days^2
  edges <- (0:max_steps) * step_days            # bin edges in days
  w <- diff(pgamma(edges, shape = shape, rate = rate))  # mass in bins s = 1..max_steps
  if (drop_first) w[1] <- 0                     # no transmission in first step
  s <- sum(w)
  if (s <= 0) stop("Degenerate GI: all mass outside support; increase max_steps.")
  w / s
}

# Convenience: build the base-case GI from a config list.
gi_from_config <- function(cfg, mean_days = NULL) {
  md <- if (is.null(mean_days)) cfg$gi$mean_days else mean_days
  sd <- md * cfg$gi$cv
  discretize_gi(mean_days = md, sd_days = sd, step_days = cfg$step_days,
                max_steps = cfg$gi$max_steps, drop_first = cfg$gi$drop_first_week)
}

# Fraction of GI probability mass that falls in the first NON-zero step.
# Used by the resolution diagnostic: when the GI mass piles into one reporting
# interval, successive transmission generations collapse and R_t degenerates.
gi_first_bin_mass <- function(mean_days, sd_days, step_days = 7, max_steps = 20) {
  shape <- (mean_days / sd_days)^2
  rate  <- mean_days / sd_days^2
  edges <- (0:max_steps) * step_days
  w <- diff(pgamma(edges, shape = shape, rate = rate))
  w[1] / sum(w)
}
