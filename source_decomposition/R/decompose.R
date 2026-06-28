# =============================================================================
# Source-decomposition counterfactual.
#
# A fraction theta of transmission is exogenous common-source (non-reproductive:
# the static proportional model is correct) and (1-theta) is propagated (the
# renewal model with herd compounding is correct). Because both the static and
# renewal counterfactuals scale LINEARLY with the incidence magnitude, splitting
# the curve into theta and (1-theta) and applying each model is exactly the
# theta-weighted average of the two whole-curve results:
#
#   averted(theta) = theta * averted_static + (1 - theta) * averted_renewal
#   pct(theta)     = theta * pct_static     + (1 - theta) * pct_renewal
#
# Endpoints: theta = 1 -> static exactly; theta = 0 -> renewal exactly.
# =============================================================================

suppressMessages({ library(dplyr); library(tidyr) })

# Per-draw static & renewal at the base scenario, blended by each outbreak's theta.
decompose_bracket <- function(scn, theta_tab, cfg) {
  base <- scn %>% filter(tau == cfg$scenario$tau_weeks, vacc_cov == cfg$scenario$coverage) %>%
    mutate(pct = ifelse(s_ch_tot > 0, 100 * s_ch_averted_tot / s_ch_tot, 0)) %>%
    select(study_id, draw, model, averted = s_ch_averted_tot, pct) %>%
    pivot_wider(names_from = model, values_from = c(averted, pct))
  d <- base %>% left_join(theta_tab[, c("study_id", "final_theta", "final_theta_lo",
                                        "final_theta_hi")], by = "study_id") %>%
    filter(is.finite(final_theta)) %>%
    mutate(
      averted_theta    = final_theta * averted_static + (1 - final_theta) * averted_renewal,
      averted_theta_lo = final_theta_hi * averted_static + (1 - final_theta_hi) * averted_renewal,
      averted_theta_hi = final_theta_lo * averted_static + (1 - final_theta_lo) * averted_renewal,
      pct_theta        = final_theta * pct_static + (1 - final_theta) * pct_renewal)
  q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
  d %>% group_by(study_id) %>%
    summarise(
      final_theta = first(final_theta),
      pct_static = median(pct_static), pct_renewal = median(pct_renewal),
      pct_theta = median(pct_theta),
      averted_static = median(averted_static), averted_renewal = median(averted_renewal),
      averted_theta = median(averted_theta),
      # theta-range sensitivity (uncertainty in the source fraction)
      averted_theta_lo = median(averted_theta_lo), averted_theta_hi = median(averted_theta_hi),
      # parameter (draw) uncertainty at the point theta
      averted_theta_p025 = q(averted_theta, 0.025), averted_theta_p975 = q(averted_theta, 0.975),
      .groups = "drop")
}

# Validation: theta=1 reproduces static, theta=0 reproduces renewal (machine precision).
validate_endpoints <- function(scn, cfg) {
  base <- scn %>% filter(tau == cfg$scenario$tau_weeks, vacc_cov == cfg$scenario$coverage) %>%
    select(study_id, draw, model, averted = s_ch_averted_tot) %>%
    pivot_wider(names_from = model, values_from = averted)
  blend <- function(th) th * base$static + (1 - th) * base$renewal
  c(err_theta1_vs_static  = max(abs(blend(1) - base$static), na.rm = TRUE),
    err_theta0_vs_renewal = max(abs(blend(0) - base$renewal), na.rm = TRUE))
}

# Cross-validate Stream A (paper) and Stream B (data) -> final theta + agreement.
cross_validate <- function(streamA, streamB) {
  a <- streamA %>% mutate(theta_paper = (theta_paper_lo + theta_paper_hi) / 2)
  m <- a %>% left_join(streamB, by = "study_id") %>%
    mutate(
      agreement = ifelse(is.finite(theta_data) & abs(theta_paper - theta_data) <= 0.2,
                         "agree", "disagree"),
      # Paper evidence (Stream A) is the primary classifier; Stream B (data) is a
      # weak cross-check given identifiability, so it only WIDENS the interval --
      # it does not move the point estimate.
      final_theta = theta_paper,
      final_theta_lo = pmin(theta_paper_lo, ifelse(is.finite(theta_data), theta_data, theta_paper_lo)),
      final_theta_hi = pmax(theta_paper_hi, ifelse(is.finite(theta_data), theta_data, theta_paper_hi)))
  m
}
