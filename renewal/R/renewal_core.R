# =============================================================================
# Renewal-equation counterfactual engine for ORI impact.
#
# Manuscript notation:
#   I_t        observed weekly incidence (the static model's `suspected_cases`,
#              which is `total_cases` relabeled in R/1-data.R)
#   w_s        discretized generation interval (sum = 1), see gi.R
#   Lambda_t   total infectiousness = sum_{s>=1} w_s I_{t-s}
#   R_t        instantaneous reproduction number
#   tau        ORI campaign timing; delta = campaign+immunity delay
#   t_eff      first protected week = tau + delta
#   pi         coverage; psi direct VE; psi_T transmission-relevant VE (<= psi)
#   eta        indirect (herd) VE   [STATIC model only]
#
# All functions are pure (no global state). Equations are commented inline so
# the vaccine-effect and propagation mechanism is auditable.
# =============================================================================

# --- Total infectiousness  Lambda_t = sum_{s>=1} w_s I_{t-s} -----------------
total_infectiousness <- function(incidence, w) {
  Tn <- length(incidence); K <- length(w)
  lam <- numeric(Tn)
  for (t in seq_len(Tn)) {
    s_max <- min(K, t - 1L)
    if (s_max >= 1L) lam[t] <- sum(w[1:s_max] * incidence[t - (1:s_max)])
  }
  lam
}

# --- Exact, self-consistent R_t reconstruction -------------------------------
# R_hat_t = I_t^obs / Lambda_t^obs. Forward-propagating with R_hat reproduces
# the observed curve exactly, so the observed outbreak IS the counterfactual
# baseline (required for a clean counterfactual; preserved by every routine here).
reconstruct_Rt <- function(incidence, w) {
  lam <- total_infectiousness(incidence, w)
  Rt  <- ifelse(lam > 0, incidence / lam, NA_real_)
  list(Rt = Rt, lambda = lam)
}

# --- Cumulative protected fraction c(t) in [0,1] -----------------------------
# step (manuscript default): c(t) = 1(t >= t_eff)
# ramp: c rises linearly 0->1 over the campaign-rollout + immunity-onset window
#       (from tau to t_eff), i.e. fully protected at t_eff.
protected_fraction <- function(Tn, tau, t_eff, shape = c("step", "ramp")) {
  shape <- match.arg(shape)
  tt <- seq_len(Tn)
  if (shape == "step") return(as.numeric(tt >= t_eff))
  width <- max(t_eff - tau, 1)
  pmin(pmax((tt - tau) / width, 0), 1)            # 0 before tau, 1 at/after t_eff
}

# --- STATIC counterfactual (manuscript Eq. 2) --------------------------------
# P = 1 - (1-eta)(1-pi*psi);  I_t^static = I_t^obs [1 - 1(t>=t_eff) P]
P_halloran <- function(pi, psi, eta) 1 - ((1 - eta) * (1 - pi * psi))

static_counterfactual <- function(incidence, t_eff, P, c_t = NULL) {
  # c_t optional: if supplied, use the same cumulative-protection shape as the
  # renewal model (so a ramp comparison is apples-to-apples). Default = step.
  if (is.null(c_t)) c_t <- as.numeric(seq_along(incidence) >= t_eff)
  incidence * (1 - c_t * P)
}

# --- RENEWAL counterfactual --------------------------------------------------
# R_t^v = [1 - c(t) pi psi_T] R_hat_t, applied from t_eff onward.
# Propagate the renewal equation using the MODIFIED history I^v so each averted
# infection lowers Lambda on later weeks -> compounding = emergent herd effect.
#
# feedback = FALSE reproduces the static mechanism: propagate using the OBSERVED
# history instead of I^v (no transmission feedback). Used by the static-identity
# unit test to document static = renewal minus herd feedback (plus fixed eta).
renewal_counterfactual <- function(incidence, w, tau, t_eff, pi, psi_T,
                                   c_shape = c("step", "ramp"),
                                   feedback = TRUE) {
  c_shape <- match.arg(c_shape)
  Tn <- length(incidence); K <- length(w)
  Rt <- reconstruct_Rt(incidence, w)$Rt
  c_t <- protected_fraction(Tn, tau, t_eff, c_shape)
  Rt_vax <- Rt * (1 - c_t * pi * psi_T)           # transmission reduction

  vax <- incidence                                # pre-intervention = observed
  for (t in seq_len(Tn)) {
    if (c_t[t] <= 0) next                          # unprotected weeks unchanged
    s_max <- min(K, t - 1L)
    hist  <- if (feedback) vax else incidence      # modified vs observed history
    lam_t <- if (s_max >= 1L) sum(w[1:s_max] * hist[t - (1:s_max)]) else 0
    vax[t] <- if (is.na(Rt_vax[t])) incidence[t] else Rt_vax[t] * lam_t
  }
  list(incidence_v = vax, Rt = Rt, Rt_vax = Rt_vax, c_t = c_t)
}

# --- ADDITIVE renewal-with-source counterfactual -----------------------------
# This is the transmission-level additive model. Observed TRUE-typhoid incidence
# T_t is decomposed INSIDE the renewal equation as
#
#   T_t = X_t + P_t,
#   P_t = R_t^P * sum_s w_s T_{t-s},
#
# where X_t is exogenous/background-source typhoid and P_t is propagated
# transmission. `source_fraction` (theta) partitions renewal-identifiable
# incidence: X_t=theta*T_t and P_t=(1-theta)*T_t when Lambda_t>0. If Lambda_t=0,
# a positive T_t cannot have arisen through the renewal term, so it is treated as
# an exogenous seed (X_t=T_t, P_t=0). This seed rule preserves the observed
# baseline exactly, including after gaps in incidence.
#
# Under vaccination q_t = 1-c(t)*coverage*psi_T. Direct protection acts on both
# acquisition routes, X_t^v=q_t*X_t and
# P_t^v=q_t*R_t^P*sum_s w_s*T_{t-s}^v. Only P_t feeds back recursively, while all
# true cases (whatever their source) contribute to later infectiousness.
additive_source_counterfactual <- function(
    incidence, w, tau, t_eff, coverage, psi_T, source_fraction,
    c_shape = c("step", "ramp")) {
  c_shape <- match.arg(c_shape)
  incidence <- as.numeric(incidence)
  Tn <- length(incidence)
  if (!Tn || any(!is.finite(incidence)) || any(incidence < 0))
    stop("incidence must be a non-negative finite vector")
  if (!is.numeric(source_fraction) ||
      !(length(source_fraction) %in% c(1L, Tn)))
    stop("source_fraction must be a scalar or have one value per incidence interval")
  theta_t <- rep(as.numeric(source_fraction), length.out = Tn)
  if (any(!is.finite(theta_t)) || any(theta_t < 0 | theta_t > 1))
    stop("source_fraction must lie in [0,1]")
  if (!is.finite(coverage) || !is.finite(psi_T) ||
      coverage < 0 || coverage > 1 || psi_T < 0 || psi_T > 1)
    stop("coverage and psi_T must lie in [0,1]")

  lambda_obs <- total_infectiousness(incidence, w)
  identifiable <- lambda_obs > 0
  source_obs <- ifelse(identifiable, theta_t * incidence, incidence)
  propagated_obs <- incidence - source_obs
  Rt_prop <- ifelse(identifiable, propagated_obs / lambda_obs, NA_real_)

  c_t <- protected_fraction(Tn, tau, t_eff, c_shape)
  q_t <- 1 - c_t * coverage * psi_T
  source_v <- source_obs
  propagated_v <- propagated_obs
  incidence_v <- incidence
  K <- length(w)

  for (t in seq_len(Tn)) {
    if (c_t[t] <= 0) next
    source_v[t] <- q_t[t] * source_obs[t]
    if (identifiable[t]) {
      s_max <- min(K, t - 1L)
      lambda_v <- if (s_max >= 1L)
        sum(w[1:s_max] * incidence_v[t - (1:s_max)]) else 0
      propagated_v[t] <- q_t[t] * Rt_prop[t] * lambda_v
    } else {
      # With no prior infectiousness the baseline policy assigns all incidence
      # to the exogenous seed, hence propagated_obs[t] is exactly zero.
      propagated_v[t] <- q_t[t] * propagated_obs[t]
    }
    incidence_v[t] <- source_v[t] + propagated_v[t]
  }

  total <- sum(incidence)
  list(
    incidence_v = incidence_v,
    source_obs = source_obs,
    propagated_obs = propagated_obs,
    source_v = source_v,
    propagated_v = propagated_v,
    Rt_prop = Rt_prop,
    lambda_obs = lambda_obs,
    identifiable = identifiable,
    seed = !identifiable & incidence > 0,
    c_t = c_t,
    q_t = q_t,
    source_fraction = theta_t,
    source_fraction_realized = if (total > 0) sum(source_obs) / total else NA_real_,
    baseline_error = max(abs(source_obs + propagated_obs - incidence)),
    counterfactual_balance_error =
      max(abs(source_v + propagated_v - incidence_v))
  )
}

# --- Impact measures ---------------------------------------------------------
# A = sum_t (I_t^obs - I_t^v);  rho = 100 * A / sum_t I_t^obs
impact_measures <- function(incidence_obs, incidence_v) {
  total <- sum(incidence_obs)
  averted <- sum(incidence_obs - incidence_v)
  list(total = total,
       averted = averted,
       pct_reduction = if (total > 0) 100 * averted / total else 0)
}

# --- Effective-IVE diagnostic (bridge to the static model) -------------------
# f_post = (sum_{t>=t_eff} I_t^obs) / (sum_t I_t^obs)
# eta_eff = 1 - [1 - rho_renewal/(100 f_post)] / (1 - pi psi_T)
# Matches TOTAL reductions; does not assume the renewal reduction factorizes.
eta_eff <- function(incidence_obs, t_eff, pct_reduction_renewal, pi, psi_T) {
  total <- sum(incidence_obs)
  f_post <- if (total > 0) sum(incidence_obs[seq_along(incidence_obs) >= t_eff]) / total else 0
  if (f_post <= 0) return(NA_real_)
  1 - (1 - (pct_reduction_renewal / 100) / f_post) / (1 - pi * psi_T)
}
