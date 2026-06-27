# =============================================================================
# Kernel + observation-model machinery for the rigorous R_t analysis.
# Reuses renewal/R/gi.R::discretize_gi. Pure functions.
# =============================================================================

# Mixture generation-interval kernel: acute symptomatic gamma + a long
# chronic-carrier tail. The carrier transmission weight omega is the carrier
# share of total infectiousness, omega = (pi_c c_c mu_carrier)/(mu_g + pi_c c_c mu_carrier)
# (prevalence x per-time infectiousness x duration). Returns a normalized weekly pmf.
discretize_gi_mix <- function(gi_mean, gi_cv, pi_c = 0, c_c = 0, mu_carrier = 180,
                              step_days = 7, max_steps = 40, drop_first = TRUE) {
  w_ac <- discretize_gi(mean_days = gi_mean, sd_days = gi_mean * gi_cv,
                        step_days = step_days, max_steps = max_steps, drop_first = drop_first)
  omega <- (pi_c * c_c * mu_carrier) / (gi_mean + pi_c * c_c * mu_carrier)
  if (omega <= 0) return(w_ac)
  w_ca <- discretize_gi(mean_days = mu_carrier, sd_days = mu_carrier * 0.5,
                        step_days = step_days, max_steps = max_steps, drop_first = drop_first)
  w <- (1 - omega) * w_ac + omega * w_ca
  w / sum(w)
}

# Time-varying PPV reshaping of the observed incidence. The true-typhoid shape is
# I_obs * (1 + delta*((t-0.5)/T - 0.5)); a constant PPV (delta = 0) leaves the
# series unchanged, so R_t and percent reduction are invariant (tested).
ppv_reshape <- function(incidence, delta = 0) {
  if (delta == 0) return(incidence)
  Tn <- length(incidence)
  factor <- 1 + delta * (((seq_len(Tn) - 0.5) / Tn) - 0.5)
  pmax(incidence * factor, 0)
}

# Transmission-relevant efficacy from the asymptomatic split:
# phi_a = fraction of transmission that is asymptomatic/subclinical (vaccine
# blocks disease, not necessarily this); psi_T = psi * (1 - phi_a).
psi_T_from_asymptomatic <- function(psi, f_a, c_a) {
  phi_a <- (f_a * c_a) / (f_a * c_a + (1 - f_a))
  list(psi_T = psi * (1 - phi_a), phi_a = phi_a)
}
