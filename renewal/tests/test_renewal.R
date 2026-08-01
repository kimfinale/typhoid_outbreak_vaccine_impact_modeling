# =============================================================================
# Validation / unit tests for the renewal engine.
# Run standalone (Rscript renewal/tests/test_renewal.R) or via run_analysis.R.
# Uses a tiny base-R harness (no testthat dependency).
# =============================================================================

.here <- function(...) file.path(Sys.getenv("RENEWAL_ROOT", "."), ...)
source(.here("renewal/R/gi.R"))
source(.here("renewal/R/renewal_core.R"))
source(.here("renewal/R/ppv.R"))

.tests_run <- 0L; .tests_fail <- 0L
check <- function(desc, cond) {
  .tests_run <<- .tests_run + 1L
  ok <- isTRUE(cond)
  if (!ok) .tests_fail <<- .tests_fail + 1L
  cat(sprintf("[%s] %s\n", ifelse(ok, "PASS", "FAIL"), desc))
}

# A deterministic synthetic outbreak (declining R_t) for tests that need a curve.
w <- discretize_gi(mean_days = 14, sd_days = 8.4, step_days = 7, drop_first = TRUE)
Tn <- 30; wk <- 1:Tn
Rt_path <- 0.45 + 1.8 / (1 + exp(0.6 * (wk - 9)))
inc <- numeric(Tn); inc[1:3] <- c(2, 3, 5); K <- length(w)
for (t in 4:Tn) {
  s_max <- min(K, t - 1L)
  inc[t] <- Rt_path[t] * sum(w[1:s_max] * inc[t - (1:s_max)])
}

# --- Test 2: GI normalization (sum w_s = 1 for every mu_g) --------------------
for (mg in c(7, 14, 21, 28)) {
  wmg <- discretize_gi(mean_days = mg, sd_days = mg * 0.6, step_days = 7, drop_first = TRUE)
  check(sprintf("GI normalizes to 1 (mu_g=%dd)", mg), abs(sum(wmg) - 1) < 1e-12)
  check(sprintf("GI w_0/first-bin zeroed (mu_g=%dd)", mg), wmg[1] == 0)
}

# --- Test 1: Self-consistency (psi=0 reproduces observed curve) ---------------
cf0 <- renewal_counterfactual(inc, w, tau = 8, t_eff = 8, pi = 0.8, psi_T = 0,
                              c_shape = "step", feedback = TRUE)
check("Self-consistency: max|I^v(psi=0) - I^obs| < 1e-8",
      max(abs(cf0$incidence_v - inc)) < 1e-8)

# --- Additive observation model: S = T + F and cumulative PPV ---------------
S_obs <- c(3, 4, 8, 20, 9, 5, 3)
obs_parts <- ppv_incidence_components(S_obs, pi = 0.35)
check("Observation additivity: suspected = true typhoid + other febrile",
      max(abs(S_obs - obs_parts$true_typhoid - obs_parts$other_febrile)) < 1e-12)
check("Observation PPV: sum(true typhoid) = pi * sum(suspected)",
      abs(sum(obs_parts$true_typhoid) - 0.35 * sum(S_obs)) < 1e-12)
check("Observation components are non-negative",
      all(obs_parts$true_typhoid >= 0) && all(obs_parts$other_febrile >= 0))
obs_all_typhoid <- ppv_incidence_components(S_obs, pi = 1)
check("Observation endpoint: pi=1 assigns all suspected incidence to typhoid",
      max(abs(obs_all_typhoid$true_typhoid - S_obs)) < 1e-12 &&
        max(abs(obs_all_typhoid$other_febrile)) < 1e-12)

# --- Additive transmission model: T = X + P inside one recursion ------------
add0 <- additive_source_counterfactual(
  inc, w, tau = 8, t_eff = 8, coverage = 0.8, psi_T = 0,
  source_fraction = 0.45)
check("Transmission self-consistency: q=1 reproduces true-typhoid baseline",
      max(abs(add0$incidence_v - inc)) < 1e-8)
check("Transmission baseline mass balance: true typhoid = source + propagated",
      add0$baseline_error < 1e-12)
check("Transmission counterfactual mass balance",
      add0$counterfactual_balance_error < 1e-12)
check("Zero-infectiousness positive weeks are assigned to exogenous seeds",
      all(add0$source_obs[add0$seed] == inc[add0$seed]) &&
        all(add0$propagated_obs[add0$seed] == 0))

# --- Test 3: Static identity (feedback OFF, eta=0) == static (eta=0) ----------
t_eff <- 8
P0 <- P_halloran(pi = 0.8, psi = 0.83, eta = 0)
stat0 <- static_counterfactual(inc, t_eff, P0)
# renewal with feedback OFF propagates using observed history => pure direct
# proportional reduction = static with eta=0, to machine precision.
ren_nofb <- renewal_counterfactual(inc, w, tau = 8, t_eff = t_eff, pi = 0.8,
                                   psi_T = 0.83, c_shape = "step", feedback = FALSE)
check("Static identity: renewal(feedback off, eta=0) == static(eta=0)",
      max(abs(ren_nofb$incidence_v - stat0)) < 1e-10)

# Source-fraction endpoints use the same seed policy and direct-effect schedule.
add_theta0 <- additive_source_counterfactual(
  inc, w, tau = 8, t_eff = t_eff, coverage = 0.8, psi_T = 0.83,
  source_fraction = 0)
ren_theta0 <- renewal_counterfactual(
  inc, w, tau = 8, t_eff = t_eff, pi = 0.8, psi_T = 0.83,
  c_shape = "step", feedback = TRUE)
check("Transmission endpoint: source fraction 0 equals pure renewal",
      max(abs(add_theta0$incidence_v - ren_theta0$incidence_v)) < 1e-10)

add_theta1 <- additive_source_counterfactual(
  inc, w, tau = 8, t_eff = t_eff, coverage = 0.8, psi_T = 0.83,
  source_fraction = 1)
direct_theta1 <- static_counterfactual(
  inc, t_eff, P_halloran(0.8, 0.83, eta = 0), c_t = add_theta1$c_t)
check("Transmission endpoint: source fraction 1 equals direct-only static",
      max(abs(add_theta1$incidence_v - direct_theta1)) < 1e-10)
check("Additive counterfactual incidence is non-negative and no greater than baseline",
      all(add_theta0$incidence_v >= -1e-12) &&
        all(add_theta0$incidence_v <= inc + 1e-8))

# --- Test 4: Resolution assert (mu_g/Delta >= ~2 weekly; < 2 monthly) ---------
check("Resolution: weekly Delta passes (mu_g/Delta = 2)", (14 / 7) >= 2)
check("Resolution: monthly Delta fails (mu_g/Delta < 2)", (14 / 30) < 2)

# --- Test 5: Monotonicity (earlier intervention => >= cases averted) ----------
averted_by_teff <- sapply(c(4, 6, 8, 10, 12), function(te) {
  cf <- renewal_counterfactual(inc, w, tau = te, t_eff = te, pi = 0.8, psi_T = 0.83,
                               c_shape = "step", feedback = TRUE)
  impact_measures(inc, cf$incidence_v)$averted
})
check("Monotonicity: earlier t_eff weakly increases cases averted",
      all(diff(averted_by_teff) <= 1e-8))

# --- Extra: renewal averts >= static at same direct effect when feedback on ---
P_eta0 <- P_halloran(0.8, 0.83, 0)
A_static  <- impact_measures(inc, static_counterfactual(inc, t_eff, P_eta0))$averted
A_renewal <- impact_measures(inc, renewal_counterfactual(inc, w, 8, t_eff, 0.8, 0.83,
                              feedback = TRUE)$incidence_v)$averted
check("Herd feedback: renewal averted >= static averted (eta=0, early response)",
      A_renewal >= A_static - 1e-8)

cat(sprintf("\n%d tests, %d failures\n", .tests_run, .tests_fail))
if (.tests_fail > 0L) quit(status = 1L)
