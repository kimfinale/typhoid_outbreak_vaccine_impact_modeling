# =============================================================================
# Guard tests for the revision outputs. Reuses the renewal engine tests for
# self-consistency / static-identity / resolution (3-5) and adds the
# revision-specific invariance (1) and impossibility (2) guards.
# Sourced by run_revision_outputs.R (objects scn_cea, cf, prep, rcfg, cfg in env).
# =============================================================================

.tn <- 0L; .tf <- 0L
chk <- function(d, ok) { .tn <<- .tn + 1L; if (!isTRUE(ok)) .tf <<- .tf + 1L
  cat(sprintf("[%s] %s\n", ifelse(isTRUE(ok), "PASS", "FAIL"), d)) }

# --- 1. Invariance: % reduction identical before/after alpha (Part D) --------
inv <- invariance_check(scn_cea, cf, cfg$confirmation$s_culture_base)
chk("Invariance: max|pct_susp - pct_true| < 1e-8", inv < 1e-8)

# --- 2. Impossibility: per-outbreak averted <= post-intervention cases -------
chk("Impossibility: averted <= post-intervention cases (all rows)",
    all(scn_cea$s_ch_averted_tot <= scn_cea$post_int_cases + 1e-6))
Pmax <- P_max_value(rcfg)
stat_med <- scn_cea %>% dplyr::filter(model == "static", vacc_cov == 0.80) %>%
  dplyr::mutate(pct = ifelse(s_ch_tot > 0, 100 * s_ch_averted_tot / s_ch_tot, 0)) %>%
  dplyr::group_by(study_id, tau) %>% dplyr::summarise(p = median(pct), .groups = "drop")
chk("Impossibility: static median % reduction <= 100*P_max per outbreak",
    all(stat_med$p <= 100 * Pmax + 1e-6))

# --- 3. Self-consistency: renewal with psi=0 reproduces observed --------------
w0 <- gi_from_config(rcfg)
sc <- max(vapply(names(prep$series), function(sid) {
  inc <- prep$series[[sid]]
  cfv <- renewal_counterfactual(inc, w0, tau = 8, t_eff = 8, pi = 0.8, psi_T = 0,
                                c_shape = rcfg$timing$c_shape, feedback = TRUE)
  max(abs(cfv$incidence_v - inc))
}, numeric(1)))
chk("Self-consistency: renewal(psi=0) reproduces observed (<1e-8)", sc < 1e-8)

# --- 4. Static identity: renewal(feedback off, eta=0) == static(eta=0) --------
inc <- prep$series[[names(prep$series)[1]]]
ren_nofb <- renewal_counterfactual(inc, w0, tau = 8, t_eff = 8, pi = 0.8, psi_T = 0.83,
                                   c_shape = "step", feedback = FALSE)$incidence_v
stat0 <- static_counterfactual(inc, 8, P_halloran(0.8, 0.83, 0))
chk("Static identity: renewal(feedback off, eta=0) == static(eta=0)",
    max(abs(ren_nofb - stat0)) < 1e-10)

# --- 5. Resolution: only mu_g/Delta >= 2 outbreaks included -------------------
res <- prep$resolution
chk("Resolution: all included have mu_g/Delta >= 2",
    all(res$mu_g_over_delta[res$included == "Y"] >= rcfg$gi$resolution_min_ratio - 1e-9))

cat(sprintf("\n%d guard tests, %d failures\n", .tn, .tf))
if (.tf > 0L) stop("guard tests failed")
