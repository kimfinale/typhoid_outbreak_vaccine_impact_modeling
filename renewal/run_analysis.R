# =============================================================================
# run_analysis.R — top-level reproducible driver for the renewal-equation
# (R_t) ORI analysis. Run from the repo root:
#     Rscript renewal/run_analysis.R
#
# PHASE 1 (this file): config, data prep, resolution table, self-consistency
# table, validation tests, sessionInfo. Figures / policy grid / Quarto report
# are added in Phase 2 once these pass.
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", ".")
setwd(ROOT)
Sys.setenv(RENEWAL_ROOT = ROOT)

suppressWarnings({
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
})

source("R/2-functions.R")            # REUSE: compute_yld/yll, add_cea_results
source("renewal/R/gi.R")
source("renewal/R/renewal_core.R")
source("renewal/R/data_prep.R")
source("renewal/R/scenario.R")
source("renewal/R/cost_daly.R")
source("renewal/R/summarise.R")
source("renewal/R/epiestim_rt.R")
source("renewal/R/figures.R")
source("renewal/R/ppv.R")            # PPV (pi) posterior propagation

cfg <- yaml::read_yaml("renewal/config.yml")
ppv_draws_override <- Sys.getenv("PPV_DRAWS_OVERRIDE", "")
tables_override <- Sys.getenv("RENEWAL_TABLES_OVERRIDE", "")
figures_override <- Sys.getenv("RENEWAL_FIGURES_OVERRIDE", "")
outputs_override <- Sys.getenv("RENEWAL_OUTPUTS_OVERRIDE", "")
if (nzchar(ppv_draws_override)) cfg$ppv$draws <- ppv_draws_override
if (nzchar(tables_override)) cfg$paths$tables <- tables_override
if (nzchar(figures_override)) cfg$paths$figures <- figures_override
if (nzchar(outputs_override)) cfg$paths$outputs <- outputs_override
set.seed(cfg$seed)
cat("=== Renewal analysis — Phase 1 ===\n")
cat("Seed:", cfg$seed, "| Sobol n:", cfg$n_sobol,
    "| GI mean:", cfg$gi$mean_days, "d (CV", cfg$gi$cv, ")",
    "| drop_first_week:", cfg$gi$drop_first_week, "\n\n")

# --- PPV (pi) posterior: suspected -> TRUE typhoid impact --------------------
pi_post <- if (isTRUE(cfg$ppv$enable)) load_pi_posterior(cfg$ppv$draws, cfg$ppv$anchor) else NULL
if (!is.null(pi_post)) {
  cat(sprintf("PPV propagation: ENABLED (%d posterior draws; anchored: %s).\n",
              pi_post$ndraw, paste(pi_post$anchor, collapse = ", ")))
  cat(sprintf("  community pi typical = %.2f [%.2f, %.2f] (median inv_logit(mu) across draws)\n\n",
              plogis(median(pi_post$mu_pi)),
              plogis(quantile(pi_post$mu_pi, .025)),
              plogis(quantile(pi_post$mu_pi, .975))))
} else cat("PPV propagation: DISABLED (raw suspected-case impact).\n\n")

dir.create(cfg$paths$tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$outputs, recursive = TRUE, showWarnings = FALSE)

# --- Data prep + resolution rule --------------------------------------------
prep <- prep_outbreaks(cfg)
write.csv(prep$resolution, file.path(cfg$paths$tables, "tab_resolution.csv"), row.names = FALSE)
n_incl <- sum(prep$resolution$included == "Y")
cat(sprintf("Eligible outbreaks: %d | included by mu_g/Delta>=%.1f: %d\n",
            nrow(prep$resolution), cfg$gi$resolution_min_ratio, n_incl))
cat("Wrote", file.path(cfg$paths$tables, "tab_resolution.csv"), "\n\n")
stopifnot(n_incl == length(prep$series))

# --- Self-consistency: psi_T = 0 must reproduce observed curve ----------------
w_base <- gi_from_config(cfg)
sc <- do.call(rbind, lapply(names(prep$series), function(sid) {
  inc <- prep$series[[sid]]
  iw  <- prep$meta[[sid]]$intervention_week
  t_eff <- if (is.na(iw)) cfg$timing$delay_base_weeks else max(round(iw), 1)
  cf <- renewal_counterfactual(inc, w_base, tau = t_eff, t_eff = t_eff,
                               pi = cfg$vaccine$coverage_base, psi_T = 0,
                               c_shape = cfg$timing$c_shape, feedback = TRUE)
  data.frame(study_id = sid, n_weeks = length(inc),
             max_abs_diff = max(abs(cf$incidence_v - inc)), stringsAsFactors = FALSE)
}))
write.csv(sc, file.path(cfg$paths$tables, "tab_selfconsistency.csv"), row.names = FALSE)
cat("Self-consistency: worst max|I^v(psi=0)-I^obs| =", format(max(sc$max_abs_diff)), "\n")
cat("Wrote", file.path(cfg$paths$tables, "tab_selfconsistency.csv"), "\n\n")
stopifnot(max(sc$max_abs_diff) < 1e-8)

# --- Validation tests --------------------------------------------------------
cat("--- Running validation tests ---\n")
test_status <- tryCatch({ source("renewal/tests/test_renewal.R", local = new.env()); 0L },
                        error = function(e) { cat("TEST ERROR:", conditionMessage(e), "\n"); 1L })

if (!identical(test_status, 0L)) { cat("Validation tests failed; stopping.\n"); quit(status = 1L) }

# =============================================================================
# PHASE 2 — Sobol scenarios, CEA, diagnostics, figures
# =============================================================================
cat("\n=== Phase 2 — scenarios + CEA + figures ===\n")
amr_props <- load_amr_props(cfg)
cost_env  <- setup_cost_env(cfg)          # NULL if cost inputs missing
have_cost <- !is.null(cost_env)
cat("Cost/DALY arm:", ifelse(have_cost, "ENABLED", "SKIPPED (missing inputs)"), "\n")
tab <- function(n) file.path(cfg$paths$tables, n)

maybe_cea <- function(df) if (have_cost) add_cost_daly(df, cost_env) else df

# --- Analysis 2: delay x coverage grid (base GI, psi_T = psi) -----------------
cat("Running delay x coverage grid (both models)...\n")
grid_raw  <- run_scenarios(prep, cfg, amr_props,
                           coverage = cfg$vaccine$coverage_grid,
                           delays   = cfg$timing$delay_grid_weeks,
                           pi_post  = pi_post)
grid_cea  <- maybe_cea(grid_raw)
summ_grid <- summarise_draws(grid_cea)
write.csv(summ_grid, tab("summary_delay_coverage.csv"), row.names = FALSE)

# --- PPV effect: base-case impact, suspected-naive vs TRUE typhoid (+ dilution) --
if (!is.null(pi_post)) {
  cat(sprintf("Computing PPV effect (regime=%s, base case, renewal)...\n", cfg$ppv$regime))
  bc <- function(df) df[df$model == "renewal" &
                        df$tau == cfg$timing$delay_base_weeks &
                        df$vacc_cov == cfg$vaccine$coverage_base, ]
  gt <- bc(grid_raw)                                               # PPV-adjusted (true typhoid)
  gs <- bc(run_scenarios(prep, cfg, amr_props, coverage = cfg$vaccine$coverage_base,
                         delays = cfg$timing$delay_base_weeks, pi_post = NULL))  # suspected-naive
  perstudy_sum <- function(df) sum(tapply(df$s_ch_averted_tot, df$study_id, median, na.rm = TRUE))
  ppv_cmp <- data.frame(
    quantity = c("total_cases_averted (per-study median, summed)",
                 "pooled_TRUE_pct_reduction_median",
                 "pooled_OBSERVED_pct_reduction_median (suspected surveillance)"),
    suspected_naive = c(round(perstudy_sum(gs)),
                        round(median(100 * gs$s_ch_averted_tot / gs$s_ch_tot, na.rm = TRUE), 1),
                        NA_real_),
    true_typhoid = c(round(perstudy_sum(gt)),
                     round(median(100 * gt$s_ch_averted_tot / gt$s_ch_tot, na.rm = TRUE), 1),
                     round(median(gt$obs_pct_reduction, na.rm = TRUE), 1)))
  write.csv(ppv_cmp, tab("tab_ppv_effect.csv"), row.names = FALSE)
  cat("Wrote", tab("tab_ppv_effect.csv"),
      "-- true typhoid cases averted (pi-scaled) + surveillance dilution (additive regime)\n")
}

# --- Analysis 1: base case (tau = 8, pi = 0.80) ------------------------------
summ_base <- summ_grid %>%
  dplyr::filter(tau == cfg$timing$delay_base_weeks, vacc_cov == cfg$vaccine$coverage_base)
pooled    <- pooled_estimates(summ_grid)
write.csv(pooled, tab("pooled_estimates.csv"), row.names = FALSE)

# --- tab_impact_summary (mirror Supp Table 2 w/ renewal column) --------------
impact_summary <- summ_grid %>%
  dplyr::filter(vacc_cov %in% cfg$vaccine$coverage_grid) %>%
  dplyr::group_by(model, tau, vacc_cov) %>%
  dplyr::summarise(
    cases_averted_median = median(s_ch_averted_median, na.rm = TRUE),
    cases_averted_iqr_low = quantile(s_ch_averted_median, 0.25, na.rm = TRUE),
    cases_averted_iqr_high = quantile(s_ch_averted_median, 0.75, na.rm = TRUE),
    pct_reduction_median = median(pct_reduction_median, na.rm = TRUE),
    cost_per_daly_median = median(cost_per_daly_averted_median, na.rm = TRUE),
    .groups = "drop")
write.csv(impact_summary, tab("tab_impact_summary.csv"), row.names = FALSE)

# --- eta_eff diagnostic ------------------------------------------------------
cat("Computing eta_eff...\n")
eta_tab <- eta_eff_table(prep, cfg, amr_props)
write.csv(eta_tab, tab("tab_eta_eff.csv"), row.names = FALSE)

# --- Analysis 3: GI sensitivity (mu_g in 7-28 d), base timing/coverage -------
cat("Running GI sensitivity...\n")
gi_sens <- do.call(rbind, lapply(cfg$gi$mean_days_sweep, function(mg) {
  r <- run_scenarios(prep, cfg, amr_props, coverage = cfg$vaccine$coverage_base,
                     delays = cfg$timing$delay_base_weeks, gi_mean_days = mg)
  s <- summarise_draws(maybe_cea(r)); s$gi_mean <- mg
  s[s$model == "renewal", ]
}))
write.csv(gi_sens, tab("sens_gi.csv"), row.names = FALSE)

# --- Analysis 4: psi_T sensitivity (0.6*psi - psi) ---------------------------
cat("Running psi_T sensitivity...\n")
psiT_fracs <- seq(cfg$vaccine$psiT_frac_lower, 1, length.out = 5)
psiT_sens <- do.call(rbind, lapply(psiT_fracs, function(fr) {
  r <- run_scenarios(prep, cfg, amr_props, coverage = cfg$vaccine$coverage_base,
                     delays = cfg$timing$delay_base_weeks, psi_T_frac = fr)
  s <- summarise_draws(maybe_cea(r)); s$psiT_frac <- fr
  s[s$model == "renewal", ]
}))
write.csv(psiT_sens, tab("sens_psiT.csv"), row.names = FALSE)

# --- Figures -----------------------------------------------------------------
cat("Building figures...\n")
save_fig(fig_Rt_panels(prep, cfg), "fig_Rt_panels", cfg, width = 10, height = 10)
save_fig(fig_Rt_with_incidence(prep, cfg), "fig_Rt_with_incidence", cfg, width = 13, height = 11)
save_fig(fig_forest_pctreduction(summ_base, pooled, cfg), "fig_forest_pctreduction", cfg, height = 7)
save_fig(fig_amplification_vs_delay(summ_grid, cfg), "fig_amplification_vs_delay", cfg)
save_fig(fig_eta_eff(eta_tab, cfg), "fig_eta_eff", cfg, height = 7)
if (have_cost) save_fig(fig_policy_grid(summ_grid, cfg), "fig_policy_grid", cfg, width = 10, height = 7)
save_fig(fig_sensitivity(gi_sens, "gi_mean", "GI mean (days)", "GI sensitivity (renewal)"),
         "fig_gi_sensitivity", cfg)
save_fig(fig_sensitivity(psiT_sens, "psiT_frac", "psi_T / psi", "Transmission-VE sensitivity (renewal)"),
         "fig_psiT_sensitivity", cfg)
# two representative outbreaks: largest (cases) and longest (weeks)
sizes <- vapply(prep$series, sum, numeric(1)); lens <- vapply(prep$series, length, numeric(1))
ex_ids <- unique(c(names(which.max(sizes)), names(which.max(lens))))
save_fig(fig_curves_examples(prep, cfg, ex_ids), "fig_curves_examples", cfg, width = 10)

# --- Reproducibility record --------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(cfg$paths$outputs, "sessionInfo.txt"))
saveRDS(list(summ_grid = summ_grid, pooled = pooled, eta_tab = eta_tab,
             gi_sens = gi_sens, psiT_sens = psiT_sens, have_cost = have_cost),
        file.path(cfg$paths$outputs, "results.rds"))
cat("\nWrote tables to", cfg$paths$tables, "and figures to", cfg$paths$figures, "\n")

# --- Render the Quarto report (if quarto is available) -----------------------
if (nzchar(Sys.which("quarto"))) {
  cat("Rendering report.qmd...\n")
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  rc <- system2("quarto", c("render", shQuote("renewal/report.qmd")),
                stdout = FALSE, stderr = FALSE)
  cat(if (rc == 0) "Wrote renewal/report.html\n" else "quarto render failed (non-fatal)\n")
} else cat("quarto not found; skipping report render.\n")

cat("=== Phase 2 complete ===\n")
