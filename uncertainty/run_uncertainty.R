# =============================================================================
# run_uncertainty.R — rigorous R_t analysis propagating GI / asymptomatic / PPV
# uncertainty through the exact reconstruction and the vaccine counterfactual.
# Run from the repo root:  Rscript uncertainty/run_uncertainty.R
# Reuses the renewal engine (renewal/R/*).
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressWarnings(if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1))
suppressMessages({ library(dplyr) })

for (f in c("gi","renewal_core","data_prep","epiestim_rt"))
  source(file.path("renewal/R", paste0(f, ".R")))
for (f in c("litparams","kernels","uncertainty","figures_unc","bayes"))
  source(file.path("uncertainty/R", paste0(f, ".R")))

rcfg <- yaml::read_yaml("renewal/config.yml")
cfg  <- yaml::read_yaml("uncertainty/config_uncertainty.yml")
set.seed(cfg$seed)
cat("=== Rigorous R_t analysis (propagated uncertainty) ===\n")
cat("MC draws:", cfg$n_mc, "| scenario tau =", cfg$scenario$tau_weeks,
    "cov =", cfg$scenario$coverage, "\n\n")
for (p in cfg$paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
tabf <- function(n) file.path(cfg$paths$tables, n)

prep <- prep_outbreaks(rcfg)
params <- sample_uncertainty_params(cfg, rcfg)
write.csv(params, tabf("mc_param_draws.csv"), row.names = FALSE)

# --- Monte Carlo propagation -------------------------------------------------
cat("Running Monte Carlo (", cfg$n_mc, "draws x", length(prep$series), "outbreaks)...\n")
mc <- run_mc(prep, rcfg, cfg, params)
ribbon <- rt_ribbon(mc)
summ <- summarise_mc(mc)
write.csv(ribbon, tabf("rt_ribbon.csv"), row.names = FALSE)
write.csv(summ$per_outbreak, tabf("tab_mc_impact.csv"), row.names = FALSE)
write.csv(summ$pooled, tabf("tab_mc_pooled.csv"), row.names = FALSE)

# --- Point-estimate baseline (base GI, no reshaping, psi_T = psi) -------------
w_base <- discretize_gi_mix(rcfg$gi$mean_days, rcfg$gi$cv, 0, 0, drop_first = rcfg$gi$drop_first_week)
te <- mc$t_eff
baseline <- do.call(rbind, lapply(names(prep$series), function(sid) {
  inc <- prep$series[[sid]]; Tn <- length(inc)
  Rt <- reconstruct_Rt(inc, w_base)$Rt
  occurred <- te <= (Tn - rcfg$tail_guard_weeks)
  rho <- if (!occurred) 0 else impact_measures(inc, renewal_counterfactual(
    inc, w_base, cfg$scenario$tau_weeks, te, cfg$scenario$coverage, rcfg$vaccine$psi_mean,
    c_shape = rcfg$timing$c_shape, feedback = TRUE)$incidence_v)$pct_reduction
  data.frame(study_id = sid, rho_base = rho, stringsAsFactors = FALSE)
}))
base_rt <- do.call(rbind, lapply(names(prep$series), function(sid) {
  inc <- prep$series[[sid]]
  data.frame(study_id = sid, week = seq_along(inc),
             rt_base = pmin(reconstruct_Rt(inc, w_base)$Rt, 6), stringsAsFactors = FALSE)
}))

# --- Global sensitivity (PRCC on pooled % reduction) -------------------------
cat("Computing PRCC global sensitivity...\n")
pd <- summ$pooled_by_draw %>% left_join(params, by = "draw")
inputs <- c("gi_mean","gi_cv","pi_c","c_c","mu_carrier","f_a","c_a",
            "ppv_level","ppv_trend","s_cult","psi")
pr <- prcc(params[match(pd$draw, params$draw), ], pd$rho, inputs = inputs)
prcc_df <- data.frame(param = names(pr), prcc = as.numeric(pr),
                      stringsAsFactors = FALSE)
prcc_df <- prcc_df[order(-abs(prcc_df$prcc)), ]
write.csv(prcc_df, tabf("tab_prcc.csv"), row.names = FALSE)

# --- Validation: PPV invariance & identifiability ----------------------------
cat("\n--- Validation ---\n")
# (1) ppv_reshape with delta=0 is the identity -> constant PPV is invariant.
inv_ppv <- max(sapply(prep$series, function(x) max(abs(ppv_reshape(x, 0) - x))))
cat(sprintf("[%s] PPV invariance: reshape(delta=0) == identity (err %.0e)\n",
            ifelse(inv_ppv < 1e-12, "PASS", "FAIL"), inv_ppv))
# (2) ppv_level and s_cult must not influence % reduction (|PRCC| ~ 0): they only
#     rescale absolute true-case numbers (non-identifiable level) — documented.
lvl_prcc <- max(abs(prcc_df$prcc[prcc_df$param %in% c("ppv_level", "s_cult")]))
cat(sprintf("[%s] Level invariance: |PRCC| of ppv_level/s_cult on %% reduction < 0.15 (max %.3f)\n",
            ifelse(lvl_prcc < 0.15, "PASS", "WARN"), lvl_prcc))

# --- Figures -----------------------------------------------------------------
cat("Building figures...\n")
save_fig_u(fig_Rt_uncertainty(prep, ribbon, base_rt, cfg), "fig_Rt_uncertainty", cfg,
           width = 13, height = 11)
save_fig_u(fig_prcc_tornado(prcc_df, cfg), "fig_prcc_tornado", cfg, width = 8, height = 5)
save_fig_u(fig_forest_uncertainty(summ$per_outbreak, baseline, cfg),
           "fig_forest_uncertainty", cfg, height = 7)

# --- Phase 6: full Bayesian renewal model (comparison) -----------------------
if (isTRUE(cfg$bayes$enable)) {
  cat("\n--- Full Bayesian renewal model (", cfg$bayes$n_chains, "chains x",
      cfg$bayes$n_iter, "iter) ---\n")
  b <- cfg$bayes
  bayes_rows <- lapply(names(prep$series), function(sid) {
    inc <- prep$series[[sid]]
    fit <- bayes_fit_outbreak(inc, cfg, rcfg, gi_lo = b$gi_lo, gi_hi = b$gi_hi,
                              n_iter = b$n_iter, n_warmup = b$n_warmup,
                              n_chains = b$n_chains, rw_sd = b$rw_sd_logRt, thin = b$thin)
    rho <- bayes_impact(inc, fit, cfg, rcfg, te)
    cat(sprintf("  %-28s rho %.1f%% (95%% CrI %.1f-%.1f); GI mean post %.1f d\n",
                sid, median(rho), quantile(rho, 0.025), quantile(rho, 0.975), median(fit$gi)))
    data.frame(study_id = sid, rho_median = median(rho),
               rho_lo = quantile(rho, 0.025), rho_hi = quantile(rho, 0.975),
               gi_mean_post = median(fit$gi), stringsAsFactors = FALSE)
  })
  bayes_per <- do.call(rbind, bayes_rows)
  write.csv(bayes_per, tabf("tab_bayes_impact.csv"), row.names = FALSE)
  save_fig_u(fig_bayes_compare(summ$per_outbreak, bayes_per, baseline, cfg),
             "fig_bayes_vs_mc", cfg, height = 7)
  cat(sprintf("Bayesian pooled %% reduction (median of per-outbreak medians): %.1f%%\n",
              median(bayes_per$rho_median)))
}

# --- Console summary ---------------------------------------------------------
cat(sprintf("\nPooled %% reduction: point estimate %.1f%% -> with full uncertainty %.1f%% (95%% CrI %.1f-%.1f)\n",
            median(baseline$rho_base), summ$pooled$rho_median,
            summ$pooled$rho_lo, summ$pooled$rho_hi))
cat("Top sensitivity drivers (|PRCC|):\n")
print(head(prcc_df, 4), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(cfg$paths$outputs, "sessionInfo.txt"))
saveRDS(list(mc = mc$impact, summ = summ, prcc = prcc_df, baseline = baseline,
             ribbon = ribbon, params = params),
        file.path(cfg$paths$outputs, "uncertainty_results.rds"))
cat("\n=== Uncertainty analysis complete ===\n")
