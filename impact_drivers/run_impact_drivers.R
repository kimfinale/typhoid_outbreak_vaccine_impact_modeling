# =============================================================================
# run_impact_drivers.R — which outbreak/environment characteristics (knowable
# AT vaccination) predict the magnitude of ORI impact.
# Run from the repo root:  Rscript impact_drivers/run_impact_drivers.R
#
# PHASE 1 (this file): assemble the per-(outbreak, delay) feature x outcome
# dataset. Phases 2-5 (empirical association, simulation sweep, interaction,
# decision tool) build on tab_feature_outcomes.csv.
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressWarnings(if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1))
suppressMessages({ library(dplyr); library(tidyr) })

source("R/2-functions.R")
for (f in c("gi","renewal_core","data_prep","scenario","cost_daly","summarise"))
  source(file.path("renewal/R", paste0(f, ".R")))
for (f in c("features","associate","simulate","decision","figures_drivers"))
  source(file.path("impact_drivers/R", paste0(f, ".R")))

rcfg <- yaml::read_yaml("renewal/config.yml")
cfg  <- yaml::read_yaml("impact_drivers/config.yml")
set.seed(cfg$seed)
cat("=== Impact drivers — Phase 1 (feature x outcome dataset) ===\n")
for (p in cfg$paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
tabf <- function(n) file.path(cfg$paths$tables, n)

prep <- prep_outbreaks(rcfg)
amr_props <- load_amr_props(rcfg)
cost_env  <- setup_cost_env(rcfg); have_cost <- !is.null(cost_env)
cat("outbreaks:", length(prep$series), "| cost arm:", ifelse(have_cost, "ON", "OFF"), "\n")

# --- Prospective predictors at each delay -----------------------------------
feat <- compute_features(prep, rcfg, cfg)

# --- Impact outcomes per (outbreak, delay): renewal, coverage 0.80 ----------
cat("Running scenarios over the delay grid...\n")
scn_raw <- run_scenarios(prep, rcfg, amr_props, coverage = cfg$coverage,
                         delays = cfg$delays_weeks)
scn_cea <- if (have_cost) add_cost_daly(scn_raw, cost_env) else scn_raw
summ <- summarise_draws(scn_cea)

# amplification ratio (renewal / static cases averted) per outbreak x delay
amp <- summ %>% select(study_id, model, tau, vacc_cov, s_ch_averted_median) %>%
  pivot_wider(names_from = model, values_from = s_ch_averted_median) %>%
  mutate(amplification = ifelse(static > 0, renewal / static, NA_real_)) %>%
  select(study_id, tau, amplification)

gdp_by <- scn_cea %>% group_by(study_id) %>%
  summarise(gdp = suppressWarnings(first(gdp[is.finite(gdp)])), .groups = "drop")

outcomes <- summ %>% filter(model == "renewal", vacc_cov == cfg$coverage) %>%
  transmute(study_id, tau,
            case_averted_per_1000 = case_averted_per_1000_median,
            pct_reduction = pct_reduction_median,
            cases_averted = s_ch_averted_median,
            cost_per_daly = cost_per_daly_averted_median,
            daly_averted = daly_averted_median,
            ori_occurred) %>%
  left_join(amp, by = c("study_id", "tau"))

# --- Join features + outcomes + GDP -----------------------------------------
dat <- feat %>%
  left_join(outcomes, by = c("study_id", "tau")) %>%
  left_join(gdp_by, by = "study_id")
write.csv(dat, tabf("tab_feature_outcomes.csv"), row.names = FALSE)

cat(sprintf("Wrote %s (%d rows)\n", tabf("tab_feature_outcomes.csv"), nrow(dat)))

# =============================================================================
# PHASE 2 — empirical association (real outbreaks; exploratory)
# =============================================================================
cat("\n=== Phase 2: empirical association ===\n")
PREDS <- c("Rt_at_tau", "r_early", "doubling_time", "Rt_mean_pre", "Rt_slope_pre",
           "cum_cases_by_tau", "cum_AR_by_tau", "incidence_at_tau", "log_pop")
OUTS  <- c("case_averted_per_1000", "pct_reduction", "cases_averted")
assoc <- associate_screen(dat, PREDS, OUTS, tau = 8)
write.csv(assoc, tabf("tab_associations.csv"), row.names = FALSE)
med <- do.call(rbind, lapply(OUTS, function(o) mediation_check(dat, o, tau = 8)))
write.csv(med, tabf("tab_mediation.csv"), row.names = FALSE)
abd <- associate_by_delay(dat, c("Rt_at_tau", "r_early", "cum_AR_by_tau"),
                          "pct_reduction", cfg$delays_weeks)
cat("Top empirical predictor per outcome (tau=8):\n")
print(assoc %>% group_by(outcome) %>% slice(1) %>%
      select(outcome, predictor, rho, ci_lo, ci_hi, n) %>% as.data.frame(), row.names = FALSE)

# =============================================================================
# PHASE 3 — simulation experiment + response surfaces
# =============================================================================
cat("\n=== Phase 3: simulation sweep ===\n")
sim <- run_simulation(rcfg, cfg, n_sim = 1500)
write.csv(sim, tabf("sim_outbreaks.csv"), row.names = FALSE)
cat("Simulated", length(unique(sim$sim_id)), "outbreaks x", length(cfg$delays_weeks), "delays =",
    nrow(sim), "rows\n")
sim_imp <- do.call(rbind, lapply(OUTS, function(o) surface_importance(sim, o)))
write.csv(sim_imp, tabf("tab_sim_importance.csv"), row.names = FALSE)
cat("Top simulation driver per outcome:\n")
print(sim_imp %>% group_by(outcome) %>% slice(1) %>% as.data.frame(), row.names = FALSE)

# =============================================================================
# PHASE 4 — characteristic x delay interaction
# =============================================================================
cat("\n=== Phase 4: delay x growth interaction ===\n")
ds <- delay_sensitivity(sim, "pct_reduction")
write.csv(ds, tabf("tab_delay_sensitivity.csv"), row.names = FALSE)
print(as.data.frame(ds), row.names = FALSE)

# =============================================================================
# PHASE 5 — impact map, decision rule, validation
# =============================================================================
cat("\n=== Phase 5: impact map + decision rule ===\n")
surf <- fit_surface_model(sim, "case_averted_per_1000")
pred_real <- predict_real_from_sim(surf, dat, "case_averted_per_1000")
write.csv(pred_real, tabf("tab_validation.csv"), row.names = FALSE)
# LOO for the primary outcome with its strongest empirical predictor, plus the
# % reduction link to cumulative burden (its strongest empirical predictor).
loo <- loo_cv_real(dat, "case_averted_per_1000", "Rt_at_tau", tau = 8)
loo_pct <- loo_cv_real(dat, "pct_reduction", "cum_AR_by_tau", tau = 8)
cat(sprintf("LOO-CV cases/1000 ~ R_t at vaccination (n=%d): R2=%.2f, cor=%.2f\n",
            loo$n, loo$loo_r2, loo$loo_cor))
cat(sprintf("LOO-CV %% reduction ~ cumulative AR by tau (n=%d): R2=%.2f, cor=%.2f\n",
            loo_pct$n, loo_pct$loo_r2, loo_pct$loo_cor))
tree <- fit_decision_tree(sim, "case_averted_per_1000")
rules <- capture.output(tree)
writeLines(rules, tabf("decision_rule.txt"))

# --- Figures -----------------------------------------------------------------
cat("Building figures...\n")
for (o in OUTS)
  save_fig_d(fig_predictor_ranking(assoc, sim_imp, o, cfg),
             paste0("fig_ranking_", o), cfg, width = 9, height = 5)
save_fig_d(fig_assoc_by_delay(abd, cfg), "fig_assoc_by_delay", cfg)
save_fig_d(fig_impact_surface(surf, dat, "case_averted_per_1000", cfg), "fig_impact_surface", cfg)
save_fig_d(fig_interaction(ds, cfg), "fig_interaction", cfg)
save_fig_d(fig_validation(pred_real, "case_averted_per_1000", cfg), "fig_validation", cfg)
if (requireNamespace("rpart.plot", quietly = TRUE)) tryCatch({
  png(file.path(cfg$paths$figures, "fig_decision_tree.png"), width = 1800, height = 1100, res = 200)
  rpart.plot::rpart.plot(tree, roundint = FALSE,
                         main = "Decision rule for cases averted / 1,000 doses")
  dev.off()
}, error = function(e) message("tree plot skipped: ", conditionMessage(e)))

writeLines(capture.output(sessionInfo()), file.path(cfg$paths$outputs, "sessionInfo.txt"))
saveRDS(list(dat = dat, assoc = assoc, med = med, abd = abd, sim_imp = sim_imp,
             ds = ds, surf = surf, pred_real = pred_real, loo = loo, loo_pct = loo_pct,
             tree = tree),
        file.path(cfg$paths$outputs, "impact_drivers.rds"))

# --- Render report -----------------------------------------------------------
if (nzchar(Sys.which("quarto"))) {
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  rc <- system2("quarto", c("render", shQuote("impact_drivers/report.qmd")),
                stdout = FALSE, stderr = FALSE)
  cat(if (rc == 0) "Wrote impact_drivers/report.html\n" else "quarto render skipped/failed\n")
}
cat("=== Impact-drivers analysis complete ===\n")
