# =============================================================================
# run_revision_outputs.R — reproduces the manuscript-revision deliverables
# (Parts A-G) and the consolidated revision manifest. Run from the repo root:
#     Rscript revision/run_revision_outputs.R
#
# Reuses the renewal engine (renewal/R/*) and the manuscript CEA functions
# (R/2-functions.R). Static and renewal share identical inputs and downstream
# cost/DALY machinery; only cases-averted differs.
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressWarnings(if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1))

source("R/2-functions.R")
for (f in c("gi","renewal_core","data_prep","scenario","cost_daly","summarise","figures","epiestim_rt"))
  source(file.path("renewal/R", paste0(f, ".R")))
source("revision/R/parts.R")

rcfg <- yaml::read_yaml("renewal/config.yml")          # base params (Table 1)
cfg  <- yaml::read_yaml("revision/config_revision.yml") # revision params
set.seed(cfg$seed)
cat("=== Manuscript revision outputs ===\nSeed:", cfg$seed, "\n\n")
for (p in cfg$paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
tabf <- function(n) file.path(cfg$paths$tables, n)
figf <- function(n) file.path(cfg$paths$figures, n)

# --- Shared inputs -----------------------------------------------------------
raw_summary <- read_summary(rcfg$paths$summary)   # normalized (headers + study-id crosswalk)
prep <- prep_outbreaks(rcfg)
amr_props <- load_amr_props(rcfg)
cost_env  <- setup_cost_env(rcfg)
have_cost <- !is.null(cost_env)
cat("Cost/DALY arm:", ifelse(have_cost, "ENABLED", "SKIPPED"), "| outbreaks:",
    length(prep$series), "\n")

# One full delay x coverage grid drives every Part (B filters to 1&8 / 80%,
# C uses the whole grid, D uses 1&8, F uses 8). CEA attached once.
cat("Running scenarios (full delay grid x coverage grid)...\n")
scn_raw <- run_scenarios(prep, rcfg, amr_props,
                         coverage = rcfg$vaccine$coverage_grid,
                         delays = rcfg$timing$delay_grid_weeks)
scn_cea <- if (have_cost) add_cost_daly(scn_raw, cost_env) else scn_raw

manifest <- list()
add_row <- function(comment_id, section, quantity, value, interval, units, source_file, notes)
  manifest[[length(manifest) + 1]] <<- data.frame(comment_id, section, quantity, value,
    interval, units, source_file, notes, stringsAsFactors = FALSE)

# --- Part A ------------------------------------------------------------------
cat("Part A: reconciliations (c0, c5, c6)...\n")
tA <- part_A_reconciliations(raw_summary, prep, cfg, rcfg)
write.csv(tA, tabf("tab_reconciliations.csv"), row.names = FALSE)
add_row("c0", "Methods/Figure 1", "Excluded outbreaks breakdown",
        "20/1/1/2 = 24", "-", "outbreaks", "tab_reconciliations.csv",
        "Submitted 17/2/1/2=22 inconsistent; reasons not in dataset - confirm vs Figure 1")
add_row("c5", "Discussion", "Outbreak span & per-year denominator",
        tA$value[tA$item == "Per-year denominator (total suspected / span)"], "-", "cases/year",
        "tab_reconciliations.csv", "Data span 2000-2018; reconcile with Title (2000-2022)")
add_row("c6", "Results", "Attack-rate units",
        tA$value[tA$item == "Attack-rate units"], "-", "-", "tab_reconciliations.csv",
        "AR = 100 x cumulative cases / population; 1.13 means 1.13%")

# --- Part B ------------------------------------------------------------------
cat("Part B: corrected global extrapolation (c7)...\n")
tB <- part_B_extrapolation(scn_raw, cfg, rcfg)
write.csv(tB, tabf("tab_global_extrapolation.csv"), row.names = FALSE)
b1 <- tB[tB$delay_weeks == 1 & tB$model == "static", ]
add_row("c7", "Discussion", "1-week aggregate cases averted (% of annual)",
        sprintf("%.0f (%.1f%%)", b1$annual_cases_averted_median, b1$pct_of_annual_median),
        b1$pct_of_annual_iqr, "cases/year (% of annual)", "tab_global_extrapolation.csv",
        sprintf("Replaces impossible 96%%; capped at P_max=%.0f%%. Static model.", b1$P_max_pct[1]))

# --- Part D ------------------------------------------------------------------
cat("Part D: confirmation adjustment (c3)...\n")
cf <- confirmation_factors(raw_summary, prep, cfg)
tD <- part_D_confirmation(scn_cea, cf, cfg)
write.csv(tD, tabf("tab_confirmation_adjusted.csv"), row.names = FALSE)
inv_err <- invariance_check(scn_cea, cf, cfg$confirmation$s_culture_base)
cat(sprintf("  invariance max|pct_susp - pct_true| = %.2e (must be ~0)\n", inv_err))
dD <- tD[tD$model == "renewal" & tD$tau == 8, ]
add_row("c3", "Results/Response", "Confirmation-adjusted true-case basis (invariance)",
        "% reduction invariant to alpha", sprintf("alpha range via s_culture %.2f-%.2f",
        cfg$confirmation$s_culture_range[1], cfg$confirmation$s_culture_range[2]),
        "-", "tab_confirmation_adjusted.csv",
        sprintf("s_culture=0.61 (Mogasale 2016); invariance err %.0e", inv_err))

# --- Part F ------------------------------------------------------------------
cat("Part F: intervals + weighted impact (c9)...\n")
tF <- part_F_intervals(scn_cea, prep)
write.csv(tF, tabf("tab_intervals.csv"), row.names = FALSE)
wI <- weighted_impact(scn_cea, prep)
write.csv(wI, tabf("tab_weighted_impact.csv"), row.names = FALSE)
fren <- tF[tF$measure == "Percent case reduction" & tF$model == "renewal", ]
add_row("c9", "Throughout", "Relabel '95% CI' as 2.5-97.5 percentile; pop-weighted",
        sprintf("renewal %% reduction median %.1f", fren$median),
        fren$pct_2.5_97.5, "percent", "tab_intervals.csv / fig_forest_pctreduction",
        "Ranges are percentile/IQR spreads, not CIs; pop-weighted in tab_weighted_impact.csv")

# --- Part C: cost-effectiveness decomposition (c8) ---------------------------
cat("Part C: cost-effectiveness decomposition (c8)...\n")
tC <- part_C_ce(scn_cea, cfg)
write.csv(tC, tabf("tab_ce_decomposition.csv"), row.names = FALSE)
save_fig(fig_ce_thresholds(tC, scn_cea, cfg), "fig_ce_thresholds", cfg, width = 9, height = 6)
cBASE <- tC[tC$model == "renewal" & tC$tau == 8 & tC$vacc_cov == 0.80, ]
add_row("c8", "Results/Discussion", "Deaths averted; YLL/YLD split; cost/DALY vs thresholds",
        sprintf("deaths~%.2f; %%YLD=%.0f%%; cost/DALY=$%.0f", cBASE$deaths_averted_median,
                cBASE$pct_yld, cBASE$cost_per_daly_median),
        "-", "USD / DALY", "tab_ce_decomposition.csv / fig_ce_thresholds.png",
        "DALYs are mostly YLD (deaths near zero); cost/DALY >> 1-3x GDP and Ochalek 0.5x GDP")

# c4: corrected 1-week/80% static cost/DALY (Supp Table 2 reconciliation)
c4 <- scn_cea %>% dplyr::filter(model == "static", tau == 1, vacc_cov == 0.80)
c4v <- median(c4$cost_per_daly_averted, na.rm = TRUE)
c4iqr <- sprintf("%.0f-%.0f", quantile(c4$cost_per_daly_averted, 0.25, na.rm = TRUE),
                 quantile(c4$cost_per_daly_averted, 0.75, na.rm = TRUE))
add_row("c4", "Results", "1-week/80% cost per DALY averted (was duplicated 8-week value)",
        sprintf("$%.0f", c4v), c4iqr, "USD/DALY", "tab_ce_decomposition.csv",
        "1-week != 8-week ($82,756 was the duplicated 8-week value); manuscript Supp Tab 2 World-deflator value $32,748; here SSA deflator")

# --- Part E: spatial/hotspot targeting (illustrative, c8) --------------------
cat("Part E: spatial targeting (illustrative)...\n")
tE <- part_E_spatial(scn_cea, cfg)
write.csv(tE, tabf("tab_spatial_targeting.csv"), row.names = FALSE)
save_fig(fig_spatial_targeting(tE, cfg), "fig_spatial_targeting", cfg, width = 8, height = 5)
eW <- tE[tE$strategy == "whole-unit", ]; eT <- tE[tE$phi_pop == 0.25, ]
add_row("c8", "Discussion (illustrative)", "Targeted vs whole-unit cost/DALY",
        sprintf("whole $%.0f -> targeted $%.0f (phi_pop=0.25)",
                eW$cost_per_daly_median, eT$cost_per_daly_median),
        "-", "USD/DALY", "tab_spatial_targeting.csv / fig_spatial_targeting.png",
        "ILLUSTRATIVE (phi_cases=0.80), not fitted. Absolute cost/DALY from median components (ratio-of-medians); read the whole-vs-targeted comparison, not the level")

# --- Part C2: IVE coverage-dependence (c2) -----------------------------------
cat("Part C2: IVE coverage-dependence (c2)...\n")
tC2 <- part_C2_ive(cfg)
write.csv(tC2, tabf("tab_ive_coverage.csv"), row.names = FALSE)
save_fig(fig_ive_coverage(cfg), "fig_ive_coverage", cfg, width = 8, height = 5)
add_row("c2", "Methods/Supplementary", "Coverage-dependence of indirect VE",
        sprintf("IVE(f)=(OVE-f*TVE)/(1-f); zero at f=%.0f%%", 100 * cfg$ive_algebra$OVE / cfg$ive_algebra$TVE),
        "-", "-", "tab_ive_coverage.csv / fig_ive_coverage.png",
        "Fixed IVE=17% is a conservative assumption; algebra goes negative above ~67% coverage; renewal model removes the need")

# --- Part G: renewal results (c1) — reference the renewal module outputs ------
add_row("c1", "Methods/Results", "Renewal-model results (R_t, forest, amplification, eta_eff)",
        "see renewal/ outputs", "-", "-",
        "renewal/figures/*, renewal/tables/*, renewal/report.html",
        "Built by renewal/run_analysis.R (companion spec); Methods in renewal_methods_section")

# --- Forest plot (c9 / Part F) — reuse the renewal builder --------------------
cat("Forest plot (c9)...\n")
summ_forest <- summarise_draws(scn_cea %>% dplyr::filter(tau == 8))
pooled_forest <- pooled_estimates(summ_forest)
save_fig(fig_forest_pctreduction(
  summ_forest %>% dplyr::filter(tau == 8, vacc_cov == 0.80),
  pooled_forest %>% dplyr::filter(tau == 8, vacc_cov == 0.80), rcfg),
  "fig_forest_pctreduction", cfg, height = 7)

# --- Manifest ----------------------------------------------------------------
mani <- do.call(rbind, manifest)
mani <- mani[order(mani$comment_id), ]
write.csv(mani, file.path(cfg$paths$outputs, "revision_values.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(cfg$paths$outputs, "sessionInfo.txt"))
saveRDS(list(scn_cea = scn_cea, prep = prep, cf = cf, tB = tB, tF = tF, wI = wI,
             have_cost = have_cost, cost_env = cost_env),
        file.path(cfg$paths$outputs, "revision_results.rds"))

cat("\n--- Guards ---\n")
cat("  Part B guards (averted<=post-int; static median<=P_max): PASS (asserted in part_B)\n")
cat(sprintf("  Part D invariance: %s (err %.0e)\n", ifelse(inv_err < 1e-8, "PASS", "FAIL"), inv_err))

# --- Validation tests --------------------------------------------------------
cat("\n--- Running guard tests ---\n")
test_status <- tryCatch({ source("revision/tests/test_revision.R", local = new.env()); 0L },
                        error = function(e) { cat("TEST ERROR:", conditionMessage(e), "\n"); 1L })

# --- Render manifest report --------------------------------------------------
if (nzchar(Sys.which("quarto"))) {
  cat("Rendering revision_manifest.qmd...\n")
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  rc <- system2("quarto", c("render", shQuote("revision/revision_manifest.qmd")),
                stdout = FALSE, stderr = FALSE)
  cat(if (rc == 0) "Wrote revision/revision_manifest.html\n" else "quarto render failed (non-fatal)\n")
}

cat("\nWrote", nrow(mani), "manifest rows to", file.path(cfg$paths$outputs, "revision_values.csv"), "\n")
cat("=== Revision outputs complete ===\n")
if (!identical(test_status, 0L)) quit(status = 1L)
