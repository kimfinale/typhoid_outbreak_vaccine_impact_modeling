# =============================================================================
# Classify transmission mode and run the two-level additive counterfactual:
#
#   suspected = true typhoid + other febrile illness
#   true typhoid = exogenous/source + propagated renewal incidence
#
# The former theta-weighted outcome calculation is retained only as a paired
# structural comparator. Run from the repository root:
#   Rscript source_decomposition/run_source_decomposition.R
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", ".")
setwd(ROOT)
suppressWarnings(
  if (requireNamespace("data.table", quietly = TRUE))
    data.table::setDTthreads(1))
suppressMessages({ library(dplyr); library(tidyr) })

for (f in c("gi", "renewal_core", "data_prep", "scenario", "cost_daly",
            "summarise", "ppv"))
  source(file.path("renewal/R", paste0(f, ".R")))
for (f in c("streamB", "decompose", "figures_sd"))
  source(file.path("source_decomposition/R", paste0(f, ".R")))

rcfg <- yaml::read_yaml("renewal/config.yml")
cfg <- yaml::read_yaml("source_decomposition/config.yml")
set.seed(cfg$seed)
cat("=== Additive source + propagated-transmission analysis ===\n")
cat("Seed:", cfg$seed, "| scenario tau =", cfg$scenario$tau_weeks,
    "cov =", cfg$scenario$coverage, "\n\n")
for (p in cfg$paths)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
tabf <- function(n) file.path(cfg$paths$tables, n)

prep <- prep_outbreaks(rcfg)
pi_post <- if (isTRUE(rcfg$ppv$enable))
  load_pi_posterior(rcfg$ppv$draws, rcfg$ppv$anchor) else NULL
if (isTRUE(rcfg$ppv$enable) && is.null(pi_post))
  stop("The configured additive observation model requires the PPV posterior")
if (!is.null(pi_post))
  cat(sprintf("Observation additivity: PPV posterior enabled (%d draws)\n",
              pi_post$ndraw))

# --- Stream A: primary paper evidence ----------------------------------------
streamA <- read.csv(
  "source_decomposition/stream_A_classification.csv",
  stringsAsFactors = FALSE, check.names = FALSE)
streamA_ren <- streamA[streamA$study_id %in% names(prep$series), ]
cat("Stream A (papers): classified", nrow(streamA), "outbreaks;",
    nrow(streamA_ren), "in the renewal set\n")

# --- Stream B: observed-curve diagnostic, not a fitted T_t decomposition -----
cat("Stream B: observed suspected-curve source diagnostics (suggestive only)...\n")
streamB <- stream_B(prep, rcfg, cfg)
write.csv(streamB, tabf("tab_theta_data.csv"), row.names = FALSE)

# --- Cross-validation: paper point, curve diagnostic widens range ------------
mode_tab <- cross_validate(streamA_ren, streamB)
write.csv(
  mode_tab[, c(
    "study_id", "mode_paper", "theta_paper_lo", "theta_paper_hi",
    "theta_data", "theta_data_lo", "theta_data_hi", "agreement",
    "final_theta", "final_theta_lo", "final_theta_hi", "confidence",
    "cs_r2", "n_waves", "n_resurg", "citation", "source_evidence")],
  tabf("tab_transmission_mode.csv"), row.names = FALSE)
cat(sprintf(
  "Cross-validation: %d agree, %d disagree (theta within 0.2)\n",
  sum(mode_tab$agreement == "agree"),
  sum(mode_tab$agreement == "disagree")))
cat("Paper modes among the renewal set:\n")
print(table(mode_tab$mode_paper))

# --- Paired primary additive model versus theta-weighted comparators ----------
cat("\nRunning paired source-formulation comparison...\n")
paired <- simulate_source_comparison(
  prep, rcfg, mode_tab, cfg, pi_post = pi_post)
cmp <- summarise_source_comparison(paired)
pooled <- summarise_pooled_source_comparison(paired)
val <- source_validation_table(paired)

write.csv(
  cmp, tabf("tab_additive_vs_theta_weighted.csv"), row.names = FALSE)
write.csv(
  pooled, tabf("tab_additive_vs_theta_weighted_pooled.csv"),
  row.names = FALSE)
write.csv(
  val, tabf("tab_additive_validation.csv"), row.names = FALSE)
saveRDS(
  paired,
  file.path(cfg$paths$outputs, "additive_vs_theta_weighted.rds"))

cat(sprintf(
  paste0(
    "Pooled true-typhoid reduction: additive %.1f%% | historical ",
    "theta-weighted %.1f%% | paired difference %+.1f percentage points\n"),
  pooled$additive_true_pct_median,
  pooled$historical_weighted_true_pct_median,
  pooled$delta_true_pp_median))
cat(sprintf("Validation: %d/%d checks passed\n",
            sum(val$pass), nrow(val)))
if (!all(val$pass)) {
  print(val[!val$pass, ], row.names = FALSE)
  stop("Additive source-model validation failed")
}

# --- Figures -----------------------------------------------------------------
cat("Building figures...\n")
save_fig_sd(
  fig_mode_distribution(mode_tab, cfg),
  "fig_mode_distribution", cfg, width = 8, height = 5)
if (requireNamespace("ggrepel", quietly = TRUE))
  save_fig_sd(
    fig_theta_crossvalidate(mode_tab, cfg),
    "fig_theta_crossvalidate", cfg, height = 7)
save_fig_sd(
  fig_additive_vs_theta(cmp, mode_tab, cfg),
  "fig_additive_vs_theta_weighted", cfg, height = 7)
save_fig_sd(
  fig_additive_difference(cmp, mode_tab, cfg),
  "fig_additive_difference", cfg, height = 6)

writeLines(
  capture.output(sessionInfo()),
  file.path(cfg$paths$outputs, "sessionInfo.txt"))
saveRDS(
  list(
    mode_tab = mode_tab, streamB = streamB, comparison = cmp,
    pooled = pooled, validation = val),
  file.path(cfg$paths$outputs, "source_decomposition.rds"))

if (nzchar(Sys.which("quarto"))) {
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  rc <- system2(
    "quarto",
    c("render", shQuote("source_decomposition/report.qmd")),
    stdout = FALSE, stderr = FALSE)
  cat(if (rc == 0)
    "Wrote source_decomposition/report.html\n"
  else
    "quarto render skipped/failed\n")
}
cat("=== Additive source analysis complete ===\n")
