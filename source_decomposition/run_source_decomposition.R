# =============================================================================
# run_source_decomposition.R — classify each outbreak's transmission mode
# (common-source vs propagated vs mixed), assign a source fraction theta, and
# run the static<->renewal decomposition bracket. Run from the repo root:
#     Rscript source_decomposition/run_source_decomposition.R
# Reuses the renewal/static engines.
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressWarnings(if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1))
suppressMessages({ library(dplyr); library(tidyr) })

for (f in c("gi","renewal_core","data_prep","scenario","cost_daly","summarise"))
  source(file.path("renewal/R", paste0(f, ".R")))
for (f in c("streamB","decompose","figures_sd"))
  source(file.path("source_decomposition/R", paste0(f, ".R")))

rcfg <- yaml::read_yaml("renewal/config.yml")
cfg  <- yaml::read_yaml("source_decomposition/config.yml")
set.seed(cfg$seed)
cat("=== Transmission-mode classification & source decomposition ===\n")
cat("Seed:", cfg$seed, "| scenario tau =", cfg$scenario$tau_weeks,
    "cov =", cfg$scenario$coverage, "\n\n")
for (p in cfg$paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
tabf <- function(n) file.path(cfg$paths$tables, n)

prep <- prep_outbreaks(rcfg)
amr_props <- load_amr_props(rcfg)

# --- Stream A (paper evidence; read from the curated classification) ---------
streamA <- read.csv("source_decomposition/stream_A_classification.csv",
                    stringsAsFactors = FALSE, check.names = FALSE)
streamA_ren <- streamA[streamA$study_id %in% names(prep$series), ]
cat("Stream A (papers): classified", nrow(streamA), "outbreaks;",
    nrow(streamA_ren), "in the renewal set\n")

# --- Stream B (data-driven theta from the incidence series) ------------------
cat("Stream B: fitting common-source model + propagation diagnostics...\n")
streamB <- stream_B(prep, rcfg, cfg)
write.csv(streamB, tabf("tab_theta_data.csv"), row.names = FALSE)

# --- Cross-validate -> final theta ------------------------------------------
mode_tab <- cross_validate(streamA_ren, streamB)
write.csv(mode_tab[, c("study_id", "mode_paper", "theta_paper_lo", "theta_paper_hi",
                       "theta_data", "theta_data_lo", "theta_data_hi", "agreement",
                       "final_theta", "final_theta_lo", "final_theta_hi", "confidence",
                       "cs_r2", "n_waves", "n_resurg", "citation", "source_evidence")],
          tabf("tab_transmission_mode.csv"), row.names = FALSE)
cat(sprintf("Cross-validation: %d agree, %d disagree (theta within 0.2)\n",
            sum(mode_tab$agreement == "agree"), sum(mode_tab$agreement == "disagree")))
cat("Paper modes among the 13:\n"); print(table(mode_tab$mode_paper))

# --- Scenarios (static + renewal) at the base scenario -----------------------
cat("\nRunning static + renewal scenarios (tau=8, 80% coverage)...\n")
scn <- run_scenarios(prep, rcfg, amr_props, coverage = cfg$scenario$coverage,
                     delays = cfg$scenario$tau_weeks)

# --- Validation: theta=1 -> static, theta=0 -> renewal -----------------------
val <- validate_endpoints(scn, cfg)
cat(sprintf("[%s] Endpoint validation: max|theta1-static|=%.0e, max|theta0-renewal|=%.0e\n",
            ifelse(max(val) < 1e-10, "PASS", "FAIL"), val[1], val[2]))

# --- Decomposition bracket at final theta ------------------------------------
dec <- decompose_bracket(scn, mode_tab, cfg)
write.csv(dec, tabf("tab_decomposition.csv"), row.names = FALSE)
cat(sprintf("\nPooled %% reduction: static (theta=1) %.1f%% | renewal (theta=0) %.1f%% | at theta %.1f%%\n",
            median(dec$pct_static), median(dec$pct_renewal), median(dec$pct_theta)))

# --- Figures -----------------------------------------------------------------
cat("Building figures...\n")
save_fig_sd(fig_mode_distribution(mode_tab, cfg), "fig_mode_distribution", cfg, width = 8, height = 5)
if (requireNamespace("ggrepel", quietly = TRUE))
  save_fig_sd(fig_theta_crossvalidate(mode_tab, cfg), "fig_theta_crossvalidate", cfg, height = 7)
save_fig_sd(fig_decomposition_bracket(dec, mode_tab, cfg), "fig_decomposition_bracket", cfg, height = 6)

writeLines(capture.output(sessionInfo()), file.path(cfg$paths$outputs, "sessionInfo.txt"))
saveRDS(list(mode_tab = mode_tab, streamB = streamB, dec = dec, val = val),
        file.path(cfg$paths$outputs, "source_decomposition.rds"))

if (nzchar(Sys.which("quarto"))) {
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  rc <- system2("quarto", c("render", shQuote("source_decomposition/report.qmd")),
                stdout = FALSE, stderr = FALSE)
  cat(if (rc == 0) "Wrote source_decomposition/report.html\n" else "quarto render skipped/failed\n")
}
cat("=== Source decomposition complete ===\n")
