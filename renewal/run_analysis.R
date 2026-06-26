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

source("renewal/R/gi.R")
source("renewal/R/renewal_core.R")
source("renewal/R/data_prep.R")

cfg <- yaml::read_yaml("renewal/config.yml")
set.seed(cfg$seed)
cat("=== Renewal analysis — Phase 1 ===\n")
cat("Seed:", cfg$seed, "| Sobol n:", cfg$n_sobol,
    "| GI mean:", cfg$gi$mean_days, "d (CV", cfg$gi$cv, ")",
    "| drop_first_week:", cfg$gi$drop_first_week, "\n\n")

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

# --- Reproducibility record --------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(cfg$paths$outputs, "sessionInfo.txt"))
cat("\nWrote", file.path(cfg$paths$outputs, "sessionInfo.txt"), "\n")
cat("=== Phase 1 complete ===\n")
if (!identical(test_status, 0L)) quit(status = 1L)
