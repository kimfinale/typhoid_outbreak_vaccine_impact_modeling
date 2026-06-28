# =============================================================================
# Phase 1 — synthetic identifiability gate (MUST pass before real-data fits).
# Simulate outbreaks across a theta gradient under the generative model, fit the
# NEUTRAL hierarchical model, and check whether theta is recoverable. PASS ->
# proceed to Phase 2. FAIL -> write IDENTIFIABILITY_FAILED.md and STOP.
#   Rscript hierarchical_theta/run_phase1_gate.R
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(cmdstanr); library(posterior); library(ggplot2) })
source("hierarchical_theta/R/hier_helpers.R")

cfg <- yaml::read_yaml("hierarchical_theta/config.yml")
set.seed(cfg$seed)
for (p in cfg$paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
tabf <- function(n) file.path(cfg$paths$tables, n)
outf <- function(n) file.path(cfg$paths$outputs, n)
cat("=== Phase 1: synthetic identifiability gate ===\n")
cat("cmdstan:", cmdstan_version(), "| R seed:", cfg$seed, "| Stan seed:", cfg$stan_seed, "\n\n")

# --- Compile -----------------------------------------------------------------
mod <- cmdstan_model(file.path(cfg$paths$stan, "renewal_source.stan"))
cat("Model compiled.\n")

# --- Simulate + fit (neutral) ------------------------------------------------
syn <- simulate_synthetic(cfg, gi_mean = 2.0)
cat("Simulated", length(syn$series), "outbreaks; true theta in {",
    paste(sort(unique(syn$true_theta)), collapse = ", "), "}\n")
sdata <- build_stan_data(syn$series, cfg, modes = syn$modes, informed = 0L)

fit <- mod$sample(
  data = sdata, seed = cfg$stan_seed,
  chains = cfg$mcmc$chains, parallel_chains = cfg$mcmc$parallel_chains,
  iter_warmup = cfg$mcmc$iter_warmup, iter_sampling = cfg$mcmc$iter_sampling,
  adapt_delta = cfg$mcmc$adapt_delta, max_treedepth = cfg$mcmc$max_treedepth,
  refresh = 200)

# --- Diagnostics -------------------------------------------------------------
diag <- fit$diagnostic_summary()
sumy <- fit$summary(c("theta", "tau_theta", "gi_mean", "phi"))
ndiv <- sum(diag$num_divergent); rhat_max <- max(sumy$rhat, na.rm = TRUE)
cat(sprintf("\nDiagnostics: divergences=%d, max Rhat=%.3f, min ESS=%.0f\n",
            ndiv, rhat_max, min(sumy$ess_bulk, na.rm = TRUE)))

# --- Recovery ----------------------------------------------------------------
theta_draws <- as_draws_matrix(fit$draws("theta"))
rec <- recovery_metrics(theta_draws, syn$true_theta)
rec$per$outbreak <- names(syn$series)
write.csv(rec$per, tabf("synthetic_recovery.csv"), row.names = FALSE)
cat(sprintf("\nRecovery: 90%% CrI coverage=%.2f | median bias=%.3f | rank cor=%.2f\n",
            rec$coverage, rec$median_bias, rec$rank_cor))

# --- Recovery figure ---------------------------------------------------------
p <- ggplot(rec$per, aes(true_theta, post_median)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_errorbar(aes(ymin = cri90_lo, ymax = cri90_hi), width = 0.01,
                colour = ifelse(rec$per$covered, "#0072B2", "#D55E00")) +
  geom_point(size = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "true theta", y = "posterior theta (median, 90% CrI)",
       title = "Synthetic identifiability gate: theta recovery",
       subtitle = sprintf("coverage=%.2f, bias=%.3f, rank cor=%.2f (orange = not covered)",
                          rec$coverage, rec$median_bias, rec$rank_cor)) +
  theme_minimal(base_size = 11)
ggsave(file.path(cfg$paths$figures, "fig_synthetic_recovery.png"), p, width = 7, height = 7, dpi = 300)

# --- GATE --------------------------------------------------------------------
g <- cfg$gate
converged <- ndiv == 0 && rhat_max < 1.05
pass <- converged &&
  rec$coverage >= g$min_cri_coverage &&
  abs(rec$median_bias) <= g$max_abs_median_bias &&
  rec$rank_cor >= g$min_rank_cor
saveRDS(list(rec = rec, diag = diag, pass = pass, syn = syn),
        outf("phase1_gate.rds"))

cat("\n=== GATE:", ifelse(pass, "PASS", "FAIL"), "===\n")
if (pass) {
  writeLines(c(sprintf("GATE PASSED %s", Sys.time()),
               sprintf("coverage=%.2f bias=%.3f rank_cor=%.2f div=%d rhat=%.3f",
                       rec$coverage, rec$median_bias, rec$rank_cor, ndiv, rhat_max)),
             outf("GATE_PASSED.txt"))
  cat("Proceed to Phase 2 (run_phase2_3.R).\n")
} else {
  writeLines(c(
    "# Identifiability gate FAILED",
    "",
    sprintf("Synthetic recovery of the source fraction theta did not meet the pre-set criteria (coverage>=%.2f, |bias|<=%.2f, rank cor>=%.2f) or the fit did not converge.",
            g$min_cri_coverage, g$max_abs_median_bias, g$min_rank_cor),
    "",
    sprintf("- 90%% CrI coverage: %.2f", rec$coverage),
    sprintf("- median bias: %.3f", rec$median_bias),
    sprintf("- rank correlation: %.2f", rec$rank_cor),
    sprintf("- divergences: %d ; max Rhat: %.3f", ndiv, rhat_max),
    "",
    "**Conclusion:** theta is weakly identified from curve shape alone. The article",
    "classification (Stream A in source_decomposition/) is PRIMARY; this Bayesian model",
    "is at most a loose consistency check with wide intervals. Real-data theta point",
    "estimates are NOT reported. See README.md."),
    "hierarchical_theta/IDENTIFIABILITY_FAILED.md")
  cat("Wrote IDENTIFIABILITY_FAILED.md. STOPPING (do not fit real data).\n")
}
