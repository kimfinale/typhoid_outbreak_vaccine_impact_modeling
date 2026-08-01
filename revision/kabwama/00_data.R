# =============================================================================
# Kabwama (Kampala 2015) comprehensive analysis — 00: assemble real inputs.
#   Rscript revision/kabwama/00_data.R
# Extracts the weekly suspected-case series and the fitted community-PPV posterior
# (per-draw pi for Kabwama, plus mu_pi/sigma_pi), and saves them for the model.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({ library(yaml) })
source("renewal/R/gi.R"); source("renewal/R/renewal_core.R")
source("renewal/R/data_prep.R"); source("renewal/R/epiestim_rt.R")
source("renewal/R/ppv.R")
OUT <- "revision/kabwama/data"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cfg  <- yaml::read_yaml("renewal/config.yml")
prep <- prep_outbreaks(cfg)
STUDY <- "Kabwama 2017"
S <- as.numeric(prep$series[[STUDY]])                 # weekly suspected-case incidence
stopifnot(length(S) > 0)
cat(sprintf("Kabwama weekly series: %d weeks, total suspected = %.0f, peak = %.0f (week %d)\n",
            length(S), sum(S), max(S), which.max(S)))
cat("weekly:", paste(round(S), collapse = " "), "\n")

# --- community-PPV posterior for Kabwama ------------------------------------
pi_post <- load_pi_posterior(cfg$ppv$draws, cfg$ppv$anchor)
pi_kab  <- as.numeric(pi_post$pi_anchor[, "Kabwama 2017"])
mu_pi   <- as.numeric(pi_post$mu_pi); sigma_pi <- as.numeric(pi_post$sigma_pi)
q <- function(v) round(quantile(v, c(.5, .05, .95)), 3)
cat(sprintf("\nKabwama PPV pi: median %.3f [90%% %.3f, %.3f]  (n=%d draws)\n",
            q(pi_kab)[1], q(pi_kab)[2], q(pi_kab)[3], length(pi_kab)))
cat(sprintf("population mu_pi -> pi_pop median %.3f [%.3f, %.3f]; sigma_pi median %.2f\n",
            q(plogis(mu_pi))[1], q(plogis(mu_pi))[2], q(plogis(mu_pi))[3], median(sigma_pi)))

saveRDS(list(study = STUDY, weekly = S, n_weeks = length(S),
             pi_draws = pi_kab, mu_pi = mu_pi, sigma_pi = sigma_pi,
             gi_mean = cfg$gi$mean_days, gi_cv = cfg$gi$cv,
             coverage = cfg$vaccine$coverage_base, psi = cfg$vaccine$psi_mean,
             psiT_frac_lower = cfg$vaccine$psiT_frac_lower,
             delay_base_weeks = cfg$timing$delay_base_weeks,
             campaign_days = cfg$timing$campaign_duration_days,
             onset_days = cfg$timing$immunity_onset_days, step_days = cfg$step_days),
        file.path(OUT, "kabwama_inputs.rds"))
write.csv(data.frame(week = seq_along(S), suspected = S),
          file.path(OUT, "kabwama_weekly.csv"), row.names = FALSE)
cat("\nSaved", file.path(OUT, "kabwama_inputs.rds"), "and kabwama_weekly.csv\n")
