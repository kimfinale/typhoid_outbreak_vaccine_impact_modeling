# Validate the additive vs multiplicative PPV regimes on one anchored outbreak.
setwd(Sys.getenv("RENEWAL_ROOT", "."))
source("R/2-functions.R"); source("renewal/R/gi.R"); source("renewal/R/renewal_core.R")
source("renewal/R/data_prep.R"); source("renewal/R/scenario.R"); source("renewal/R/cost_daly.R")
source("renewal/R/summarise.R"); source("renewal/R/ppv.R")
cfg <- yaml::read_yaml("renewal/config.yml"); set.seed(cfg$seed)
prep <- prep_outbreaks(cfg); amr <- load_amr_props(cfg)
pi_post <- load_pi_posterior(cfg$ppv$draws, cfg$ppv$anchor)
stopifnot(!is.null(pi_post))

run1 <- function(reg) run_scenarios(prep, cfg, amr,
  coverage = cfg$vaccine$coverage_base, delays = cfg$timing$delay_base_weeks,
  pi_post = pi_post, ppv_regime = reg)

for (sid in c("Kabwama 2017", "Neil 2012")) {
  cat("\n==============", sid, "(base case tau=8, cov=0.8, renewal) ==============\n")
  a <- run1("additive"); a <- a[a$study_id == sid & a$model == "renewal", ]
  m <- run1("multiplicative"); m <- m[m$study_id == sid & m$model == "renewal", ]
  cat(sprintf("pi_ppv median            : %.3f\n", median(a$pi_ppv)))
  cat(sprintf("suspected total          : %.0f\n", a$suspected_tot[1]))
  cat(sprintf("[ADD] true total (med)   : %.0f   (true/susp = %.3f ~ pi)\n",
              median(a$s_ch_tot), median(a$s_ch_tot) / a$suspected_tot[1]))
  cat(sprintf("[ADD] true cases averted : %.0f\n", median(a$s_ch_averted_tot)))
  cat(sprintf("[ADD] TRUE %%reduction     : %.1f%%\n", median(100 * a$s_ch_averted_tot / a$s_ch_tot)))
  cat(sprintf("[ADD] OBSERVED %%reduction : %.1f%%  (diluted by background)\n", median(a$obs_pct_reduction, na.rm = TRUE)))
  cat(sprintf("[MUL] true cases averted : %.0f\n", median(m$s_ch_averted_tot)))
  cat(sprintf("[MUL] %%reduction (=obs)   : %.1f%%\n", median(100 * m$s_ch_averted_tot / m$s_ch_tot)))
}
