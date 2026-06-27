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
source("impact_drivers/R/features.R")

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

writeLines(capture.output(sessionInfo()), file.path(cfg$paths$outputs, "sessionInfo.txt"))
saveRDS(list(dat = dat, prep = prep), file.path(cfg$paths$outputs, "impact_drivers_phase1.rds"))

cat(sprintf("\nWrote %s (%d rows = %d outbreaks x %d delays)\n",
            tabf("tab_feature_outcomes.csv"), nrow(dat),
            length(unique(dat$study_id)), length(cfg$delays_weeks)))
cat("\nFeatures:", paste(setdiff(names(feat), c("study_id","tau")), collapse = ", "), "\n")
cat("Outcomes: case_averted_per_1000, pct_reduction, cases_averted, cost_per_daly, amplification\n")
cat("=== Phase 1 complete ===\n")
