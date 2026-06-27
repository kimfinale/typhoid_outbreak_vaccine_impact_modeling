# =============================================================================
# Manuscript-revision deliverables, Parts A-G. Each function returns a data
# frame (and writes nothing); the orchestrator writes CSVs and figures.
#
# Reuses the renewal engine (renewal/R/*) and the manuscript CEA functions
# (R/2-functions.R via cost_daly.R). Static and renewal results share identical
# inputs and downstream cost/DALY machinery — only cases-averted differs.
# =============================================================================

suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

q025 <- function(x) as.numeric(stats::quantile(x, 0.025, na.rm = TRUE))
q975 <- function(x) as.numeric(stats::quantile(x, 0.975, na.rm = TRUE))
iqr_lo <- function(x) as.numeric(stats::quantile(x, 0.25, na.rm = TRUE))
iqr_hi <- function(x) as.numeric(stats::quantile(x, 0.75, na.rm = TRUE))
fmt_ci <- function(lo, hi, d = 1) sprintf("%.*f-%.*f", d, lo, d, hi)

# ============================================================================
# Part A — Data reconciliations (c0, c5, c6)
# ============================================================================
part_A_reconciliations <- function(raw_summary, prep, cfg, rcfg) {
  d <- raw_summary[!raw_summary$study_id %in% RENEWAL_EXCLUSIONS, ]   # 19 eligible
  num <- function(x) suppressWarnings(as.numeric(x))

  # c5: span from the data. Denominator basis = the analysis set's observed
  # case total (sum of the time-series series), which matches the manuscript's
  # 26,298 (-> /18 = 1,461/yr), NOT the all-19 tot_suspected sum.
  yr_min <- min(num(d$outbreak_start_year), na.rm = TRUE)
  yr_max <- max(num(d$outbreak_end_year),  na.rm = TRUE)
  span_excl <- yr_max - yr_min                  # difference convention (manuscript "18 years")
  span_incl <- span_excl + 1                    # inclusive count
  analysis_total <- sum(vapply(prep$series, sum, numeric(1)))
  denom_per_year <- analysis_total / cfg$span_years_manuscript

  # c6: attack-rate units — verify AR = 100 * cumulative cases / population
  ar_check <- d %>% mutate(ar_data = num(attack_rate),
                           ar_recomputed = 100 * num(tot_suspected_cases) / num(population)) %>%
    filter(is.finite(ar_data), is.finite(ar_recomputed))
  ar_ratio <- median(ar_check$ar_data / ar_check$ar_recomputed, na.rm = TRUE)
  ar_units <- if (abs(ar_ratio - 1) < 0.25) "percent (per 100)" else
              if (abs(ar_ratio - 10) < 2.5) "per 1,000" else "indeterminate"
  ar_median <- median(num(d$attack_rate), na.rm = TRUE)

  data.frame(
    comment_id = c("c0", "c0", "c5", "c5", "c5", "c5", "c6", "c6"),
    item = c("Excluded outbreaks (total)", "Excluded breakdown (Gavi/prior/short/nested)",
             "Min outbreak year", "Max outbreak year", "Span (years)",
             "Per-year denominator (total suspected / span)",
             "Attack-rate units", "Median attack rate"),
    value = c(length(RENEWAL_EXCLUSIONS), "20 / 1 / 1 / 2 = 24 (Figure 1)",
              yr_min, yr_max, sprintf("%d (diff) or %d (inclusive)", span_excl, span_incl),
              sprintf("%.0f / %d = %.1f per year", analysis_total, cfg$span_years_manuscript, denom_per_year),
              ar_units, sprintf("%.2f", ar_median)),
    source_field = c("R/1-data.R exclusion list", "Figure 1 (not in dataset)",
                     "outbreak_start_year", "outbreak_end_year", "outbreak_start/end_year",
                     "tot_suspected_cases", "attack_rate vs 100*tot_suspected/population",
                     "attack_rate"),
    recommendation = c(
      "19 retained for analysis; 24 excluded",
      "Use 20/1/1/2 (=24); submitted 17/2/1/2 (=22) is inconsistent - confirm vs Figure 1",
      "Reconcile Title/Methods (2000-2022) with data", "Data max is 2018, not 2020/2022",
      sprintf("Use %d years if inclusive; manuscript uses %d", span_incl, cfg$span_years_manuscript),
      "Denominator basis for the annual extrapolation (Part B)",
      "State explicitly as percent; 1.13 means 1.13%", "median 1.13% (per 100)"),
    stringsAsFactors = FALSE)
}

# ============================================================================
# Part B — Corrected global extrapolation (c7)
# Per-draw aggregate across outbreaks, annualized by the reconciled span.
# GUARDS: per-outbreak averted <= post-intervention cases; aggregate % <= P_max.
# ============================================================================
P_max_value <- function(rcfg) {
  P_halloran(rcfg$vaccine$coverage_base, rcfg$vaccine$psi_mean, rcfg$vaccine$eta_static)
}

part_B_extrapolation <- function(scn_raw, cfg, rcfg, span_years = cfg$span_years_manuscript) {
  scn_raw <- scn_raw %>% filter(vacc_cov == cfg$extrapolation$coverage)  # 80% coverage
  stopifnot(all(scn_raw$s_ch_averted_tot <= scn_raw$post_int_cases + 1e-6))  # GUARD 1
  Pmax <- P_max_value(rcfg)
  total_observed <- scn_raw %>% filter(model == "renewal", tau == tau[1], draw == 1) %>%
    summarise(t = sum(s_ch_tot)) %>% pull(t)
  annual_observed <- total_observed / span_years

  agg <- scn_raw %>%
    group_by(model, tau, draw) %>%
    summarise(averted = sum(s_ch_averted_tot), .groups = "drop") %>%
    mutate(annual_averted = averted / span_years,
           pct_of_annual = 100 * annual_averted / annual_observed)
  if (any(agg$pct_of_annual > 100 + 1e-6))                                   # GUARD: <= observed
    stop("Aggregate % averted exceeds 100% of observed.")

  out <- agg %>% group_by(delay_weeks = tau, model) %>%
    summarise(annual_cases_averted_median = median(annual_averted),
              annual_averted_iqr = fmt_ci(iqr_lo(annual_averted), iqr_hi(annual_averted), 0),
              annual_averted_pct95 = fmt_ci(q025(annual_averted), q975(annual_averted), 0),
              pct_of_annual_median = median(pct_of_annual),
              pct_of_annual_iqr = fmt_ci(iqr_lo(pct_of_annual), iqr_hi(pct_of_annual)),
              annual_observed = round(annual_observed),
              P_max_pct = round(100 * Pmax, 1), .groups = "drop")
  # GUARD: the reported STATIC point estimate (median) cannot exceed P_max
  # (the 96% impossibility). The renewal model legitimately exceeds it (herd
  # compounding) and individual high-psi draws may too; only the static median
  # is bounded.
  st <- out[out$model == "static", ]
  if (any(st$pct_of_annual_median > 100 * Pmax + 1e-6))
    stop(sprintf("Static median %% averted exceeds P_max (%.1f%%)", 100 * Pmax))
  out
}

# ============================================================================
# Part D — Confirmation adjustment for suspected cases (c3)
# alpha_o = c_o / s_culture, capped at 1.  c_o = confirmed/suspected per outbreak.
# Only true-case-linked quantities scale; vaccination cost is FIXED:
#   cases_true   = alpha * cases_averted
#   DALY_true    = alpha * DALY_averted
#   tr_cost_true = alpha * coi_averted
#   net_cost     = vacc_cost - tr_cost_true
#   cost/DALY    = (V - alpha*Tr) / (alpha*D)
# INVARIANCE: percent case reduction is unchanged (num & denom both scale by alpha).
# ============================================================================
confirmation_factors <- function(raw_summary, prep, cfg) {
  num <- function(x) suppressWarnings(as.numeric(x))
  d <- raw_summary
  # c_o = culture positivity = confirmed-by-culture (AI) / tested-by-culture (AG).
  # Correct denominator (the TESTED, not all suspected); proportion of clinically
  # suspected cases that are true typhoid follows as PPV = c_o / s_culture.
  co <- num(d$lab_confirmed_cases) / num(d$lab_tested_cases)
  co[is.finite(co)] <- pmin(co[is.finite(co)], 1)
  cf <- data.frame(study_id = d$study_id, c_o = co, stringsAsFactors = FALSE)
  cf <- cf[cf$study_id %in% names(prep$series), ]
  sub <- d[d$study_id %in% names(prep$series), ]
  pooled_co <- sum(num(sub$lab_confirmed_cases), na.rm = TRUE) /
               sum(num(sub$lab_tested_cases), na.rm = TRUE)
  cf$c_o_source <- ifelse(is.finite(cf$c_o), "per-outbreak (AI/AG)", "pooled fallback")
  cf$c_o[!is.finite(cf$c_o)] <- pooled_co               # fallback for missing AG/AI
  attr(cf, "pooled_co") <- pooled_co
  cf
}

part_D_confirmation <- function(scn_cea, cf, cfg, s_culture = cfg$confirmation$s_culture_base) {
  base <- scn_cea %>%
    filter(tau == 1 | tau == rcfg_base_tau(), vacc_cov == 0.80) %>%
    left_join(cf[, c("study_id", "c_o")], by = "study_id") %>%
    mutate(alpha = pmin(c_o / s_culture, 1))            # cap at 1
  per <- base %>%
    group_by(study_id, model, tau) %>%
    summarise(
      alpha = first(alpha),
      pct_reduction = median(ifelse(s_ch_tot > 0, 100 * s_ch_averted_tot / s_ch_tot, 0)),
      cases_averted_suspected = median(s_ch_averted_tot),
      cases_averted_true = median(alpha * s_ch_averted_tot),
      daly_true = median(alpha * daly_averted, na.rm = TRUE),
      vacc_cost = median(vacc_cost, na.rm = TRUE),
      tr_cost_true = median(alpha * coi_averted, na.rm = TRUE),
      cost_per_daly_true = median(ifelse(alpha * daly_averted > 0,
                                  (vacc_cost - alpha * coi_averted) / (alpha * daly_averted), NA),
                                  na.rm = TRUE),
      .groups = "drop")
  per
}

# Invariance check: % reduction identical before/after alpha (per outbreak/draw).
invariance_check <- function(scn_cea, cf, s_culture) {
  d <- scn_cea %>% left_join(cf[, c("study_id", "c_o")], by = "study_id") %>%
    mutate(alpha = pmin(c_o / s_culture, 1),
           pct_susp = ifelse(s_ch_tot > 0, 100 * s_ch_averted_tot / s_ch_tot, 0),
           pct_true = ifelse(alpha * s_ch_tot > 0,
                             100 * (alpha * s_ch_averted_tot) / (alpha * s_ch_tot), 0))
  max(abs(d$pct_susp - d$pct_true), na.rm = TRUE)
}

# ============================================================================
# Part F — Statistical reporting (c9): IQR + 2.5-97.5, pop/cases-weighted
# ============================================================================
part_F_intervals <- function(scn_cea, prep) {
  measures <- c(pct_reduction = "Percent case reduction",
                s_ch_averted_tot = "Cases averted",
                case_averted_per_1000_OCV = "Cases averted per 1,000 doses",
                cost_per_daly_averted = "Cost per DALY averted")
  d <- scn_cea %>%
    filter(vacc_cov == 0.80, tau == 8) %>%
    mutate(pct_reduction = ifelse(s_ch_tot > 0, 100 * s_ch_averted_tot / s_ch_tot, 0))
  out <- lapply(names(measures), function(v) {
    x <- d[[v]]
    data.frame(measure = measures[v], model = NA, median = median(x, na.rm = TRUE),
               iqr = fmt_ci(iqr_lo(x), iqr_hi(x), 2),
               pct_2.5_97.5 = fmt_ci(q025(x), q975(x), 2), stringsAsFactors = FALSE)
  })
  # split by model
  out <- lapply(names(measures), function(v) {
    d %>% group_by(model) %>%
      summarise(measure = measures[v], median = median(.data[[v]], na.rm = TRUE),
                iqr = fmt_ci(iqr_lo(.data[[v]]), iqr_hi(.data[[v]]), 2),
                pct_2.5_97.5 = fmt_ci(q025(.data[[v]]), q975(.data[[v]]), 2), .groups = "drop")
  })
  bind_rows(out) %>% select(measure, model, median, iqr, pct_2.5_97.5)
}

# Population- and cases-weighted pooled % reduction (alongside the median).
weighted_impact <- function(scn_cea, prep) {
  permo <- scn_cea %>% filter(vacc_cov == 0.80, tau == 8) %>%
    mutate(pct = ifelse(s_ch_tot > 0, 100 * s_ch_averted_tot / s_ch_tot, 0)) %>%
    group_by(study_id, model) %>%
    summarise(pct = median(pct), pop = first(population), cases = first(s_ch_tot), .groups = "drop")
  permo %>% group_by(model) %>%
    summarise(median_pct = median(pct, na.rm = TRUE),
              pop_weighted_pct = {ok <- is.finite(pop) & is.finite(pct)
                                  sum(pop[ok] * pct[ok]) / sum(pop[ok])},
              cases_weighted_pct = {ok <- is.finite(cases) & is.finite(pct)
                                    sum(cases[ok] * pct[ok]) / sum(cases[ok])}, .groups = "drop")
}

rcfg_base_tau <- function() 8L   # set by orchestrator if needed

# ============================================================================
# Part C — Cost-effectiveness decomposition (c8)
# Deaths averted, YLL/YLD split (%YLD), cost/DALY vs 0.5x (Ochalek), 1x, 3x GDP.
# ============================================================================
part_C_ce <- function(scn_cea, cfg) {
  permo <- scn_cea %>%
    group_by(study_id, model, tau, vacc_cov) %>%
    summarise(deaths_averted = median(death_averted_tot, na.rm = TRUE),
              yll = median(yll_averted, na.rm = TRUE),
              yld = median(yld_averted, na.rm = TRUE),
              daly = median(daly_averted, na.rm = TRUE),
              cost_per_daly = median(cost_per_daly_averted, na.rm = TRUE),
              gdp = first(gdp), .groups = "drop")
  mult <- as.numeric(unlist(cfg$thresholds$gdp_multipliers))   # yaml -> numeric
  permo %>% group_by(model, tau, vacc_cov) %>%
    summarise(
      deaths_averted_median = median(deaths_averted, na.rm = TRUE),
      yll_median = median(yll, na.rm = TRUE),
      yld_median = median(yld, na.rm = TRUE),
      pct_yld = 100 * median(yld, na.rm = TRUE) /
                (median(yld, na.rm = TRUE) + median(yll, na.rm = TRUE)),
      cost_per_daly_median = median(cost_per_daly, na.rm = TRUE),
      frac_ce_0.5xGDP = mean(cost_per_daly < mult[1] * gdp, na.rm = TRUE),
      frac_ce_1xGDP   = mean(cost_per_daly < mult[2] * gdp, na.rm = TRUE),
      frac_ce_3xGDP   = mean(cost_per_daly < mult[3] * gdp, na.rm = TRUE),
      .groups = "drop")
}

fig_ce_thresholds <- function(ce_tab, scn_cea, cfg) {
  d <- ce_tab %>% filter(vacc_cov == 0.80)
  gdp_med <- median(scn_cea$gdp, na.rm = TRUE)
  thr <- data.frame(label = c("0.5x GDP (Ochalek)", "1x GDP", "3x GDP"),
                    y = as.numeric(unlist(cfg$thresholds$gdp_multipliers)) * gdp_med)
  ggplot(d, aes(tau, cost_per_daly_median, colour = model)) +
    geom_line() + geom_point() +
    geom_hline(data = thr, aes(yintercept = y, linetype = label), colour = "grey40") +
    scale_y_log10() +
    scale_colour_manual(values = c(renewal = "#0072B2", static = "#D55E00")) +
    labs(x = "ORI delay (weeks)", y = "Cost per DALY averted (log, USD)",
         colour = "Model", linetype = "Threshold (median GDP)",
         title = "Cost-effectiveness vs thresholds (80% coverage)",
         caption = sprintf("Median country GDP/capita = $%.0f; Ochalek opportunity cost = 0.5x GDP",
                           gdp_med)) +
    theme_minimal(base_size = 10)
}

# ============================================================================
# Part E — Spatial / hotspot targeting (ILLUSTRATIVE, c8)
# Vaccinating the sub-area holding phi_cases of cases needs only phi_pop of the
# population. Doses (hence vacc cost) scale by phi_pop; averted cases / DALYs /
# treatment cost scale by phi_cases. NOT a fitted estimate.
# ============================================================================
part_E_spatial <- function(scn_cea, cfg, tau0 = 8, cov0 = 0.80) {
  base <- scn_cea %>% filter(model == "renewal", tau == tau0, vacc_cov == cov0) %>%
    group_by(study_id) %>%
    summarise(averted = median(s_ch_averted_tot, na.rm = TRUE),
              daly = median(daly_averted, na.rm = TRUE),
              coi = median(coi_averted, na.rm = TRUE),
              cod = median(cod_averted, na.rm = TRUE),
              vacc_cost = median(vacc_cost, na.rm = TRUE),
              vacc_dose = median(vacc_dose, na.rm = TRUE), .groups = "drop")
  phi_cases <- cfg$spatial$phi_cases
  rows <- lapply(c(1.0, as.numeric(unlist(cfg$spatial$phi_pop_sweep))), function(phi_pop) {
    targeted <- phi_pop < 1.0
    fc <- if (targeted) phi_cases else 1.0
    d <- base %>% mutate(
      vacc_cost_s = vacc_cost * phi_pop,
      daly_s = daly * fc,
      net_cost_s = vacc_cost_s - (coi + cod) * fc,
      cost_per_daly_s = ifelse(daly_s > 0, net_cost_s / daly_s, NA),
      per1000_s = ifelse(vacc_dose * phi_pop > 0, 1000 * averted * fc / (vacc_dose * phi_pop), NA))
    data.frame(strategy = if (targeted) sprintf("targeted (phi_pop=%.2f)", phi_pop) else "whole-unit",
               phi_pop = phi_pop, phi_cases = if (targeted) phi_cases else 1.0,
               cost_per_daly_median = median(d$cost_per_daly_s, na.rm = TRUE),
               cases_per_1000_median = median(d$per1000_s, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

fig_spatial_targeting <- function(sp_tab, cfg) {
  ggplot(sp_tab, aes(reorder(strategy, phi_pop), cost_per_daly_median,
                     fill = phi_pop < 1.0)) +
    geom_col() +
    scale_fill_manual(values = c("FALSE" = "#D55E00", "TRUE" = "#0072B2"),
                      labels = c("Whole-unit", "Targeted"), name = NULL) +
    labs(x = NULL, y = "Cost per DALY averted (USD, median)",
         title = "Illustrative: whole-unit vs hotspot-targeted delivery",
         caption = sprintf("ILLUSTRATIVE: %.0f%% of cases captured in phi_pop of population (Poncin 2022 Harare anchor)",
                           100 * cfg$spatial$phi_cases)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ============================================================================
# Part C2 — IVE coverage-dependence (c2)
# IVE(f) = (OVE - f*TVE)/(1-f); zero-crossing at f* = OVE/TVE.
# Shows the submitted fixed IVE=0.17 is inconsistent with the algebra above f*.
# ============================================================================
part_C2_ive <- function(cfg) {
  a <- cfg$ive_algebra
  cov_grid <- as.numeric(unlist(a$coverage_grid))
  ive_f <- function(f) (a$OVE - f * a$TVE) / (1 - f)
  fstar <- a$OVE / a$TVE
  data.frame(
    coverage = cov_grid,
    ive_algebraic = round(sapply(cov_grid, ive_f), 3),
    ive_fixed_submitted = a$ive_fixed,
    zero_crossing_coverage = round(fstar, 3),
    note = ifelse(cov_grid > fstar, "algebraic IVE < 0 (implausible)", "algebraic IVE > 0"),
    stringsAsFactors = FALSE)
}

fig_ive_coverage <- function(cfg) {
  a <- cfg$ive_algebra
  cov_grid <- as.numeric(unlist(a$coverage_grid))
  fstar <- a$OVE / a$TVE
  grid <- data.frame(f = seq(0, 0.95, by = 0.01))
  grid$ive <- (a$OVE - grid$f * a$TVE) / (1 - grid$f)
  pts <- data.frame(f = cov_grid,
                    ive = (a$OVE - cov_grid * a$TVE) / (1 - cov_grid))
  ggplot(grid, aes(f, ive)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_hline(yintercept = a$ive_fixed, linetype = "dashed", colour = "#D55E00") +
    geom_vline(xintercept = fstar, linetype = "dotted", colour = "grey40") +
    geom_line(colour = "#0072B2", linewidth = 1) +
    geom_point(data = pts, size = 2, colour = "#0072B2") +
    scale_x_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-1, 0.6)) +
    annotate("text", x = fstar, y = 0.5, hjust = -0.05,
             label = sprintf("zero at f=%.0f%%", 100 * fstar), size = 3) +
    annotate("text", x = 0.05, y = a$ive_fixed + 0.04, hjust = 0,
             label = sprintf("submitted fixed IVE = %.0f%%", 100 * a$ive_fixed),
             colour = "#D55E00", size = 3) +
    labs(x = "Vaccine coverage f", y = "Algebraic IVE(f)",
         title = "Indirect VE implied by the algebra is coverage-dependent",
         subtitle = "IVE(f) = (OVE - f*TVE)/(1-f), Khanam Vi-TT: OVE=0.57, TVE=0.85") +
    theme_minimal(base_size = 10)
}
