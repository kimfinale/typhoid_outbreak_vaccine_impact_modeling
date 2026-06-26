# =============================================================================
# Figures and tables for the renewal analysis (ggplot2). Each builder returns a
# ggplot; save_fig() writes PNG + PDF at >= 300 dpi.
# =============================================================================

suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

save_fig <- function(plot, name, cfg, width = 9, height = 6, dpi = 300) {
  base <- file.path(cfg$paths$figures, name)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  invisible(base)
}

COL_MODEL <- c(renewal = "#0072B2", static = "#D55E00")

# --- fig_Rt_panels: exact R_t + EpiEstim ribbon, small multiples -------------
fig_Rt_panels <- function(prep, cfg) {
  w <- gi_from_config(cfg)
  rows <- lapply(names(prep$series), function(sid) {
    inc <- prep$series[[sid]]
    exact <- reconstruct_Rt(inc, w)$Rt
    ee <- epiestim_rt(inc, cfg)
    d <- data.frame(study_id = sid, week = seq_along(inc), exact = exact)
    if (!is.null(ee)) d <- left_join(d, ee, by = "week")
    iw <- prep$meta[[sid]]$intervention_week
    d$intervention_week <- if (is.na(iw)) NA_real_ else iw
    d
  })
  df <- bind_rows(rows)
  ggplot(df, aes(week)) +
    {if ("lo" %in% names(df)) geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#0072B2", alpha = 0.2)} +
    {if ("mean_r" %in% names(df)) geom_line(aes(y = mean_r), colour = "#0072B2", linewidth = 0.4)} +
    geom_line(aes(y = exact), colour = "black", linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
    geom_vline(aes(xintercept = intervention_week), linetype = "dashed",
               colour = "red", linewidth = 0.3, na.rm = TRUE) +
    facet_wrap(~study_id, scales = "free", ncol = 3) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(x = "Week", y = expression(R[t]),
         title = "Reconstructed R_t (black = exact, blue = EpiEstim 95% CrI)") +
    theme_minimal(base_size = 9)
}

# --- fig_forest_pctreduction: static vs renewal per outbreak -----------------
fig_forest_pctreduction <- function(summ_base, pooled, cfg) {
  d <- summ_base %>% arrange(model, pct_reduction_median) %>%
    mutate(study_id = factor(study_id, levels = unique(study_id[model == "renewal"][order(
      summ_base$pct_reduction_median[summ_base$model == "renewal"])])))
  ren_pool <- pooled %>% filter(model == "renewal")
  ggplot(d, aes(pct_reduction_median, study_id, colour = model)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbarh(aes(xmin = pct_reduction_low, xmax = pct_reduction_high),
                   position = position_dodge(width = 0.5), height = 0.3) +
    geom_vline(data = ren_pool, aes(xintercept = pooled_pct_median),
               linetype = "dashed", colour = COL_MODEL["renewal"]) +
    geom_vline(data = ren_pool, aes(xintercept = popwt_pct),
               linetype = "dotted", colour = COL_MODEL["renewal"]) +
    scale_colour_manual(values = COL_MODEL) +
    labs(x = "% cases averted (median, 95% UI)", y = NULL, colour = "Model",
         title = "Per-outbreak % reduction: static vs renewal",
         subtitle = "dashed = renewal pooled median; dotted = population-weighted") +
    theme_minimal(base_size = 10)
}

# --- fig_amplification_vs_delay: renewal/static averted ratio vs delay --------
fig_amplification_vs_delay <- function(summ_grid, cfg) {
  amp <- summ_grid %>%
    select(study_id, model, tau, vacc_cov, s_ch_averted_median) %>%
    pivot_wider(names_from = model, values_from = s_ch_averted_median) %>%
    filter(vacc_cov == cfg$vaccine$coverage_base, static > 0) %>%
    mutate(ratio = renewal / static)
  pooled <- amp %>% group_by(tau) %>%
    summarise(ratio = median(ratio, na.rm = TRUE), .groups = "drop")
  ggplot(amp, aes(tau, ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_jitter(width = 0.2, alpha = 0.4, colour = "grey50", size = 1) +
    geom_line(data = pooled, colour = "#009E73", linewidth = 1) +
    geom_point(data = pooled, colour = "#009E73", size = 2) +
    labs(x = "ORI delay tau (weeks)", y = "Renewal / static cases averted",
         title = "Amplification of impact by response timing",
         subtitle = ">1: static underestimates (emergent herd effect)") +
    theme_minimal(base_size = 10)
}

# --- fig_eta_eff: per-outbreak eta_eff vs assumed 17% ------------------------
fig_eta_eff <- function(eta_tab, cfg) {
  d <- eta_tab %>% mutate(study_id = reorder(study_id, eta_eff_median))
  ggplot(d, aes(eta_eff_median, study_id)) +
    geom_vline(xintercept = cfg$vaccine$eta_static, linetype = "dashed", colour = "red") +
    geom_point(size = 2, colour = "#0072B2") +
    geom_errorbarh(aes(xmin = eta_eff_low, xmax = eta_eff_high), height = 0.3, colour = "#0072B2") +
    scale_x_continuous(labels = scales::percent) +
    labs(x = expression(eta[eff]), y = NULL,
         title = "Effective indirect VE implied by the renewal model",
         subtitle = paste0("red = assumed eta = ", 100 * cfg$vaccine$eta_static, "%")) +
    theme_minimal(base_size = 10)
}

# --- fig_policy_grid: 4 metrics by delay x coverage, renewal vs static --------
fig_policy_grid <- function(summ_grid, cfg) {
  metrics <- c(pct_reduction_median = "% cases averted",
               s_ch_averted_median = "Cases averted",
               cost_per_daly_averted_median = "Cost per DALY averted",
               case_averted_per_1000_median = "Cases averted / 1000 doses")
  d <- summ_grid %>%
    select(study_id, model, tau, vacc_cov, names(metrics)) %>%
    pivot_longer(names(metrics), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metrics[metric], levels = metrics)) %>%
    group_by(model, tau, vacc_cov, metric) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = "drop")
  ggplot(d, aes(tau, value, colour = model, linetype = factor(vacc_cov))) +
    geom_line() + geom_point(size = 1) +
    facet_wrap(~metric, scales = "free_y") +
    scale_colour_manual(values = COL_MODEL) +
    labs(x = "ORI delay tau (weeks)", colour = "Model", linetype = "Coverage",
         y = NULL, title = "Policy grid: impact by delay x coverage") +
    theme_minimal(base_size = 9)
}

# --- fig_gi_sensitivity / fig_psiT_sensitivity -------------------------------
fig_sensitivity <- function(sens_df, xvar, xlab, title) {
  ggplot(sens_df, aes(.data[[xvar]], pct_reduction_median, group = study_id)) +
    geom_line(alpha = 0.3, colour = "grey60") +
    stat_summary(aes(group = 1), fun = median, geom = "line",
                 colour = "#0072B2", linewidth = 1) +
    labs(x = xlab, y = "% cases averted (median)", title = title) +
    theme_minimal(base_size = 10)
}

# --- fig_curves_examples: observed / static CF / renewal CF for 2 outbreaks ---
fig_curves_examples <- function(prep, cfg, ids) {
  w <- gi_from_config(cfg)
  tau <- cfg$timing$delay_base_weeks; cov <- cfg$vaccine$coverage_base
  delta_w <- cfg$timing$immunity_onset_days / cfg$step_days +
    round(cfg$timing$campaign_duration_days / cfg$step_days / 2)
  te <- max(round(tau + delta_w), 1L)
  rows <- lapply(ids, function(sid) {
    inc <- prep$series[[sid]]
    psi <- cfg$vaccine$psi_mean
    ren <- renewal_counterfactual(inc, w, tau, te, cov, psi,
                                  c_shape = cfg$timing$c_shape, feedback = TRUE)$incidence_v
    P <- P_halloran(cov, psi, cfg$vaccine$eta_static)
    sta <- static_counterfactual(inc, te, P)
    data.frame(study_id = sid, week = seq_along(inc),
               Observed = inc, `Static CF` = sta, `Renewal CF` = ren,
               t_eff = te, check.names = FALSE)
  })
  df <- bind_rows(rows) %>%
    pivot_longer(c("Observed", "Static CF", "Renewal CF"),
                 names_to = "series", values_to = "cases")
  ggplot(df, aes(week, cases, colour = series)) +
    geom_line(linewidth = 0.7) +
    geom_vline(aes(xintercept = t_eff), linetype = "dashed", colour = "grey50") +
    facet_wrap(~study_id, scales = "free") +
    scale_colour_manual(values = c("Observed" = "black",
                                   "Static CF" = "#D55E00", "Renewal CF" = "#0072B2")) +
    labs(x = "Week", y = "Weekly cases", colour = NULL,
         title = "Observed vs counterfactual curves (base case)") +
    theme_minimal(base_size = 10)
}
