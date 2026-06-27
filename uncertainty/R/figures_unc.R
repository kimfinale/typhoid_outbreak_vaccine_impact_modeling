# =============================================================================
# Figures for the rigorous (propagated-uncertainty) R_t analysis.
# =============================================================================
suppressMessages({ library(ggplot2); library(dplyr) })

save_fig_u <- function(plot, name, cfg, width = 9, height = 6, dpi = 300) {
  base <- file.path(cfg$paths$figures, name)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = dpi)
  ok <- tryCatch({ ggsave(paste0(base, ".pdf"), plot, width = width, height = height,
                          device = cairo_pdf); TRUE }, error = function(e) FALSE)
  if (!ok) ok <- tryCatch({ ggsave(paste0(base, ".pdf"), plot, width = width,
                                   height = height); TRUE }, error = function(e) FALSE)
  if (!ok) warning("Could not write ", base, ".pdf; PNG written.")
}

# Observed incidence + R_t with the full MC uncertainty ribbon, per outbreak.
fig_Rt_uncertainty <- function(prep, ribbon, base_rt, cfg, ncol = 4) {
  stopifnot(requireNamespace("patchwork", quietly = TRUE))
  units <- lapply(names(prep$series), function(sid) {
    inc <- prep$series[[sid]]; wk <- seq_along(inc)
    rb <- ribbon[ribbon$study_id == sid, ]
    br <- base_rt[base_rt$study_id == sid, ]
    p_inc <- ggplot(data.frame(week = wk, cases = inc), aes(week, cases)) +
      geom_col(fill = "grey55", width = 0.9) +
      labs(title = sid, x = NULL, y = "cases") +
      theme_minimal(base_size = 7) +
      theme(plot.title = element_text(size = 7.5, face = "bold"),
            axis.text.x = element_blank(), plot.margin = margin(2, 4, 0, 2))
    p_rt <- ggplot(rb, aes(week, rt_median)) +
      geom_ribbon(aes(ymin = rt_lo, ymax = rt_hi), fill = "#0072B2", alpha = 0.2) +
      geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
      geom_line(colour = "#0072B2", linewidth = 0.5) +
      geom_line(data = br, aes(week, rt_base), colour = "black", linetype = "dashed",
                linewidth = 0.4, inherit.aes = FALSE) +
      coord_cartesian(ylim = c(0, 6)) +
      labs(x = "week", y = expression(R[t])) +
      theme_minimal(base_size = 7) + theme(plot.margin = margin(0, 4, 2, 2))
    patchwork::wrap_plots(p_inc, p_rt, ncol = 1, heights = c(1, 1))
  })
  patchwork::wrap_plots(units, ncol = ncol) +
    patchwork::plot_annotation(
      title = "R_t with propagated GI / asymptomatic / PPV uncertainty",
      subtitle = "Blue = MC median + 95% credible band; black dashed = point-estimate (base GI). Grey bars = observed incidence.",
      theme = theme(plot.title = element_text(size = 12, face = "bold"),
                    plot.subtitle = element_text(size = 9)))
}

# Tornado of partial rank correlation coefficients (global sensitivity).
fig_prcc_tornado <- function(prcc_df, cfg) {
  d <- prcc_df %>% mutate(param = reorder(param, abs(prcc)))
  ggplot(d, aes(prcc, param, fill = prcc > 0)) +
    geom_col() +
    geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = c("FALSE" = "#D55E00", "TRUE" = "#0072B2"), guide = "none") +
    labs(x = "Partial rank correlation with pooled % reduction", y = NULL,
         title = "Global sensitivity: what drives the impact estimate",
         subtitle = "Asymptomatic-transmission fraction (f_a, c_a) and psi dominate; GI shape is secondary; ppv_level/s_cult ~ 0 (level-invariant)") +
    theme_minimal(base_size = 10)
}

# Method comparison: exact-MC vs full-Bayesian % reduction per outbreak.
fig_bayes_compare <- function(summ_per, bayes_per, baseline, cfg) {
  d <- summ_per %>%
    transmute(study_id, method = "Exact + MC", med = rho_median, lo = rho_lo, hi = rho_hi) %>%
    bind_rows(bayes_per %>%
      transmute(study_id, method = "Full Bayesian", med = rho_median, lo = rho_lo, hi = rho_hi)) %>%
    left_join(baseline, by = "study_id") %>%
    mutate(study_id = reorder(study_id, med))
  ggplot(d, aes(med, study_id, colour = method)) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.4, position = position_dodge(width = 0.6)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    geom_point(aes(x = rho_base, y = study_id), colour = "black", shape = 4, size = 2,
               inherit.aes = FALSE) +
    scale_colour_manual(values = c("Exact + MC" = "#0072B2", "Full Bayesian" = "#009E73")) +
    labs(x = "% cases averted (median, 95% CrI)", y = NULL, colour = "Method",
         title = "Exact-reconstruction + MC vs full Bayesian renewal model",
         subtitle = "black x = point estimate (base GI). Agreement = conclusions robust to inference approach.") +
    theme_minimal(base_size = 10)
}

# Before/after forest: point-estimate vs full-uncertainty % reduction per outbreak.
fig_forest_uncertainty <- function(summ_per, baseline, cfg) {
  d <- summ_per %>% left_join(baseline, by = "study_id") %>%
    mutate(study_id = reorder(study_id, rho_median))
  ggplot(d, aes(rho_median, study_id)) +
    geom_errorbarh(aes(xmin = rho_lo, xmax = rho_hi), height = 0.3, colour = "#0072B2") +
    geom_point(colour = "#0072B2", size = 2) +
    geom_point(aes(x = rho_base), colour = "black", shape = 4, size = 2) +
    labs(x = "% cases averted (median, 95% credible interval)", y = NULL,
         title = "Per-outbreak % reduction with full uncertainty",
         subtitle = "Blue = MC median + 95% CrI (GI+asymptomatic+PPV); black x = point-estimate (base GI)") +
    theme_minimal(base_size = 10)
}
