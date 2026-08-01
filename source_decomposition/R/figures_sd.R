# =============================================================================
# Figures for the transmission-mode classification & source decomposition.
# =============================================================================
suppressMessages({ library(ggplot2); library(dplyr) })

save_fig_sd <- function(plot, name, cfg, width = 9, height = 6, dpi = 300) {
  base <- file.path(cfg$paths$figures, name)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = dpi)
  ok <- tryCatch({ ggsave(paste0(base, ".pdf"), plot, width = width, height = height,
                          device = cairo_pdf); TRUE }, error = function(e) FALSE)
  if (!ok) tryCatch(ggsave(paste0(base, ".pdf"), plot, width = width, height = height),
                    error = function(e) warning("no pdf for ", name))
}

MODE_COL <- c(A_common_source = "#D55E00", C_mixed = "#9467BD", B_propagated = "#0072B2",
              A = "#D55E00", C = "#9467BD", B = "#0072B2")

# Distribution of paper-based modes across the renewal outbreaks.
fig_mode_distribution <- function(mode_tab, cfg) {
  d <- mode_tab %>% mutate(mode = factor(mode_paper, levels = c("A", "C", "B")))
  ggplot(d, aes(mode, fill = mode)) +
    geom_bar() + scale_fill_manual(values = MODE_COL, guide = "none") +
    scale_x_discrete(labels = c(A = "A: common-source", C = "C: mixed", B = "B: propagated")) +
    labs(x = NULL, y = "outbreaks",
         title = "Transmission mode of the 13 renewal-analysed outbreaks (paper evidence)",
         subtitle = "Most are mixed, motivating source and propagated incidence within one additive recursion") +
    theme_minimal(base_size = 11)
}

# Stream A (paper) vs Stream B (data) source fraction theta.
fig_theta_crossvalidate <- function(mode_tab, cfg) {
  d <- mode_tab %>% mutate(theta_paper = (theta_paper_lo + theta_paper_hi) / 2)
  ggplot(d, aes(theta_paper, theta_data, colour = mode_paper)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_abline(slope = 1, intercept = 0.2, linetype = "dotted", colour = "grey75") +
    geom_abline(slope = 1, intercept = -0.2, linetype = "dotted", colour = "grey75") +
    geom_errorbar(aes(ymin = theta_data_lo, ymax = theta_data_hi), width = 0.01, alpha = 0.5) +
    geom_errorbar(aes(xmin = theta_paper_lo, xmax = theta_paper_hi),
                  orientation = "y", width = 0.01, alpha = 0.5) +
    geom_point(size = 2.5) +
    ggrepel::geom_text_repel(aes(label = study_id), size = 2.6, max.overlaps = 20) +
    scale_colour_manual(values = MODE_COL, name = "Paper mode") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "theta (paper, Stream A)", y = "theta (data, Stream B)",
         title = "Cross-validation of the source fraction theta",
         subtitle = "Dotted band = agreement within 0.2. Identifiability limits Stream B (data alone).") +
    theme_minimal(base_size = 10)
}

# Frozen historical static<->renewal outcome-weighted bracket.
fig_decomposition_bracket <- function(dec, mode_tab, cfg) {
  d <- dec %>% left_join(mode_tab[, c("study_id", "mode_paper")], by = "study_id") %>%
    mutate(study_id = reorder(study_id, pct_theta))
  ggplot(d, aes(y = study_id)) +
    geom_segment(aes(x = pct_renewal, xend = pct_static, yend = study_id),
                 colour = "grey70", linewidth = 1.5) +
    geom_point(aes(x = pct_renewal), colour = "#0072B2", size = 2) +
    geom_point(aes(x = pct_static), colour = "#D55E00", size = 2) +
    geom_point(aes(x = pct_theta, colour = mode_paper), size = 3, shape = 18) +
    scale_colour_manual(values = MODE_COL, name = "Mode") +
    labs(x = "% cases averted (tau=8, 80% coverage)", y = NULL,
         title = "Legacy comparator: theta-weighted static-renewal outcomes",
         subtitle = "Retained for audit only; theta is not used this way in the primary additive recursion") +
    theme_minimal(base_size = 10)
}

# Primary additive recursion versus the exact historical theta-weighted result.
fig_additive_vs_theta <- function(cmp, mode_tab, cfg) {
  d <- cmp %>%
    left_join(mode_tab[, c("study_id", "mode_paper")], by = "study_id") %>%
    mutate(study_id = reorder(study_id, pct_additive_median))
  ggplot(d, aes(y = study_id)) +
    geom_segment(
      aes(x = pct_theta_weighted_historical_median,
          xend = pct_additive_median, yend = study_id),
      colour = "grey70", linewidth = 1.2) +
    geom_errorbar(
      aes(xmin = pct_additive_structural_lo,
          xmax = pct_additive_structural_hi),
      orientation = "y", width = 0.24,
      colour = "#0072B2", linewidth = 0.5) +
    geom_point(
      aes(x = pct_theta_weighted_historical_median),
      colour = "#D55E00", shape = 17, size = 2.5) +
    geom_point(
      aes(x = pct_additive_median, colour = mode_paper),
      size = 2.8) +
    scale_colour_manual(values = MODE_COL, name = "Paper mode") +
    labs(
      x = "% true-typhoid cases averted (tau=8, 80% coverage)",
      y = NULL,
      title = "Additive renewal-with-source model versus theta weighting",
      subtitle = paste(
        "Circles = primary additive recursion; triangles = historical",
        "theta-weighted outcome; blue bars = structural theta range")) +
    theme_minimal(base_size = 10)
}

# Paired percentage-point difference against the historical implementation.
fig_additive_difference <- function(cmp, mode_tab, cfg) {
  d <- cmp %>%
    left_join(mode_tab[, c("study_id", "mode_paper")], by = "study_id")
  p <- ggplot(
    d,
    aes(theta, delta_true_pp_vs_historical_median,
        colour = mode_paper)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_errorbar(
      aes(ymin = delta_true_pp_vs_historical_p025,
          ymax = delta_true_pp_vs_historical_p975),
      width = 0.015, alpha = 0.65) +
    geom_point(size = 2.7) +
    scale_colour_manual(values = MODE_COL, name = "Paper mode") +
    labs(
      x = "Source fraction theta (used inside the additive recursion)",
      y = "Additive minus historical theta-weighted reduction (percentage points)",
      title = "Paired structural difference between source formulations",
      subtitle = "Same PPV, vaccine, timing, coverage, GI, and theta draw in both calculations") +
    theme_minimal(base_size = 10)
  if (requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(
      aes(label = study_id), size = 2.6, max.overlaps = 20)
  p
}
