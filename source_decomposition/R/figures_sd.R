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
         subtitle = "Most are MIXED -- the static<->renewal bracket, not pure static, is the right framing") +
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
    geom_errorbarh(aes(xmin = theta_paper_lo, xmax = theta_paper_hi), height = 0.01, alpha = 0.5) +
    geom_point(size = 2.5) +
    ggrepel::geom_text_repel(aes(label = study_id), size = 2.6, max.overlaps = 20) +
    scale_colour_manual(values = MODE_COL, name = "Paper mode") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "theta (paper, Stream A)", y = "theta (data, Stream B)",
         title = "Cross-validation of the source fraction theta",
         subtitle = "Dotted band = agreement within 0.2. Identifiability limits Stream B (data alone).") +
    theme_minimal(base_size = 10)
}

# Per-outbreak static<->renewal bracket with the theta point estimate.
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
         title = "Source decomposition: static<->renewal bracket with theta point estimate",
         subtitle = "Blue = renewal (theta=0); orange = static (theta=1); diamond = theta point") +
    theme_minimal(base_size = 10)
}
