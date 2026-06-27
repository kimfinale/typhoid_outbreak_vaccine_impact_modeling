# =============================================================================
# Figures for the impact-drivers analysis.
# =============================================================================
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

save_fig_d <- function(plot, name, cfg, width = 9, height = 6, dpi = 300) {
  base <- file.path(cfg$paths$figures, name)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = dpi)
  ok <- tryCatch({ ggsave(paste0(base, ".pdf"), plot, width = width, height = height,
                          device = cairo_pdf); TRUE }, error = function(e) FALSE)
  if (!ok) tryCatch(ggsave(paste0(base, ".pdf"), plot, width = width, height = height),
                    error = function(e) warning("no pdf for ", name))
}

# Predictor importance: empirical |rho| (real) next to simulation deviance-explained.
fig_predictor_ranking <- function(assoc, sim_imp, outcome, cfg) {
  emp <- assoc %>% filter(outcome == !!outcome) %>%
    transmute(predictor, value = abs(rho), source = "Empirical |rho| (n=13)")
  sim <- sim_imp %>% filter(outcome == !!outcome) %>%
    transmute(predictor, value = dev_explained, source = "Simulation dev. explained")
  d <- bind_rows(emp, sim) %>% filter(predictor %in% unique(sim$predictor))
  ggplot(d, aes(reorder(predictor, value), value, fill = source)) +
    geom_col(position = "dodge") + coord_flip() +
    scale_fill_manual(values = c("Empirical |rho| (n=13)" = "#D55E00",
                                 "Simulation dev. explained" = "#0072B2")) +
    labs(x = NULL, y = "Importance", fill = NULL,
         title = paste0("Drivers of ", outcome),
         subtitle = "Prospective predictors ranked: real-data |Spearman rho| vs simulation deviance explained") +
    theme_minimal(base_size = 10)
}

# How the R_t / growth association shifts across vaccination delay.
fig_assoc_by_delay <- function(abd, cfg) {
  ggplot(abd, aes(tau, rho, colour = predictor)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_line(linewidth = 0.8) + geom_point() +
    labs(x = "Vaccination delay tau (weeks)", y = "Spearman rho with % reduction",
         colour = NULL, title = "Predictive value shifts with delay",
         subtitle = "R_t-at-vaccination tracks impact at the decision time; early growth fades as outbreaks burn out") +
    theme_minimal(base_size = 10)
}

# Impact map: expected primary outcome over (R_t-at-vaccination x delay) from the
# simulation surface, with the real outbreaks overlaid.
fig_impact_surface <- function(surface_model, real_dat, outcome, cfg) {
  grid <- expand.grid(Rt_at_tau = seq(0.3, 3.5, by = 0.05),
                      tau = seq(1, 12, by = 1),
                      log_pop = median(real_dat$log_pop, na.rm = TRUE))
  grid$z <- as.numeric(predict(surface_model, newdata = grid))
  ymax <- 3.5
  # winsorize real R_t to the surface range (a few early-week reconstructions
  # spike on tiny counts); mark winsorized points with a caret.
  real <- real_dat %>% filter(is.finite(.data[[outcome]]), is.finite(Rt_at_tau)) %>%
    mutate(Rt_disp = pmin(pmax(Rt_at_tau, 0.3), ymax),
           winsor = Rt_at_tau > ymax)
  ggplot(grid, aes(tau, Rt_at_tau)) +
    geom_raster(aes(fill = z), interpolate = TRUE) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "white") +
    geom_point(data = real, aes(tau, Rt_disp, size = .data[[outcome]], shape = winsor),
               fill = "white", colour = "black", alpha = 0.85) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
    scale_fill_viridis_c(option = "C") +
    coord_cartesian(ylim = c(0.3, ymax)) +
    labs(x = "Vaccination delay tau (weeks)", y = "R_t at vaccination",
         fill = outcome, size = "observed",
         title = paste0("Impact map: ", outcome),
         subtitle = "Simulation surface (fill); real outbreaks overlaid (triangles = R_t winsorized to 3.5). Dashed = R_t=1.") +
    theme_minimal(base_size = 10)
}

# Delay-sensitivity by R_t stratum (the interaction: every week counts more when growing).
fig_interaction <- function(ds, cfg) {
  ggplot(ds, aes(rt_stratum, slope_per_week, fill = rt_stratum)) +
    geom_col() + geom_hline(yintercept = 0, colour = "grey50") +
    scale_fill_viridis_d(option = "C", guide = "none") +
    labs(x = "R_t at vaccination (stratum)",
         y = "d(% reduction)/d(delay)  [pts per week]",
         title = "Delay matters more for still-growing outbreaks",
         subtitle = "Slope of impact vs delay within R_t strata (simulation)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# External validation: sim-trained prediction vs observed impact (real outbreaks).
fig_validation <- function(pred_df, outcome, cfg) {
  ggplot(pred_df, aes(predicted, actual)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(colour = "#0072B2", alpha = 0.7) +
    labs(x = "Predicted (simulation surface)", y = "Observed (real outbreak)",
         title = paste0("Do real outbreaks lie on the simulated surface? (", outcome, ")"),
         subtitle = "Each point = one outbreak x delay; dashed = 1:1") +
    theme_minimal(base_size = 10)
}
