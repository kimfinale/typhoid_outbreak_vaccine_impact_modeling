#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

dir.create("revision/figures", recursive = TRUE, showWarnings = FALSE)

pi_est <- read.csv("latent_class_ppv/tables/final_pi_community.csv", check.names = FALSE)
sel <- read.csv("latent_class_ppv/tables/final_selection_outbreak.csv", check.names = FALSE)
audit <- read.csv("latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv", check.names = FALSE)
impact <- read.csv("renewal/tables/tab_ppv_selection_sensitivity.csv", check.names = FALSE)
effect <- read.csv("renewal/tables/tab_ppv_effect.csv", check.names = FALSE)
source_cmp <- read.csv(
  "source_decomposition/tables/tab_additive_vs_theta_weighted_pooled.csv",
  check.names = FALSE)

dm <- posterior::as_draws_matrix(readRDS("latent_class_ppv/outputs/final_draws.rds"))
tested <- do.call(rbind, lapply(seq_len(nrow(pi_est)), function(i) {
  p <- dm[, paste0("pi[", i, "]")]
  q0 <- dm[, paste0("q0[", i, "]")]
  q1 <- dm[, paste0("q1[", i, "]")]
  pt <- p * q1 / (p * q1 + (1 - p) * q0)
  data.frame(study = pi_est$study[i], pi_tested_med = median(pt),
             pi_tested_lo = unname(quantile(pt, .05)),
             pi_tested_hi = unname(quantile(pt, .95)))
}))

grade <- audit %>%
  filter(include_for_main_ppv_model) %>%
  select(study, verification_adjustment_grade)

e <- pi_est %>%
  left_join(tested, by = "study") %>%
  left_join(grade, by = "study")

ord <- e %>% arrange(pi_med) %>% pull(study)
e$study <- factor(e$study, levels = ord)

long <- bind_rows(
  e %>% transmute(study, grade = verification_adjustment_grade,
                  estimand = "Observed blood-culture positivity",
                  med = positivity, lo = positivity, hi = positivity),
  e %>% transmute(study, grade = verification_adjustment_grade,
                  estimand = "PPV among cultured suspects",
                  med = pi_tested_med, lo = pi_tested_lo, hi = pi_tested_hi),
  e %>% transmute(study, grade = verification_adjustment_grade,
                  estimand = "PPV among all suspected cases",
                  med = pi_med, lo = pi_lo, hi = pi_hi)
)
long$estimand <- factor(long$estimand, levels = c(
  "Observed blood-culture positivity", "PPV among cultured suspects",
  "PPV among all suspected cases"))

connect <- long %>%
  group_by(study) %>% summarise(x = min(med), xend = max(med), .groups = "drop")

cols <- c("Observed blood-culture positivity" = "#8A8F98",
          "PPV among cultured suspects" = "#D17B0F",
          "PPV among all suspected cases" = "#007C91")
shapes <- c("Observed blood-culture positivity" = 16,
            "PPV among cultured suspects" = 17,
            "PPV among all suspected cases" = 18)

p1 <- ggplot(long, aes(med, study)) +
  geom_segment(data = connect, aes(x = x, xend = xend, y = study, yend = study),
               inherit.aes = FALSE, colour = "#D9DDE2", linewidth = 0.8) +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = estimand), orientation = "y",
                width = 0.12, linewidth = 0.65, alpha = 0.9) +
  geom_point(aes(colour = estimand, shape = estimand), size = 3) +
  scale_colour_manual(values = cols) + scale_shape_manual(values = shapes) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Fraction of suspected cases", y = NULL,
       title = "Culture positivity is not the PPV of a suspected-case definition",
       subtitle = "Points are medians; horizontal intervals are 90% credible intervals") +
  guides(colour = guide_legend(title = NULL, nrow = 3, byrow = TRUE),
         shape = guide_legend(title = NULL, nrow = 3, byrow = TRUE)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "#4E5968"))

ggsave("revision/figures/fig_main_ppv_evidence.png", p1, width = 8.2, height = 5.7, dpi = 400)
ggsave("revision/figures/fig_main_ppv_evidence.pdf", p1, width = 8.2, height = 5.7)

cmp_cases <- data.frame(
  model = c("Additive source + propagated", "Historical theta-weighted"),
  med = c(source_cmp$additive_cases_median,
          source_cmp$historical_weighted_cases_median),
  lo = c(source_cmp$additive_cases_p025,
         source_cmp$historical_weighted_cases_p025),
  hi = c(source_cmp$additive_cases_p975,
         source_cmp$historical_weighted_cases_p975))
cmp_cases$model <- factor(cmp_cases$model, levels = rev(cmp_cases$model))
p2a <- ggplot(cmp_cases, aes(med, model, colour = model)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.12, linewidth = 0.8) +
  geom_point(size = 3.4) +
  geom_text(aes(label = comma(round(med))), hjust = -0.35, size = 3.1) +
  scale_colour_manual(values = c("Additive source + propagated" = "#007C91",
                                 "Historical theta-weighted" = "#75808F"), guide = "none") +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = "True typhoid cases averted", y = NULL,
       title = "A Paired structural comparison",
       subtitle = "Median and 95% simulation interval") +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

red <- data.frame(
  model = rep(c("Additive", "Historical weighted"), each = 2),
  measure = rep(c("True typhoid", "Suspected-case surveillance"), 2),
  value = c(source_cmp$additive_true_pct_median,
            source_cmp$additive_observed_pct_median,
            source_cmp$historical_weighted_true_pct_median,
            source_cmp$historical_weighted_observed_pct_median))

p2b <- ggplot(red, aes(value, model, colour = measure, shape = measure)) +
  geom_line(aes(group = model), colour = "#D9DDE2", linewidth = 0.8) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c("True typhoid" = "#007C91",
                                 "Suspected-case surveillance" = "#D17B0F")) +
  scale_shape_manual(values = c("True typhoid" = 18,
                                "Suspected-case surveillance" = 16)) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 30),
                     breaks = seq(0, 30, 5)) +
  labs(x = "Percentage reduction", y = NULL,
       title = "B True versus surveillance-visible impact") +
  guides(colour = guide_legend(title = NULL), shape = guide_legend(title = NULL)) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

p2 <- p2a + p2b + plot_annotation(
  title = "Two-level additive ORI model versus the historical theta-weighted calculation",
  subtitle = sprintf("Pooled paired difference: %+.1f percentage points (%.1f to %.1f)",
                     source_cmp$delta_true_pp_median,
                     source_cmp$delta_true_pp_p025,
                     source_cmp$delta_true_pp_p975),
  tag_levels = "A") & theme(plot.title = element_text(face = "bold"))

ggsave("revision/figures/fig_main_ppv_impact.png", p2, width = 11.5, height = 5.7, dpi = 400)
ggsave("revision/figures/fig_main_ppv_impact.pdf", p2, width = 11.5, height = 5.7)

cat("Wrote main-manuscript PPV figures to revision/figures\n")
