# =============================================================================
# Phase 2 + 3 — real-data fits (M1 neutral, M2 article-informed) and the
# convergent-validity comparison. ONLY runs if the Phase-1 gate passed.
#   Rscript hierarchical_theta/run_phase2_3.R
# =============================================================================

ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(cmdstanr); library(posterior); library(ggplot2); library(dplyr) })
source("hierarchical_theta/R/hier_helpers.R")

cfg <- yaml::read_yaml("hierarchical_theta/config.yml")
if (!file.exists(file.path(cfg$paths$outputs, "GATE_PASSED.txt")))
  stop("Phase-1 gate has not passed (no GATE_PASSED.txt). Run run_phase1_gate.R first; ",
       "if it failed, the article classification is primary and real-data theta is not reported.")
rcfg <- yaml::read_yaml("renewal/config.yml")
set.seed(cfg$seed)
tabf <- function(n) file.path(cfg$paths$tables, n); outf <- function(n) file.path(cfg$paths$outputs, n)
cat("=== Phase 2/3: real-data hierarchical theta ===\n")

for (f in c("gi","renewal_core","data_prep","scenario","summarise"))
  source(file.path("renewal/R", paste0(f, ".R")))
prep <- prep_outbreaks(rcfg)
series <- prep$series; sids <- names(series)

# Article-derived mode/theta (Stream A) — covariate + validation target.
mt <- read.csv(cfg$mode_table, stringsAsFactors = FALSE, check.names = FALSE)
mt <- mt[match(sids, mt$study_id), ]
mode_int <- c(A = 1L, B = 2L, C = 3L)[mt$mode_paper]
article_theta <- (mt$theta_paper_lo + mt$theta_paper_hi) / 2

mod <- cmdstan_model(file.path(cfg$paths$stan, "renewal_source.stan"))
fit_one <- function(informed) {
  sdata <- build_stan_data(series, cfg, modes = mode_int, informed = informed)
  mod$sample(data = sdata, seed = cfg$stan_seed,
             chains = cfg$mcmc$chains, parallel_chains = cfg$mcmc$parallel_chains,
             iter_warmup = cfg$mcmc$iter_warmup, iter_sampling = cfg$mcmc$iter_sampling,
             adapt_delta = cfg$mcmc$adapt_delta, max_treedepth = cfg$mcmc$max_treedepth,
             refresh = 0)
}
cat("Fitting M1 (neutral)...\n"); M1 <- fit_one(0L)
cat("Fitting M2 (article-informed)...\n"); M2 <- fit_one(1L)

# --- Diagnostics -------------------------------------------------------------
diag_row <- function(fit, label) {
  s <- fit$summary(c("theta","tau_theta","gi_mean","phi")); d <- fit$diagnostic_summary()
  data.frame(model = label, max_rhat = max(s$rhat, na.rm = TRUE),
             min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
             divergences = sum(d$num_divergent),
             tau_theta = fit$summary("tau_theta")$median,
             gi_mean = fit$summary("gi_mean")$median)
}
diagtab <- rbind(diag_row(M1, "M1_neutral"), diag_row(M2, "M2_informed"))
write.csv(diagtab, tabf("tab_diagnostics.csv"), row.names = FALSE)
cat("Diagnostics:\n"); print(diagtab, row.names = FALSE)

# --- Posterior theta ---------------------------------------------------------
th <- function(fit) {
  m <- as_draws_matrix(fit$draws("theta"))
  data.frame(study_id = sids,
             med = apply(m, 2, median), lo = apply(m, 2, quantile, 0.05),
             hi = apply(m, 2, quantile, 0.95))
}
t1 <- th(M1); t2 <- th(M2)
rcols <- function(fit) as_draws_matrix(fit$draws("theta_realized")) |> apply(2, median)

out <- data.frame(
  study_id = sids, mode_paper = mt$mode_paper,
  article_theta_lo = mt$theta_paper_lo, article_theta_hi = mt$theta_paper_hi, article_theta,
  M1_theta_med = t1$med, M1_theta_lo = t1$lo, M1_theta_hi = t1$hi,
  M2_theta_med = t2$med, M2_theta_lo = t2$lo, M2_theta_hi = t2$hi,
  theta_realized_M1 = rcols(M1),
  agree_M1 = abs(t1$med - article_theta) <= 0.2)
write.csv(out, tabf("tab_theta_hierarchical.csv"), row.names = FALSE)

rank_cor <- suppressWarnings(cor(t1$med, article_theta, method = "spearman"))
cat(sprintf("\nConvergent validity (M1 neutral vs article): rank cor = %.2f\n", rank_cor))
cat(sprintf("tau_theta (M1) = %.2f (pooling strength)\n", diagtab$tau_theta[1]))

# --- Figures -----------------------------------------------------------------
MODE_COL <- c(A = "#D55E00", C = "#9467BD", B = "#0072B2")
cv <- ggplot(out, aes(article_theta, M1_theta_med, colour = mode_paper)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60") +
  geom_errorbar(aes(ymin = M1_theta_lo, ymax = M1_theta_hi), width = 0.01, alpha = 0.6) +
  geom_errorbar(aes(xmin = article_theta_lo, xmax = article_theta_hi),
                orientation = "y", width = 0.01, alpha = 0.6) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(aes(label = study_id), size = 2.6, max.overlaps = 20) +
  scale_colour_manual(values = MODE_COL, name = "Article mode") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "article theta (Stream A band)", y = "M1 neutral posterior theta (median, 90% CrI)",
       title = "Convergent validity: data-driven (neutral) vs article-derived theta",
       subtitle = sprintf("Spearman rank cor = %.2f. Agreement here = curve shape recovers field mode (no article prior).", rank_cor)) +
  theme_minimal(base_size = 10)
ggsave(file.path(cfg$paths$figures, "fig_theta_convergent_validity.png"), cv, width = 8, height = 7, dpi = 300)

mu_grand <- plogis(M1$summary("mu_theta")$median)
shr <- ggplot(out %>% mutate(study_id = reorder(study_id, M1_theta_med)), aes(M1_theta_med, study_id)) +
  geom_vline(xintercept = mu_grand, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = M1_theta_lo, xmax = M1_theta_hi),
                orientation = "y", width = 0.3, colour = "#0072B2") +
  geom_point(colour = "#0072B2", size = 2) +
  geom_point(aes(x = article_theta), colour = "black", shape = 4, size = 2) +
  labs(x = "theta", y = NULL, title = "Partial pooling of theta (M1 neutral)",
       subtitle = "blue = pooled posterior (median, 90% CrI); x = article theta; dashed = grand mean") +
  theme_minimal(base_size = 10)
ggsave(file.path(cfg$paths$figures, "fig_theta_shrinkage.png"), shr, width = 8, height = 6, dpi = 300)

# --- Posterior-predictive checks (a few outbreaks) ---------------------------
yrep1 <- as_draws_matrix(M1$draws("y_rep"))
for (o in seq_len(min(4, length(sids)))) {
  idx <- with(build_stan_data(series, cfg, mode_int, 0L), start[o]:(start[o] + Tn[o] - 1))
  qs <- apply(yrep1[, idx, drop = FALSE], 2, quantile, c(0.05, 0.5, 0.95))
  dpp <- data.frame(week = seq_along(idx), obs = series[[o]],
                    med = qs[2, ], lo = qs[1, ], hi = qs[3, ])
  pp <- ggplot(dpp, aes(week)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#0072B2", alpha = 0.25) +
    geom_line(aes(y = med), colour = "#0072B2") +
    geom_point(aes(y = obs), colour = "black", size = 1) +
    labs(title = paste0("PPC: ", sids[o]), x = "week", y = "incidence") +
    theme_minimal(base_size = 10)
  ggsave(file.path(cfg$paths$figures, paste0("fig_ppc_", gsub("[^A-Za-z0-9]+", "_", sids[o]), ".png")),
         pp, width = 6, height = 4, dpi = 300)
}

# Do not convert theta into a post hoc weighted static-plus-renewal impact.
# If the identifiability gate ever passes, posterior theta draws must be passed
# to the additive source counterfactual in renewal/R/renewal_core.R, where source
# and propagated incidence share one renewal recursion. The historical weighting
# is retained only in source_decomposition/ as a labeled comparator.

M1$save_object(outf("fit_M1.rds")); M2$save_object(outf("fit_M2.rds"))
writeLines(capture.output(sessionInfo()), outf("sessionInfo.txt"))
if (nzchar(Sys.which("quarto"))) {
  Sys.setenv(RENEWAL_ROOT = normalizePath(ROOT))
  system2("quarto", c("render", shQuote("hierarchical_theta/report_hierarchical_theta.qmd")),
          stdout = FALSE, stderr = FALSE)
}
cat("=== Phase 2/3 complete ===\n")
