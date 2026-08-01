# =============================================================================
# High-impact main-text figures for the PPV adjustment in the ORI analysis.
#   Rscript revision/make_main_text_figures.R
# Reads the fitted latent-class draws/tables (9-study paired set, Vallenas
# excluded) and the renewal PPV-propagation outputs. Colorblind-safe Okabe-Ito
# palette (validated); static print figures at 300 dpi.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({ library(ggplot2); library(posterior) })
outdir <- "revision/main_text_figures"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
il <- function(x) 1/(1+exp(-x))

# --- palette (Okabe-Ito; validated CVD-safe) ---------------------------------
BLUE <- "#0072B2"   # true typhoid / community
VERM <- "#D55E00"   # suspected / surveillance-observed
GREY <- "#7A7A7A"
grade_cols <- c(`negligible/low` = "#9ECAE1", moderate = "#4292C6", major = "#08519C")  # sequential
base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10.5),
        axis.title = element_text(colour = "grey25"),
        legend.position = "top")

draws <- as_draws_matrix(readRDS("latent_class_ppv/outputs/final_draws.rds"))
pars  <- read.csv("latent_class_ppv/tables/final_parameters.csv", check.names = FALSE)
pic   <- read.csv("latent_class_ppv/tables/final_pi_community.csv", check.names = FALSE, stringsAsFactors = FALSE)
eff   <- read.csv("renewal/tables/tab_ppv_effect.csv", check.names = FALSE, stringsAsFactors = FALSE)
ssens <- read.csv("renewal/tables/tab_ppv_selection_sensitivity.csv", check.names = FALSE, stringsAsFactors = FALSE)
src   <- read.csv("source_decomposition/tables/tab_additive_vs_theta_weighted_pooled.csv",
                  check.names = FALSE, stringsAsFactors = FALSE)
grade_from_bias <- c(`0` = "negligible/low", `1` = "moderate", `2` = "major")

# =============================================================================
# FIGURE 1 — Blood-culture sensitivity increases with cultured volume and is
# far below the conventional 0.5 assumption.
# =============================================================================
a0 <- as.numeric(draws[, "alpha0"]); a1 <- as.numeric(draws[, "alpha1"])
gv <- seq(1.5, 15, 0.1)
curve <- data.frame(vol = gv,
  med = sapply(gv, function(v) median(il(a0 + a1*log(v)))),
  lo  = sapply(gv, function(v) quantile(il(a0 + a1*log(v)), .05)),
  hi  = sapply(gv, function(v) quantile(il(a0 + a1*log(v)), .95)))
H <- read.csv("latent_class_ppv/data/mogasale2016_paired_bc_bmc_seed.csv", check.names = FALSE)
vh <- H$blood_vol_mL; vh[is.na(vh)] <- median(vh, na.rm = TRUE)
pts <- data.frame(vol = vh,
                  se_raw = (H$a_BCpos_BMpos + H$b_BCpos_BMneg) /
                           (H$a_BCpos_BMpos + H$b_BCpos_BMneg + H$c_BCneg_BMpos))
se5 <- pars$med[pars$param == "Se_BC_5mL"]

fig1 <- ggplot(curve, aes(vol, med)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = GREY) +
  annotate("text", x = 14.5, y = 0.47, label = "conventional 0.5 assumption",
           hjust = 1, size = 3.2, colour = GREY) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = BLUE, alpha = 0.18) +
  geom_line(colour = BLUE, linewidth = 1.1) +
  geom_point(data = pts, aes(vol, se_raw), inherit.aes = FALSE,
             colour = "grey30", size = 2.2, alpha = 0.8) +
  annotate("point", x = 5, y = se5, colour = VERM, size = 2.6) +
  annotate("text", x = 5.25, y = se5 + 0.045,
           label = sprintf("%.2f at 5 mL", se5), hjust = 0, size = 3.4, colour = VERM) +
  coord_cartesian(ylim = c(0, 1), xlim = c(1.5, 15)) +
  scale_x_continuous(breaks = c(2,5,8,10,15)) +
  labs(x = "Cultured blood volume (mL)", y = expression("Blood-culture sensitivity  " * Se[BC]),
       title = "Blood culture detects only about two-thirds of true typhoid",
       subtitle = "Posterior median and 90% CrI (line/band); points are historic paired-study detection proportions") +
  base_theme
ggsave(file.path(outdir, "fig1_se_bc_volume.png"), fig1, width = 7, height = 4.6, dpi = 300)

# =============================================================================
# FIGURE 2 — Community PPV is heterogeneous and mostly low; culture positivity
# is not PPV. Observed positivity (hollow) vs fitted PPV among all suspected.
# =============================================================================
d2 <- data.frame(study = pic$study, positivity = pic$positivity,
                 pi_med = pic$pi_med, pi_lo = pic$pi_lo, pi_hi = pic$pi_hi,
                 grade = grade_from_bias[as.character(pic$bias_score)],
                 stringsAsFactors = FALSE)
d2$grade <- factor(d2$grade, levels = names(grade_cols))
d2$study <- factor(d2$study, levels = d2$study[order(d2$pi_med)])
pop_med <- pars$med[pars$param == "pi_population_median"]
pop_lo  <- pars$lo.5.[pars$param == "pi_population_median"]
pop_hi  <- pars$hi.95.[pars$param == "pi_population_median"]

fig2 <- ggplot(d2, aes(y = study)) +
  annotate("rect", xmin = pop_lo, xmax = pop_hi, ymin = -Inf, ymax = Inf,
           fill = BLUE, alpha = 0.07) +
  annotate("segment", x = pop_med, xend = pop_med, y = -Inf, yend = Inf,
           colour = BLUE, linetype = "dotted") +
  geom_segment(aes(x = positivity, xend = pi_med, y = study, yend = study),
               colour = "grey75", linewidth = 0.5) +
  geom_point(aes(x = positivity), shape = 21, fill = "white", colour = "grey45", size = 2.6) +
  geom_linerange(aes(xmin = pi_lo, xmax = pi_hi, colour = grade), linewidth = 0.8) +
  geom_point(aes(x = pi_med, colour = grade), size = 3) +
  scale_colour_manual(values = grade_cols, name = "Verification-adjustment grade", drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = pop_med, y = 0.55, label = "population PPV", colour = BLUE,
           size = 3.1, hjust = -0.05) +
  labs(x = "PPV of the suspected case definition (among all suspected cases)", y = NULL,
       title = "True-typhoid PPV is low and varies widely across outbreaks",
       subtitle = "Filled = fitted PPV (90% CrI); hollow = raw culture positivity k/n. Positivity is not PPV.") +
  base_theme + theme(panel.grid.major.y = element_blank())
ggsave(file.path(outdir, "fig2_community_ppv.png"), fig2, width = 8.2, height = 4.8, dpi = 300)

# =============================================================================
# FIGURE 3 — Translating suspected-case impact into true typhoid burden.
# 3 panels: (A) cases averted, (B) % reduction true vs observable,
#           (C) robustness of the absolute burden to the selection prior.
# =============================================================================
tt <- function(q) as.numeric(eff$true_typhoid[eff$quantity == q |
        grepl(substr(q,1,20), eff$quantity, fixed = TRUE)][1])
avert_true <- src$additive_cases_median
avert_susp <- src$historical_weighted_cases_median
red_true   <- src$additive_true_pct_median
red_obs    <- src$additive_observed_pct_median

pA <- data.frame(lab = c("Historical\ntheta-weighted", "Additive source\n+ propagated"),
                 val = c(avert_susp, avert_true),
                 col = c(GREY, BLUE))
pA$lab <- factor(pA$lab, levels = pA$lab)
gA <- ggplot(pA, aes(lab, val, fill = col)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = format(round(val), big.mark = ",")), vjust = -0.4, size = 3.8, colour = "grey20") +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(x = NULL, y = "Cases averted", title = "A  Paired structural comparison",
       subtitle = "Pooled median true typhoid cases") +
  base_theme + theme(legend.position = "none")

pB <- data.frame(lab = c("True typhoid", "Visible in suspected\nsurveillance"),
                 val = c(red_true, red_obs), col = c(BLUE, VERM))
pB$lab <- factor(pB$lab, levels = pB$lab)
gB <- ggplot(pB, aes(lab, val, fill = col)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", val)), vjust = -0.4, size = 3.8, colour = "grey20") +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(x = NULL, y = "Reduction (%)", title = "B  Percent reduction",
       subtitle = "Surveillance dilutes the visible effect") +
  base_theme + theme(legend.position = "none")

scen_lab <- c(no_enrichment = "No enrichment", no_grade = "No grade effect",
              skeptical = "Skeptical priors", binary_grade = "Binary grade",
              base = "Base model", permissive = "Permissive priors")
ssens$label <- ifelse(ssens$scenario %in% names(scen_lab), scen_lab[ssens$scenario], ssens$scenario)
ssens$label <- factor(ssens$label, levels = ssens$label[order(ssens$true_cases_averted)])
ssens$scenario <- ssens$label
rng <- range(ssens$true_cases_averted)
gC <- ggplot(ssens, aes(true_cases_averted, scenario)) +
  annotate("rect", xmin = rng[1], xmax = rng[2], ymin = -Inf, ymax = Inf,
           fill = BLUE, alpha = 0.07) +
  geom_point(colour = BLUE, size = 2.8) +
  geom_text(aes(label = format(round(true_cases_averted), big.mark = ",")),
            hjust = -0.25, size = 3, colour = "grey30") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.18))) +
  labs(x = "True typhoid cases averted", y = NULL,
       title = "C  Observation-layer sensitivity",
       subtitle = sprintf("Conditional pure-renewal workflow: %s-%s cases; true reduction %.1f-%.1f%%",
                          format(rng[1], big.mark=","), format(rng[2], big.mark=","),
                          min(ssens$true_pct_reduction), max(ssens$true_pct_reduction))) +
  base_theme

ok <- requireNamespace("patchwork", quietly = TRUE)
if (ok) {
  library(patchwork)
  fig3 <- (gA | gB) / gC +
    patchwork::plot_annotation(
      title = "Additive source-plus-propagated ORI model and structural comparisons",
      theme = theme(plot.title = element_text(face = "bold", size = 14)))
  ggsave(file.path(outdir, "fig3_burden_translation.png"), fig3, width = 9, height = 8, dpi = 300)
} else {
  ggsave(file.path(outdir, "fig3A_cases_averted.png"), gA, width = 4.5, height = 4.2, dpi = 300)
  ggsave(file.path(outdir, "fig3B_pct_reduction.png"), gB, width = 4.5, height = 4.2, dpi = 300)
  ggsave(file.path(outdir, "fig3C_robustness.png"), gC, width = 7, height = 4, dpi = 300)
}

cat("Figures written to", outdir, "\n")
cat(sprintf("Se_BC(5mL)=%.2f  pop PPV=%.2f  true averted=%s  true red=%.1f%%  obs red=%.1f%%\n",
            se5, pop_med, format(round(avert_true), big.mark=","), red_true, red_obs))
cat("patchwork combined:", ok, "\n")
