# ============================================================================
# Summary figures + tables for the PPV (latent-class) and ORI (renewal) work.
#   Rscript manuscript/make_summary_figs.R
# Outputs -> manuscript/figures/*.png  and  manuscript/tables/*.csv
# ============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
figdir <- "manuscript/figures"; tabdir <- "manuscript/tables"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)
th <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13))
sav <- function(p, f, w = 8, h = 5) ggsave(file.path(figdir, f), p, width = w, height = h, dpi = 130)

par <- read.csv("latent_class_ppv/tables/final_parameters.csv", check.names = FALSE)
g <- function(nm) par[par$param == nm, c("med", "lo.5.", "hi.95.")]

# ---- FIG 1A: culture Se/Sp -- blood-culture Se vs volume + BM Se + Sp --------
a0 <- g("alpha0")$med; a1 <- g("alpha1")$med
vol <- seq(1, 15, 0.1); curve <- data.frame(vol, Se = plogis(a0 + a1 * log(vol)))
anch <- data.frame(vol = c(2, 5, 10),
                   Se = c(g("Se_BC_2mL")$med, g("Se_BC_5mL")$med, g("Se_BC_10mL")$med),
                   lo = c(g("Se_BC_2mL")$lo, g("Se_BC_5mL")$lo, g("Se_BC_10mL")$lo),
                   hi = c(g("Se_BC_2mL")$hi, g("Se_BC_5mL")$hi, g("Se_BC_10mL")$hi))
ant <- read.csv("latent_class_ppv/data/antillon2018_volume_sensitivity.csv")
sebm <- g("Se_BM"); spbc <- g("Sp_BC")
p1a <- ggplot() +
  annotate("rect", xmin = 1, xmax = 15, ymin = sebm$lo, ymax = sebm$hi, alpha = .15, fill = "#2C7FB8") +
  annotate("segment", x = 1, xend = 15, y = sebm$med, yend = sebm$med, colour = "#2C7FB8", linewidth = 1) +
  annotate("text", x = 12.5, y = sebm$med + .03, label = "Bone marrow Se", colour = "#2C7FB8", size = 3.5) +
  annotate("rect", xmin = 1, xmax = 15, ymin = spbc$lo, ymax = spbc$hi, alpha = .12, fill = "grey30") +
  annotate("text", x = 12.5, y = spbc$med - .035, label = "Blood-culture Sp", colour = "grey30", size = 3.5) +
  geom_point(data = ant, aes(blood_vol_mL, Se_BC), colour = "grey55", size = 1.8, alpha = .8) +
  geom_line(data = curve, aes(vol, Se), colour = "#D95F02", linewidth = 1.3) +
  geom_errorbar(data = anch, aes(vol, ymin = lo, ymax = hi), width = .3, colour = "#D95F02") +
  geom_point(data = anch, aes(vol, Se), colour = "#D95F02", size = 3) +
  annotate("text", x = 10.4, y = anch$Se[3] - .05, label = "Blood-culture Se", colour = "#D95F02", size = 3.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(title = "A. Diagnostic accuracy of blood vs bone-marrow culture",
       subtitle = "Latent-class estimates; grey = Antillon 2018 study-level Se points",
       x = "Blood volume cultured (mL)", y = "Sensitivity / Specificity") + th
sav(p1a, "fig1a_culture_accuracy.png", 8, 5)

# ---- FIG 1B: PPV by setting -- hospital phi + community pi forest ------------
phi <- read.csv("latent_class_ppv/tables/final_phi_hospital.csv")
pic <- read.csv("latent_class_ppv/tables/final_pi_community.csv")
fh <- data.frame(setting = "Hospital (phi)", study = paste0(phi$study, " (", phi$vol_mL, " mL)"),
                 med = phi$phi_med, lo = phi$phi_lo, hi = phi$phi_hi)
fc <- data.frame(setting = "Community/outbreak (pi)", study = pic$study,
                 med = pic$pi_med, lo = pic$pi_lo, hi = pic$pi_hi)
fdf <- rbind(fh, fc); fdf <- fdf[order(fdf$setting, fdf$med), ]; fdf$study <- factor(fdf$study, levels = fdf$study)
p1b <- ggplot(fdf, aes(med, study, colour = setting)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .3) + geom_point(size = 2.6) +
  scale_colour_manual(values = c("Hospital (phi)" = "#1B9E77", "Community/outbreak (pi)" = "#7570B3")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(title = "B. PPV of the clinical case definition by setting",
       subtitle = "Hospital-referred (phi) vs community/outbreak surveillance (pi)",
       x = "PPV  =  P(true typhoid | suspected)", y = NULL, colour = NULL) +
  th + theme(legend.position = "top")
sav(p1b, "fig1b_ppv_by_setting.png", 8, 6)

# ---- FIG 1C: PPV spectrum across case definitions x prevalence --------------
spec <- read.csv("latent_class_ppv/tables/ppv_spectrum_table.csv", check.names = FALSE)
sl <- spec %>% pivot_longer(c("endemic 0.03", "moderate 0.10", "outbreak 0.30", "high 0.50"),
                            names_to = "prev", values_to = "PPV")
sl$prev <- factor(sl$prev, levels = c("endemic 0.03", "moderate 0.10", "outbreak 0.30", "high 0.50"))
sl$variant <- factor(sl$variant, levels = spec$variant)
p1c <- ggplot(sl, aes(prev, PPV, group = variant, colour = variant)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(title = "C. PPV surface across case definitions and prevalence",
       subtitle = "Syndromic definitions cap specificity; serologic tiers lift PPV",
       x = "True-typhoid prevalence among febrile (setting)", y = "PPV", colour = "Case definition") +
  th + theme(legend.position = "right", legend.text = element_text(size = 8))
sav(p1c, "fig1c_ppv_spectrum.png", 9, 5)

# ============================================================================
# ORI figures -- reuse the renewal pipeline data
# ============================================================================
source("R/2-functions.R"); source("renewal/R/gi.R"); source("renewal/R/renewal_core.R")
source("renewal/R/data_prep.R"); source("renewal/R/scenario.R"); source("renewal/R/ppv.R")
cfg <- yaml::read_yaml("renewal/config.yml"); set.seed(cfg$seed)
prep <- prep_outbreaks(cfg)
pim_med <- setNames(pic$pi_med, pic$study)

# ---- FIG 2A: observed suspected S(t) vs model true typhoid I(t) = S - B ------
reps <- intersect(c("Kabwama 2017", "Neil 2012", "Aye 2004"), names(prep$series))
dec <- do.call(rbind, lapply(reps, function(sid) {
  S <- prep$series[[sid]]
  pv <- if (sid %in% names(pim_med)) pim_med[[sid]] else plogis(median(load_pi_posterior(cfg$ppv$draws, cfg$ppv$anchor)$mu))
  I <- ppv_true_incidence(S, pv)
  rbind(data.frame(sid = sprintf("%s (pi=%.2f)", sid, pv), week = seq_along(S), val = S, series = "Suspected (observed)"),
        data.frame(sid = sprintf("%s (pi=%.2f)", sid, pv), week = seq_along(S), val = I, series = "True typhoid I(t) (model)"))
}))
p2a <- ggplot(dec, aes(week, val, fill = series)) +
  geom_area(data = subset(dec, series == "Suspected (observed)"), alpha = .35, fill = "#9E9E9E") +
  geom_area(data = subset(dec, series == "True typhoid I(t) (model)"), alpha = .75, fill = "#D7301F") +
  facet_wrap(~sid, scales = "free", ncol = 3) +
  labs(title = "A. Observed suspected cases vs model-predicted true typhoid (S = I + B)",
       subtitle = "Grey = suspected surveillance; red = de-backgrounded true typhoid I(t); gap = non-typhoid background B(t)",
       x = "Epidemiological week", y = "Cases", fill = NULL) + th
sav(p2a, "fig2a_suspected_vs_true.png", 11, 4)

# ---- FIG 2B: ORI impact across coverage x timing ----------------------------
g2 <- read.csv("renewal/tables/summary_delay_coverage.csv")
g2 <- g2[g2$model == "renewal" & g2$ori_occurred != FALSE, ]
pool <- g2 %>% group_by(tau, vacc_cov) %>%
  summarise(pct = median(pct_reduction_median, na.rm = TRUE),
            averted = sum(s_ch_averted_median, na.rm = TRUE), .groups = "drop")
pool$vacc_cov <- factor(pool$vacc_cov)
p2b1 <- ggplot(pool, aes(tau, pct, colour = vacc_cov, group = vacc_cov)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "B. ORI impact: true-typhoid % reduction by campaign timing x coverage",
       subtitle = "Pooled median across developing-country outbreaks (additive PPV regime)",
       x = "Campaign delay tau (weeks to protection)", y = "True typhoid cases averted (%)", colour = "Coverage") +
  th + theme(legend.position = "top")
sav(p2b1, "fig2b_ori_impact_grid.png", 8, 5)
p2b2 <- ggplot(pool, aes(tau, averted, colour = vacc_cov, group = vacc_cov)) +
  geom_line(linewidth = 1) + geom_point(size = 2) + scale_colour_brewer(palette = "Set1") +
  labs(title = "C. True typhoid cases averted by timing x coverage",
       subtitle = "Summed across outbreaks", x = "Campaign delay tau (weeks)",
       y = "True typhoid cases averted", colour = "Coverage") + th + theme(legend.position = "top")
sav(p2b2, "fig2c_cases_averted_grid.png", 8, 5)

# ---- Tables ----------------------------------------------------------------
t_ppv <- par; write.csv(t_ppv, file.path(tabdir, "tab_ppv_parameters.csv"), row.names = FALSE)
t_ori <- pool %>% arrange(vacc_cov, tau)
write.csv(t_ori, file.path(tabdir, "tab_ori_impact_grid.csv"), row.names = FALSE)
cat("Wrote figures to", figdir, "and tables to", tabdir, "\n")
cat("PPV: Se_BC 0.52/0.62/0.68 @2/5/10mL; Se_BM 0.90; Sp 0.99\n")
cat("ORI base (tau=8,cov=0.8): pooled true %reduction",
    round(pool$pct[pool$tau == 8 & pool$vacc_cov == 0.8], 1), "%\n")
