# =============================================================================
# Per-outbreak allocation brackets across the renewal set.
#   Rscript revision/run_outbreak_allocation_brackets.R
# Each outbreak is bracketed over ONLY the allocations that are epidemiologically
# plausible for it (from transmission mode, case definition, and documented
# reporting artefacts), not the uniform five. Background-fixed counterfactual
# throughout, so suspected_reduction = pi * true_reduction exactly.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
source("revision/allocation_scenarios.R")   # ALLOCS, vaccinate(), renewal helpers
outdir <- "revision/main_text_figures"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
dir.create("revision/tables", showWarnings=FALSE, recursive=TRUE)

cfg  <- yaml::read_yaml("renewal/config.yml")
prep <- prep_outbreaks(cfg); w <- gi_from_config(cfg)
coverage <- cfg$vaccine$coverage_base
psi_T    <- 0.8*cfg$vaccine$psi_mean
delta    <- round((cfg$timing$campaign_duration_days+cfg$timing$immunity_onset_days)/cfg$step_days)
REL_TEFF <- 0.35
# TIMING: "production" = config base ORI delay (tau) + campaign/immunity delay, absolute weeks
#         "relative"   = protection at 35% of each outbreak's duration (illustrative)
TIMING <- Sys.getenv("TIMING", "production")
teff_of <- function(Tn) if (identical(TIMING,"production"))
  cfg$timing$delay_base_weeks + delta else max(2L, round(REL_TEFF*Tn))
tag <- if (identical(TIMING,"production")) "production" else "rel35"
cat(sprintf("TIMING = %s (protection week %s)\n", TIMING,
            if (identical(TIMING,"production")) as.character(cfg$timing$delay_base_weeks+delta) else "35% of duration"))

pic  <- read.csv("latent_class_ppv/tables/final_pi_community.csv", check.names=FALSE, stringsAsFactors=FALSE)
pars <- read.csv("latent_class_ppv/tables/final_parameters.csv", check.names=FALSE)
PI_POP <- pars$med[pars$param=="pi_population_median"]      # fitted population PPV for unanchored

# Plausible allocations per outbreak (first = primary); see transmission-mode review.
PLAUSIBLE <- list(
  "Lewis 2005"                  = c("Additive","Multiplicative"),
  "Davis 2018"                  = c("Reporting-shock","Additive"),
  "N'Cho 2019"                  = c("Additive","Time-varying PPV","Reporting-shock"),
  "Muyembe-Tamfum 2009"         = c("Multiplicative","Additive"),
  "Neil 2012"                   = c("Renewal-consistent","Multiplicative","Additive"),
  "Qamar 2018"                  = c("Renewal-consistent","Multiplicative"),
  "Yousafzai 2019"              = c("Renewal-consistent","Multiplicative"),
  "Muti 2014"                   = c("Renewal-consistent","Multiplicative","Time-varying PPV"),
  "Polonsky 2014, Dzivaresekwa" = c("Renewal-consistent","Multiplicative","Time-varying PPV"),
  "Polonsky 2014_Kuwadzana"     = c("Renewal-consistent","Multiplicative","Time-varying PPV"),
  "Ali 2017"                    = c("Time-varying PPV","Additive","Renewal-consistent"),
  "Aye 2004"                    = c("Additive","Multiplicative"),
  "Kabwama 2017"                = c("Additive","Multiplicative","Renewal-consistent"))

rows <- list()
for (ob in names(prep$series)) {
  S <- prep$series[[ob]]; Tn <- length(S)
  pi_o <- pic$pi_med[match(ob, pic$study)]
  pi_src <- if (is.na(pi_o)) "population" else "anchored"
  if (is.na(pi_o)) pi_o <- PI_POP
  t_eff <- teff_of(Tn)
  plaus <- PLAUSIBLE[[ob]]; if (is.null(plaus)) plaus <- names(ALLOCS)
  for (nm in names(ALLOCS)) {
    I <- if (nm=="Time-varying PPV") ALLOCS[[nm]](S, pi_o, 2) else ALLOCS[[nm]](S, pi_o)
    Bkg <- pmax(S - I, 0)
    vc <- vaccinate(I, Bkg, w, t_eff, coverage, psi_T, delta, cfg$timing$c_shape)
    rows[[length(rows)+1L]] <- data.frame(
      outbreak=ob, pi=pi_o, pi_src=pi_src, Tn=Tn, t_eff=t_eff, missed=(t_eff>=Tn), alloc=nm,
      plausible = nm %in% plaus, primary = identical(nm, plaus[1]),
      susp_total=sum(S), true_total=sum(I),
      true_averted=sum(I)-sum(vc$Iv),
      true_red=100*(sum(I)-sum(vc$Iv))/sum(I),
      susp_red=100*(sum(S)-sum(vc$Sv))/sum(S), stringsAsFactors=FALSE)
  }
}
res <- do.call(rbind, rows)
write.csv(res, sprintf("revision/tables/outbreak_allocation_brackets_%s.csv", tag), row.names=FALSE)

miss <- unique(res$outbreak[res$missed])
if (length(miss)) cat(sprintf("Campaign lands at/after outbreak end (no impact) for %d outbreak(s): %s\n",
                             length(miss), paste(miss, collapse=", ")))

pl <- res[res$plausible, ]
per <- do.call(rbind, lapply(split(pl, pl$outbreak), function(d) data.frame(
  outbreak=d$outbreak[1], pi=d$pi[1], pi_src=d$pi_src[1], missed=d$missed[1],
  susp_total=d$susp_total[1], true_total=d$true_total[1],
  primary_alloc=d$alloc[d$primary][1], primary_true_red=d$true_red[d$primary][1],
  true_red_lo=min(d$true_red), true_red_hi=max(d$true_red),
  averted_lo=min(d$true_averted), averted_hi=max(d$true_averted), stringsAsFactors=FALSE)))
per <- per[order(-per$primary_true_red), ]

cat("=== per-outbreak plausible-allocation brackets (protection at 35% of duration) ===\n")
print(data.frame(outbreak=per$outbreak, PPV=round(per$pi,2), src=per$pi_src,
                 primary=per$primary_alloc,
                 true_red=sprintf("%.0f%% [%.0f-%.0f]", per$primary_true_red, per$true_red_lo, per$true_red_hi),
                 averted=sprintf("%.0f [%.0f-%.0f]", (per$averted_lo+per$averted_hi)/2, per$averted_lo, per$averted_hi)),
      row.names=FALSE)

agg_susp <- sum(per$susp_total); agg_true <- sum(per$true_total)
agg_lo <- sum(per$averted_lo);   agg_hi <- sum(per$averted_hi)
cat(sprintf("\n=== AGGREGATE across %d renewal outbreaks ===\n", nrow(per)))
cat(sprintf("suspected cases %s | true typhoid %s (mean PPV %.2f)\n",
            format(round(agg_susp), big.mark=","), format(round(agg_true), big.mark=","), agg_true/agg_susp))
cat(sprintf("true typhoid averted envelope: %s - %s\n", format(round(agg_lo), big.mark=","), format(round(agg_hi), big.mark=",")))
cat(sprintf("aggregate TRUE reduction: %.1f%% - %.1f%%\n", 100*agg_lo/agg_true, 100*agg_hi/agg_true))
cat(sprintf("aggregate SUSPECTED (surveillance-visible) reduction: %.1f%% - %.1f%%\n",
            100*agg_lo/agg_susp, 100*agg_hi/agg_susp))

# ---- figure: per-outbreak bracket ------------------------------------------
cols <- c(Multiplicative="#0072B2", Additive="#D55E00", `Time-varying PPV`="#009E73",
          `Renewal-consistent`="#CC79A7", `Reporting-shock`="#E69F00")
pl$outbreak <- factor(pl$outbreak, levels=rev(per$outbreak))
seg <- per; seg$outbreak <- factor(seg$outbreak, levels=rev(per$outbreak))
g <- ggplot() +
  geom_segment(data=seg, aes(y=outbreak, yend=outbreak, x=true_red_lo, xend=true_red_hi),
               colour="grey78", linewidth=1.2) +
  geom_point(data=pl, aes(y=outbreak, x=true_red, colour=alloc), size=3) +
  geom_point(data=pl[pl$primary,], aes(y=outbreak, x=true_red), shape=21, fill=NA, colour="grey20", size=5) +
  scale_colour_manual(values=cols) +
  scale_x_continuous(labels=function(x)paste0(x,"%")) +
  labs(x=sprintf("True typhoid averted (%%) - %s", if (identical(TIMING,"production"))
         sprintf("production ORI timing (protection at week %d)", cfg$timing$delay_base_weeks+delta)
         else "protection at 35%% of outbreak duration"), y=NULL,
       title="Per-outbreak impact bracket over epidemiologically plausible allocations",
       subtitle=sprintf("Ringed point = primary allocation for that outbreak. Aggregate: %s-%s true cases averted (%.0f-%.0f%% of true typhoid; only %.0f-%.0f%% visible in suspected surveillance).",
                        format(round(agg_lo), big.mark=","), format(round(agg_hi), big.mark=","),
                        100*agg_lo/agg_true, 100*agg_hi/agg_true, 100*agg_lo/agg_susp, 100*agg_hi/agg_susp)) +
  theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        legend.position="top", legend.title=element_blank(),
        plot.title=element_text(face="bold", size=12.5), plot.subtitle=element_text(colour="grey40", size=8.8))
ggsave(file.path(outdir, sprintf("fig_outbreak_allocation_brackets_%s.png", tag)), g, width=9.5, height=6, dpi=300)
cat(sprintf("\nWrote fig_outbreak_allocation_brackets_%s.png + outbreak_allocation_brackets_%s.csv\n", tag, tag))