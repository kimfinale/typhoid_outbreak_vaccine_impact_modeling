# =============================================================================
# Multi-outbreak vaccine-impact figures (three design approaches), single
# vaccination scenario. Uses the four outbreaks with both a renewal time series
# and a fitted local community PPV: Aye 2004, Neil 2012, Kabwama 2017, N'Cho 2019.
#   Rscript revision/make_multi_outbreak_figures.R
# Two additive levels: suspected = true typhoid + other febrile illness, and
# true typhoid = source/background + propagated renewal incidence.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({ library(ggplot2); library(patchwork) })
for (f in c("data_prep.R","gi.R","renewal_core.R","ppv.R","epiestim_rt.R")) source(file.path("renewal/R", f))
cfg <- yaml::read_yaml("renewal/config.yml")
outdir <- "revision/main_text_figures"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

coverage <- cfg$vaccine$coverage_base            # 0.80
psi_T    <- 0.8 * cfg$vaccine$psi_mean            # ~0.66 transmission VE
REL_TEFF <- 0.35                                  # single scenario: protection onset at 35% of outbreak duration
OUTBREAKS <- c("Aye 2004","Neil 2012","Kabwama 2017","N'Cho 2019")
ob_cols <- c("Aye 2004"="#E69F00","Neil 2012"="#0072B2","Kabwama 2017"="#D55E00","N'Cho 2019"="#009E73")
BLUE<-"#0072B2"; GREY<-"#BDBDBD"
th <- theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(colour="grey93"), legend.position="top", legend.title=element_blank(),
        plot.title=element_text(face="bold", size=12.5), plot.subtitle=element_text(colour="grey40", size=9.5),
        strip.text=element_text(face="bold", size=10), axis.title=element_text(colour="grey25"))

prep <- prep_outbreaks(cfg); w <- gi_from_config(cfg)
pic <- read.csv("latent_class_ppv/tables/final_pi_community.csv", check.names=FALSE, stringsAsFactors=FALSE)
mode <- read.csv("source_decomposition/tables/tab_transmission_mode.csv",
                 check.names=FALSE, stringsAsFactors=FALSE)

per <- lapply(OUTBREAKS, function(ob) {
  S <- prep$series[[ob]]; Tn <- length(S); wk <- seq_len(Tn)
  pi_o <- pic$pi_med[match(ob, pic$study)]
  theta_o <- mode$final_theta[match(ob, mode$study_id)]
  if (!is.finite(theta_o)) stop("No source fraction for ", ob)
  obs <- ppv_incidence_components(S, pi_o)
  I <- obs$true_typhoid; Bkg <- obs$other_febrile
  t_eff <- max(2L, round(REL_TEFF*Tn))
  tau <- max(1,t_eff-round((cfg$timing$campaign_duration_days+cfg$timing$immunity_onset_days)/cfg$step_days))
  add <- additive_source_counterfactual(
    I, w, tau=tau, t_eff=t_eff, coverage=coverage, psi_T=psi_T,
    source_fraction=theta_o, c_shape=cfg$timing$c_shape)
  Iv <- add$incidence_v
  Sv <- Iv + Bkg
  list(ob=ob, Tn=Tn, wk=wk, pi=pi_o, theta=theta_o, S=S, I=I, Bkg=Bkg,
       Iv=Iv, Sv=Sv, t_eff=t_eff, Rt_prop=add$Rt_prop,
       true_red=100*(sum(I)-sum(Iv))/sum(I), susp_red=100*(sum(S)-sum(Sv))/sum(S))
})
names(per) <- OUTBREAKS
lab <- setNames(sprintf("%s  (PPV %.2f; theta %.2f)", OUTBREAKS,
                        sapply(per,function(p)p$pi), sapply(per,function(p)p$theta)), OUTBREAKS)
for (p in per) cat(sprintf("%-14s PPV %.2f theta %.2f protect wk %d/%d true -%.0f%% suspected -%.0f%%\n",
                           p$ob, p$pi, p$theta, p$t_eff, p$Tn, p$true_red, p$susp_red))

# =============================================================================
# APPROACH A - Small multiples: per-outbreak case decomposition + vaccine
# =============================================================================
stk <- do.call(rbind, lapply(per, function(p) rbind(
  data.frame(ob=unname(lab[p$ob]), week=p$wk, val=p$Bkg, part="Non-typhoid febrile background"),
  data.frame(ob=unname(lab[p$ob]), week=p$wk, val=p$I,   part="True typhoid"))))
stk$part <- factor(stk$part, levels=c("Non-typhoid febrile background","True typhoid"))
vline <- do.call(rbind, lapply(per, function(p) data.frame(ob=unname(lab[p$ob]), week=p$wk, Iv=p$Iv)))
ann <- do.call(rbind, lapply(per, function(p) data.frame(ob=unname(lab[p$ob]),
              txt=sprintf("true -%.0f%%\nsuspected -%.0f%%", p$true_red, p$susp_red))))
gA <- ggplot() +
  geom_area(data=stk, aes(week, val, fill=part), alpha=0.9) +
  geom_line(data=vline, aes(week, Iv), colour="#111111", linewidth=0.8, linetype="42") +
  geom_text(data=ann, aes(x=Inf, y=Inf, label=txt), hjust=1.05, vjust=1.15, size=3, colour="grey25") +
  facet_wrap(~ob, scales="free", ncol=2) +
  scale_fill_manual(values=c("True typhoid"=BLUE,"Non-typhoid febrile background"=GREY)) +
  labs(x="Week from outbreak onset", y="Cases / week",
       title="Approach A - Small multiples: true typhoid vs background, per outbreak",
       subtitle="Stacked = reported suspected; dashed = true typhoid after a single vaccination campaign (protection at ~35% of the outbreak)") + th
ggsave(file.path(outdir,"fig_multi_A_smallmultiples.png"), gA, width=9, height=7, dpi=300)

# =============================================================================
# APPROACH B - Normalized overlay: all outbreaks on common progress axis
# =============================================================================
ovR <- do.call(rbind, lapply(per, function(p)
  data.frame(ob=p$ob, prog=p$wk/p$Tn, Rt=p$Rt_prop)))
ovR <- ovR[is.finite(ovR$Rt), ]
ovC <- do.call(rbind, lapply(per, function(p) rbind(
  data.frame(ob=p$ob, prog=p$wk/p$Tn, val=p$S/max(p$S), series="Observed"),
  data.frame(ob=p$ob, prog=p$wk/p$Tn, val=p$Sv/max(p$S), series="Vaccinated"))))
gB1 <- ggplot(ovR, aes(prog, Rt, colour=ob)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey55") +
  geom_line(linewidth=0.9) + scale_colour_manual(values=ob_cols) +
  coord_cartesian(ylim=c(0, 4)) +
  labs(x=NULL, y=expression(R[t]), title="Approach B - Overlaid on a common outbreak-progress axis",
       subtitle="Propagated-channel coefficient conditional on each paper-informed source fraction") + th
gB2 <- ggplot(ovC, aes(prog, val, colour=ob, linetype=series)) +
  geom_line(linewidth=0.9) + scale_colour_manual(values=ob_cols) +
  scale_linetype_manual(values=c(Observed="solid", Vaccinated="42")) +
  labs(x="Outbreak progress (week / duration)", y="Suspected cases (/ peak)",
       subtitle="Reported suspected cases, normalized to each outbreak's peak: observed (solid) vs vaccinated (dashed)") + th
ggsave(file.path(outdir,"fig_multi_B_overlay.png"), gB1/gB2, width=8.5, height=8, dpi=300)

# =============================================================================
# APPROACH C - Dumbbell summary: true vs surveillance-visible reduction
# =============================================================================
dd <- do.call(rbind, lapply(per, function(p) data.frame(ob=p$ob, pi=p$pi, theta=p$theta,
        true=p$true_red, susp=p$susp_red)))
dd$lab <- sprintf("%s (PPV %.2f; theta %.2f)", dd$ob, dd$pi, dd$theta)
dd$lab <- factor(dd$lab, levels=dd$lab[order(dd$true)])
gC <- ggplot(dd) +
  geom_segment(aes(y=lab, yend=lab, x=susp, xend=true), colour="grey70", linewidth=1.1) +
  geom_point(aes(y=lab, x=true), colour=BLUE, size=4) +
  geom_point(aes(y=lab, x=susp), colour="#8A8A8A", size=4) +
  geom_text(aes(y=lab, x=true, label=sprintf("%.0f%%", true)), colour=BLUE, vjust=-1.1, size=3.2) +
  geom_text(aes(y=lab, x=susp, label=sprintf("%.0f%%", susp)), colour="#6A6A6A", vjust=-1.1, size=3.2) +
  scale_x_continuous(limits=c(0, max(dd$true)*1.1), labels=function(x)paste0(x,"%")) +
  labs(x="Reduction after one vaccination campaign", y=NULL,
       title="Approach C - Summary: true typhoid vs surveillance-visible reduction",
       subtitle="Blue = true typhoid averted; grey = reduction visible in suspected-case surveillance under the two-level additive model.") +
  th + theme(legend.position="none", panel.grid.major.y=element_blank())
ggsave(file.path(outdir,"fig_multi_C_dumbbell.png"), gC, width=8, height=4.2, dpi=300)

cat("\nWrote 3 figures to", outdir, "\n")
