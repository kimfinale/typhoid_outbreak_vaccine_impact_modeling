# =============================================================================
# Illustrative multipanel outbreak figure: R_t, reported suspected cases, and
# the true-typhoid / background-febrile decomposition, with the vaccine
# counterfactual (R_t and case reductions) at representative campaign times.
#   Rscript revision/make_outbreak_mechanism_figure.R
# Two additive levels: suspected = true typhoid + other febrile illness, then
# true typhoid = source/background + propagated renewal incidence.
# Self-contained; reads committed renewal code + the fitted community PPV.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({ library(ggplot2) })
for (f in c("data_prep.R","gi.R","renewal_core.R","ppv.R","epiestim_rt.R"))
  source(file.path("renewal/R", f))
cfg <- yaml::read_yaml("renewal/config.yml")
outdir <- "revision/main_text_figures"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- choices ----------------------------------------------------------------
OUTBREAK  <- Sys.getenv("OUTBREAK", "Kabwama 2017")   # anchored + 24 weeks
coverage  <- cfg$vaccine$coverage_base                 # 0.80
psi_T     <- 0.8 * cfg$vaccine$psi_mean                 # ~0.66: midpoint of transmission-VE range [0.6,1]*psi
delta_wk  <- (cfg$timing$campaign_duration_days + cfg$timing$immunity_onset_days)/cfg$step_days  # ~5 wk
campaigns <- data.frame(label=c("Early response","Later response"),
                        tau=c(2, 8), col=c("#D55E00","#009E73"), stringsAsFactors=FALSE)  # weeks from onset

BLUE<-"#0072B2"; GREY<-"#BDBDBD"; INK<-"#333333"
theme_ts <- theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(colour="grey93"),
        plot.title=element_text(face="bold", size=11.5), plot.subtitle=element_text(colour="grey40", size=9.5),
        legend.position="top", legend.title=element_blank(), axis.title=element_text(colour="grey25"))

# ---- data + GI + PPV --------------------------------------------------------
prep <- prep_outbreaks(cfg)
S <- prep$series[[OUTBREAK]]
if (is.null(S)) stop("Outbreak not in renewal set: ", OUTBREAK)
Tn <- length(S); wk <- seq_len(Tn)
w  <- gi_from_config(cfg)
pic <- read.csv("latent_class_ppv/tables/final_pi_community.csv", check.names=FALSE, stringsAsFactors=FALSE)
pi_o <- pic$pi_med[match(OUTBREAK, pic$study)]
if (is.na(pi_o)) stop("No fitted PPV for ", OUTBREAK)

# Primary additive observation allocation.
obs_components <- ppv_incidence_components(S, pi_o)
I <- obs_components$true_typhoid
Bkg <- obs_components$other_febrile
mode <- read.csv("source_decomposition/tables/tab_transmission_mode.csv",
                 check.names=FALSE, stringsAsFactors=FALSE)
theta_o <- mode$final_theta[match(OUTBREAK, mode$study_id)]
if (!is.finite(theta_o)) stop("No source fraction for ", OUTBREAK)

# Baseline propagated-channel coefficient for display.
baseline <- additive_source_counterfactual(
  I, w, tau=1, t_eff=Tn + 1, coverage=0, psi_T=psi_T,
  source_fraction=theta_o, c_shape="step")
Rt_prop <- baseline$Rt_prop

# ---- additive vaccine counterfactual per campaign ---------------------------
propagate_true <- function(I, w, tau, t_eff) {
  add <- additive_source_counterfactual(
    I, w, tau=tau, t_eff=t_eff, coverage=coverage, psi_T=psi_T,
    source_fraction=theta_o, c_shape=cfg$timing$c_shape)
  list(Iv=add$incidence_v, Sv=add$incidence_v + Bkg,
       c_t=add$c_t, q_t=add$q_t)
}
scen <- lapply(seq_len(nrow(campaigns)), function(i) {
  tau <- campaigns$tau[i]; t_eff <- tau + delta_wk
  cf <- propagate_true(I, w, tau, t_eff)
  true_red <- 100*(sum(I)-sum(cf$Iv))/sum(I)
  susp_red <- 100*(sum(S)-sum(cf$Sv))/sum(S)
  list(i=i, label=campaigns$label[i], col=campaigns$col[i], tau=tau, t_eff=t_eff,
       Iv=cf$Iv, Sv=cf$Sv, true_red=true_red, susp_red=susp_red)
})
cat(sprintf("%s: pi=%.2f, theta=%.2f, true-typhoid=%d/%d suspected\n",
            OUTBREAK, pi_o, theta_o, round(sum(I)), round(sum(S))))
for (s in scen) cat(sprintf("  %-15s campaign wk %d (protected wk %.0f): true reduction %.0f%%, suspected reduction %.0f%%\n",
                            s$label, s$tau, s$t_eff, s$true_red, s$susp_red))

# =============================================================================
# PANEL A — reproduction number
# =============================================================================
rtdf <- data.frame(week=wk, mean_r=Rt_prop)
# Vaccinated propagated coefficient, conditional on the source allocation.
vax_rt <- do.call(rbind, lapply(scen, function(s) {
  keep <- wk >= s$t_eff & is.finite(Rt_prop)
  data.frame(week=wk[keep], r=Rt_prop[keep]*(1 - coverage*psi_T), label=s$label)
}))
gA <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", colour="grey55") +
  geom_vline(xintercept=sapply(scen,function(s)s$t_eff),
             linetype="dotted", colour=campaigns$col) +
  geom_line(data=rtdf, aes(week, mean_r), colour=INK, linewidth=1) +
  geom_line(data=vax_rt, aes(week, r, colour=label), linewidth=1) +
  scale_colour_manual(values=setNames(campaigns$col, campaigns$label)) +
  coord_cartesian(ylim=c(0, min(4, max(rtdf$mean_r, na.rm=TRUE)*1.05))) +
  labs(x=NULL, y=expression(R[t]^P), title="A  Propagated-channel coefficient",
       subtitle="Conditional on the paper-informed source allocation; coloured = vaccinated from protection onset; R=1 dashed") +
  theme_ts

# =============================================================================
# PANEL B — reported suspected cases + vaccinated suspected
# =============================================================================
obs <- data.frame(week=wk, S=S)
vax_s <- do.call(rbind, lapply(scen, function(s) data.frame(week=wk, S=s$Sv, label=s$label)))
subB <- paste(sapply(scen, function(s) sprintf("%s: -%.0f%%", s$label, s$susp_red)), collapse="   ")
gB <- ggplot() +
  geom_col(data=obs, aes(week, S), fill=GREY, width=0.85) +
  geom_line(data=vax_s, aes(week, S, colour=label), linewidth=0.9) +
  scale_colour_manual(values=setNames(campaigns$col, campaigns$label)) +
  labs(x=NULL, y="Suspected cases / week",
       title="B  Reported suspected cases",
       subtitle=paste0("Observed (bars) vs vaccinated. Surveillance-visible reduction is diluted:  ", subB)) +
  theme_ts

# =============================================================================
# PANEL C — true typhoid vs background febrile, with true-typhoid reduction
# =============================================================================
stackdf <- rbind(
  data.frame(week=wk, val=Bkg, part="Non-typhoid febrile background"),
  data.frame(week=wk, val=I,   part="True typhoid"))
stackdf$part <- factor(stackdf$part, levels=c("Non-typhoid febrile background","True typhoid"))
vax_true <- do.call(rbind, lapply(scen, function(s) data.frame(week=wk, Iv=s$Iv, label=s$label)))
gC <- ggplot() +
  geom_area(data=stackdf, aes(week, val, fill=part), alpha=0.9) +
  geom_line(data=vax_true, aes(week, Iv, colour=label), linewidth=1, linetype="42") +
  scale_fill_manual(values=c("True typhoid"=BLUE, "Non-typhoid febrile background"=GREY)) +
  scale_colour_manual(values=setNames(campaigns$col, campaigns$label)) +
  labs(x="Week from outbreak onset", y="Cases / week",
       title="C  Inferred true typhoid vs background febrile",
       subtitle=sprintf("Stacked = suspected (PPV %.2f; source fraction %.2f); dashed = true typhoid under vaccination (true reduction: early -%.0f%%, later -%.0f%%)",
                        pi_o, theta_o, scen[[1]]$true_red, scen[[2]]$true_red)) +
  theme_ts

ok <- requireNamespace("patchwork", quietly=TRUE)
fig <- if (ok) { library(patchwork); (gA / gB / gC) +
  patchwork::plot_annotation(
    title=sprintf("Vaccine impact on transmission, reported cases, and true typhoid - %s", OUTBREAK),
    subtitle="The vaccine cuts true typhoid transmission; non-typhoid febrile illness continues, so reported (suspected) cases fall far less than true typhoid",
    theme=theme(plot.title=element_text(face="bold", size=13),
                plot.subtitle=element_text(colour="grey40", size=10))) } else gA
ggsave(file.path(outdir, "fig_outbreak_mechanism.png"), fig, width=8.5, height=10, dpi=300)
cat("Wrote", file.path(outdir, "fig_outbreak_mechanism.png"), "\n")
