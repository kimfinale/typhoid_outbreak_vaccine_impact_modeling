# =============================================================================
# Kabwama (Kampala 2015) comprehensive analysis — 02: run models & scenarios.
#   Rscript revision/kabwama/02_analysis.R
#
# Kabwama is a MIXED outbreak, common-source dominant (theta=0.60; waterborne
# spring + vended beverages). The primary model adds common-source and propagated
# typhoid incidence inside one renewal recursion. The former theta-weighted
# static-plus-renewal outcome is retained only as a paired historical comparator.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
source("renewal/R/renewal_core.R")
source("revision/kabwama/01_model.R")
OUT <- "revision/kabwama"; TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
inp <- readRDS(file.path(OUT, "data/kabwama_inputs.rds"))
S <- inp$weekly; Tw <- length(S); wk <- seq_len(Tw)

## ---- fixed parameters -------------------------------------------------------
w      <- gi_weekly(inp$gi_mean, inp$gi_cv)
# ONE vaccine efficacy in BOTH arms: the static direct effect and the renewal
# first-step R_t reduction must be identical (they differ only in whether indirect
# effects compound). psi = 0.83. psiT (transmission VE in [0.6psi,psi]) is retained
# ONLY as a separate transmission-VE sensitivity, not a base-case static/renewal gap.
psi    <- inp$psi; VE <- psi; psiT_sens <- 0.8 * psi
delta  <- round((inp$campaign_days + inp$onset_days) / inp$step_days)
t_eff_base <- inp$delay_base_weeks + delta; kap_base <- inp$coverage
INVEST_WK <- 6; p_symp_base <- 0.40; phi_NTS <- 0.15
THETA_CS <- 0.60                                        # common-source fraction (source_decomposition)
par <- list(incub_d = 10, infper_d = 14, p_symp = p_symp_base, xi = 0.5,
            theta = 0.029, gc = 0.10, muC = 1/(3*365), N = 1.5e6, psi = psi, seed = 8)
pi_med <- median(inp$pi_draws); tot_susp <- sum(S)
q3 <- function(v) quantile(v, c(.5,.05,.95), na.rm=TRUE)

# Additive primary impact plus the paired historical weighted comparator.
impact_one <- function(I, R, t_eff, kappa) {
  vr <- vaccinate_renewal(S, I, w, R, t_eff, kappa, VE)   # same VE as static
  sa <- static_averted(I, t_eff, kappa, VE)
  ra <- vr$typ_averted
  wa <- THETA_CS * sa + (1 - THETA_CS) * ra
  add <- additive_source_counterfactual(
    I, w, tau = t_eff, t_eff = t_eff, coverage = kappa, psi_T = VE,
    source_fraction = THETA_CS, c_shape = "step")
  aa <- sum(I - add$incidence_v)
  tot <- sum(I)
  list(static = 100*sa/tot, renewal = 100*ra/tot,
       additive = 100*aa/tot, weighted = 100*wa/tot,
       additive_averted = aa, weighted_averted = wa,
       obs_additive = 100*aa/tot_susp,
       incidence_v = add$incidence_v)
}

## ===========================================================================
## A. Outbreak decomposition (median PPV)
## ===========================================================================
I_med <- debackground(S, pi_med); B_med <- S - I_med
draws_typ <- inp$pi_draws * tot_susp
decomp <- data.frame(
  quantity = c("suspected cases (observed)", "true typhoid symptomatic (PPV)",
               "non-typhoid background (all)", "  non-typhoid BACTERIAL (~NTS/Shigella)",
               "  non-bacterial (malaria/viral)", "total typhoid infections (incl asymptomatic)",
               "  asymptomatic typhoid", "chronic carriers seeded (theta x infections)"),
  value = c(round(tot_susp),
            sprintf("%.0f [%.0f, %.0f]", q3(draws_typ)[1], q3(draws_typ)[2], q3(draws_typ)[3]),
            round(tot_susp*(1-pi_med)), round(tot_susp*(1-pi_med)*phi_NTS),
            round(tot_susp*(1-pi_med)*(1-phi_NTS)),
            sprintf("%.0f [%.0f, %.0f]", q3(draws_typ)[1]/p_symp_base, q3(draws_typ)[2]/p_symp_base, q3(draws_typ)[3]/p_symp_base),
            round(pi_med*tot_susp*(1-p_symp_base)/p_symp_base),
            round(pi_med*tot_susp/p_symp_base*par$theta)))
write.csv(decomp, file.path(TAB, "tab_decomposition.csv"), row.names = FALSE)
cat("=== A. Outbreak decomposition (median PPV pi =", round(pi_med,3), ") ===\n"); print(decomp, row.names = FALSE)

## ===========================================================================
## B. R_t + primary additive impact (production timing)
## ===========================================================================
R_raw <- reconstruct_R(I_med, w)                        # self-consistent (for counterfactual)
R_med <- cori_R(I_med, w)                               # smoothed (for display only)
cat(sprintf("\n=== B. Typhoid R_t: raw peak %.1f, smoothed peak %.1f (wk %d) ===\n",
            max(R_raw,na.rm=TRUE), max(R_med,na.rm=TRUE), which.max(R_med)))
cat(sprintf("  NOTE: peak R_t >> typical typhoid (2-4) -> common-source exposures misread as propagation; renewal is an UPPER bracket.\n"))
ib <- impact_one(I_med, R_raw, t_eff_base, kap_base)
cat(sprintf("Production timing (protection wk %d, coverage %.2f):\n", t_eff_base, kap_base))
cat(sprintf("  STATIC (direct)   true reduction %.1f%%\n", ib$static))
cat(sprintf("  RENEWAL (upper)   true reduction %.1f%%\n", ib$renewal))
cat(sprintf("  ADDITIVE source+propagated (theta=%.2f, PRIMARY) true reduction %.1f%%  -> OBSERVED suspected %.1f%%\n",
            THETA_CS, ib$additive, ib$obs_additive))
cat(sprintf("  HISTORICAL theta-weighted comparator true reduction %.1f%% (difference %+.1f pp)\n",
            ib$weighted, ib$additive - ib$weighted))
cat(sprintf("  additive typhoid symptomatic averted %.0f; total infections averted %.0f\n",
            ib$additive_averted, ib$additive_averted/p_symp_base))

## ===========================================================================
## C1. timing x coverage grid (additive primary and weighted comparator)
## ===========================================================================
delays <- c(1,2,4,6,8,10,12); covs <- c(0.3,0.5,0.65,0.8,0.9)
set.seed(1); idx <- sample(length(inp$pi_draws), 300)
Ipre <- lapply(idx, function(j) { I <- debackground(S, inp$pi_draws[j]); list(I=I, R=reconstruct_R(I, w)) })
grid <- expand.grid(delay = delays, coverage = covs)
gridres <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  d <- grid$delay[i]; kp <- grid$coverage[i]; te <- d + delta
  m <- vapply(Ipre, function(z) { r <- impact_one(z$I, z$R, te, kp)
    c(r$static, r$renewal, r$additive, r$weighted) }, numeric(4))
  data.frame(delay=d, protection_wk=te, coverage=kp,
             static_med=median(m[1,]), renewal_med=median(m[2,]),
             additive_med=median(m[3,]), additive_lo=quantile(m[3,],.05),
             additive_hi=quantile(m[3,],.95), weighted_med=median(m[4,]),
             delta_pp_med=median(m[3,]-m[4,]), observed_med=median(m[3,])*pi_med)
}))
write.csv(gridres, file.path(TAB, "tab_impact_grid.csv"), row.names = FALSE)
cat("\n=== C1. Impact grid: true % reduction (static | renewal | ADDITIVE[90% interval] | weighted comparator) ===\n")
print(within(gridres, { additive <- sprintf("%.0f [%.0f-%.0f]", additive_med, additive_lo, additive_hi)
  static <- round(static_med); renewal <- round(renewal_med); observed <- round(observed_med,1)
  weighted <- round(weighted_med)
})[, c("delay","protection_wk","coverage","static","renewal","additive","weighted","observed")], row.names=FALSE)

## C2. time-varying reporting (surveillance ramp at investigation wk 6) --------
rep_scn <- do.call(rbind, lapply(c(1.0,0.7,0.5), function(r_pre) {
  Sadj <- S; Sadj[wk < INVEST_WK] <- S[wk < INVEST_WK] / r_pre
  Ia <- debackground(Sadj, pi_med); Ra <- reconstruct_R(Ia, w)
  add <- additive_source_counterfactual(
    Ia, w, tau=t_eff_base, t_eff=t_eff_base, coverage=kap_base,
    psi_T=VE, source_fraction=THETA_CS, c_shape="step")
  aa <- sum(Ia - add$incidence_v)
  data.frame(pre_ramp_reporting=r_pre, Rt_peak_smoothed=round(max(cori_R(Ia,w),na.rm=TRUE),1),
             additive_true_red=round(100*aa/sum(Ia),1))
}))
write.csv(rep_scn, file.path(TAB, "tab_reporting_scenario.csv"), row.names = FALSE)
cat("\n=== C2. Time-varying reporting (ramp at wk 6) ===\n"); print(rep_scn, row.names = FALSE)

## C3. asymptomatic fraction sweep --------------------------------------------
asy <- do.call(rbind, lapply(c(0.30,0.40,0.50,0.60), function(ps) data.frame(
  p_symptomatic=ps, total_infections=round(pi_med*tot_susp/ps),
  asymptomatic=round(pi_med*tot_susp*(1-ps)/ps),
  tot_infections_averted_wk13=round(ib$additive_averted/ps))))
write.csv(asy, file.path(TAB, "tab_asymptomatic_sweep.csv"), row.names = FALSE)
cat("\n=== C3. Asymptomatic fraction sweep ===\n"); print(asy, row.names = FALSE)

## C4. observation model: additive (primary) vs multiplicative ----------------
mult <- vaccinate_multiplicative(S, w, reconstruct_R(S, w), t_eff_base, kap_base, VE)
obsmod <- data.frame(model=c("additive (primary): background fixed",
                             "multiplicative: I=pi*S (no dilution)"),
  true_pct_reduction=c(round(ib$additive,1), round(100*mult$red,1)),
  observed_pct_reduction=c(round(ib$obs_additive,1), round(100*mult$red,1)))
write.csv(obsmod, file.path(TAB, "tab_obs_model.csv"), row.names = FALSE)
cat("\n=== C4. Observation model ===\n"); print(obsmod, row.names = FALSE)

## ===========================================================================
## D. Comprehensive burden overlay: asymptomatic + carriers (analytic)
## Asymptomatic transmission is already inside R_t (which is reconstructed from
## symptomatic cases generated by the whole transmission network). The overlay
## partitions the true infection burden: total = symptomatic / p_symp.
## ===========================================================================
p <- p_symp_base
T_week <- I_med / p                                      # total typhoid infections / week
A_week <- T_week - I_med                                 # asymptomatic / week
Tv_week <- ib$incidence_v / p                            # with-ORI total under additive model
tot_inf <- sum(T_week); asy_inf <- sum(A_week); carr <- par$theta * tot_inf
tot_averted <- ib$additive_averted / p; carr_averted <- par$theta * tot_averted
cat(sprintf("\n=== D. Comprehensive burden (p_symp %.2f) ===\n", p))
cat(sprintf("total typhoid infections %.0f | symptomatic %.0f | asymptomatic %.0f | chronic carriers formed %.0f\n",
            tot_inf, sum(I_med), asy_inf, carr))
cat(sprintf("production ORI: total infections averted %.0f | carriers averted %.0f (additive source model)\n",
            tot_averted, carr_averted))

saveRDS(list(I_med=I_med,B_med=B_med,R_raw=R_raw,R_med=R_med,ib=ib,gridres=gridres,
             T_week=T_week,A_week=A_week,Tv_week=Tv_week,decomp=decomp),
        file.path(OUT, "kabwama_results.rds"))

## ===========================================================================
## Figures
## ===========================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  ## Fig 1: decomposition
  dc <- rbind(data.frame(week=wk, layer="typhoid (symptomatic)", cases=I_med),
              data.frame(week=wk, layer="non-typhoid bacterial (~NTS)", cases=B_med*phi_NTS),
              data.frame(week=wk, layer="other febrile (malaria/viral)", cases=B_med*(1-phi_NTS)))
  dc$layer <- factor(dc$layer, levels=c("other febrile (malaria/viral)","non-typhoid bacterial (~NTS)","typhoid (symptomatic)"))
  g1 <- ggplot(dc, aes(week, cases, fill=layer)) + geom_area(alpha=.92) +
    geom_vline(xintercept=INVEST_WK, linetype="dashed", colour="grey30") +
    scale_fill_manual(values=c("other febrile (malaria/viral)"="#56B4E9","non-typhoid bacterial (~NTS)"="#E69F00","typhoid (symptomatic)"="#D55E00")) +
    labs(x="Epidemic week", y="Suspected cases / week", fill=NULL,
         title="Kabwama (Kampala 2015): decomposition of the suspected-case curve",
         subtitle=sprintf("PPV pi=%.2f -> ~%.0f%% of suspected are typhoid; the rest is a TCV-invariant non-typhoid background. Dashed = investigation/treatment-centre rollout.", pi_med, 100*pi_med)) +
    theme_minimal(base_size=11) + theme(legend.position="top", plot.title=element_text(face="bold",size=12))
  ggsave(file.path(FIG,"fig1_decomposition.png"), g1, width=9, height=4.6, dpi=300)
  ## Fig 2: R_t
  g2 <- ggplot(data.frame(week=wk, Rt=R_med), aes(week, Rt)) +
    geom_hline(yintercept=1, linetype="dotted", colour="grey50") +
    geom_line(colour="#0072B2", linewidth=1) + geom_point(colour="#0072B2", size=1.5) +
    geom_vline(xintercept=t_eff_base, linetype="dashed", colour="#D55E00") +
    labs(x="Epidemic week", y=expression(R[t]~"(typhoid)"),
         title="Reconstructed typhoid R_t (Cori) — inflated by common-source exposures",
         subtitle=sprintf("Peak R_t=%.1f is far above typical typhoid (2-4): the waterborne common source is misread as propagation. Dashed=production protection wk %d.", max(R_med,na.rm=TRUE), t_eff_base)) +
    theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold",size=12))
  ggsave(file.path(FIG,"fig2_Rt.png"), g2, width=8, height=4.2, dpi=300)
  ## Fig 3: impact grid (additive primary)
  g3 <- ggplot(gridres, aes(protection_wk, additive_med, colour=factor(coverage))) +
    geom_ribbon(aes(ymin=additive_lo, ymax=additive_hi, fill=factor(coverage)), alpha=.12, colour=NA) +
    geom_line(linewidth=1) + geom_point(size=1.5) + geom_vline(xintercept=t_eff_base, linetype="dashed", colour="grey40") +
    labs(x="Protection week", y="True typhoid reduction (%) — additive source model", colour="coverage", fill="coverage",
         title="Vaccine impact vs response speed (additive source + propagated transmission)",
         subtitle=sprintf("theta=%.2f allocates baseline incidence inside the recursion. Dashed=production timing (wk %d).", THETA_CS, t_eff_base)) +
    theme_minimal(base_size=11) + theme(legend.position="top", plot.title=element_text(face="bold",size=12))
  ggsave(file.path(FIG,"fig3_impact_grid.png"), g3, width=9, height=4.6, dpi=300)
  ## Fig 4: comprehensive burden (symptomatic + asymptomatic) + with-ORI
  bd <- rbind(data.frame(week=wk, layer="symptomatic (observed via PPV)", cases=I_med),
              data.frame(week=wk, layer="asymptomatic (unobserved)", cases=A_week))
  bd$layer <- factor(bd$layer, levels=c("asymptomatic (unobserved)","symptomatic (observed via PPV)"))
  g4 <- ggplot(bd, aes(week, cases, fill=layer)) + geom_area(alpha=.9) +
    geom_line(data=data.frame(week=wk, cases=Tv_week), aes(week, cases), inherit.aes=FALSE, linetype="dashed", linewidth=.9, colour="grey20") +
    scale_fill_manual(values=c("symptomatic (observed via PPV)"="#D55E00","asymptomatic (unobserved)"="#F0C86E")) +
    labs(x="Epidemic week", y="Typhoid infections / week", fill=NULL,
         title="Comprehensive typhoid burden: symptomatic + asymptomatic",
          subtitle=sprintf("Dashed = total infections under production ORI (additive source model). p_symp=%.0f%% -> ~%.0f%% asymptomatic.", 100*p_symp_base, 100*(1-p_symp_base))) +
    theme_minimal(base_size=11) + theme(legend.position="top", plot.title=element_text(face="bold",size=12))
  ggsave(file.path(FIG,"fig4_burden.png"), g4, width=9, height=4.6, dpi=300)
  cat("\nWrote figures to", FIG, "\n")
}
cat("\nDONE. Tables in", TAB, "\n")
