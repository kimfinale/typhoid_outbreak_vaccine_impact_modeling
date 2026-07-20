# =============================================================================
# Full dynamic SIRCV model vs the renewal-equation model, on the 2015 Kampala
# typhoid outbreak (our Kabwama 2017 anchor).
#
#   Rscript revision/seir_comparison/kampala_seir_vs_renewal.R
#
# The SIRCV model is a faithful re-implementation of Lee et al., PLOS NTD 2025
# (doi:10.1371/journal.pntd.0013566; repo modeling-computation/Typhoid_modeling):
#   dS = bN - lambda S - phi + omega_v V + omega R - mu S
#   dI = lambda S - delta I - mu I
#   dR = delta(1-theta-alpha) I - omega R - mu R
#   dC = delta theta I - mu C
#   dV = phi - omega_v V - mu V,   lambda = beta(t)(I + gamma C)/N
# Table-1 parameters and Table-4 fitted betas (beta1..4 = 0.19,0.22,0.14,0.01),
# reporting incidence(t)=f I(t), f=0.2, N0=1.4e6, R0=2.49. Vaccination moves
# phi = kappa*v*N0/T from S to V over the campaign window (S->V).
#
# The renewal model is applied to the SAME outbreak: R_t is reconstructed from the
# SIRCV baseline new-infection flow with the model's intrinsic acute generation
# interval, then reduced multiplicatively by (1 - c(t) kappa v) over the identical
# campaign window (c(t) ramps 0->1 over the campaign, then holds), and propagated
# forward. Both models therefore see the same baseline and the same intervention;
# they differ only in transmission structure (nonlinear compartmental depletion +
# carrier reservoir vs linear renewal feedback on a data-reconstructed R_t).
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
OUT <- "revision/seir_comparison"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
PUB <- file.path(OUT, "published_seir_sc1.csv")   # provenance: xlsx export

## ---- parameters (Table 1 + Table 4) ----------------------------------------
yr   <- 365
b    <- 0.039/yr; mu <- 0.007/yr
delta<- 1/11.8; theta <- 0.029; alpha <- 0.025
omega<- 1/(2*yr); omegav <- 1/(4*yr)
gamma<- 0.0009; f <- 0.2; N0 <- 1.4e6; R0 <- 2.49
beta_p <- c(0.19, 0.22, 0.14, 0.01)               # P1..P4 (Table 4)
# day index: Jan1=0 .. Jun12=162 (163 days). Feb10=40, Feb24=54, Mar17=75.
TMAX <- 162; D0_P2 <- 40; D0_P3 <- 54; D0_P4 <- 75
beta_t <- function(t) if (t < D0_P2) beta_p[1] else if (t < D0_P3) beta_p[2] else if (t < D0_P4) beta_p[3] else beta_p[4]
WINDOWS <- list(`Scenario 1 (early)` = c(D0_P2, D0_P3),  # P2, 14 d
                `Scenario 2 (late)`  = c(D0_P3, D0_P4),  # P3, 21 d
                `Scenario 3 (comb.)` = c(D0_P2, D0_P4))  # P2+P3, 35 d
COVER <- c(0.70, 0.50, 0.30); VE <- 0.8

## ---- SIRCV integrator (RK4, sub-daily) -------------------------------------
# Returns per-day: reported prevalence f*I (paper's "cases"), and daily NEW
# infections J = integral of lambda*S over the day (renewal input).
sircv <- function(kappa=0, window=NULL) {
  Tcamp <- if (is.null(window)) Inf else (window[2]-window[1])
  phi_rate <- if (is.null(window) || kappa==0) 0 else kappa*VE*N0/Tcamp   # people/day
  phi_t <- function(t) if (!is.null(window) && t>=window[1] && t<window[2]) phi_rate else 0
  deriv <- function(t, y) {
    S<-y[1]; I<-y[2]; R<-y[3]; C<-y[4]; V<-y[5]; N<-S+I+R+C+V
    lam <- beta_t(t)*(I + gamma*C)/N
    phi <- min(phi_t(t), max(S,0))                 # cannot vaccinate more than S
    c(b*N - lam*S - phi + omegav*V + omega*R - mu*S,
      lam*S - delta*I - mu*I,
      delta*(1-theta-alpha)*I - omega*R - mu*R,
      delta*theta*I - mu*C,
      phi - omegav*V - mu*V,
      lam*S)                                       # 6th slot: cumulative new-inf flow
  }
  dt <- 0.05; y <- c(N0-1, 1, 0, 0, 0, 0)
  fI <- numeric(TMAX+1); J <- numeric(TMAX+1); fI[1] <- f*y[2]
  cumJ_prev <- 0
  for (day in 1:TMAX) {
    for (k in 1:round(1/dt)) {
      t <- (day-1) + (k-1)*dt
      k1<-deriv(t,y); k2<-deriv(t+dt/2,y+dt/2*k1); k3<-deriv(t+dt/2,y+dt/2*k2); k4<-deriv(t+dt,y+dt*k3)
      y <- y + dt/6*(k1+2*k2+2*k3+k4)
    }
    fI[day+1] <- f*y[2]
    J[day+1]  <- y[6]-cumJ_prev; cumJ_prev <- y[6]   # new infections during this day
  }
  list(fI=fI, J=J, cumfI=sum(fI), cumJ=sum(J))
}

base <- sircv(kappa=0)
cat(sprintf("SIRCV baseline: cumulative f*I = %.0f  (published xlsx no_vacc = 8741)\n", base$cumfI))
pub <- read.csv(PUB, stringsAsFactors=FALSE)
cat(sprintf("baseline day-1..3 f*I: %.4f %.4f %.4f  | xlsx: %.4f %.4f %.4f\n",
            base$fI[1],base$fI[2],base$fI[3], pub$no_vacc[1],pub$no_vacc[2],pub$no_vacc[3]))

## ---- renewal machinery on the SIRCV baseline new-infection flow ------------
# intrinsic acute GI of this SIR: exponential, mean 1/delta = 11.8 d (CV=1).
gi_disc <- function(mean, cv, maxs=60) {
  shape <- 1/cv^2; rate <- shape/mean
  G <- function(x) pgamma(x, shape, rate)
  w <- sapply(1:maxs, function(s) G(s)-G(s-1)); w/sum(w)
}
w <- gi_disc(11.8, 1.0)
Lambda <- function(J, t, w) { s <- seq_len(min(t-1, length(w))); if(!length(s)) 0 else sum(w[s]*J[t-s+1]) }
# reconstruct R_t (index t corresponds to day t-1; J[1]=day0)
Rt_hat <- rep(NA_real_, TMAX+1)
for (t in 2:(TMAX+1)) { L <- Lambda(base$J, t, w); Rt_hat[t] <- if (L>0) base$J[t]/L else NA }
# period-mean reconstructed R_t vs SEIR E(Rt)
per_idx <- list(P2=(D0_P2+1):(D0_P3), P3=(D0_P3+1):(D0_P4), P4=(D0_P4+1):(TMAX+1))
cat("\nRenewal-reconstructed mean R_t vs SEIR E(R_t) [Table 4/5]:\n")
for (nm in names(per_idx)) cat(sprintf("  %s: renewal %.2f  vs SEIR %s\n", nm,
    mean(Rt_hat[per_idx[[nm]]], na.rm=TRUE), c(P2="2.87",P3="1.86",P4="0.15")[nm]))

# renewal counterfactual: reduce R_t by (1 - c(t) kappa VE), c ramps over campaign
renewal_vacc <- function(kappa, window) {
  kv <- kappa*VE; t0<-window[1]; t1<-window[2]
  cfun <- function(day) if (day<t0) 0 else if (day>=t1) 1 else (day-t0)/(t1-t0)
  Jv <- base$J
  for (t in 2:(TMAX+1)) {
    day <- t-1
    if (day < t0) { Jv[t] <- base$J[t]; next }      # pre-campaign identical
    Rv <- Rt_hat[t]*(1 - cfun(day)*kv)
    L  <- Lambda(Jv, t, w)
    Jv[t] <- if (is.na(Rv)||is.na(L)) base$J[t] else Rv*L
  }
  1 - sum(Jv)/sum(base$J)                            # reduction rate on new infections
}

## ---- run both models over 3 scenarios x 3 coverages ------------------------
rows <- list()
pubcols <- list(`Scenario 1 (early)`="sc1", `Scenario 2 (late)`="sc2", `Scenario 3 (comb.)`="sc3")
pubfile <- c(sc1=file.path(OUT,"published_seir_sc1.csv"), sc2=file.path(OUT,"published_seir_sc2.csv"),
             sc3=file.path(OUT,"published_seir_sc3.csv"))
for (sc in names(WINDOWS)) {
  win <- WINDOWS[[sc]]; pc <- pubcols[[sc]]; pubd <- read.csv(pubfile[pc], stringsAsFactors=FALSE)
  base_cum <- sum(pubd$no_vacc)
  for (kap in COVER) {
    # --- SEIR (my re-implementation): stock AND flow reduction (verified equal) ---
    sv <- sircv(kappa=kap, window=win)
    seir_red      <- 1 - sv$cumfI/base$cumfI          # on reported stock f*I
    seir_flow_red <- 1 - sv$cumJ/base$cumJ            # on new-infection flow (matches renewal scale)
    seir_av  <- base$cumfI - sv$cumfI
    # --- published SEIR (xlsx) for reference (NOT the like-for-like comparator) ---
    pcol <- paste0(pc, "_", round(kap*100)); pub_cum <- sum(pubd[[pcol]])
    pub_red <- 1 - pub_cum/base_cum; pub_av <- base_cum - pub_cum
    # --- renewal (matched intervention), reduction on new-infection flow ---
    ren_red <- renewal_vacc(kap, win); ren_av <- ren_red*base$cumfI
    rows[[length(rows)+1L]] <- data.frame(
      scenario=sc, coverage=kap,
      seir_reduction=100*seir_red, seir_flow_reduction=100*seir_flow_red, seir_averted=seir_av,
      published_reduction=100*pub_red, published_averted=pub_av,
      renewal_reduction=100*ren_red, renewal_averted=ren_av,
      stringsAsFactors=FALSE)
  }
}
res <- do.call(rbind, rows)
# Like-for-like gap: re-implementation SEIR (flow) vs renewal (flow). Using the PUBLISHED reduction
# instead would fold in the baseline-reproduction discrepancy (see report section 2).
res$gap_points <- res$seir_flow_reduction - res$renewal_reduction
BASE_PUB <- 8741                                    # published no-vacc cumulative (xlsx)
# NB: averted-on-common-baseline are DEFINITIONAL rescalings, not independent validation:
# published_reduction/100*BASE_PUB is identically the published averted. Kept only for reference.
res$seir_averted_common    <- res$published_reduction/100*BASE_PUB
res$renewal_averted_common <- res$renewal_reduction/100*BASE_PUB
write.csv(res, file.path(OUT,"comparison_results.csv"), row.names=FALSE)

cat("\n================ SEIR (re-impl / published) vs RENEWAL ================\n")
cat("reduction rate (%) of cumulative cases vs no-vaccination baseline\n\n")
fmt <- function(x) sprintf("%5.1f", x)
cat(sprintf("%-22s %5s | %8s %8s %8s %6s | %8s %8s\n","scenario","cov",
            "SEIR*","publ.","renewal","gap","SEIRavt","renAvt"))
cat("(* = my re-implementation; publ. = paper xlsx; averted on common baseline 8741)\n")
for (i in 1:nrow(res)) with(res[i,], cat(sprintf(
  "%-22s %4.0f%% | %8s %8s %8s %5.1f | %8.0f %8.0f\n",
  scenario, coverage*100, fmt(seir_reduction), fmt(published_reduction),
  fmt(renewal_reduction), gap_points, seir_averted_common, renewal_averted_common)))
cat("\nWrote", file.path(OUT,"comparison_results.csv"), "\n")

## ---- figure ----------------------------------------------------------------
ok <- requireNamespace("ggplot2", quietly=TRUE)
if (ok) {
  library(ggplot2)
  long <- rbind(
    data.frame(scenario=res$scenario, coverage=res$coverage, model="Full dynamic (SIRCV)", reduction=res$seir_reduction),
    data.frame(scenario=res$scenario, coverage=res$coverage, model="Renewal equation",     reduction=res$renewal_reduction))
  long$cov <- factor(paste0(long$coverage*100,"%"), levels=c("30%","50%","70%"))
  cols <- c("Full dynamic (SIRCV)"="#0072B2","Renewal equation"="#D55E00")
  g <- ggplot(long, aes(cov, reduction, fill=model)) +
    geom_col(position=position_dodge(0.7), width=0.62) +
    geom_text(aes(label=sprintf("%.0f",reduction)), position=position_dodge(0.7), vjust=-0.3, size=3) +
    facet_wrap(~scenario) + scale_fill_manual(values=cols) +
    scale_y_continuous(limits=c(0,100), expand=expansion(mult=c(0,0.08))) +
    labs(x="Vaccine coverage", y="Reduction in cumulative cases (%)", fill=NULL,
         title="Vaccine impact on the 2015 Kampala typhoid outbreak: full dynamic vs renewal model",
         subtitle="Matched intervention (efficacy 0.8, identical campaign windows). Both reduce transmission by coverage x efficacy.") +
    theme_minimal(base_size=11) +
    theme(legend.position="top", panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", size=12), plot.subtitle=element_text(colour="grey40", size=9))
  ggsave(file.path(OUT,"fig_seir_vs_renewal.png"), g, width=10, height=4.2, dpi=300)
  cat("Wrote", file.path(OUT,"fig_seir_vs_renewal.png"), "\n")
}
