# =============================================================================
# How close can the renewal R_t get to the dynamic-model R_t by tuning the
# generation interval?  Reconstruct R_hat_t from the SIRCV baseline new-infection
# flow under a grid of generation intervals (mean, CV, and a carrier-tail
# mixture) and compare to the SIRCV's exact daily R_t(t).
#   Rscript revision/seir_comparison/gi_tuning.R
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
OUT <- "revision/seir_comparison"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

## ---- SIRCV parameters (Table 1 + Table 4) ----------------------------------
yr<-365; b<-0.039/yr; mu<-0.007/yr; delta<-1/11.8; theta<-0.029; alpha<-0.025
omega<-1/(2*yr); omegav<-1/(4*yr); gamma<-0.0009; f<-0.2; N0<-1.4e6
beta_p<-c(0.19,0.22,0.14,0.01); TMAX<-162; D0_P2<-40; D0_P3<-54; D0_P4<-75
beta_t<-function(t) if(t<D0_P2)beta_p[1] else if(t<D0_P3)beta_p[2] else if(t<D0_P4)beta_p[3] else beta_p[4]
Rt_coef <- (mu+delta*theta*gamma)/(mu*(delta+mu))     # E(Rt)=beta*S/N*Rt_coef

## ---- baseline SIRCV: new-infection flow J_t and exact daily R_t -------------
sim <- function() {
  deriv<-function(t,y){S<-y[1];I<-y[2];R<-y[3];C<-y[4];V<-y[5];N<-S+I+R+C+V
    lam<-beta_t(t)*(I+gamma*C)/N
    c(b*N-lam*S+omegav*V+omega*R-mu*S, lam*S-delta*I-mu*I,
      delta*(1-theta-alpha)*I-omega*R-mu*R, delta*theta*I-mu*C, -omegav*V-mu*V, lam*S)}
  dt<-0.05; y<-c(N0-1,1,0,0,0,0); J<-numeric(TMAX+1); S_day<-numeric(TMAX+1)
  S_day[1]<-y[1]; cj<-0
  for(day in 1:TMAX){for(k in 1:round(1/dt)){t<-(day-1)+(k-1)*dt
    k1<-deriv(t,y);k2<-deriv(t+dt/2,y+dt/2*k1);k3<-deriv(t+dt/2,y+dt/2*k2);k4<-deriv(t+dt,y+dt*k3)
    y<-y+dt/6*(k1+2*k2+2*k3+k4)}
    J[day+1]<-y[6]-cj; cj<-y[6]; S_day[day+1]<-y[1]}
  N_day <- N0  # ~constant
  list(J=J, S=S_day, Rt_seir = sapply(0:TMAX, function(d) beta_t(d)*S_day[d+1]/N_day*Rt_coef))
}
B <- sim()

## ---- generation intervals & reconstruction ---------------------------------
gi_gamma <- function(mean, cv, maxs=90){ sh<-1/cv^2; rt<-sh/mean
  w<-sapply(1:maxs, function(s) pgamma(s,sh,rt)-pgamma(s-1,sh,rt)); w/sum(w) }
gi_mix <- function(mean_acute, cv_acute, w_carrier, mean_carrier, maxs=90){
  wa<-gi_gamma(mean_acute,cv_acute,maxs); wc<-gi_gamma(mean_carrier,1.0,maxs)
  w<-(1-w_carrier)*wa + w_carrier*wc; w/sum(w) }
recon <- function(J,w){ R<-rep(NA_real_,length(J))
  for(t in 2:length(J)){ s<-seq_len(min(t-1,length(w))); L<-sum(w[s]*J[t-s+1]); R[t]<-if(L>0)J[t]/L else NA }
  R }
pidx <- list(P2=(D0_P2+1):D0_P3, P3=(D0_P3+1):D0_P4, P4=(D0_P4+1):(TMAX+1))  # day+1 indices
pmean <- function(R,idx) sapply(pidx, function(ix) mean(R[ix], na.rm=TRUE))
# incidence-weighted daily RMSE of R_hat vs SEIR R_t over the informative window
score <- function(R){ ix<-21:141; wgt<-B$J[ix]; d<-(R[ix]-B$Rt_seir[ix])
  sqrt(sum(wgt*d^2, na.rm=TRUE)/sum(wgt, na.rm=TRUE)) }

seir_pm <- pmean(B$Rt_seir, pidx)   # SEIR period means (≈ E(Rt) Table 5)
cat(sprintf("SIRCV E(Rt) period means: P2=%.2f P3=%.2f P4=%.2f  (paper 2.87/1.86/0.15)\n\n",
            seir_pm[1],seir_pm[2],seir_pm[3]))

## ---- (1) sweep GI mean at CV=1, and CV at mean=11.8 -------------------------
cat("=== GI SWEEP: reconstructed R_hat period means vs SEIR ===\n")
cat(sprintf("%-28s %6s %6s %6s | %6s\n","generation interval","P2","P3","P4","wRMSE"))
report <- function(lab,w){ R<-recon(B$J,w); pm<-pmean(R,pidx)
  cat(sprintf("%-28s %6.2f %6.2f %6.2f | %6.3f\n",lab,pm[1],pm[2],pm[3],score(R))); invisible(R) }
cat("-- SEIR target --\n"); cat(sprintf("%-28s %6.2f %6.2f %6.2f | %6.3f\n","SIRCV E(Rt)",seir_pm[1],seir_pm[2],seir_pm[3],0))
cat("-- vary mean (CV=1, exponential) --\n")
for(m in c(11.8,13,14,15,16,18)) report(sprintf("gamma mean=%.1f CV=1.0",m), gi_gamma(m,1.0))
cat("-- vary CV (mean=11.8) --\n")
for(cv in c(1.0,0.85,0.71,0.6,0.5)) report(sprintf("gamma mean=11.8 CV=%.2f",cv), gi_gamma(11.8,cv))
cat("-- carrier mixture (acute Exp 11.8 + 10% long tail) --\n")
for(mc in c(30,60,90)) report(sprintf("90%%Exp11.8 + 10%%Exp%d",mc), gi_mix(11.8,1.0,0.10,mc))

## ---- (2) best-fit (mean, CV) by minimizing incidence-weighted daily RMSE ----
obj <- function(p){ m<-p[1]; cv<-p[2]; if(m<5||m>30||cv<0.3||cv>1.3) return(1e6); score(recon(B$J, gi_gamma(m,cv))) }
best_all <- optim(c(14,0.8), obj, method="Nelder-Mead")
# growth-phase-only fit (P2-P3 window, days 41-74)
score_g <- function(R){ ix<-41:74; wgt<-B$J[ix]; d<-(R[ix]-B$Rt_seir[ix]); sqrt(sum(wgt*d^2,na.rm=TRUE)/sum(wgt,na.rm=TRUE)) }
obj_g <- function(p){ m<-p[1]; cv<-p[2]; if(m<5||m>30||cv<0.3||cv>1.3) return(1e6); score_g(recon(B$J, gi_gamma(m,cv))) }
best_g <- optim(c(15,0.7), obj_g, method="Nelder-Mead")
cat(sprintf("\nBest-fit GI (all-phase wRMSE): mean=%.1f d, CV=%.2f  -> wRMSE=%.3f\n", best_all$par[1],best_all$par[2],best_all$value))
Rbest <- recon(B$J, gi_gamma(best_all$par[1],best_all$par[2])); pmb<-pmean(Rbest,pidx)
cat(sprintf("   period means: P2=%.2f P3=%.2f P4=%.2f  (SEIR 2.87/1.86/0.15)\n", pmb[1],pmb[2],pmb[3]))
cat(sprintf("Best-fit GI (growth-phase only): mean=%.1f d, CV=%.2f  -> growth wRMSE=%.3f\n", best_g$par[1],best_g$par[2],best_g$value))
Rbg <- recon(B$J, gi_gamma(best_g$par[1],best_g$par[2])); pmg<-pmean(Rbg,pidx)
cat(sprintf("   period means: P2=%.2f P3=%.2f P4=%.2f\n", pmg[1],pmg[2],pmg[3]))

## ---- (3) figure: daily R_t trajectories ------------------------------------
if(requireNamespace("ggplot2",quietly=TRUE)){ library(ggplot2)
  Racute <- recon(B$J, gi_gamma(11.8,1.0))
  df <- rbind(
    data.frame(day=0:TMAX, Rt=B$Rt_seir,            series="SIRCV dynamic model  E(R_t)"),
    data.frame(day=0:TMAX, Rt=Racute,               series="Renewal  acute GI (mean 11.8, CV 1)"),
    data.frame(day=0:TMAX, Rt=Rbest,                series=sprintf("Renewal  best-fit GI (mean %.0f, CV %.2f)",best_all$par[1],best_all$par[2])))
  df <- df[df$day>=15 & df$day<=150,]
  cols <- c("SIRCV dynamic model  E(R_t)"="black")
  cols[grep("acute",names(table(df$series)),value=TRUE)]<-"#D55E00"
  cols[grep("best-fit",names(table(df$series)),value=TRUE)]<-"#0072B2"
  g <- ggplot(df, aes(day, Rt, colour=series, linetype=series)) +
    geom_vline(xintercept=c(D0_P2,D0_P3,D0_P4), colour="grey80", linewidth=0.4) +
    geom_hline(yintercept=1, colour="grey60", linetype="dotted") +
    geom_line(linewidth=0.9) +
    scale_colour_manual(values=c("black","#D55E00","#0072B2")) +
    scale_linetype_manual(values=c("solid","solid","solid")) +
    annotate("text", x=c(47,64,110), y=3.4, label=c("P2","P3","P4"), colour="grey50", size=3.5) +
    labs(x="Day (0 = Jan 1, 2015)", y=expression(R[t]), colour=NULL, linetype=NULL,
         title="Renewal R_t vs dynamic-model R_t on the 2015 Kampala outbreak: effect of the generation interval",
         subtitle="Tuning the GI aligns the renewal R_t with the dynamic model in the growth phases; the decline phase cannot be matched by any single GI.") +
    theme_minimal(base_size=11) + theme(legend.position="top", panel.grid.minor=element_blank(),
         plot.title=element_text(face="bold",size=11.5), plot.subtitle=element_text(colour="grey40",size=9))
  ggsave(file.path(OUT,"fig_gi_tuning.png"), g, width=10, height=4.6, dpi=300)
  cat("\nWrote", file.path(OUT,"fig_gi_tuning.png"), "\n")
}
