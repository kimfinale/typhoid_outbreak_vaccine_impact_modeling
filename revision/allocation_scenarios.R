# =============================================================================
# Weekly-allocation scenario framework for the PPV decomposition.
#   source("revision/allocation_scenarios.R")  # provides run_allocation_scenarios()
# The cumulative PPV pins sum(I) = pi*sum(S) but NOT the temporal allocation.
# Five allocations of true typhoid I(t) (each satisfies the cumulative constraint),
# under the mechanistically REQUIRED background-fixed counterfactual (a typhoid
# vaccine cannot reduce non-typhoid fever, so B(t)=S(t)-I(t) is held fixed).
# Consequence: suspected_reduction = pi * true_reduction EXACTLY for every
# allocation; the allocation only changes the magnitude of the true reduction
# (via how much true typhoid the campaign catches).
# =============================================================================
suppressMessages({ library(ggplot2) })
for (f in c("data_prep.R","gi.R","renewal_core.R","ppv.R","epiestim_rt.R"))
  source(file.path("renewal/R", f))

il <- function(x) 1/(1+exp(-x))
.norm_to <- function(x, target) { x <- pmax(x,0); s <- sum(x); if (s>0) x*target/s else x }

# ---- five allocations: each returns I(t) with sum(I) = pi*sum(S) --------------
alloc_multiplicative <- function(S, pi) pi*S
alloc_additive <- function(S, pi) ppv_true_incidence(S, pi)          # excess over flat endemic baseline
alloc_timevarying <- function(S, pi, gamma=2) {                       # declining PPV -> front-loaded typhoid
  z <- (seq_along(S)-mean(seq_along(S)))/sd(seq_along(S))
  g <- function(lam) sum(il(lam - gamma*z)*S) - pi*sum(S)
  lam <- tryCatch(uniroot(g, c(-12,12))$root, error=function(e) qlogis(pi))
  il(lam - gamma*z)*S
}
alloc_renewal <- function(S, pi) {                                    # smooth transmission-consistent typhoid
  sm <- as.numeric(stats::filter(S, c(1,2,3,2,1)/9, sides=2)); sm[is.na(sm)] <- S[is.na(sm)]
  .norm_to(sm, pi*sum(S))
}
alloc_reporting <- function(S, pi) {                                  # sharpest spike attributed to background
  d <- c(0, diff(S)); t_spike <- which.max(d)
  Sds <- S; if (t_spike>1 && t_spike<length(S)) Sds[t_spike] <- mean(S[c(t_spike-1,t_spike+1)])
  x <- pmax(Sds - (1-pi)*mean(Sds), 0); .norm_to(x, pi*sum(S))
}
ALLOCS <- list(Multiplicative=alloc_multiplicative, Additive=alloc_additive,
               `Time-varying PPV`=alloc_timevarying, `Renewal-consistent`=alloc_renewal,
               `Reporting-shock`=alloc_reporting)

# ---- background-fixed vaccine counterfactual on a given I(t) ------------------
vaccinate <- function(I, Bkg, w, t_eff, coverage, psi_T, delta, cshape="step") {
  Tn <- length(I); K <- length(w)
  Rt <- reconstruct_Rt(I, w)$Rt
  c_t <- protected_fraction(Tn, max(1,t_eff-delta), t_eff, cshape)
  Rt_vax <- Rt*(1 - c_t*coverage*psi_T)
  Iv <- I
  for (t in seq_len(Tn)) { if (c_t[t]<=0) next; sm <- min(K,t-1L)
    lam <- if (sm>=1L) sum(w[1:sm]*Iv[t-(1:sm)]) else 0
    Iv[t] <- if (is.na(Rt_vax[t])) I[t] else Rt_vax[t]*lam }
  list(Iv=Iv, Sv=Iv+Bkg)
}

run_allocation_scenarios <- function(outbreak, cfg=NULL, coverage=0.80, psi_T=NULL,
                                     rel_teff=0.35, gamma=2) {
  if (is.null(cfg)) cfg <- yaml::read_yaml("renewal/config.yml")
  if (is.null(psi_T)) psi_T <- 0.8*cfg$vaccine$psi_mean
  prep <- prep_outbreaks(cfg); w <- gi_from_config(cfg)
  S <- prep$series[[outbreak]]; if (is.null(S)) stop("not a renewal outbreak: ", outbreak)
  Tn <- length(S)
  pic <- read.csv("latent_class_ppv/tables/final_pi_community.csv", check.names=FALSE, stringsAsFactors=FALSE)
  pi_o <- pic$pi_med[match(outbreak, pic$study)]
  if (is.na(pi_o)) { pi_o <- il(median(qlogis(pmin(pmax(pic$pi_med,1e-3),1-1e-3)))); attr(pi_o,"source")<-"population" }
  t_eff <- max(2L, round(rel_teff*Tn))
  delta <- round((cfg$timing$campaign_duration_days+cfg$timing$immunity_onset_days)/cfg$step_days)
  res <- lapply(names(ALLOCS), function(nm) {
    I <- if (nm=="Time-varying PPV") ALLOCS[[nm]](S, pi_o, gamma) else ALLOCS[[nm]](S, pi_o)
    Bkg <- pmax(S - I, 0)
    vc <- vaccinate(I, Bkg, w, t_eff, coverage, psi_T, delta, cfg$timing$c_shape)
    data.frame(alloc=nm, week=seq_len(Tn), S=S, I=I, Bkg=Bkg, Iv=vc$Iv, Sv=vc$Sv,
               true_red=100*(sum(I)-sum(vc$Iv))/sum(I),
               susp_red=100*(sum(S)-sum(vc$Sv))/sum(S), stringsAsFactors=FALSE)
  })
  df <- do.call(rbind, res)
  df$alloc <- factor(df$alloc, levels=names(ALLOCS))
  list(df=df, pi=pi_o, t_eff=t_eff, Tn=Tn, outbreak=outbreak,
       summ=unique(df[,c("alloc","true_red","susp_red")]))
}
