# =============================================================================
# Kabwama (Kampala 2015) comprehensive analysis — 01: model functions.
# Sourced by 02_analysis.R. Pure functions; no side effects.
#
# Observation model:  suspected S_t = typhoid symptomatic I_t + non-typhoid
# background B_t (NTS + other bacterial + non-bacterial). PPV pi pins the
# cumulative split (sum I = pi * sum S). The vaccine acts ONLY on typhoid, so
# B_t is invariant to TCV regardless of its bacterial composition (this is why
# the observed suspected reduction is diluted to pi x the true reduction).
# Asymptomatic typhoid is carried as an overlay (total infections = I / p_symp)
# and explicitly in the SEIAR-C transmission model.
# =============================================================================

## ---- generation interval (weekly bins) -------------------------------------
gi_weekly <- function(mean_days = 14, cv = 0.6, step_days = 7, maxw = 12) {
  shape <- 1 / cv^2; rate <- shape / mean_days
  w <- vapply(1:maxw, function(s) pgamma(s * step_days, shape, rate) -
                pgamma((s - 1) * step_days, shape, rate), numeric(1))
  w / sum(w)
}

## ---- additive de-backgrounding (matches renewal/R/ppv.R::ppv_true_incidence)
# B = (1-pi)*mean(S) held constant; I = max(S-B,0) renormalised so sum(I)=pi*sum(S).
debackground <- function(S, pi) {
  if (pi >= 1) return(S)
  B <- (1 - pi) * mean(S)
  I <- pmax(S - B, 0)
  s <- sum(I)
  if (s > 0) I * (pi * sum(S) / s) else I
}

## ---- renewal reconstruction + background-fixed vaccine counterfactual ------
# The primary source-plus-propagated counterfactual is implemented centrally as
# renewal/R/renewal_core.R::additive_source_counterfactual. The local pure-
# renewal and static functions below are retained for the paired historical
# theta-weighted comparator and structural sensitivities.
Lambda_t <- function(I, t, w) {
  smax <- min(t - 1L, length(w)); if (smax < 1L) return(0)
  sum(w[1:smax] * I[t - (1:smax)])
}
reconstruct_R <- function(I, w) {
  Tn <- length(I); R <- rep(NA_real_, Tn)
  for (t in 2:Tn) { L <- Lambda_t(I, t, w); R[t] <- if (L > 0) I[t] / L else NA_real_ }
  R
}
# Windowed Cori posterior-mean R_t (smoothed; gamma(a,b) prior, mean a/b). Stable
# at low counts and at the de-backgrounded onset, unlike the raw ratio.
cori_R <- function(I, w, tau = 3, a = 1, b = 0.2) {
  Tn <- length(I); Lam <- vapply(seq_len(Tn), function(t) Lambda_t(I, t, w), numeric(1))
  R <- rep(NA_real_, Tn)
  for (t in seq_len(Tn)) { win <- max(1, t - tau + 1):t
    sL <- sum(Lam[win]); if (sL > 0) R[t] <- (a + sum(I[win])) / (b + sL) }
  R
}

## ---- static direct-effect model (Halloran; common-source arm) --------------
# Post-protection incidence reduced by the direct VE only (no emergent indirect
# effect). Appropriate for the common-source component of a mixed outbreak with
# ongoing exposure. Returns cases averted.
static_averted <- function(I, t_eff, kappa, psi) {
  post <- I; if (t_eff > 1) post[seq_len(min(t_eff - 1L, length(I)))] <- 0
  kappa * psi * sum(post)
}
# Reduce typhoid R_t by (1 - c(t)*kappa*psiT) from protection week t_eff (step,
# or linear ramp over ramp_wk), propagate typhoid incidence forward. Background
# B_t = S_t - I_t held FIXED. Returns list(I, Iv, Sv, ...).
vaccinate_renewal <- function(S, I, w, R, t_eff, kappa, psiT, ramp_wk = 1) {
  Tn <- length(I); Iv <- I; B <- S - I
  cfun <- function(t) if (t < t_eff) 0 else min(1, (t - t_eff + 1) / ramp_wk)
  for (t in 2:Tn) {
    if (t < t_eff) { Iv[t] <- I[t]; next }
    Rv <- R[t] * (1 - cfun(t) * kappa * psiT)
    L  <- Lambda_t(Iv, t, w)
    Iv[t] <- if (is.na(Rv) || is.na(L)) I[t] else Rv * L
  }
  Sv <- Iv + B
  list(I = I, Iv = Iv, Sv = Sv, B = B,
       typ_averted = sum(I) - sum(Iv),
       true_red = 1 - sum(Iv) / sum(I),
       obs_red  = 1 - sum(Sv) / sum(S))         # = pi * true_red by construction
}

## ---- multiplicative observation model (alternative): I = pi*S --------------
vaccinate_multiplicative <- function(S, w, R, t_eff, kappa, psiT, ramp_wk = 1) {
  # reconstruct on suspected directly; % reduction is pi-invariant
  Tn <- length(S); Sv <- S
  cfun <- function(t) if (t < t_eff) 0 else min(1, (t - t_eff + 1) / ramp_wk)
  for (t in 2:Tn) {
    if (t < t_eff) { Sv[t] <- S[t]; next }
    Rv <- R[t] * (1 - cfun(t) * kappa * psiT)
    L  <- Lambda_t(Sv, t, w)
    Sv[t] <- if (is.na(Rv) || is.na(L)) S[t] else Rv * L
  }
  list(Sv = Sv, red = 1 - sum(Sv) / sum(S))     # true_red == obs_red here
}

## ---- SEIAR-C mechanistic model (asymptomatic + carrier explicit) -----------
# Susceptible->Exposed->(symptomatic Is + asymptomatic Ia)->Recovered, fraction
# theta -> chronic Carrier. Waterborne/direct transmission folded into a
# time-varying beta(t) taken from the reconstructed typhoid R_t, so the model
# reproduces the observed typhoid curve while resolving the unobserved
# compartments (asymptomatic infections, carriers). Vaccination moves kappa*psi
# of S to V over a 14-day campaign starting at protection week t_eff.
# Time unit = days (RK4, dt=0.5); returns weekly symptomatic/total incidence.
seiar_c <- function(R_week, par, vacc = NULL) {
  with(par, {
    Tw <- length(R_week); Td <- Tw * 7L
    Rday <- approx((seq_len(Tw) - 0.5) * 7, ifelse(is.na(R_week), 1, R_week),
                   xout = 0:(Td - 1), rule = 2)$y
    sig <- 1 / incub_d; gam <- 1 / infper_d
    Deff <- (p_symp + xi * (1 - p_symp)) / gam
    y <- c(S = N - seed / p_symp, E = seed / p_symp * 0.5, Is = seed,
           Ia = seed * (1 - p_symp) / p_symp, C = 0, R = 0, V = 0)
    vcov <- if (is.null(vacc)) 0 else vacc$kappa * psi
    t_eff_day <- if (is.null(vacc)) Inf else vacc$t_eff * 7
    dt <- 0.5; inc_s_day <- numeric(Td); inc_tot_day <- numeric(Td)
    rates <- function(y, betat, phi) {
      Nnow <- sum(y[1:7]); lam <- betat * (y["Is"] + xi * y["Ia"] + gc * y["C"]) / Nnow
      c(S = -lam * y["S"] - phi,
        E =  lam * y["S"] - sig * y["E"],
        Is = p_symp * sig * y["E"] - gam * y["Is"],
        Ia = (1 - p_symp) * sig * y["E"] - gam * y["Ia"],
        C =  theta * gam * (y["Is"] + y["Ia"]) - muC * y["C"],
        R = (1 - theta) * gam * (y["Is"] + y["Ia"]) - muC * y["R"],
        V =  phi - muC * y["V"])
    }
    for (day in 1:Td) {
      betat <- Rday[day] / Deff
      phi <- if (day >= t_eff_day && day < t_eff_day + 14 && vcov > 0) vcov * N / 14 else 0
      for (k in 1:(1/dt)) {
        k1 <- rates(y, betat, phi); k2 <- rates(y + dt/2 * k1, betat, phi)
        k3 <- rates(y + dt/2 * k2, betat, phi); k4 <- rates(y + dt * k3, betat, phi)
        y <- y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        inc_s_day[day]   <- inc_s_day[day]   + p_symp * sig * y["E"] * dt
        inc_tot_day[day] <- inc_tot_day[day] + sig * y["E"] * dt
      }
    }
    wk <- rep(seq_len(Tw), each = 7)
    list(inc_s_week  = as.numeric(tapply(inc_s_day, wk, sum)),
         inc_tot_week = as.numeric(tapply(inc_tot_day, wk, sum)),
         sympt_total = sum(inc_s_day), tot_infections = sum(inc_tot_day),
         carriers_end = as.numeric(y["C"]))
  })
}
