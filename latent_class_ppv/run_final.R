# =============================================================================
# FINAL fit (real data, beta-free): consolidate the historic paired blood/bone-
# marrow studies + the community outbreak positivity, estimate Se_BC(volume),
# Se_BM, hospital PPV phi_s and community PPV pi_o, and draw inferences.
#   Rscript latent_class_ppv/run_final.R
# =============================================================================
ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(cmdstanr); library(posterior); library(ggplot2); library(dplyr) })
for (d in c("figures","tables","outputs")) dir.create(file.path("latent_class_ppv", d), showWarnings = FALSE)
il <- function(x) 1/(1+exp(-x))
OUT_VOL <- 5    # explicit fallback only when study blood volume is unreported

# ---- 1. consolidate HISTORIC paired blood/bone-marrow data -------------------
H <- read.csv("latent_class_ppv/data/mogasale2016_paired_bc_bmc_seed.csv",
              check.names = FALSE, stringsAsFactors = FALSE)
cells <- as.matrix(H[, c("a_BCpos_BMpos","b_BCpos_BMneg","c_BCneg_BMpos","d_BCneg_BMneg")])
vol_h <- H$blood_vol_mL
vol_h[is.na(vol_h)] <- median(vol_h, na.rm = TRUE)     # impute Hirsowitz (no volume) at the median
cat("=== Historic paired studies (n=", nrow(H), ") ===\n", sep = "")
print(data.frame(study = H$study, vol_mL = round(vol_h,1),
                 BCpp = cells[,1], BCpn = cells[,2], BCnp = cells[,3], BCnn = cells[,4],
                 N = rowSums(cells), raw_BC_sens_vs_any = round((cells[,1]+cells[,2])/rowSums(cells[,1:3]),2)),
      row.names = FALSE)

# ---- 2. consolidate COMMUNITY OUTBREAK positivity from authoritative audit ----
# The merged audit carries specimen-scope corrections, source-overlap warnings,
# and explicit model-inclusion decisions. Do not rebuild this set from the raw
# outbreak summary or the older hand-maintained community CSV.
D <- read.csv("latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv",
              check.names = FALSE, stringsAsFactors = FALSE)
include_main <- D$include_for_main_ppv_model %in% c(TRUE,"TRUE","True")
audited_frame_n <- suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*$", "\\1",
                                                    D$audited_suspected_testing_frame)))
Vv <- data.frame(
  study = D$study[include_main],
  suspected = audited_frame_n[include_main],
  tested = suppressWarnings(as.numeric(D$audited_blood_tested[include_main])),
  confirmed_cx = suppressWarnings(as.numeric(D$audited_blood_positive[include_main])),
  blood_volume_ml = suppressWarnings(as.numeric(D$blood_volume_ml[include_main])),
  independence_group = D$independence_group[include_main],
  verification_adjustment_grade = D$verification_adjustment_grade[include_main],
  stringsAsFactors = FALSE
)
Vv$positivity <- Vv$confirmed_cx / Vv$tested
if (anyNA(Vv$suspected) || anyNA(Vv$tested) || anyNA(Vv$confirmed_cx) ||
    any(Vv$confirmed_cx > Vv$tested) || any(Vv$tested > Vv$suspected) ||
    any(Vv$tested <= 0)) {
  stop("Invalid audited blood-culture counts in main PPV set")
}
if (anyDuplicated(Vv$independence_group[nzchar(Vv$independence_group)])) {
  stop("Main PPV set contains duplicate independence groups; resolve overlap in audit")
}
if (!"blood_volume_ml" %in% names(Vv)) Vv$blood_volume_ml <- NA_real_
Vv$blood_volume_ml <- suppressWarnings(as.numeric(Vv$blood_volume_ml))
Vv$blood_volume_imputed <- is.na(Vv$blood_volume_ml) | Vv$blood_volume_ml <= 0
Vv$blood_volume_ml[Vv$blood_volume_imputed] <- OUT_VOL
Vv$bias_score <- unname(c(`negligible/low`=0, moderate=1, major=2)[Vv$verification_adjustment_grade])
if (anyNA(Vv$bias_score)) stop("Unmapped verification adjustment grade in main PPV set")
Vv$logit_tested_fraction <- qlogis((Vv$tested + 0.5) / (Vv$suspected + 1))
cat("\n=== Community outbreak points (tight; study-specific blood volume) ===\n")
print(data.frame(study = Vv$study, suspected=Vv$suspected, tested = Vv$tested, conf = Vv$confirmed_cx,
                 positivity = round(Vv$positivity,3), blood_volume_ml=Vv$blood_volume_ml,
                 volume_imputed=Vv$blood_volume_imputed, bias_score=Vv$bias_score), row.names = FALSE)

standata <- function(out_vol, bias_score=Vv$bias_score,
                     mu_delta_sd=0.5, beta_bias_mu=0.5,
                     beta_bias_sd=0.4, sigma_delta_sd=0.5) list(
  S = nrow(H), N_hist = as.integer(rowSums(cells)), cells = matrix(as.integer(cells), ncol = 4),
  logv_hist = log(vol_h),
  O = nrow(Vv), N_suspected=as.integer(Vv$suspected),
  n_out = as.integer(Vv$tested), k_out = as.integer(Vv$confirmed_cx),
  logv_out = log(ifelse(Vv$blood_volume_imputed, out_vol, Vv$blood_volume_ml)),
  bias_score=bias_score, logit_tested_fraction=Vv$logit_tested_fraction,
  mu_delta_sd=mu_delta_sd, beta_bias_mu=beta_bias_mu,
  beta_bias_sd=beta_bias_sd, sigma_delta_sd=sigma_delta_sd,
  se_bm_a = 12, se_bm_b = 2,          # Se_BM ~ Beta(12,2): mean 0.857, favours high (bone marrow)
  sp_a = 200, sp_b = 1,               # Sp ~ 0.995
  alpha1_mu = 0.36, alpha1_sd = 0.15) # volume slope from Antillon (0.51@2mL -> 0.65@10mL)

# ---- 3. fit -----------------------------------------------------------------
cat("\n=== compile + fit (real data) ===\n")
mod <- cmdstan_model("latent_class_ppv/stan/model_final.stan")
fit <- mod$sample(data = standata(OUT_VOL), seed = 2026, chains = 4, parallel_chains = 4,
                  iter_warmup = 1500, iter_sampling = 1500, adapt_delta = 0.98,
                  max_treedepth = 12, init = 0.5, refresh = 0, show_messages = FALSE)
dg <- fit$diagnostic_summary()
fit_summary <- fit$summary()
fit_diagnostics <- data.frame(
  chains = 4L,
  post_warmup_draws = 4L * 1500L,
  divergences = sum(dg$num_divergent),
  max_rhat = max(fit_summary$rhat, na.rm = TRUE),
  min_ess_bulk = min(fit_summary$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(fit_summary$ess_tail, na.rm = TRUE))
cat(sprintf("Diagnostics: divergences=%d, max Rhat=%.3f, min ESS=%.0f\n",
            fit_diagnostics$divergences, fit_diagnostics$max_rhat,
            fit_diagnostics$min_ess_bulk))
write.csv(fit_diagnostics, "latent_class_ppv/tables/fit_diagnostics.csv", row.names = FALSE)
qs <- function(nm){ v<-as.numeric(fit$draws(nm, format="draws_matrix")); c(med=median(v), lo=quantile(v,.05), hi=quantile(v,.95)) }

# ---- 4. inferences ----------------------------------------------------------
cat("\n=== KEY ESTIMATES (median [90% CrI]) ===\n")
for (nm in c("Se_BC_2mL","Se_BC_5mL","Se_BC_10mL","Se_BM","alpha1","tau","Sp_BC",
             "pi_population_median","sigma_pi","mu_delta","beta_bias","sigma_delta")) {
  q <- qs(nm); cat(sprintf("  %-11s %.2f [%.2f, %.2f]\n", nm, q[1], q[2], q[3]))
}
phi <- fit$summary("phi"); pi <- fit$summary("pi")
cat("\n-- hospital PPV phi_s (historic studies) --\n")
print(data.frame(study=H$study, phi_med=round(phi$median,2), lo=round(phi$q5,2), hi=round(phi$q95,2)), row.names=FALSE)
cat("\n-- community surveillance PPV pi_o (outbreaks) --\n")
print(data.frame(study=Vv$study, pi_med=round(pi$median,2), lo=round(pi$q5,2), hi=round(pi$q95,2)), row.names=FALSE)
cat(sprintf("\nSanity: median hospital phi = %.2f vs median community pi = %.2f (expect community <= hospital-ish)\n",
            median(phi$median), median(pi$median)))

# outbreak-volume sensitivity (the one remaining assumption)
cat("\n=== community pi_o sensitivity to assumed outbreak blood volume ===\n")
cat(sprintf("%-8s %-10s %-24s\n","vol_mL","Se_BC","median pi_o (per outbreak)"))
volume_sensitivity <- list()
for (vv in c(3,5,8,10)) {
  f2 <- mod$sample(data = standata(vv), seed=2026, chains=2, iter_warmup=800, iter_sampling=800,
                   adapt_delta=0.98, init=0.5, refresh=0, show_messages=FALSE)
  a0_v <- as.numeric(f2$draws("alpha0", format="draws_matrix"))
  a1_v <- as.numeric(f2$draws("alpha1", format="draws_matrix"))
  se_draw <- il(a0_v + a1_v * log(vv))
  pim <- f2$summary("pi")$median
  se <- median(se_draw)
  cat(sprintf("%-8d %-10.2f %s\n", vv, se,
      paste(sprintf("%s=%.2f", Vv$study, pim), collapse=", ")))
  volume_sensitivity[[length(volume_sensitivity) + 1L]] <- data.frame(
    assumed_blood_volume_mL = vv,
    Se_BC_median = se,
    Se_BC_lo = as.numeric(quantile(se_draw, 0.05)),
    Se_BC_hi = as.numeric(quantile(se_draw, 0.95)),
    study = Vv$study,
    pi_median = pim,
    stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, volume_sensitivity),
          "latent_class_ppv/tables/volume_sensitivity.csv", row.names = FALSE)

# Selection layer is structurally weakly identified. Pre-specified sensitivity
# scenarios vary both enrichment priors and audit-grade coding. These fits are
# reported as structural sensitivity, not model-selection competitors.
selection_scenarios <- list(
  no_enrichment = list(score=rep(0, nrow(Vv)), mu_sd=0.02, beta_mu=0,
                       beta_sd=0.02, sigma_sd=0.02),
  skeptical = list(score=Vv$bias_score, mu_sd=0.25, beta_mu=0,
                   beta_sd=0.25, sigma_sd=0.25),
  base = list(score=Vv$bias_score, mu_sd=0.5, beta_mu=0.5,
              beta_sd=0.4, sigma_sd=0.5),
  permissive = list(score=Vv$bias_score, mu_sd=1.0, beta_mu=0.75,
                    beta_sd=0.75, sigma_sd=1.0),
  binary_grade = list(score=as.numeric(Vv$bias_score > 0), mu_sd=0.5,
                      beta_mu=0.5, beta_sd=0.4, sigma_sd=0.5),
  no_grade = list(score=rep(0, nrow(Vv)), mu_sd=0.5, beta_mu=0.5,
                  beta_sd=0.4, sigma_sd=0.5)
)
selection_sensitivity <- list()
set.seed(20260716)
for (scenario_name in names(selection_scenarios)) {
  sc <- selection_scenarios[[scenario_name]]
  fs <- mod$sample(
    data=standata(OUT_VOL, bias_score=sc$score, mu_delta_sd=sc$mu_sd,
                  beta_bias_mu=sc$beta_mu, beta_bias_sd=sc$beta_sd,
                  sigma_delta_sd=sc$sigma_sd),
    seed=2026, chains=2, parallel_chains=2, iter_warmup=800,
    iter_sampling=800, adapt_delta=0.98, max_treedepth=12,
    init=0.5, refresh=0, show_messages=FALSE)
  saveRDS(fs$draws(), file.path("latent_class_ppv", "outputs",
                               paste0("selection_sensitivity_", scenario_name, ".rds")))
  dm_s <- posterior::as_draws_matrix(fs$draws(c("pi", "q0", "q1", "mu_pi", "sigma_pi",
                                                "mu_delta", "beta_bias", "sigma_delta")))
  pi_s <- as.matrix(dm_s[, paste0("pi[", seq_len(nrow(Vv)), "]"), drop=FALSE])
  q0_s <- as.matrix(dm_s[, paste0("q0[", seq_len(nrow(Vv)), "]"), drop=FALSE])
  q1_s <- as.matrix(dm_s[, paste0("q1[", seq_len(nrow(Vv)), "]"), drop=FALSE])
  pi_tested_s <- pi_s*q1_s/(pi_s*q1_s + (1-pi_s)*q0_s)
  deflation_s <- pi_s/pi_tested_s
  pred_s <- plogis(rnorm(nrow(dm_s), dm_s[,"mu_pi"], dm_s[,"sigma_pi"]))
  dg_s <- fs$diagnostic_summary(); sm_s <- fs$summary()
  for (i in seq_len(nrow(Vv))) {
    selection_sensitivity[[length(selection_sensitivity)+1L]] <- data.frame(
      scenario=scenario_name, study=Vv$study[i], bias_score=sc$score[i],
      pi_all_med=median(pi_s[,i]), pi_all_lo=quantile(pi_s[,i],.05), pi_all_hi=quantile(pi_s[,i],.95),
      pi_tested_med=median(pi_tested_s[,i]), pi_tested_lo=quantile(pi_tested_s[,i],.05), pi_tested_hi=quantile(pi_tested_s[,i],.95),
      deflation_med=median(deflation_s[,i]), deflation_lo=quantile(deflation_s[,i],.05), deflation_hi=quantile(deflation_s[,i],.95),
      pi_population_med=median(plogis(dm_s[,"mu_pi"])),
      sigma_pi_med=median(dm_s[,"sigma_pi"]),
      pi_new_med=median(pred_s), pi_new_lo=quantile(pred_s,.05), pi_new_hi=quantile(pred_s,.95),
      mu_delta_med=median(dm_s[,"mu_delta"]), beta_bias_med=median(dm_s[,"beta_bias"]),
      sigma_delta_med=median(dm_s[,"sigma_delta"]),
      divergences=sum(dg_s$num_divergent), max_rhat=max(sm_s$rhat,na.rm=TRUE),
      stringsAsFactors=FALSE)
  }
}
write.csv(do.call(rbind, selection_sensitivity),
          "latent_class_ppv/tables/selection_prior_sensitivity.csv", row.names=FALSE)

# ---- 5. figures + tables ----------------------------------------------------
# Se_BC vs volume curve with historic raw points
gv <- seq(1.5, 15, 0.25); a0 <- as.numeric(fit$draws("alpha0",format="draws_matrix")); a1 <- as.numeric(fit$draws("alpha1",format="draws_matrix"))
curve <- data.frame(vol=gv, med=sapply(gv, function(v) median(il(a0+a1*log(v)))),
                    lo=sapply(gv, function(v) quantile(il(a0+a1*log(v)),.05)),
                    hi=sapply(gv, function(v) quantile(il(a0+a1*log(v)),.95)))
hist_pts <- data.frame(vol=vol_h, se_raw=(cells[,1]+cells[,2])/rowSums(cells[,1:3]))
gg1 <- ggplot(curve, aes(vol, med)) + geom_ribbon(aes(ymin=lo, ymax=hi), fill="#0072B2", alpha=0.2) +
  geom_line(colour="#0072B2", linewidth=1) +
  geom_point(data=hist_pts, aes(vol, se_raw), inherit.aes=FALSE, colour="grey30") +
  coord_cartesian(ylim=c(0,1)) +
  labs(x="blood volume (mL)", y="Se_BC", title="Blood-culture sensitivity vs volume (final model)",
       subtitle="line+band = model (90% CrI); points = historic raw BC detection") + theme_minimal(base_size=11)
ggsave("latent_class_ppv/figures/fig_final_se_volume.png", gg1, width=7, height=4.5, dpi=300)

ppv <- rbind(
  data.frame(study=H$study, ppv=phi$median, lo=phi$q5, hi=phi$q95, kind="hospital (phi_s)"),
  data.frame(study=Vv$study, ppv=pi$median, lo=pi$q5, hi=pi$q95, kind="community (pi_o)"))
gg2 <- ggplot(ppv, aes(reorder(study,ppv), ppv, colour=kind)) +
  geom_pointrange(aes(ymin=lo, ymax=hi)) + coord_flip(ylim=c(0,1)) +
  scale_colour_manual(values=c("hospital (phi_s)"="#D55E00","community (pi_o)"="#0072B2"), name=NULL) +
  labs(x=NULL, y="PPV of the clinical case definition", title="PPV: hospital (phi_s) vs community (pi_o)",
       subtitle="hospital PPVs generally exceed community; both are local, not pooled") + theme_minimal(base_size=11)
ggsave("latent_class_ppv/figures/fig_final_ppv.png", gg2, width=7.5, height=5.5, dpi=300)

parameter_names <- c("Se_BC_2mL","Se_BC_5mL","Se_BC_10mL","Se_BM","alpha0","alpha1","tau","Sp_BC",
                     "mu_pi","sigma_pi","pi_population_median","mu_delta","beta_bias","sigma_delta")
est <- rbind(data.frame(param=parameter_names, t(sapply(parameter_names, qs))))
write.csv(est, "latent_class_ppv/tables/final_parameters.csv", row.names = FALSE)
write.csv(data.frame(study=H$study, vol_mL=round(vol_h,1), phi_med=phi$median, phi_lo=phi$q5, phi_hi=phi$q95),
          "latent_class_ppv/tables/final_phi_hospital.csv", row.names = FALSE)
write.csv(data.frame(study=Vv$study, suspected=Vv$suspected, tested=Vv$tested,
                     conf=Vv$confirmed_cx, positivity=Vv$positivity, bias_score=Vv$bias_score,
                     pi_med=pi$median, pi_lo=pi$q5, pi_hi=pi$q95),
          "latent_class_ppv/tables/final_pi_community.csv", row.names = FALSE)
delta <- fit$summary("delta"); q0 <- fit$summary("q0"); q1 <- fit$summary("q1")
write.csv(data.frame(
  study=Vv$study,
  delta_med=delta$median, delta_lo=delta$q5, delta_hi=delta$q95,
  q0_med=q0$median, q0_lo=q0$q5, q0_hi=q0$q95,
  q1_med=q1$median, q1_lo=q1$q5, q1_hi=q1$q95),
  "latent_class_ppv/tables/final_selection_outbreak.csv", row.names=FALSE)
saveRDS(fit$draws(), "latent_class_ppv/outputs/final_draws.rds")
cat("\n=== final fit complete ===\n")
