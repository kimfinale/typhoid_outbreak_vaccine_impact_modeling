# =============================================================================
# SENSITIVITY: does case-definition STRICTNESS explain community-PPV variation?
#   Rscript revision/run_strictness_sensitivity.R
#
# Self-contained; does NOT modify model_final.stan or run_final.R. Fits the
# strictness-covariate variant on the 8 community anchors, reports the strictness
# slope gamma, the residual sigma_pi (vs the base 1.35), and predicts PPV for the
# 9 unanchored renewal outbreaks under pooled vs strictness-conditioned draws.
#
# Strictness is coded on four transparent axes from each suspected-case
# definition (higher = more specific / more likely true typhoid):
#   fever_dur : fever-duration threshold      0 none | 1 >=3d | 2 >=7d
#   n_sym     : required additional symptoms   0..3 (capped)
#   exclude   : "malaria-negative"/"no other cause" clause   0/1
#   gate      : severity or epi-link gate      0 | 1 partial | 2 contact-of-confirmed | 3 peritonitis
#   strictness = fever_dur + n_sym + exclude + gate
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({ library(cmdstanr); library(posterior) })
outdir <- "revision/tables/strictness_sensitivity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
il <- function(x) 1/(1+exp(-x))

# ---- transparent strictness coding (fever_dur, n_sym, exclude, gate) ---------
code <- function(fever_dur, n_sym, exclude, gate) fever_dur + n_sym + exclude + gate
strict_tbl <- rbind(
  # 8 community anchors (fit set)
  data.frame(study="Dhadwal 2008",          set="anchor", fever_dur=0,n_sym=0,exclude=0,gate=0),
  data.frame(study="Makungo 2020",          set="anchor", fever_dur=0,n_sym=1,exclude=0,gate=0),
  data.frame(study="N'Cho 2019",            set="anchor", fever_dur=0,n_sym=1,exclude=0,gate=0),
  data.frame(study="Neil 2012",             set="anchor", fever_dur=0,n_sym=1,exclude=0,gate=1),
  data.frame(study="Aye 2004",              set="anchor", fever_dur=2,n_sym=1,exclude=0,gate=0),
  data.frame(study="Kabwama 2017",          set="anchor", fever_dur=1,n_sym=2,exclude=1,gate=0),
  data.frame(study="Hu 2022",               set="anchor", fever_dur=1,n_sym=2,exclude=1,gate=0),
  data.frame(study="Swaddiwudhipong 2001",  set="anchor", fever_dur=1,n_sym=3,exclude=1,gate=0),
  # 9 unanchored renewal outbreaks (predict only)
  data.frame(study="Lewis 2005",            set="unanchored", fever_dur=0,n_sym=0,exclude=0,gate=0),
  data.frame(study="Polonsky 2014, Dzivaresekwa", set="unanchored", fever_dur=0,n_sym=1,exclude=0,gate=0),
  data.frame(study="Polonsky 2014_Kuwadzana",     set="unanchored", fever_dur=0,n_sym=1,exclude=0,gate=0),
  data.frame(study="Muti 2014",             set="unanchored", fever_dur=1,n_sym=1,exclude=0,gate=0),
  data.frame(study="Qamar 2018",            set="unanchored", fever_dur=1,n_sym=0,exclude=1,gate=0),
  data.frame(study="Davis 2018",            set="unanchored", fever_dur=0,n_sym=3,exclude=0,gate=0),
  data.frame(study="Ali 2017",              set="unanchored", fever_dur=1,n_sym=1,exclude=1,gate=0),
  data.frame(study="Muyembe-Tamfum 2009",   set="unanchored", fever_dur=0,n_sym=0,exclude=0,gate=3),
  data.frame(study="Yousafzai 2019",        set="unanchored", fever_dur=1,n_sym=0,exclude=1,gate=2))
strict_tbl$strictness <- with(strict_tbl, code(fever_dur,n_sym,exclude,gate))
write.csv(strict_tbl, "revision/case_definition_strictness.csv", row.names=FALSE)

# ---- historic paired data (identical to run_final.R) ------------------------
H <- read.csv("latent_class_ppv/data/mogasale2016_paired_bc_bmc_seed.csv", check.names=FALSE, stringsAsFactors=FALSE)
cells <- as.matrix(H[, c("a_BCpos_BMpos","b_BCpos_BMneg","c_BCneg_BMpos","d_BCneg_BMneg")])
vol_h <- H$blood_vol_mL; vol_h[is.na(vol_h)] <- median(vol_h, na.rm=TRUE)

# ---- community anchors (identical build to run_final.R) ----------------------
D <- read.csv("latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv", check.names=FALSE, stringsAsFactors=FALSE)
inc <- D$include_for_main_ppv_model %in% c(TRUE,"TRUE","True")
frame_n <- suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*$","\\1", D$audited_suspected_testing_frame)))
Vv <- data.frame(study=D$study[inc], suspected=frame_n[inc],
                 tested=suppressWarnings(as.numeric(D$audited_blood_tested[inc])),
                 confirmed_cx=suppressWarnings(as.numeric(D$audited_blood_positive[inc])),
                 grade=D$verification_adjustment_grade[inc], stringsAsFactors=FALSE)
Vv$bias_score <- unname(c(`negligible/low`=0,moderate=1,major=2)[Vv$grade])
Vv$logit_tested_fraction <- qlogis((Vv$tested+0.5)/(Vv$suspected+1))
# strictness for the 8 anchors, centered/scaled on the anchor set
sa <- strict_tbl[strict_tbl$set=="anchor", ]
Vv$strictness_raw <- sa$strictness[match(Vv$study, sa$study)]
str_mean <- mean(Vv$strictness_raw); str_sd <- sd(Vv$strictness_raw)
Vv$strictness_z <- (Vv$strictness_raw - str_mean)/str_sd
stopifnot(!anyNA(Vv$strictness_raw))

standata <- list(
  S=nrow(H), N_hist=as.integer(rowSums(cells)), cells=matrix(as.integer(cells),ncol=4),
  logv_hist=log(vol_h), O=nrow(Vv), N_suspected=as.integer(Vv$suspected),
  n_out=as.integer(Vv$tested), k_out=as.integer(Vv$confirmed_cx),
  logv_out=log(rep(5, nrow(Vv))), bias_score=Vv$bias_score, strictness=Vv$strictness_z,
  logit_tested_fraction=Vv$logit_tested_fraction,
  mu_delta_sd=0.5, beta_bias_mu=0.5, beta_bias_sd=0.4, sigma_delta_sd=0.5,
  se_bm_a=12, se_bm_b=2, sp_a=200, sp_b=1, alpha1_mu=0.36, alpha1_sd=0.15, gamma_sd=1.0)

mod <- cmdstan_model("latent_class_ppv/stan/model_strictness_sensitivity.stan")
fit <- mod$sample(data=standata, seed=2026, chains=4, parallel_chains=4,
                  iter_warmup=1500, iter_sampling=1500, adapt_delta=0.98,
                  max_treedepth=12, init=0.5, refresh=0, show_messages=FALSE)
dg <- fit$diagnostic_summary()
cat(sprintf("Diagnostics: divergences=%d  max Rhat=%.3f\n",
            sum(dg$num_divergent), max(fit$summary()$rhat, na.rm=TRUE)))

dm <- as_draws_matrix(fit$draws(c("gamma","sigma_pi","mu_pi","pi_population_median")))
q <- function(v) c(med=median(v), lo=quantile(v,.05), hi=quantile(v,.95))
g <- as.numeric(dm[,"gamma"]); sp <- as.numeric(dm[,"sigma_pi"]); mp <- as.numeric(dm[,"mu_pi"])
cat("\n=== STRICTNESS COVARIATE FIT (8 anchors) ===\n")
cat(sprintf("gamma (logit PPV per +1 SD strictness): %.2f [%.2f, %.2f]\n", q(g)[1],q(g)[2],q(g)[3]))
cat(sprintf("  prior N(0, 1);  posterior mass P(gamma>0) = %.2f\n", mean(g>0)))
cat(sprintf("sigma_pi (residual): %.2f [%.2f, %.2f]   (base model without covariate: 1.35)\n", q(sp)[1],q(sp)[2],q(sp)[3]))
cat(sprintf("pi at mean strictness (pop median): %.2f [%.2f, %.2f]\n",
            q(il(mp))[1], q(il(mp))[2], q(il(mp))[3]))

# per-anchor fitted pi
pia <- fit$summary("pi")
cat("\n=== per-anchor fitted PPV (with strictness covariate) ===\n")
print(data.frame(study=Vv$study, strictness=Vv$strictness_raw,
                 pi_med=round(pia$median,2), lo=round(pia$q5,2), hi=round(pia$q95,2)), row.names=FALSE)

# transportability: predicted PPV for the 9 unanchored outbreaks,
# pooled (strictness held at anchor mean) vs strictness-conditioned
un <- strict_tbl[strict_tbl$set=="unanchored", ]
un$str_z <- (un$strictness - str_mean)/str_sd
pred <- function(str_z) {
  # marginal over z_pi ~ N(0, sigma_pi): median = inv_logit(mu_pi + gamma*str_z)
  sapply(str_z, function(sz) median(il(mp + g*sz)))
}
un$pi_pooled <- pred(0)                 # ignoring strictness (anchor-mean)
un$pi_conditioned <- pred(un$str_z)     # conditioned on own definition strictness
un$ratio <- round(un$pi_conditioned/un$pi_pooled, 2)
cat("\n=== unanchored renewal outbreaks: pooled vs strictness-conditioned PPV ===\n")
print(data.frame(study=un$study, strictness=un$strictness,
                 pi_pooled=round(un$pi_pooled,3), pi_conditioned=round(un$pi_conditioned,3),
                 ratio=un$ratio), row.names=FALSE)

write.csv(data.frame(study=Vv$study, strictness=Vv$strictness_raw, pi_med=pia$median,
                     pi_lo=pia$q5, pi_hi=pia$q95),
          file.path(outdir,"anchor_pi_with_strictness.csv"), row.names=FALSE)
write.csv(un[,c("study","strictness","pi_pooled","pi_conditioned","ratio")],
          file.path(outdir,"unanchored_transportability.csv"), row.names=FALSE)
write.csv(data.frame(gamma_med=q(g)[1], gamma_lo=q(g)[2], gamma_hi=q(g)[3], p_gamma_gt0=mean(g>0),
                     sigma_pi_med=q(sp)[1], sigma_pi_base=1.35,
                     pop_pi_med=q(il(mp))[1]),
          file.path(outdir,"strictness_summary.csv"), row.names=FALSE)
cat("\nWrote strictness coding + sensitivity outputs.\n")
