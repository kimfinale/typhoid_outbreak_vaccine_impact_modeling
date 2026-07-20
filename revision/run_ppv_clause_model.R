# =============================================================================
# Exclusion-clause covariate model for community PPV (recommended sensitivity).
#   Rscript revision/run_ppv_clause_model.R
# Fits logit(pi_o) = mu_pi + delta*clause_o + sigma_pi z_o, delta>=0, alongside a
# like-for-like pooled model (use_covariate=0). Reports the clause effect, the
# two group PPVs, per-anchor and unanchored predictions, a prior-sensitivity
# sweep, and a leave-one-anchor-out (PSIS-LOO) comparison vs pooling.
# Self-contained; does not modify model_final.stan/run_final.R.
# =============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({ library(cmdstanr); library(posterior) })
outdir <- "revision/tables/clause_model"; dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
il <- function(x) 1/(1+exp(-x)); q <- function(v) c(med=median(v), lo=quantile(v,.05), hi=quantile(v,.95))

# ---- historic paired data (as run_final.R) ----------------------------------
H <- read.csv("latent_class_ppv/data/mogasale2016_paired_bc_bmc_seed.csv", check.names=FALSE, stringsAsFactors=FALSE)
cells <- as.matrix(H[, c("a_BCpos_BMpos","b_BCpos_BMneg","c_BCneg_BMpos","d_BCneg_BMneg")])
vol_h <- H$blood_vol_mL; vol_h[is.na(vol_h)] <- median(vol_h, na.rm=TRUE)

# ---- community anchors (as run_final.R) + clause covariate -------------------
D <- read.csv("latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv", check.names=FALSE, stringsAsFactors=FALSE)
inc <- D$include_for_main_ppv_model %in% c(TRUE,"TRUE","True")
frame_n <- suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*$","\\1", D$audited_suspected_testing_frame)))
Vv <- data.frame(study=D$study[inc], suspected=frame_n[inc],
                 tested=suppressWarnings(as.numeric(D$audited_blood_tested[inc])),
                 confirmed_cx=suppressWarnings(as.numeric(D$audited_blood_positive[inc])),
                 grade=D$verification_adjustment_grade[inc], stringsAsFactors=FALSE)
Vv$bias_score <- unname(c(`negligible/low`=0,moderate=1,major=2)[Vv$grade])
Vv$logit_tested_fraction <- qlogis((Vv$tested+0.5)/(Vv$suspected+1))
strict <- read.csv("revision/case_definition_strictness.csv", check.names=FALSE, stringsAsFactors=FALSE)
Vv$clause <- strict$exclude[match(Vv$study, strict$study)]      # objective exclusion-clause indicator
stopifnot(!anyNA(Vv$clause))
cat("Anchors and exclusion clause:\n")
print(data.frame(study=Vv$study, tested=Vv$tested, pos=Vv$confirmed_cx,
                 positivity=round(Vv$confirmed_cx/Vv$tested,2), clause=Vv$clause), row.names=FALSE)

standata <- function(delta_sd=0.5, use_cov=1L) list(
  S=nrow(H), N_hist=as.integer(rowSums(cells)), cells=matrix(as.integer(cells),ncol=4),
  logv_hist=log(vol_h), O=nrow(Vv), N_suspected=as.integer(Vv$suspected),
  n_out=as.integer(Vv$tested), k_out=as.integer(Vv$confirmed_cx),
  logv_out=log(rep(5,nrow(Vv))), bias_score=Vv$bias_score, clause=as.numeric(Vv$clause),
  logit_tested_fraction=Vv$logit_tested_fraction,
  mu_delta_sd=0.5, beta_bias_mu=0.5, beta_bias_sd=0.4, sigma_delta_sd=0.5,
  se_bm_a=12, se_bm_b=2, sp_a=200, sp_b=1, alpha1_mu=0.36, alpha1_sd=0.15,
  delta_sd=delta_sd, use_covariate=use_cov)

CSDIR <- file.path(Sys.getenv("TMP", tempdir()), "clause_cmdstan"); dir.create(CSDIR, showWarnings=FALSE)
mod <- cmdstan_model("latent_class_ppv/stan/model_ppv_clause.stan")
fit_one <- function(delta_sd, use_cov) mod$sample(data=standata(delta_sd,use_cov), seed=2026,
  chains=4, parallel_chains=4, iter_warmup=1200, iter_sampling=1200, adapt_delta=0.99,
  max_treedepth=12, init=0.5, refresh=0, show_messages=FALSE, output_dir=CSDIR)

cat("\n=== PRIMARY: clause covariate, delta ~ half-Normal(0, 0.5) ===\n")
fit <- fit_one(0.5, 1L)
dg <- fit$diagnostic_summary(); cat(sprintf("divergences=%d  max Rhat=%.3f\n", sum(dg$num_divergent), max(fit$summary()$rhat,na.rm=TRUE)))
dm <- as_draws_matrix(fit$draws(c("delta","sigma_pi","mu_pi","pi_pooled","pi_clause","pi_new_pooled","pi_new_clause")))
d <- as.numeric(dm[,"delta"])
cat(sprintf("delta (logit-PPV clause effect): %.2f [%.2f, %.2f]   P(delta>0)=%.2f (>0 by construction)\n",
            q(d)[1],q(d)[2],q(d)[3], mean(d>1e-6)))
cat(sprintf("clause effect as PPV odds ratio: %.2f [%.2f, %.2f]\n", q(exp(d))[1],q(exp(d))[2],q(exp(d))[3]))
cat(sprintf("sigma_pi (residual): %.2f [%.2f, %.2f]   (pooled base ~1.35)\n", q(as.numeric(dm[,"sigma_pi"]))[1],q(as.numeric(dm[,"sigma_pi"]))[2],q(as.numeric(dm[,"sigma_pi"]))[3]))
cat(sprintf("group population PPV: no-clause %.2f [%.2f, %.2f] | clause %.2f [%.2f, %.2f]\n",
            q(as.numeric(dm[,"pi_pooled"]))[1],q(as.numeric(dm[,"pi_pooled"]))[2],q(as.numeric(dm[,"pi_pooled"]))[3],
            q(as.numeric(dm[,"pi_clause"]))[1],q(as.numeric(dm[,"pi_clause"]))[2],q(as.numeric(dm[,"pi_clause"]))[3]))
cat(sprintf("posterior-predictive PPV for a NEW outbreak: no-clause %.2f [%.2f, %.2f] | clause %.2f [%.2f, %.2f]\n",
            q(as.numeric(dm[,"pi_new_pooled"]))[1],q(as.numeric(dm[,"pi_new_pooled"]))[2],q(as.numeric(dm[,"pi_new_pooled"]))[3],
            q(as.numeric(dm[,"pi_new_clause"]))[1],q(as.numeric(dm[,"pi_new_clause"]))[2],q(as.numeric(dm[,"pi_new_clause"]))[3]))

pia <- fit$summary("pi")
cat("\n-- per-anchor PPV (clause model) --\n")
print(data.frame(study=Vv$study, clause=Vv$clause, pi_med=round(pia$median,2), lo=round(pia$q5,2), hi=round(pia$q95,2)), row.names=FALSE)

# ---- leave-one-anchor-out (PSIS-LOO) vs pooled ------------------------------
cat("\n=== MODEL COMPARISON on the 8 community observations (WAIC-based ELPD) ===\n")
fit_pool <- fit_one(0.5, 0L)
elpd_pw <- function(f) { ll <- as_draws_matrix(f$draws("log_lik"))
  lppd <- apply(ll, 2, function(cc) log(mean(exp(cc)))); pw <- apply(ll, 2, var); lppd - pw }
ec <- elpd_pw(fit); ep <- elpd_pw(fit_pool); dif <- ec - ep
cat(sprintf("elpd_waic: clause %.2f | pooled %.2f\n", sum(ec), sum(ep)))
cat(sprintf("elpd_diff (clause - pooled) = %+.2f, SE %.2f  (positive favours clause; |diff| < 2*SE = indistinguishable)\n",
            sum(dif), sqrt(length(dif))*sd(dif)))
if (requireNamespace("loo", quietly=TRUE)) {
  cat("\nPSIS-LOO:\n")
  print(loo::loo_compare(list(clause=loo::loo(as_draws_matrix(fit$draws("log_lik"))),
                              pooled=loo::loo(as_draws_matrix(fit_pool$draws("log_lik"))))))
}

# ---- unanchored predictions by clause status --------------------------------
un <- strict[strict$set=="unanchored", c("study","exclude")]
un$group <- ifelse(un$exclude==1, "clause", "no-clause")
pp <- c(no_clause=median(as.numeric(dm[,"pi_pooled"])), clause=median(as.numeric(dm[,"pi_clause"])))
un$assigned_pi <- ifelse(un$exclude==1, pp["clause"], pp["no_clause"])
cat("\n=== assigned PPV for unanchored outbreaks (clause model, median) ===\n")
print(un[order(un$exclude), c("study","exclude","assigned_pi")], row.names=FALSE)
cat(sprintf("(pooled model would assign %.2f to all)\n", median(as.numeric(as_draws_matrix(fit_pool$draws("pi_pooled"))[,1]))))

write.csv(data.frame(study=Vv$study, clause=Vv$clause, positivity=Vv$confirmed_cx/Vv$tested,
                     pi_med=pia$median, pi_lo=pia$q5, pi_hi=pia$q95),
          file.path(outdir,"anchor_pi.csv"), row.names=FALSE)
cat("\nCore outputs written to", outdir, "\n")

# ---- prior sensitivity (guarded; secondary) ---------------------------------
cat("\n=== PRIOR SENSITIVITY on delta_sd ===\n")
tryCatch({
  psens <- do.call(rbind, lapply(c(0.25,0.5,1.0), function(sd_) {
    f <- fit_one(sd_, 1L); dd <- as.numeric(as_draws_matrix(f$draws("delta"))[,1])
    sp <- as.numeric(as_draws_matrix(f$draws("sigma_pi"))[,1])
    data.frame(delta_sd=sd_, delta_med=median(dd), delta_lo=quantile(dd,.05), delta_hi=quantile(dd,.95),
               OR_med=median(exp(dd)), sigma_pi_med=median(sp))
  }))
  print(psens, row.names=FALSE)
  write.csv(psens, file.path(outdir,"prior_sensitivity.csv"), row.names=FALSE)
}, error=function(e) cat("prior sweep failed (non-fatal):", conditionMessage(e), "\n"))
