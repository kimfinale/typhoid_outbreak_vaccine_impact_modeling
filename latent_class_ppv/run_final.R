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
OUT_VOL <- 5    # assumed community-outbreak blood volume (mL); swept below

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

# ---- 2. consolidate COMMUNITY OUTBREAK positivity (tight syndromic set) -------
D <- read.csv("latent_class_ppv/data/community_surveillance_ppv.csv", check.names = FALSE, stringsAsFactors = FALSE)
Vv <- D[D$community_syndromic %in% c(TRUE,"TRUE","True"), ]
cat("\n=== Community outbreak points (tight; assumed vol =", OUT_VOL, "mL) ===\n")
print(data.frame(study = Vv$study, tested = Vv$tested, conf = Vv$confirmed_cx,
                 positivity = round(Vv$positivity,3)), row.names = FALSE)

standata <- function(out_vol) list(
  S = nrow(H), N_hist = as.integer(rowSums(cells)), cells = matrix(as.integer(cells), ncol = 4),
  logv_hist = log(vol_h),
  O = nrow(Vv), n_out = as.integer(Vv$tested), k_out = as.integer(Vv$confirmed_cx),
  logv_out = rep(log(out_vol), nrow(Vv)),
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
cat(sprintf("Diagnostics: divergences=%d, max Rhat=%.3f, min ESS=%.0f\n",
            sum(dg$num_divergent),
            max(fit$summary()$rhat, na.rm = TRUE), min(fit$summary()$ess_bulk, na.rm = TRUE)))
qs <- function(nm){ v<-as.numeric(fit$draws(nm, format="draws_matrix")); c(med=median(v), lo=quantile(v,.05), hi=quantile(v,.95)) }

# ---- 4. inferences ----------------------------------------------------------
cat("\n=== KEY ESTIMATES (median [90% CrI]) ===\n")
for (nm in c("Se_BC_2mL","Se_BC_5mL","Se_BC_10mL","Se_BM","alpha1","tau","Sp_BC")) {
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
for (vv in c(3,5,8,10)) {
  f2 <- mod$sample(data = standata(vv), seed=2026, chains=2, iter_warmup=800, iter_sampling=800,
                   adapt_delta=0.98, init=0.5, refresh=0, show_messages=FALSE)
  se <- median(as.numeric(f2$draws(sprintf("Se_BC_%dmL", ifelse(vv%in%c(2,5,10), vv, 5)), format="draws_matrix")))
  pim <- f2$summary("pi")$median
  cat(sprintf("%-8d %-10.2f %s\n", vv, il(median(as.numeric(f2$draws("alpha0",format="draws_matrix")))+
      median(as.numeric(f2$draws("alpha1",format="draws_matrix")))*log(vv)),
      paste(sprintf("%s=%.2f", Vv$study, pim), collapse=", ")))
}

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

est <- rbind(
  data.frame(param=c("Se_BC_2mL","Se_BC_5mL","Se_BC_10mL","Se_BM","alpha0","alpha1","tau","Sp_BC"),
             t(sapply(c("Se_BC_2mL","Se_BC_5mL","Se_BC_10mL","Se_BM","alpha0","alpha1","tau","Sp_BC"), qs))))
write.csv(est, "latent_class_ppv/tables/final_parameters.csv", row.names = FALSE)
write.csv(data.frame(study=H$study, vol_mL=round(vol_h,1), phi_med=phi$median, phi_lo=phi$q5, phi_hi=phi$q95),
          "latent_class_ppv/tables/final_phi_hospital.csv", row.names = FALSE)
write.csv(data.frame(study=Vv$study, tested=Vv$tested, conf=Vv$confirmed_cx, positivity=Vv$positivity,
                     pi_med=pi$median, pi_lo=pi$q5, pi_hi=pi$q95),
          "latent_class_ppv/tables/final_pi_community.csv", row.names = FALSE)
saveRDS(fit$draws(), "latent_class_ppv/outputs/final_draws.rds")
cat("\n=== final fit complete ===\n")
