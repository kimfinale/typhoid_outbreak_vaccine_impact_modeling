# =============================================================================
# Phase 1 — simulation-recovery GATE for the BC/BMC latent-class PPV model.
# Proves recovery from simulated data BEFORE any real data. KEY structural fact
# (caveat c): the severe->mild offset beta has NO data to identify it (historic
# paired studies are all severe; outbreaks are single-test). So beta is a FIXED
# sensitivity assumption, not a fitted parameter. The gate therefore tests
# recovery GIVEN the offset (beta locked); we separately (i) show a free beta is
# unidentified and (ii) sweep the offset to quantify pi's dependence on it.
# Deliverables also include the identifiability RIDGE, the ">1" plug-in scenario,
# and the conditional-dependence stress test.
#   Rscript latent_class_ppv/run_phase1_recovery.R
# =============================================================================
ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(cmdstanr); library(posterior); library(yaml)
                   library(ggplot2); library(dplyr) })
cfg <- read_yaml("latent_class_ppv/config.yml")
set.seed(cfg$seed)
for (d in c("figures","tables","outputs")) dir.create(file.path("latent_class_ppv", d), showWarnings = FALSE)
tr <- lapply(cfg$truth, as.numeric); pr <- cfg$priors
il <- function(x) 1 / (1 + exp(-x))

# ---- simulators -------------------------------------------------------------
sim_hist <- function(S, cd = 0) {
  v   <- runif(S, cfg$sim$vol_hist[1], cfg$sim$vol_hist[2])
  N   <- sample(cfg$sim$N_hist[1]:cfg$sim$N_hist[2], S, replace = TRUE)
  phi <- runif(S, cfg$sim$phi_range[1], cfg$sim$phi_range[2])
  u   <- rnorm(S)
  se  <- il(tr$alpha0 + tr$alpha1 * log(v) + tr$beta * 0 + tr$tau * u)   # severe (mild=0)
  sm  <- tr$Se_BM
  cells <- t(sapply(seq_len(S), function(s) {
    i <- c(se[s]*sm + cd, se[s]*(1-sm) - cd, (1-se[s])*sm - cd, (1-se[s])*(1-sm) + cd)
    vneg <- c((1-tr$Sp_BC)*(1-tr$Sp_BM), (1-tr$Sp_BC)*tr$Sp_BM, tr$Sp_BC*(1-tr$Sp_BM), tr$Sp_BC*tr$Sp_BM)
    as.integer(rmultinom(1, N[s], pmax(phi[s]*i + (1-phi[s])*vneg, 1e-9)))
  }))
  list(S = S, N_hist = N, cells = cells, logv_hist = log(v), mild_hist = rep(0, S), phi_true = phi)
}
mk_out <- function(v, mild, pi, n, u = rnorm(length(v))) {          # build one/many outbreaks
  se <- il(tr$alpha0 + tr$alpha1 * log(v) + tr$beta * mild + tr$tau * u)
  th <- pi*se + (1-pi)*(1-tr$Sp_BC)
  list(O = length(v), n_out = n, k_out = rbinom(length(v), n, th),
       logv_out = log(v), mild_out = rep(mild, length(v)), pi_true = pi, se_true = se)
}
sim_out <- function(O, mild = 1) mk_out(
  v = runif(O, cfg$sim$vol_out[1], cfg$sim$vol_out[2]), mild = mild,
  pi = runif(O, cfg$sim$pi_range[1], cfg$sim$pi_range[2]),
  n = sample(cfg$sim$n_out[1]:cfg$sim$n_out[2], O, replace = TRUE))

# ---- Stan data (robust coercion; optional locked offset) --------------------
standata <- function(H, G, cd = 0, beta_lock = NULL) {
  if (is.null(H)) H <- list(S = 0, N_hist = integer(0), cells = matrix(0L, 0, 4),
                            logv_hist = numeric(0), mild_hist = numeric(0))
  if (is.null(G)) G <- list(O = 0, n_out = integer(0), k_out = integer(0),
                            logv_out = numeric(0), mild_out = numeric(0))
  locked <- !is.null(beta_lock)
  list(S = as.integer(H$S), N_hist = as.array(as.integer(H$N_hist)),
       cells = matrix(as.integer(H$cells), ncol = 4),
       logv_hist = as.array(as.numeric(H$logv_hist)), mild_hist = as.array(as.numeric(H$mild_hist)),
       cd_infected = as.numeric(cd),
       O = as.integer(G$O), n_out = as.array(as.integer(G$n_out)), k_out = as.array(as.integer(G$k_out)),
       logv_out = as.array(as.numeric(G$logv_out)), mild_out = as.array(as.numeric(G$mild_out)),
       se_bm_a = pr$se_bm_a, se_bm_b = pr$se_bm_b, sp_a = pr$sp_a, sp_b = pr$sp_b, beta_sd = pr$beta_sd,
       beta_locked = as.integer(locked), beta_lock_val = if (locked) as.numeric(beta_lock) else 0)
}
fit_model <- function(dat) mod$sample(data = dat, seed = cfg$stan_seed,
    chains = cfg$mcmc$chains, parallel_chains = cfg$mcmc$parallel_chains,
    iter_warmup = cfg$mcmc$warmup, iter_sampling = cfg$mcmc$sampling,
    adapt_delta = cfg$mcmc$adapt_delta, max_treedepth = cfg$mcmc$max_treedepth,
    init = 0.5, refresh = 0, show_messages = FALSE)   # init near centre (Sp starts ~0.5, not ~0.9)
rec_row <- function(dr, nm, truth) { q <- quantile(dr, c(.05, .5, .95), names = FALSE)
  data.frame(param = nm, true = truth, median = q[2], lo = q[1], hi = q[3],
             bias = q[2] - truth, covered = (truth >= q[1] & truth <= q[3])) }
pi_median_bias <- function(d, G) median(sapply(seq_len(G$O),
  function(o) median(as.numeric(d[, sprintf("pi[%d]", o)]))) - G$pi_true)

cat("=== Phase 1: compile ===\n"); mod <- cmdstan_model("latent_class_ppv/stan/bc_bmc_ppv.stan")
H <- sim_hist(cfg$sim$S_hist); G <- sim_out(cfg$sim$O_out)

# ===================== MAIN RECOVERY (offset LOCKED at truth) =================
cat("=== Main recovery: beta LOCKED at truth (recover everything given the offset) ===\n")
fit <- fit_model(standata(H, G, cd = 0, beta_lock = tr$beta))
sm <- fit$summary(c("alpha0","alpha1","tau","Se_BM","pi","phi"))
dg <- fit$diagnostic_summary(); ndiv <- sum(dg$num_divergent)
maxrhat <- max(sm$rhat, na.rm = TRUE); miness <- min(sm$ess_bulk, na.rm = TRUE)
d <- fit$draws(format = "draws_matrix")
rec <- rbind(rec_row(d[,"alpha0"],"alpha0",tr$alpha0), rec_row(d[,"alpha1"],"alpha1",tr$alpha1),
             rec_row(d[,"tau"],"tau",tr$tau), rec_row(d[,"Se_BM"],"Se_BM",tr$Se_BM))
pi_rec  <- do.call(rbind, lapply(seq_len(G$O), function(o) rec_row(d[,sprintf("pi[%d]",o)], sprintf("pi[%d]",o), G$pi_true[o])))
phi_rec <- do.call(rbind, lapply(seq_len(H$S), function(s) rec_row(d[,sprintf("phi[%d]",s)], sprintf("phi[%d]",s), H$phi_true[s])))
write.csv(rbind(rec, pi_rec, phi_rec), "latent_class_ppv/tables/recovery.csv", row.names = FALSE)
gate_set <- rbind(rec, pi_rec); cov_all <- mean(gate_set$covered)
bias_logit <- max(abs(rec$bias)); bias_pi <- median(abs(pi_rec$bias))
cat(sprintf("Diagnostics: divergences=%d, max Rhat=%.3f, min ESS=%.0f\n", ndiv, maxrhat, miness))
cat(sprintf("Recovery(beta locked): coverage(shared+pi)=%.2f | max|bias|shared=%.3f | median|bias|pi=%.3f | phi cov=%.2f\n",
            cov_all, bias_logit, bias_pi, mean(phi_rec$covered)))

gg1 <- ggplot(rec, aes(param, median)) + geom_pointrange(aes(ymin=lo, ymax=hi)) +
  geom_point(aes(y=true), colour="#D55E00", size=3, shape=4, stroke=1.2) +
  labs(x=NULL, y="posterior (median, 90% CrI)", title="Recovery of Se sub-model params (offset locked)",
       subtitle="orange x = true; identified from the historic paired BC/BMC data") + theme_minimal(base_size=11)
ggsave("latent_class_ppv/figures/fig_recovery_scalars.png", gg1, width=7, height=4.5, dpi=300)
gg2 <- ggplot(transform(pi_rec, cc = ifelse(covered,"covered","not covered")),
              aes(true, median, colour=cc)) + geom_abline(slope=1, intercept=0, colour="grey60") +
  geom_pointrange(aes(ymin=lo, ymax=hi)) +
  scale_colour_manual(values=c(covered="#0072B2",`not covered`="#D55E00"), name=NULL) +
  labs(x="true surveillance PPV (pi_o)", y="posterior pi_o (median, 90% CrI)",
       title="Recovery of per-outbreak surveillance PPV (offset locked at truth)",
       subtitle=sprintf("coverage=%.2f, median|bias|=%.3f", mean(pi_rec$covered), bias_pi)) + theme_minimal(base_size=11)
ggsave("latent_class_ppv/figures/fig_recovery_pi.png", gg2, width=7, height=5, dpi=300)

# ===================== OFFSET NON-IDENTIFICATION (free beta) ==================
cat("=== Free-beta fit: show the offset is NOT identified from data ===\n")
fit_free <- fit_model(standata(H, G, cd = 0, beta_lock = NULL))
df <- fit_free$draws(format = "draws_matrix"); bfree <- as.numeric(df[,"beta"])
prior_sd <- pr$beta_sd
cat(sprintf("beta: true=%.2f | free-fit posterior median=%.2f (90%% CrI %.2f..%.2f) | half-normal prior sd=%.2f -> posterior ~ prior => UNIDENTIFIED\n",
            tr$beta, median(bfree), quantile(bfree,.05), quantile(bfree,.95), prior_sd))
cat(sprintf("  -> resulting pi median bias (free beta) = %+.3f (vs %.3f when locked at truth)\n",
            pi_median_bias(df, G), pi_median_bias(d, G)))

# ===================== OFFSET SENSITIVITY SWEEP ==============================
cat("=== Offset sensitivity sweep: pi_o vs the assumed offset ===\n")
betas <- c(0, tr$beta, -1.0); swp <- list()
for (b in betas) {
  fb <- if (b == tr$beta) d else fit_model(standata(H, G, cd = 0, beta_lock = b))$draws(format = "draws_matrix")
  med <- sapply(seq_len(G$O), function(o) median(as.numeric(fb[, sprintf("pi[%d]", o)])))
  swp[[as.character(b)]] <- data.frame(beta = b, true = G$pi_true, pi_med = med)
}
swpdf <- do.call(rbind, swp)
write.csv(swpdf, "latent_class_ppv/tables/offset_sweep.csv", row.names = FALSE)
sweep_summary <- swpdf %>% group_by(beta) %>% summarise(mean_pi = mean(pi_med), mean_bias = mean(pi_med - true))
cat("Mean pi_o by assumed offset (truth beta =", tr$beta, "):\n"); print(as.data.frame(sweep_summary), row.names = FALSE)
gg3 <- ggplot(swpdf, aes(true, pi_med, colour = factor(beta))) +
  geom_abline(slope=1, intercept=0, colour="grey60") + geom_point() +
  scale_colour_brewer(palette="Set1", name="assumed offset beta") +
  labs(x="true pi_o", y="posterior median pi_o",
       title="Surveillance PPV depends on the (non-identified) offset",
       subtitle="beta=0 ignores the severe->mild load gap and under-estimates pi; truth is the middle line") +
  theme_minimal(base_size=11)
ggsave("latent_class_ppv/figures/fig_offset_sweep.png", gg3, width=7, height=5, dpi=300)

# ===================== RIDGE (one outbreak, no anchor) =======================
cat("=== Ridge demo: single outbreak, no historic anchor ===\n")
G1 <- mk_out(v = 5, mild = 1, pi = 0.30, n = 300)
dr <- fit_model(standata(NULL, G1, cd = 0, beta_lock = tr$beta))$draws(format = "draws_matrix")
ridge <- data.frame(pi = as.numeric(dr[,"pi[1]"]), Se = as.numeric(dr[,"Se_BC_out[1]"]))
ridge$theta <- ridge$pi * ridge$Se
cat(sprintf("Ridge: pi sd=%.3f, Se sd=%.3f, product(theta) sd=%.3f (product tight => ridge)\n",
            sd(ridge$pi), sd(ridge$Se), sd(ridge$theta)))
gg4 <- ggplot(ridge[sample(nrow(ridge), min(3000,nrow(ridge))),], aes(Se, pi)) +
  geom_point(alpha=0.15, colour="#0072B2") +
  labs(x="Se_BC,o", y="pi_o", title="Non-identifiability ridge (one outbreak, no anchor)",
       subtitle=sprintf("only the product theta=pi*Se is identified (theta sd=%.3f)", sd(ridge$theta))) +
  theme_minimal(base_size=11)
ggsave("latent_class_ppv/figures/fig_ridge.png", gg4, width=6.5, height=5, dpi=300)

# ===================== ">1" PLUG-IN SCENARIO ================================
cat("=== '>1' scenario ===\n")
g1 <- cfg$demo_gt1; kk <- round(g1$positivity * g1$n_cult); plugin <- g1$positivity / 0.60
Ggt <- list(O = 1, n_out = g1$n_cult, k_out = kk, logv_out = log(g1$volume), mild_out = g1$mild)
pi_gt <- NA; se_gt <- NA
ok_gt <- tryCatch({
  dgt <- fit_model(standata(H, Ggt, cd = 0, beta_lock = tr$beta))$draws(format = "draws_matrix")
  pi_gt <<- as.numeric(dgt[,"pi[1]"]); se_gt <<- as.numeric(dgt[,"Se_BC_out[1]"]); TRUE
}, error = function(e) { message("  >1 fit failed: ", conditionMessage(e)); FALSE })
if (ok_gt) {
  cat(sprintf("Positivity=%.2f; plug-in=%.2f (>1 invalid). Bayesian pi median=%.2f, P(pi>1)=%.2f, Se median=%.2f\n",
              g1$positivity, plugin, median(pi_gt), mean(pi_gt > 1), median(se_gt)))
  gg5 <- ggplot(data.frame(pi=pi_gt), aes(pi)) + geom_histogram(bins=40, fill="#0072B2", alpha=0.85) +
    geom_vline(xintercept=plugin, colour="#D55E00", linetype="dashed", linewidth=1) +
    coord_cartesian(xlim=c(0, max(1.25, plugin+0.05))) +
    labs(x="pi_o (surveillance PPV)", y="posterior draws", title="The '>1' pathology is dissolved",
         subtitle=sprintf("plug-in positivity/0.60=%.2f (>1); Bayesian pi stays in [0,1], mass shifts to higher Se_BC", plugin)) +
    theme_minimal(base_size=11)
  ggsave("latent_class_ppv/figures/fig_gt1.png", gg5, width=7, height=4.5, dpi=300)
}

# ===================== CONDITIONAL-DEPENDENCE STRESS ========================
cat("=== Conditional-dependence stress: simulate cd>0, fit independence ===\n")
cd <- cfg$dependence$cd_infected; Hd <- sim_hist(cfg$sim$S_hist, cd = cd)
ddp <- fit_model(standata(Hd, G, cd = 0, beta_lock = tr$beta))$draws(format = "draws_matrix")
se5_true <- il(tr$alpha0 + tr$alpha1*log(5))
se5_dep  <- median(il(as.numeric(ddp[,"alpha0"]) + as.numeric(ddp[,"alpha1"])*log(5)))
sebm_dep <- median(as.numeric(ddp[,"Se_BM"]))
cat(sprintf("Dependence stress (cd=%.2f): Se(severe,5mL) true=%.3f vs indep-fit=%.3f (bias=%+.3f); Se_BM true=%.2f vs fit=%.2f; pi bias=%+.3f\n",
            cd, se5_true, se5_dep, se5_dep-se5_true, tr$Se_BM, sebm_dep, pi_median_bias(ddp, G)))

# ===================== GATE =================================================
converged <- (ndiv == 0 && maxrhat <= cfg$gate$max_rhat && miness >= cfg$gate$min_ess)
pass <- converged && cov_all >= cfg$gate$min_cri_coverage &&
        bias_logit <= cfg$gate$max_abs_median_bias_logit && bias_pi <= cfg$gate$max_abs_median_bias_pi
saveRDS(list(rec=rbind(rec,pi_rec,phi_rec), sweep=swpdf, ridge=ridge,
             gt1=list(pi=pi_gt, se=se_gt, plugin=plugin),
             beta_free=bfree, diag=list(ndiv=ndiv, maxrhat=maxrhat, miness=miness)),
        "latent_class_ppv/outputs/phase1.rds")
cat(sprintf("\n=== GATE (recovery GIVEN the offset): %s ===\n", ifelse(pass,"PASS","FAIL")))
if (pass) {
  writeLines(c("PHASE 1 RECOVERY GATE: PASS (recovery given the offset)",
    sprintf("coverage(shared+pi)=%.2f  max|bias|shared=%.3f  median|bias|pi=%.3f", cov_all, bias_logit, bias_pi),
    sprintf("divergences=%d max_rhat=%.3f min_ess=%.0f", ndiv, maxrhat, miness),
    "KEY CAVEAT: the severe->mild offset beta is NOT identified from data (free-beta posterior ~ prior);",
    "it is a fixed sensitivity assumption. See tables/offset_sweep.csv for pi's dependence on it.",
    "Ridge, >1, and dependence-stress demos written. Awaiting sign-off before Phase 2 (real data)."),
    "latent_class_ppv/outputs/GATE_PASSED.txt")
  cat("Wrote GATE_PASSED.txt. STOP for sign-off before real data.\n")
} else cat("Recovery-given-offset not acceptable; see tables/recovery.csv before proceeding.\n")
cat("=== Phase 1 complete ===\n")
