# ============================================================================
# POC: ORI (outbreak-response immunization) impact on a DEVELOPING-COUNTRY
# typhoid outbreak, integrating the clinical-case-definition PPV (pi).
#
# Outbreak : Kabwama 2017 -- Kampala, Uganda 2015 (mixed, source-dominant
#            outbreak; DAILY SUSPECTED-case surveillance, ~10,152 suspected).
# PPV      : community pi = 0.23 [0.155, 0.418]  (fitted in latent_class_ppv;
#            same Kampala outbreak) -> only ~23% of suspected are true typhoid.
#
# Uses two additive levels: suspected = true typhoid + other febrile illness;
# true typhoid = background/common-source + propagated incidence. The former
# multiplicative PPV calculation is superseded by this demonstration.
#   Rscript manuscript/poc_ori_kabwama.R
# ============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
source("renewal/R/gi.R")
source("renewal/R/renewal_core.R")
source("renewal/R/ppv.R")

# --- Kabwama daily suspected -> weekly (7-day bins from first date) ----------
ts <- read.csv("data/Typhoid_Outbreak_Time_Series_2000_2022_Timeseries.csv", stringsAsFactors = FALSE)
names(ts)[1:6] <- c("study", "start", "end", "n_pat", "suspected", "confirmed")
k <- ts[trimws(ts$study) == "Kabwama 2017", ]
k$susp <- as.numeric(k$suspected); k$susp[is.na(k$susp)] <- as.numeric(k$n_pat)[is.na(k$susp)]
d0 <- as.Date(k$start)
wk <- as.integer(floor(as.numeric(d0 - min(d0)) / 7)) + 1L
weekly <- tapply(k$susp, wk, sum)
weekly <- as.numeric(weekly[order(as.integer(names(weekly)))])
Tn <- length(weekly)
cat(sprintf("Kabwama weekly: %d weeks, total suspected = %.0f, peak = %.0f (wk %d)\n\n",
            Tn, sum(weekly), max(weekly), which.max(weekly)))

# --- typhoid GI (weekly) -----------------------------------------------------
w <- discretize_gi(mean_days = 14, sd_days = 8.4, step_days = 7)
true_med <- ppv_true_incidence(weekly, 0.23)
Rt <- reconstruct_Rt(true_med, w)$Rt
cat("R_t (weekly, typhoid GI): peak =", sprintf("%.2f", max(Rt, na.rm = TRUE)),
    " at wk", which.max(Rt), "\n")

# --- ORI scenario grid: campaign timing x PPV -------------------------------
# detection = first week suspected>100 (clear outbreak signal); campaign starts
# tau weeks after start; immunity delay 2 wk -> t_eff. coverage 0.6, psi_T 0.8.
detect_wk <- which(weekly > 100)[1]
cat("outbreak-signal week (suspected>100):", detect_wk, "( peak wk", which.max(weekly), ")\n\n")
coverage <- 0.6; psi_T <- 0.8; immun_delay <- 2; theta_source <- 0.60
pi_ppv <- c(med = 0.23, lo = 0.155, hi = 0.418)

taus <- c(detect_wk, detect_wk + 2, detect_wk + 4, detect_wk + 6)
cat("=== ORI impact by campaign start (coverage=0.6, psi_T=0.8, +2wk immunity) ===\n")
cat(sprintf("%-16s %8s %12s %14s %14s\n", "campaign start", "t_eff", "true %reduc", "observed %", "avert_TRUE"))
rows <- list()
for (tau in taus) {
  te <- tau + immun_delay
  if (te > Tn) next
  add <- additive_source_counterfactual(
    true_med, w, tau = tau, t_eff = te, coverage = coverage,
    psi_T = psi_T, source_fraction = theta_source, c_shape = "step")
  im <- impact_measures(true_med, add$incidence_v)
  observed_pct <- 100 * im$averted / sum(weekly)
  cat(sprintf("wk %-13d %8d %11.1f%% %13.1f%% %14.0f\n",
              tau, te, im$pct_reduction, observed_pct, im$averted))
  rows[[length(rows) + 1]] <- data.frame(
    tau = tau, t_eff = te, true_pct = im$pct_reduction,
    observed_pct = observed_pct, avert_true = im$averted)
}
res <- do.call(rbind, rows)

# --- PPV sensitivity at the earliest campaign -------------------------------
tau1 <- taus[1]; te1 <- tau1 + immun_delay
add1 <- additive_source_counterfactual(
  true_med, w, tau = tau1, t_eff = te1, coverage = coverage,
  psi_T = psi_T, source_fraction = theta_source, c_shape = "step")
cat(sprintf("\n=== PPV sensitivity (earliest campaign wk %d; additive observation model) ===\n", tau1))
for (nm in names(pi_ppv)) {
  Tpi <- ppv_true_incidence(weekly, pi_ppv[nm])
  cf <- additive_source_counterfactual(
    Tpi, w, tau = tau1, t_eff = te1, coverage = coverage,
    psi_T = psi_T, source_fraction = theta_source, c_shape = "step")
  im <- impact_measures(Tpi, cf$incidence_v)
  cat(sprintf("  pi=%.3f (%s): %.0f true cases averted; %.1f%% true; %.1f%% observed\n",
              pi_ppv[nm], nm, im$averted, im$pct_reduction,
              100 * im$averted / sum(weekly)))
}
cat("\nInterpretation: PPV does not cancel under the additive observation model because\n",
    "other febrile illness is held fixed. Theta allocates source and propagated true\n",
    "typhoid inside one renewal recursion, rather than weighting completed outcomes.\n")

# --- plot -------------------------------------------------------------------
png("manuscript/data/poc_ori_kabwama.png", width = 1000, height = 450)
par(mfrow = c(1, 2))
plot(weekly, type = "h", lwd = 6, col = "#B0B0B0", xlab = "week", ylab = "suspected cases/week",
     main = "Kabwama 2017 (Kampala): ORI counterfactual")
lines(add1$incidence_v + (weekly - true_med), type = "h", lwd = 6, col = "#4C78A8")
abline(v = te1, lty = 2, col = "red"); legend("topright", c("observed suspected", "with ORI"),
       fill = c("#B0B0B0", "#4C78A8"), bty = "n")
plot(Rt, type = "b", pch = 19, xlab = "week", ylab = "R_t", main = "Reconstructed R_t (typhoid GI)")
abline(h = 1, lty = 2, col = "red")
dev.off()
cat("\nwrote manuscript/data/poc_ori_kabwama.png\n")
write.csv(res, "manuscript/data/poc_ori_kabwama_results.csv", row.names = FALSE)
