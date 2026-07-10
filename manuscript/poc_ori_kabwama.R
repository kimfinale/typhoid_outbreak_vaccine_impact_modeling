# ============================================================================
# POC: ORI (outbreak-response immunization) impact on a DEVELOPING-COUNTRY
# typhoid outbreak, integrating the clinical-case-definition PPV (pi).
#
# Outbreak : Kabwama 2017 -- Kampala, Uganda 2015 (large urban propagated
#            outbreak; DAILY SUSPECTED-case surveillance, ~10,152 suspected).
# PPV      : community pi = 0.23 [0.155, 0.418]  (fitted in latent_class_ppv;
#            same Kampala outbreak) -> only ~23% of suspected are true typhoid.
#
# Reuses the renewal ORI engine (renewal/R/renewal_core.R): reconstruct R_t from
# the observed curve, then a vaccination counterfactual with transmission
# feedback (emergent herd effect). Averted counts come out in SUSPECTED units;
# TRUE typhoid averted = pi * averted_suspected (constant-PPV multiplicative
# regime -> %reduction is pi-invariant, but absolute true burden scales by pi).
#   Rscript manuscript/poc_ori_kabwama.R
# ============================================================================
setwd(Sys.getenv("RENEWAL_ROOT", "."))
source("renewal/R/gi.R")
source("renewal/R/renewal_core.R")

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

# --- typhoid GI (weekly) + R_t reconstruction -------------------------------
w <- discretize_gi(mean_days = 14, sd_days = 8.4, step_days = 7)
Rt <- reconstruct_Rt(weekly, w)$Rt
cat("R_t (weekly, typhoid GI): peak =", sprintf("%.2f", max(Rt, na.rm = TRUE)),
    " at wk", which.max(Rt), "\n")

# --- ORI scenario grid: campaign timing x PPV -------------------------------
# detection = first week suspected>100 (clear outbreak signal); campaign starts
# tau weeks after start; immunity delay 2 wk -> t_eff. coverage 0.6, psi_T 0.8.
detect_wk <- which(weekly > 100)[1]
cat("outbreak-signal week (suspected>100):", detect_wk, "( peak wk", which.max(weekly), ")\n\n")
coverage <- 0.6; psi_T <- 0.8; immun_delay <- 2
pi_ppv <- c(med = 0.23, lo = 0.155, hi = 0.418)

taus <- c(detect_wk, detect_wk + 2, detect_wk + 4, detect_wk + 6)
cat("=== ORI impact by campaign start (coverage=0.6, psi_T=0.8, +2wk immunity) ===\n")
cat(sprintf("%-16s %8s %10s %12s %14s\n", "campaign start", "t_eff", "%reduc", "avert_susp", "avert_TRUE(pi=.23)"))
rows <- list()
for (tau in taus) {
  te <- tau + immun_delay
  if (te > Tn) next
  ren <- renewal_counterfactual(weekly, w, tau = tau, t_eff = te, pi = coverage,
                                psi_T = psi_T, c_shape = "step", feedback = TRUE)
  im <- impact_measures(weekly, ren$incidence_v)
  avert_true <- pi_ppv["med"] * im$averted
  cat(sprintf("wk %-13d %8d %9.1f%% %12.0f %14.0f\n",
              tau, te, im$pct_reduction, im$averted, avert_true))
  rows[[length(rows) + 1]] <- data.frame(tau = tau, t_eff = te, pct = im$pct_reduction,
                                         avert_susp = im$averted, avert_true = avert_true)
}
res <- do.call(rbind, rows)

# --- PPV sensitivity at the earliest campaign -------------------------------
tau1 <- taus[1]; te1 <- tau1 + immun_delay
ren1 <- renewal_counterfactual(weekly, w, tau = tau1, t_eff = te1, pi = coverage, psi_T = psi_T)
av_susp1 <- impact_measures(weekly, ren1$incidence_v)$averted
cat(sprintf("\n=== PPV propagation (earliest campaign wk %d): TRUE typhoid averted = pi * %.0f suspected ===\n",
            tau1, av_susp1))
for (nm in names(pi_ppv))
  cat(sprintf("  pi=%.3f (%s): %.0f true cases averted\n", pi_ppv[nm], nm, pi_ppv[nm] * av_susp1))
cat(sprintf("\nInterpretation: %% reduction (%.1f%%) is PPV-invariant (constant-pi multiplicative\n", res$pct[1]),
    "regime). Absolute TRUE burden averted -- the CEA/DALY numerator -- scales by pi=0.23,\n",
    "so ~4x fewer true cases are averted than the suspected count implies. Campaign is sized\n",
    "on suspected counts but only ~23% are vaccine-preventable typhoid.\n")

# --- plot -------------------------------------------------------------------
png("manuscript/data/poc_ori_kabwama.png", width = 1000, height = 450)
par(mfrow = c(1, 2))
plot(weekly, type = "h", lwd = 6, col = "#B0B0B0", xlab = "week", ylab = "suspected cases/week",
     main = "Kabwama 2017 (Kampala): ORI counterfactual")
lines(ren1$incidence_v, type = "h", lwd = 6, col = "#4C78A8")
abline(v = te1, lty = 2, col = "red"); legend("topright", c("observed suspected", "with ORI"),
       fill = c("#B0B0B0", "#4C78A8"), bty = "n")
plot(Rt, type = "b", pch = 19, xlab = "week", ylab = "R_t", main = "Reconstructed R_t (typhoid GI)")
abline(h = 1, lty = 2, col = "red")
dev.off()
cat("\nwrote manuscript/data/poc_ori_kabwama.png\n")
write.csv(res, "manuscript/data/poc_ori_kabwama_results.csv", row.names = FALSE)
