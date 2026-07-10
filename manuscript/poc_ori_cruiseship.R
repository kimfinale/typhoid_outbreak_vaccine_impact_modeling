# POC: estimate weekly R_t for the digitized cruise-ship outbreak (Ooms 2024) using
# the typhoid generation interval (mean 14 d, sd 8.4 d) from renewal/R/gi.R, and
# illustrate the PPV (pi) correction for a suspected-case observation regime.
# Honest caveats: (1) counts are hand-digitized from Fig1 (approximate; sum ~40 vs 52
# reported symptomatic); (2) this is a COMMON-SOURCE waterborne outbreak whose R_t
# reflects source/exposure dynamics + control (bottled water 06 Apr, evacuation 15 Apr),
# not sustained person-to-person transmission -- so it is a pipeline demo, not an ORI
# effectiveness estimate. A propagated outbreak is needed for the latter.
#   Rscript manuscript/poc_ori_cruiseship.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages(library(EpiEstim))

d <- read.csv("manuscript/data/ts_poc_cruiseship.csv", stringsAsFactors = FALSE)
inc <- d$confirmed
cat("weeks:", nrow(d), " total confirmed (digitized):", sum(inc), " (reported symptomatic: 52)\n\n")

# --- typhoid GI discretized to WEEKLY steps (matches renewal/R/gi.R defaults) ---
mean_d <- 14; sd_d <- 8.4; step <- 7
shape <- (mean_d / sd_d)^2; rate <- mean_d / sd_d^2
maxw <- 8
edges <- seq(0, maxw * step, by = step)
w <- diff(pgamma(edges, shape = shape, rate = rate)); w <- w / sum(w)
si <- c(0, w)                                  # SI mass at lag 0,1,2,... weeks (0 at lag0)
si <- si / sum(si)

# --- R_t via EpiEstim on the weekly series (sliding weekly windows) ---
Tn <- length(inc)
res <- estimate_R(incid = inc,
                  method = "non_parametric_si",
                  config = make_config(list(si_distr = si,
                                            t_start = 2:(Tn - 1), t_end = 3:Tn)))
rt <- res$R
cat("=== weekly R_t (EpiEstim, typhoid GI) ===\n")
for (i in seq_len(nrow(rt))) {
  wk <- d$week_start[rt$t_end[i]]
  cat(sprintf("  wk %s : R_t = %.2f [%.2f, %.2f]\n", wk,
              rt$`Mean(R)`[i], rt$`Quantile.0.025(R)`[i], rt$`Quantile.0.975(R)`[i]))
}
cat(sprintf("\npeak R_t = %.2f ; final-week R_t = %.2f\n",
            max(rt$`Mean(R)`), rt$`Mean(R)`[nrow(rt)]))

# --- PPV (pi) correction illustration (suspected-case observation regime) ---
# If surveillance had reported SUSPECTED cases, observed = confirmed / (pi * Se_BC-ish);
# constant-pi multiplicative regime: R_t is INVARIANT to a constant pi (ratio-based),
# but the OBSERVED outbreak SIZE (hence ORI resource sizing) scales by 1/pi.
pis <- c(Aye = 0.26, Kabwama = 0.23, Neil = 0.65)
cat("\n=== suspected-case inflation of outbreak size by community pi (size = confirmed/pi) ===\n")
for (nm in names(pis))
  cat(sprintf("  pi=%.2f (%s): implied suspected total ~ %.0f (x%.1f)\n",
              pis[nm], nm, sum(inc) / pis[nm], 1 / pis[nm]))
cat("\nNote: under constant pi the % reduction from ORI is invariant, but suspected-case\n",
    "counts (what triggers/размеры a campaign) are inflated ~", sprintf("%.1f-%.1fx", 1/max(pis), 1/min(pis)),
    "; additive background (B) would instead DILUTE the observed % reduction.\n")

png("manuscript/data/ts_poc_cruiseship_Rt.png", width = 900, height = 500)
par(mfrow = c(1, 2))
barplot(inc, names.arg = substr(d$week_start, 6, 10), las = 2, col = "#4C78A8",
        main = "Cruise-ship outbreak (digitized)", ylab = "confirmed cases/week")
plot(rt$t_end, rt$`Mean(R)`, type = "b", pch = 19, ylim = c(0, max(rt$`Quantile.0.975(R)`)),
     xlab = "week index", ylab = "R_t", main = "Weekly R_t (typhoid GI)")
abline(h = 1, lty = 2, col = "red")
arrows(rt$t_end, rt$`Quantile.0.025(R)`, rt$t_end, rt$`Quantile.0.975(R)`,
       angle = 90, code = 3, length = 0.03, col = "grey")
dev.off()
cat("\nwrote manuscript/data/ts_poc_cruiseship_Rt.png\n")
