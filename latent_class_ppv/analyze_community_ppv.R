# =============================================================================
# Community-surveillance PPV from the outbreak Summary sheet.
# Builds a dataset of suspected-case culture positivity (+ attack rate) and
# estimates a PARTIALLY-POOLED community-surveillance PPV distribution.
# Key questions: is PPV a transferable constant? does attack rate predict it?
#   Rscript latent_class_ppv/analyze_community_ppv.R
# =============================================================================
ROOT <- Sys.getenv("RENEWAL_ROOT", "."); setwd(ROOT)
suppressMessages({ library(cmdstanr); library(posterior); library(ggplot2) })
dir.create("latent_class_ppv/figures", showWarnings = FALSE)
dir.create("latent_class_ppv/tables",  showWarnings = FALSE)
num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))

d <- read.csv("data/Typhoid_Outbreak_Time_Series_2000_2022_Summary.csv",
              check.names = FALSE, stringsAsFactors = FALSE)
D <- data.frame(study = d[["Study ID"]], country = d[["Country"]], year = d[["Outbreak start year"]],
  suspected = num(d[["Total suspected cases"]]),
  tested = num(d[["Lab tested cases by blood/bone marrow culture"]]),
  confirmed_cx = num(d[["Confirmed cases by blood/ bone marrow culture"]]),
  attack_pct = num(d[["Attack rate (%)"]]), population = num(d[["Population"]]),
  stringsAsFactors = FALSE)
D <- D[!is.na(D$study) & D$study != "" & !is.na(D$tested) & D$tested > 0 & !is.na(D$confirmed_cx), ]
D$positivity <- D$confirmed_cx / D$tested
# require an attack rate (the user's community-surveillance subset)
D <- D[!is.na(D$attack_pct) & D$attack_pct > 0, ]

# EXCLUDE confirmed-only / selectively-tested points: positivity >= 0.8 implies
# pi = positivity/Se_BC > 1 (impossible); tested==suspected means only confirmed
# cases were counted. These are not syndromic community surveillance.
D$exclude_reason <- ifelse(D$positivity >= 0.80,
  "positivity>=0.8 (confirmed-only / selective testing; pi would exceed 1)",
  ifelse(!is.na(D$suspected) & D$tested >= D$suspected,
         "tested>=suspected (confirmed-case series, not syndromic)", ""))
D$usable <- D$exclude_reason == ""
# TIGHTEN to genuinely-comparable community syndromic surveillance with EXPLICIT criteria.
# Drop hospital-presentation / no-explicit-definition series (curated from the case defs):
#   Lewis 2005 = hospital admission, "admitting physician suspects typhoid" (no criteria)
#   Roy 2016   = no case definition; presented with fever to a medical college / nursing homes
tight_drop <- c("Lewis 2005", "Roy 2016")
D$community_syndromic <- D$usable & !(D$study %in% tight_drop)
write.csv(D, "latent_class_ppv/data/community_surveillance_ppv.csv", row.names = FALSE)
cat("=== community-surveillance points (from AN + AI/AG) ===\n")
print(D[, c("study","country","year","tested","confirmed_cx","positivity","attack_pct","usable","community_syndromic")], row.names = FALSE)

Vfull <- D[D$usable, ]; V <- D[D$community_syndromic, ]     # tightened primary set
cat(sprintf("\nFull syndromic-ish: %d | TIGHT community-syndromic: %d (%s)\n",
            nrow(Vfull), nrow(V), paste(V$study, collapse = ", ")))
cat(sprintf("tight positivity: %s\n", paste(sprintf("%s=%.3f", V$study, V$positivity), collapse = ", ")))

# ---- partial-pooling on the TIGHT set (community Se_BC informative ~0.55) ----
mod <- cmdstan_model("latent_class_ppv/stan/community_ppv.stan")
fitp <- function(V, se_a, se_b, chains = 4, warm = 1000, samp = 1000)
  mod$sample(data = list(O = nrow(V), n = as.integer(V$tested), k = as.integer(V$confirmed_cx),
             se_a = se_a, se_b = se_b), seed = 1, chains = chains, parallel_chains = chains,
             iter_warmup = warm, iter_sampling = samp, adapt_delta = 0.97, refresh = 0, show_messages = FALSE)
fit <- fitp(V, 12, 10)
q <- function(f, nm) { v <- f$draws(nm, format = "draws_matrix"); quantile(as.numeric(v), c(.05, .5, .95)) }
cat("\n=== TIGHT partial-pooled community PPV (Se_BC ~ Beta(12,10), mean 0.545) ===\n")
cat(sprintf("pi_mean: %.2f (90%% CrI %.2f..%.2f)\n", q(fit,"pi_mean")[2], q(fit,"pi_mean")[1], q(fit,"pi_mean")[3]))
cat(sprintf("sigma  : %.2f (90%% CrI %.2f..%.2f)  [data-driven spread]\n", q(fit,"sigma")[2], q(fit,"sigma")[1], q(fit,"sigma")[3]))
cat(sprintf("pi_new : %.2f (90%% CrI %.2f..%.2f)  [transferable prior]\n", q(fit,"pi_new")[2], q(fit,"pi_new")[1], q(fit,"pi_new")[3]))

# ---- does the SPREAD move with Se_BC? sweep the prior mean; report pi_mean AND sigma ----
cat(sprintf("\n=== level vs spread across the Se_BC prior (c_max = %.2f) ===\n", max(V$positivity)))
cat(sprintf("%-9s %-9s %-9s\n", "Se_BC~", "pi_mean", "sigma"))
for (m in c(0.45, 0.50, 0.55, 0.60, 0.66, 0.75)) {
  a <- m*22; b <- 22-a; f2 <- fitp(V, a, b, chains = 2, warm = 800, samp = 800)
  cat(sprintf("%-9.2f %-9.2f %-9.2f\n", m,
      median(as.numeric(f2$draws("pi_mean", format = "draws_matrix"))),
      median(as.numeric(f2$draws("sigma",   format = "draws_matrix")))))
}
per <- fit$summary("pi")
V$pi_med <- per$median; V$pi_lo <- per$q5; V$pi_hi <- per$q95
write.csv(V[,c("study","country","year","tested","confirmed_cx","positivity","attack_pct","pi_med","pi_lo","pi_hi")],
          "latent_class_ppv/tables/community_ppv_estimates.csv", row.names = FALSE)

gg <- ggplot(V, aes(reorder(study, pi_med), pi_med)) +
  geom_pointrange(aes(ymin = pi_lo, ymax = pi_hi), colour = "#0072B2") +
  geom_hline(yintercept = q(fit, "pi_mean")[2], colour = "#D55E00", linetype = "dashed") +
  coord_flip(ylim = c(0,1)) +
  labs(x = NULL, y = "community-surveillance PPV (pi_o, 90% CrI)",
       title = "Community-surveillance PPV (tightened set: explicit syndromic criteria)",
       subtitle = "dashed = pooled mean; level set by the Se_BC anchor (~0.55), spread is data-driven") +
  theme_minimal(base_size = 11)
ggsave("latent_class_ppv/figures/fig_community_ppv.png", gg, width = 7.5, height = 4.5, dpi = 300)
cat("\n=== community PPV analysis complete ===\n")
