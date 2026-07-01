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
write.csv(D, "latent_class_ppv/data/community_surveillance_ppv.csv", row.names = FALSE)
cat("=== community-surveillance points (from AN + AI/AG) ===\n")
print(D[, c("study","country","year","suspected","tested","confirmed_cx","positivity","attack_pct","usable")], row.names = FALSE)

V <- D[D$usable, ]
cat(sprintf("\nUsable (clean syndromic) points: %d\n", nrow(V)))
cat(sprintf("positivity: median=%.3f, range %.3f..%.3f (%.1f-fold)\n",
            median(V$positivity), min(V$positivity), max(V$positivity), max(V$positivity)/min(V$positivity)))
cat(sprintf("Spearman cor(positivity, log attack_pct) = %.2f  (attack-rate denominator is inconsistent across studies -> unreliable covariate)\n",
            suppressWarnings(cor(V$positivity, log(V$attack_pct), method = "spearman"))))

# ---- partial-pooling estimate (community Se_BC informative ~0.55) ------------
mod <- cmdstan_model("latent_class_ppv/stan/community_ppv.stan")
dat <- list(O = nrow(V), n = as.integer(V$tested), k = as.integer(V$confirmed_cx),
            se_a = 12, se_b = 10)   # Se_BC ~ Beta(12,10): mean 0.545, sd ~0.10 (community/mild)
fit <- mod$sample(data = dat, seed = 1, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
                  refresh = 0, show_messages = FALSE)
s <- fit$summary(c("pi_mean","sigma","Se_BC","pi_new","pi"))
q <- function(nm) { v <- fit$draws(nm, format = "draws_matrix"); quantile(as.numeric(v), c(.05,.5,.95)) }
cat("\n=== partial-pooled community-surveillance PPV (Se_BC ~ Beta(12,10), mean 0.545) ===\n")
cat(sprintf("pi_mean (typical community PPV): %.2f (90%% CrI %.2f..%.2f)\n", q("pi_mean")[2], q("pi_mean")[1], q("pi_mean")[3]))
cat(sprintf("sigma (between-setting logit SD): %.2f (90%% CrI %.2f..%.2f)  [data-driven]\n", q("sigma")[2], q("sigma")[1], q("sigma")[3]))
cat(sprintf("pi_new (predictive PPV, NEW community outbreak = transferable prior): %.2f (90%% CrI %.2f..%.2f)\n",
            q("pi_new")[2], q("pi_new")[1], q("pi_new")[3]))
# level sensitivity to the Se_BC prior
for (m in c(0.50, 0.60, 0.66)) {
  a <- m*22; b <- 22-a
  f2 <- mod$sample(data = c(dat[c("O","n","k")], list(se_a=a, se_b=b)), seed=1, chains=2,
                   iter_warmup=800, iter_sampling=800, adapt_delta=0.95, refresh=0, show_messages=FALSE)
  pm <- median(as.numeric(f2$draws("pi_mean", format="draws_matrix")))
  cat(sprintf("  [sensitivity] if community Se_BC~%.2f -> pi_mean ~ %.2f\n", m, pm))
}
per <- fit$summary("pi")
V$pi_med <- per$median; V$pi_lo <- per$q5; V$pi_hi <- per$q95
write.csv(V[,c("study","country","year","tested","confirmed_cx","positivity","attack_pct","pi_med","pi_lo","pi_hi")],
          "latent_class_ppv/tables/community_ppv_estimates.csv", row.names = FALSE)

gg <- ggplot(V, aes(reorder(study, pi_med), pi_med)) +
  geom_pointrange(aes(ymin = pi_lo, ymax = pi_hi), colour = "#0072B2") +
  geom_hline(yintercept = q("pi_mean")[2], colour = "#D55E00", linetype = "dashed") +
  coord_flip(ylim = c(0,1)) +
  labs(x = NULL, y = "community-surveillance PPV (pi_o, 90% CrI)",
       title = "Community-surveillance PPV varies ~3-fold across settings",
       subtitle = "dashed = pooled mean; PPV is a transferable DISTRIBUTION, not a constant (Se_BC~0.55)") +
  theme_minimal(base_size = 11)
ggsave("latent_class_ppv/figures/fig_community_ppv.png", gg, width = 7.5, height = 4.5, dpi = 300)
cat("\n=== community PPV analysis complete ===\n")
