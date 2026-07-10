# Fit the two-stratum (outbreak vs endemic-surveillance) community PPV model on the
# expanded 40-PDF positivity extraction. Characterises the endemicity gradient in pi.
#   Rscript manuscript/fit_ppv_strata.R
setwd(Sys.getenv("RENEWAL_ROOT", "."))
suppressMessages({ library(cmdstanr); library(posterior) })

d <- read.csv("manuscript/data/ppv_positivity_table.csv", stringsAsFactors = FALSE)
d <- d[grepl("outbreak|surveillance", d$stratum) & !grepl("Misra", d$study) &
       !is.na(d$tested_bc) & !is.na(d$confirmed_bc), ]        # drop Misra (confirmed-enriched hospital MDR)
# Assign strata by KNOWN study TYPE (community outbreak investigation vs endemic
# surveillance program) -- the tested/suspected ratio heuristic misfiles borderline
# surveillance studies (e.g. Voysey, ratio 0.68) into 'outbreak'.
outbreak_pat <- "Aye 2004|Neil 2012|Kabwama 2017|Cho 2019"
d$outbreak <- as.integer(grepl(outbreak_pat, d$study))
d$stratum <- ifelse(d$outbreak == 1, "outbreak (selective testing)", "endemic surveillance (all-tested)")
d <- d[order(-d$outbreak, d$study), ]
cat("studies:", nrow(d), " | outbreak:", sum(d$outbreak), " surveillance:", sum(1 - d$outbreak), "\n")
print(d[, c("study", "tested_bc", "confirmed_bc", "outbreak")], row.names = FALSE)

sd <- list(O = nrow(d), n = as.integer(d$tested_bc), k = as.integer(d$confirmed_bc),
           outbreak = as.numeric(d$outbreak), se_a = 24, se_b = 15)   # Se_BC ~ Beta(24,15): mean 0.615
mod <- cmdstan_model("latent_class_ppv/stan/community_ppv_strata.stan")
fit <- mod$sample(data = sd, seed = 2026, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000, refresh = 0)

dm <- fit$draws(format = "draws_matrix")
qs <- function(v) round(quantile(as.numeric(v), c(.5, .025, .975)), 3)
gq <- c("pi_surv", "pi_out", "OR_out_vs_surv", "Se_BC", "sigma", "pi_new_out", "pi_new_surv")
tab <- data.frame(param = gq, t(sapply(gq, function(g) qs(dm[, g]))))
names(tab) <- c("param", "med", "lo", "hi")
write.csv(tab, "latent_class_ppv/tables/community_ppv_strata_params.csv", row.names = FALSE)

pistud <- sapply(seq_len(nrow(d)), function(i) qs(dm[, paste0("pi[", i, "]")]))
res <- data.frame(study = d$study, stratum = d$stratum, tested = d$tested_bc, conf = d$confirmed_bc,
                  positivity = round(d$confirmed_bc / d$tested_bc, 3),
                  pi_med = pistud[1, ], pi_lo = pistud[2, ], pi_hi = pistud[3, ])
write.csv(res, "latent_class_ppv/tables/community_ppv_strata.csv", row.names = FALSE)
saveRDS(fit$draws(c("pi_surv", "pi_out", "OR_out_vs_surv", "pi_new_out", "pi_new_surv",
                    "mu", "beta", "sigma", "Se_BC")), "latent_class_ppv/outputs/community_strata_draws.rds")

cat("\n=== stratum PPV (median [95% CrI]) ===\n")
cat(sprintf("  endemic surveillance pi = %.3f [%.3f, %.3f]\n", tab$med[1], tab$lo[1], tab$hi[1]))
cat(sprintf("  outbreak             pi = %.3f [%.3f, %.3f]\n", tab$med[2], tab$lo[2], tab$hi[2]))
cat(sprintf("  OR (outbreak vs surv)   = %.1f [%.1f, %.1f]\n", tab$med[3], tab$lo[3], tab$hi[3]))
cat("wrote community_ppv_strata.csv + _params.csv + community_strata_draws.rds\n")
