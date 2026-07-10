# =============================================================================
# PPV (pi) posterior propagation into the ORI impact analysis.
#
# The outbreak curves are SUSPECTED cases. pi = P(true typhoid | suspected) is the
# community PPV of the clinical case definition, fitted in latent_class_ppv. Under a
# constant-pi (multiplicative) regime the reconstructed R_t and the % reduction are
# pi-INVARIANT, but the absolute TRUE-typhoid burden = pi * (suspected quantity).
# So we scale the CASE columns by a per-(study, draw) pi and leave death_averted_tot
# unchanged (deaths are already true-typhoid events -> YLL is pi-invariant; YLD, COI,
# cases/dose, and true cases scale by pi).
#
# Only 3 community outbreaks were fitted (Aye 2004, Neil 2012, Kabwama 2017 =
# pi[1..3]); the hierarchical hyperparameters mu/sigma were not saved, so we
# RECONSTRUCT them per posterior draw from the 3 anchored logit-pi (mu_d = mean,
# sigma_d = sd), then draw an unanchored outbreak's pi from inv_logit(mu_d+sigma_d*z).
# This propagates both hyperparameter posterior uncertainty and between-outbreak
# heterogeneity.
# =============================================================================

# Load the community-PPV posterior. Returns NULL if the draws file is absent.
load_pi_posterior <- function(path = "latent_class_ppv/outputs/final_draws.rds",
                              anchor = c("Aye 2004", "Neil 2012", "Kabwama 2017")) {
  if (!file.exists(path)) { warning("PPV draws not found: ", path); return(NULL) }
  suppressMessages(requireNamespace("posterior", quietly = TRUE))
  dm <- posterior::as_draws_matrix(readRDS(path))
  cols <- paste0("pi[", seq_along(anchor), "]")
  if (!all(cols %in% colnames(dm))) stop("pi[] columns not found in PPV draws")
  pim <- as.matrix(dm[, cols]); colnames(pim) <- anchor
  lg <- qlogis(pmin(pmax(pim, 1e-6), 1 - 1e-6))
  list(pi_anchor = pim, mu = rowMeans(lg), sigma = apply(lg, 1, stats::sd),
       anchor = anchor, ndraw = nrow(pim))
}

# Per-(study, scenario-draw) pi matrix [ndraw x length(study_ids)].
# Anchored studies use their own posterior (index-recycled to ndraw); others draw
# from the reconstructed community hyper-posterior. Deterministic given `seed`.
build_pi_matrix <- function(post, study_ids, ndraw, seed = 20240201L) {
  idx <- ((seq_len(ndraw) - 1L) %% post$ndraw) + 1L
  M <- matrix(NA_real_, nrow = ndraw, ncol = length(study_ids),
              dimnames = list(NULL, study_ids))
  set.seed(seed)
  for (s in study_ids) {
    if (s %in% post$anchor) {
      M[, s] <- post$pi_anchor[idx, s]
    } else {
      z <- stats::rnorm(ndraw)
      M[, s] <- plogis(post$mu[idx] + post$sigma[idx] * z)
    }
  }
  M
}

# Case columns to scale from SUSPECTED -> TRUE typhoid (death_averted_tot excluded).
.ppv_case_cols <- c("s_ch_tot", "post_int_cases", "s_ch_averted_tot",
                    "s_ch_averted_amr_tot", "s_ch_averted_resistant_tot",
                    "s_ch_averted_mdr_tot", "s_ch_averted_fqns_tot", "s_ch_averted_xdr_tot")

# Scale a run_scenarios() data frame's case columns by pi in place; add `pi_ppv`.
# (MULTIPLICATIVE / constant-pi regime.)
apply_ppv_scaling <- function(df, post, study_ids, ndraw, seed = 20240201L) {
  M <- build_pi_matrix(post, study_ids, ndraw, seed)
  pv <- M[cbind(df$draw, match(df$study_id, colnames(M)))]
  cc <- intersect(.ppv_case_cols, names(df))
  df[cc] <- df[cc] * pv
  df$pi_ppv <- pv
  df
}

# --- ADDITIVE regime: S(t) = I(t) + B(t) ------------------------------------
# Non-typhoid febrile background B(t) is treated as CONSTANT over the outbreak
# (independent of typhoid transmission): B = (1 - pi) * mean(S), then I = max(S-B,0)
# renormalised so sum(I) = pi * sum(S) (self-consistent with the cumulative PPV).
# R_t is then reconstructed from the de-backgrounded I(t) (peakier than S), and the
# vaccine acts on I(t) only -> the OBSERVED suspected %reduction is diluted by B.
# (Adjustable: swap in a baseline-derived B if pre-outbreak endemic data exist.)
ppv_true_incidence <- function(S, pi) {
  if (pi >= 1) return(S)
  B <- (1 - pi) * mean(S)
  I <- pmax(S - B, 0)
  s <- sum(I)
  if (s > 0) I * (pi * sum(S) / s) else I
}
