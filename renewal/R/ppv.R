# =============================================================================
# PPV (pi) posterior propagation into the ORI impact analysis.
#
# The outbreak curves are SUSPECTED cases. pi = P(true typhoid | suspected) is
# the community PPV of the clinical case definition, fitted in latent_class_ppv.
# The primary observation model is additive: suspected S_t = true typhoid T_t
# + other febrile F_t. Vaccination acts on T_t and F_t is held fixed, so PPV can
# affect both absolute burden and percentage reductions. The multiplicative
# alternative T_t = pi*S_t is retained as a structural sensitivity; only under
# that alternative does pi cancel from R_t and percentage reductions.
#
# Deaths are already treated as true-typhoid events and are not PPV-scaled.
#
# Community outbreaks and their pi[] order are read from the fitted output and
# cross-checked against the warning-bearing merged audit. model_final.stan fits
# these pi values within a fitted logit-normal population distribution. The
# posterior draws of mu_pi and sigma_pi define a genuine posterior predictive
# distribution for unanchored outbreaks.
# =============================================================================

# Load the community-PPV posterior. Returns NULL if the draws file is absent.
load_pi_posterior <- function(path = "latent_class_ppv/outputs/final_draws.rds",
                              anchor = NULL,
                              pi_table = "latent_class_ppv/tables/final_pi_community.csv",
                              audit_path = "latent_class_ppv/data/merged_outbreak_ppv_modeling_audit.csv") {
  if (!file.exists(path)) { warning("PPV draws not found: ", path); return(NULL) }
  if (is.null(anchor)) {
    if (!file.exists(pi_table)) stop("PPV study-order table not found: ", pi_table)
    anchor <- read.csv(pi_table, stringsAsFactors = FALSE)$study
  }
  if (file.exists(audit_path)) {
    audit <- read.csv(audit_path, check.names = FALSE, stringsAsFactors = FALSE)
    included <- audit$study[audit$include_for_main_ppv_model %in% c(TRUE, "TRUE", "True")]
    if (!identical(anchor, included)) {
      stop("Fitted PPV study order does not match merged audit inclusion order")
    }
    groups <- audit$independence_group[match(anchor, audit$study)]
    if (anyDuplicated(groups[nzchar(groups)])) {
      stop("Fitted PPV anchors contain duplicate independence groups")
    }
  }
  suppressMessages(requireNamespace("posterior", quietly = TRUE))
  dm <- posterior::as_draws_matrix(readRDS(path))
  cols <- paste0("pi[", seq_along(anchor), "]")
  if (!all(cols %in% colnames(dm))) stop("pi[] columns not found in PPV draws")
  if (!all(c("mu_pi", "sigma_pi") %in% colnames(dm))) {
    stop("PPV draws lack hierarchical mu_pi/sigma_pi; unanchored prediction is invalid")
  }
  pim <- as.matrix(dm[, cols]); colnames(pim) <- anchor
  mu_pi <- as.numeric(dm[, "mu_pi"])
  sigma_pi <- as.numeric(dm[, "sigma_pi"])
  list(pi_anchor = pim, mu_pi = mu_pi, sigma_pi = sigma_pi,
       # Backward-compatible aliases for scripts written before the fitted
       # hyperparameters were exposed explicitly.
       mu = mu_pi, sigma = sigma_pi,
       anchor = anchor, ndraw = nrow(pim), supports_unanchored_prediction = TRUE)
}

# Per-(study, scenario-draw) pi matrix [ndraw x length(study_ids)].
# Anchored studies use their own posterior (index-recycled to ndraw). Unanchored
# studies receive independent posterior-predictive draws conditional on each
# draw's fitted mu_pi and sigma_pi.
build_pi_matrix <- function(post, study_ids, ndraw, seed = 20240201L) {
  idx <- ((seq_len(ndraw) - 1L) %% post$ndraw) + 1L
  M <- matrix(NA_real_, nrow = ndraw, ncol = length(study_ids),
              dimnames = list(NULL, study_ids))
  set.seed(seed)
  for (s in study_ids) {
    if (s %in% post$anchor) {
      M[, s] <- post$pi_anchor[idx, s]
    } else {
      M[, s] <- plogis(stats::rnorm(ndraw, post$mu_pi[idx], post$sigma_pi[idx]))
    }
  }
  M
}

# Case columns used by the multiplicative sensitivity scaler
# (death_averted_tot excluded).
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

# --- ADDITIVE observation model: S_t = T_t + F_t -----------------------------
# S_t is suspected febrile incidence, T_t latent true typhoid, and F_t other
# febrile illness. A flat candidate background level, (1-pi)*mean(S), defines
# the temporal shape of outbreak excess. After truncation at zero, T is
# renormalised to satisfy the cumulative PPV exactly:
#
#   sum(T) = pi*sum(S),  F = S-T.
#
# Consequently F is fixed interval-by-interval under the vaccine
# counterfactual, but the final residual F_t is not generally perfectly flat
# after truncation and renormalisation. Returning both components makes that
# mass balance explicit and prevents prose from overstating the flatness.
ppv_incidence_components <- function(S, pi) {
  S <- as.numeric(S)
  if (!length(S) || any(!is.finite(S)) || any(S < 0))
    stop("S must be a non-negative finite vector")
  if (length(pi) != 1L || !is.finite(pi) || pi < 0 || pi > 1)
    stop("pi must be a scalar in [0,1]")

  target <- pi * sum(S)
  flat_candidate <- (1 - pi) * mean(S)
  excess <- pmax(S - flat_candidate, 0)
  if (target == 0 || sum(excess) == 0) {
    true_typhoid <- numeric(length(S))
  } else {
    true_typhoid <- excess * (target / sum(excess))
  }
  other_febrile <- S - true_typhoid
  tol <- 1e-10 * max(1, max(S))
  if (any(true_typhoid < -tol) || any(other_febrile < -tol))
    stop("PPV allocation produced a negative observation component")
  true_typhoid[abs(true_typhoid) < tol] <- 0
  other_febrile[abs(other_febrile) < tol] <- 0

  list(
    suspected = S,
    true_typhoid = true_typhoid,
    other_febrile = other_febrile,
    ppv = pi,
    flat_background_candidate = flat_candidate,
    mass_balance_error = max(abs(S - true_typhoid - other_febrile)),
    ppv_error = abs(sum(true_typhoid) - target)
  )
}

# Compatibility wrapper used by existing scenario and figure scripts.
ppv_true_incidence <- function(S, pi)
  ppv_incidence_components(S, pi)$true_typhoid
