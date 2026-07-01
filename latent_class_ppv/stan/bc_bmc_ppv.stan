// =============================================================================
// Joint estimation of blood-culture sensitivity (Se_BC) and clinical-case-
// definition PPV (pi) for typhoid, with a bacterial-LOAD (volume + severity)
// sub-model and an explicit severe->mild transfer offset.
//
// HISTORIC two-test studies (blood + bone marrow, same patients): a
// conditional-independence latent-class 2x2 that estimates Se_BC AND Se_BM
// (bone marrow treated as IMPERFECT) and the study-level hospital PPV phi_s.
//   cells = (a,b,c,d) = (BC+BM+, BC+BM-, BC-BM+, BC-BM-) ~ Multinomial(N, p)
//
// OUTBREAK single-test studies (blood culture only): the anchored Se_BC
// identifies each outbreak's surveillance PPV pi_o.
//   k_o ~ Binomial(n_o, pi_o*Se_BC_o + (1-pi_o)*(1-Sp_BC))
//
// Se_BC sub-model (load axis):
//   logit(Se_BC) = alpha0 + alpha1*log(volume) + beta*mild + tau*u
//   mild = 1 for a lower-load (mild/surveillance) population, 0 for the
//   severe hospital paired-study baseline. beta <= 0 (mild cannot raise Se).
//
// Transportable (alpha0, alpha1, beta, tau, Se_BM) vs LOCAL (phi_s, pi_o):
//   only the Se sub-model + Se_BM are shared; every prevalence/PPV is local.
//
// cd_infected: optional FIXED conditional-dependence covariance among the
//   truly infected (0 = independence). Used to stress-test the independence
//   assumption (simulate with cd>0, fit with cd=0).
// =============================================================================
data {
  // ---- historic two-test studies ----
  int<lower=0> S;
  array[S] int<lower=1> N_hist;              // enrolled & tested by both
  array[S, 4] int<lower=0> cells;            // a,b,c,d (outbreak-major order above)
  vector[S] logv_hist;                       // log blood volume (mL)
  vector<lower=0, upper=1>[S] mild_hist;     // severity/load indicator (0 = severe baseline)
  real cd_infected;                          // fixed dependence covariance (0 = independence)
  // ---- outbreak single-test studies ----
  int<lower=0> O;
  array[O] int<lower=1> n_out;
  array[O] int<lower=0> k_out;
  vector[O] logv_out;
  vector<lower=0, upper=1>[O] mild_out;
  // ---- prior hyperparameters ----
  real se_bm_a; real se_bm_b;                // Se_BM ~ Beta (informative, high)
  real sp_a;    real sp_b;                    // Sp_BC, Sp_BM ~ Beta (near 1)
  real beta_sd;                              // half-normal scale for the mild offset
  int<lower=0, upper=1> beta_locked;         // 1 = fix the offset (sensitivity mode)
  real beta_lock_val;                        // fixed offset value when locked (<=0)
}
parameters {
  real alpha0;                               // Se volume-model intercept (logit)
  real alpha1;                               // volume slope (per log mL)
  real<upper=0> beta;                        // mild (lower-load) offset, sign-constrained
  real<lower=0> tau;                         // between-study SD (logit Se)
  vector[S] u_hist;                          // non-centred study effects
  vector[O] u_out;                           // non-centred outbreak effects
  real<lower=0, upper=1> Se_BM;              // bone-marrow sensitivity (imperfect)
  real<lower=0, upper=1> Sp_BC;
  real<lower=0, upper=1> Sp_BM;
  vector<lower=0, upper=1>[S] phi;           // hospital PPV per historic study (LOCAL)
  vector<lower=0, upper=1>[O] pi;            // surveillance PPV per outbreak (LOCAL)
}
transformed parameters {
  real beta_use = (beta_locked == 1) ? beta_lock_val : beta;   // fixed offset in sensitivity mode
  vector<lower=0, upper=1>[S] Se_BC_hist;
  vector<lower=0, upper=1>[O] Se_BC_out;
  for (s in 1:S)
    Se_BC_hist[s] = inv_logit(alpha0 + alpha1 * logv_hist[s] + beta_use * mild_hist[s] + tau * u_hist[s]);
  for (o in 1:O)
    Se_BC_out[o]  = inv_logit(alpha0 + alpha1 * logv_out[o]  + beta_use * mild_out[o]  + tau * u_out[o]);
}
model {
  // ---- priors (transportable Se sub-model + test properties) ----
  alpha0 ~ normal(0, 1.5);
  alpha1 ~ normal(0, 1);                     // volume effect, weakly-informative
  beta   ~ normal(0, beta_sd);               // half-normal (beta<=0): weakly-informative offset
  tau    ~ normal(0, 1);                     // half-normal
  u_hist ~ std_normal();
  u_out  ~ std_normal();
  Se_BM  ~ beta(se_bm_a, se_bm_b);
  Sp_BC  ~ beta(sp_a, sp_b);
  Sp_BM  ~ beta(sp_a, sp_b);
  // ---- LOCAL prevalences: independent, NOT pooled ----
  phi ~ beta(1, 1);
  pi  ~ beta(1, 1);
  // ---- historic likelihood: conditional-independence LCM (+ optional dependence) ----
  for (s in 1:S) {
    real se = Se_BC_hist[s];
    real sm = Se_BM;
    // infected-class joint cell probs (cd_infected shifts mass onto agreement cells)
    real i_pp = se * sm + cd_infected;
    real i_pn = se * (1 - sm) - cd_infected;
    real i_np = (1 - se) * sm - cd_infected;
    real i_nn = (1 - se) * (1 - sm) + cd_infected;
    // uninfected-class (independence; specificities near 1)
    real v_pp = (1 - Sp_BC) * (1 - Sp_BM);
    real v_pn = (1 - Sp_BC) * Sp_BM;
    real v_np = Sp_BC * (1 - Sp_BM);
    real v_nn = Sp_BC * Sp_BM;
    vector[4] p;
    p[1] = phi[s] * i_pp + (1 - phi[s]) * v_pp;
    p[2] = phi[s] * i_pn + (1 - phi[s]) * v_pn;
    p[3] = phi[s] * i_np + (1 - phi[s]) * v_np;
    p[4] = phi[s] * i_nn + (1 - phi[s]) * v_nn;
    target += multinomial_lpmf(cells[s] | p);
  }
  // ---- outbreak likelihood: single-test binomial ----
  for (o in 1:O) {
    real theta = pi[o] * Se_BC_out[o] + (1 - pi[o]) * (1 - Sp_BC);
    k_out[o] ~ binomial(n_out[o], theta);
  }
}
generated quantities {
  // Predictive Se_BC at reference volumes for severe (mild=0) and mild (mild=1),
  // marginal over between-study heterogeneity (a fresh u ~ N(0,tau)).
  real u_new = normal_rng(0, tau);
  real Se_severe_5mL = inv_logit(alpha0 + alpha1 * log(5)  + u_new);
  real Se_mild_5mL   = inv_logit(alpha0 + alpha1 * log(5)  + beta_use + u_new);
  real Se_mild_10mL  = inv_logit(alpha0 + alpha1 * log(10) + beta_use + u_new);
}
