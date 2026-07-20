// =============================================================================
// SENSITIVITY variant of model_final.stan: augments the community-PPV mean with
// a SINGLE objective binary covariate (alternative-diagnosis exclusion clause)
// and a NON-NEGATIVE increment:
//     logit(pi_o) = mu_pi + delta * clause_o + sigma_pi z_o,   delta >= 0.
// delta ~ half-Normal(0, delta_sd) encodes the epidemiologically certain sign
// (an exclusion clause cannot lower PPV) and shrinks toward pooling as the data
// are silent. use_covariate=0 recovers the pooled model (delta forced out) for a
// like-for-like leave-one-anchor-out comparison. log_lik is the per-outbreak
// community multinomial, for PSIS-LOO. NOT a production model.
// =============================================================================
data {
  int<lower=0> S;
  array[S] int<lower=1> N_hist;
  array[S, 4] int<lower=0> cells;
  vector[S] logv_hist;
  int<lower=0> O;
  array[O] int<lower=1> N_suspected;
  array[O] int<lower=1> n_out;
  array[O] int<lower=0> k_out;
  vector[O] logv_out;
  vector[O] bias_score;
  vector<lower=0, upper=1>[O] clause;            // exclusion-clause indicator
  vector[O] logit_tested_fraction;
  real<lower=0> mu_delta_sd;
  real beta_bias_mu;
  real<lower=0> beta_bias_sd;
  real<lower=0> sigma_delta_sd;
  real se_bm_a; real se_bm_b;
  real sp_a;    real sp_b;
  real alpha1_mu; real alpha1_sd;
  real<lower=0> delta_sd;                         // half-Normal scale for the clause effect
  int<lower=0, upper=1> use_covariate;            // 0 = pooled, 1 = clause covariate
}
parameters {
  real alpha0;
  real alpha1;
  real<lower=0> tau;
  vector[S] u_hist;
  vector[O] u_out;
  real<lower=0, upper=1> Se_BM;
  real<lower=0, upper=1> Sp_BC;
  real<lower=0, upper=1> Sp_BM;
  vector<lower=0, upper=1>[S] phi;
  real mu_pi;
  real<lower=0> sigma_pi;
  vector[O] z_pi;
  real<lower=0> delta;                            // clause effect (logit-PPV), non-negative
  real mu_delta;
  real<lower=0> beta_bias;
  real<lower=0> sigma_delta;
  vector[O] z_delta;
  vector[O] logit_q0;
}
transformed parameters {
  vector<lower=0, upper=1>[S] Se_BC_hist;
  vector<lower=0, upper=1>[O] Se_BC_out;
  vector<lower=0, upper=1>[O] pi;
  vector[O] delta_sel;
  vector<lower=0, upper=1>[O] q0;
  vector<lower=0, upper=1>[O] q1;
  for (s in 1:S) Se_BC_hist[s] = inv_logit(alpha0 + alpha1 * logv_hist[s] + tau * u_hist[s]);
  for (o in 1:O) Se_BC_out[o]  = inv_logit(alpha0 + alpha1 * logv_out[o]  + tau * u_out[o]);
  pi = inv_logit(mu_pi + use_covariate * delta * clause + sigma_pi * z_pi);
  delta_sel = mu_delta + beta_bias * bias_score + sigma_delta * z_delta;
  q0 = inv_logit(logit_q0);
  q1 = inv_logit(logit_q0 + delta_sel);
}
model {
  alpha0 ~ normal(0, 1.5);
  alpha1 ~ normal(alpha1_mu, alpha1_sd);
  tau    ~ normal(0, 1);
  u_hist ~ std_normal();
  u_out  ~ std_normal();
  Se_BM  ~ beta(se_bm_a, se_bm_b);
  Sp_BC  ~ beta(sp_a, sp_b);
  Sp_BM  ~ beta(sp_a, sp_b);
  phi ~ beta(1, 1);
  mu_pi ~ normal(logit(0.25), 1.0);
  sigma_pi ~ normal(0, 0.8);
  z_pi ~ std_normal();
  delta ~ normal(0, delta_sd);                    // half-Normal (lower=0)
  mu_delta ~ normal(0, mu_delta_sd);
  beta_bias ~ normal(beta_bias_mu, beta_bias_sd);
  sigma_delta ~ normal(0, sigma_delta_sd);
  z_delta ~ std_normal();
  logit_q0 ~ normal(logit_tested_fraction, 1.0);
  for (s in 1:S) {
    real se = Se_BC_hist[s];
    real sm = Se_BM;
    vector[4] p;
    p[1] = phi[s] * se * sm            + (1 - phi[s]) * (1 - Sp_BC) * (1 - Sp_BM);
    p[2] = phi[s] * se * (1 - sm)      + (1 - phi[s]) * (1 - Sp_BC) * Sp_BM;
    p[3] = phi[s] * (1 - se) * sm      + (1 - phi[s]) * Sp_BC * (1 - Sp_BM);
    p[4] = phi[s] * (1 - se) * (1 - sm)+ (1 - phi[s]) * Sp_BC * Sp_BM;
    target += multinomial_lpmf(cells[s] | p);
  }
  for (o in 1:O) {
    vector[3] p;
    array[3] int y;
    p[1] = pi[o] * q1[o] * Se_BC_out[o] + (1 - pi[o]) * q0[o] * (1 - Sp_BC);
    p[2] = pi[o] * q1[o] * (1 - Se_BC_out[o]) + (1 - pi[o]) * q0[o] * Sp_BC;
    p[3] = pi[o] * (1 - q1[o]) + (1 - pi[o]) * (1 - q0[o]);
    y[1] = k_out[o]; y[2] = n_out[o] - k_out[o]; y[3] = N_suspected[o] - n_out[o];
    target += multinomial_lpmf(y | p);
  }
}
generated quantities {
  real Se_BC_5mL = inv_logit(alpha0 + alpha1 * log(5));
  real pi_pooled  = inv_logit(mu_pi);                          // no clause
  real pi_clause  = inv_logit(mu_pi + use_covariate * delta);  // with clause
  real pi_new_pooled = inv_logit(normal_rng(mu_pi, sigma_pi));
  real pi_new_clause = inv_logit(normal_rng(mu_pi + use_covariate * delta, sigma_pi));
  vector[O] log_lik;
  for (o in 1:O) {
    vector[3] p;
    array[3] int y;
    p[1] = pi[o] * q1[o] * Se_BC_out[o] + (1 - pi[o]) * q0[o] * (1 - Sp_BC);
    p[2] = pi[o] * q1[o] * (1 - Se_BC_out[o]) + (1 - pi[o]) * q0[o] * Sp_BC;
    p[3] = pi[o] * (1 - q1[o]) + (1 - pi[o]) * (1 - q0[o]);
    y[1] = k_out[o]; y[2] = n_out[o] - k_out[o]; y[3] = N_suspected[o] - n_out[o];
    log_lik[o] = multinomial_lpmf(y | p);
  }
}
