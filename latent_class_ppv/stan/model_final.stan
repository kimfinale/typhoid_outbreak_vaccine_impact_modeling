// =============================================================================
// FINAL model (beta-free): joint estimation of blood-culture sensitivity Se_BC
// (as a function of blood VOLUME), bone-marrow sensitivity Se_BM, and the
// clinical-case-definition PPV under hospital (phi_s) and community (pi_o)
// ascertainment. The severe->mild offset was dropped (literature verdict:
// severity is not a supported, identifiable determinant of Se_BC; VOLUME is).
//
// HISTORIC paired blood+bone-marrow studies -> conditional-independence latent-
//   class 2x2 (Hui-Walter): identifies Se_BC level, Se_BM, and hospital PPV phi_s.
//   cells = (BC+BM+, BC+BM-, BC-BM+, BC-BM-) ~ Multinomial(N, p)
// OUTBREAK studies -> hierarchical PPV among ALL suspected cases plus a
// selection model for entry into the blood-cultured subset. Observed counts are
// culture-positive, culture-negative, and untested among N_suspected.
//
// Se sub-model:  logit(Se_BC) = alpha0 + alpha1*log(volume) + tau*u
//   alpha1 = volume slope (Antillon-informed prior; historic volumes are noisy)
//   alpha0 = Se level, tau = between-study heterogeneity.
// Transportable: alpha0, alpha1, tau, Se_BM and the population PPV
// distribution (mu_pi, sigma_pi). Local: phi_s, pi_o and selection intensity.
// =============================================================================
data {
  int<lower=0> S;                              // historic paired studies
  array[S] int<lower=1> N_hist;
  array[S, 4] int<lower=0> cells;              // BC+BM+, BC+BM-, BC-BM+, BC-BM-
  vector[S] logv_hist;                         // log blood volume (mL)
  int<lower=0> O;                              // outbreak single-test studies
  array[O] int<lower=1> N_suspected;
  array[O] int<lower=1> n_out;
  array[O] int<lower=0> k_out;
  vector[O] logv_out;                          // log blood volume (assumed) for outbreaks
  vector[O] bias_score;                        // 0=low, 1=moderate, 2=major selection concern
  vector[O] logit_tested_fraction;              // weak prior center for non-case testing probability
  real<lower=0> mu_delta_sd;
  real beta_bias_mu;
  real<lower=0> beta_bias_sd;
  real<lower=0> sigma_delta_sd;
  real se_bm_a; real se_bm_b;                  // Se_BM ~ Beta (favours high; breaks LCM symmetry)
  real sp_a;    real sp_b;                      // Sp ~ Beta (near 1)
  real alpha1_mu; real alpha1_sd;              // volume-slope prior (Antillon)
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
  vector<lower=0, upper=1>[S] phi;             // hospital PPV per historic study (LOCAL)
  real mu_pi;                                  // population mean on logit scale
  real<lower=0> sigma_pi;                      // genuine between-outbreak heterogeneity
  vector[O] z_pi;
  real mu_delta;                               // baseline log-odds selection enrichment
  real<lower=0> beta_bias;                     // larger audit bias grade -> stronger enrichment
  real<lower=0> sigma_delta;                   // residual selection heterogeneity
  vector[O] z_delta;
  vector[O] logit_q0;                          // P(cultured | not true typhoid)
}
transformed parameters {
  vector<lower=0, upper=1>[S] Se_BC_hist;
  vector<lower=0, upper=1>[O] Se_BC_out;
  vector<lower=0, upper=1>[O] pi;              // PPV among ALL suspected cases
  vector[O] delta;
  vector<lower=0, upper=1>[O] q0;
  vector<lower=0, upper=1>[O] q1;
  for (s in 1:S) Se_BC_hist[s] = inv_logit(alpha0 + alpha1 * logv_hist[s] + tau * u_hist[s]);
  for (o in 1:O) Se_BC_out[o]  = inv_logit(alpha0 + alpha1 * logv_out[o]  + tau * u_out[o]);
  pi = inv_logit(mu_pi + sigma_pi * z_pi);
  delta = mu_delta + beta_bias * bias_score + sigma_delta * z_delta;
  q0 = inv_logit(logit_q0);
  q1 = inv_logit(logit_q0 + delta);
}
model {
  alpha0 ~ normal(0, 1.5);
  alpha1 ~ normal(alpha1_mu, alpha1_sd);       // Antillon-informed volume slope
  tau    ~ normal(0, 1);                        // half-normal
  u_hist ~ std_normal();
  u_out  ~ std_normal();
  Se_BM  ~ beta(se_bm_a, se_bm_b);
  Sp_BC  ~ beta(sp_a, sp_b);
  Sp_BM  ~ beta(sp_a, sp_b);
  phi ~ beta(1, 1);
  mu_pi ~ normal(logit(0.25), 1.0);
  sigma_pi ~ normal(0, 0.8);
  z_pi ~ std_normal();
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
    p[1] = pi[o] * q1[o] * Se_BC_out[o] +
           (1 - pi[o]) * q0[o] * (1 - Sp_BC);
    p[2] = pi[o] * q1[o] * (1 - Se_BC_out[o]) +
           (1 - pi[o]) * q0[o] * Sp_BC;
    p[3] = pi[o] * (1 - q1[o]) + (1 - pi[o]) * (1 - q0[o]);
    y[1] = k_out[o];
    y[2] = n_out[o] - k_out[o];
    y[3] = N_suspected[o] - n_out[o];
    target += multinomial_lpmf(y | p);
  }
}
generated quantities {
  real Se_BC_2mL  = inv_logit(alpha0 + alpha1 * log(2));   // marginal-of-heterogeneity mean curve
  real Se_BC_5mL  = inv_logit(alpha0 + alpha1 * log(5));
  real Se_BC_10mL = inv_logit(alpha0 + alpha1 * log(10));
  real ProCase_BC = inv_logit(alpha0 + alpha1 * log(5));   // ~"proportion detected" at ref volume
  real pi_population_median = inv_logit(mu_pi);
  real pi_new = inv_logit(normal_rng(mu_pi, sigma_pi));
}
