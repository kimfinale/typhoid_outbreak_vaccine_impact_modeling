// =============================================================================
// Hierarchical Bayesian renewal-with-source model for the source fraction theta.
//
//   I[o,t] ~ NegBinomial2( mu[o,t], phi )
//   mu[o,t] = s[o,t] + R_tr[o] * Lambda[o,t]
//   Lambda[o,t] = sum_{k>=1} w[k] * I[o, t-k]      (OBSERVED lagged incidence)
//
// s = exogenous common-source term: a single smooth gamma pulse in time x amplitude
//     (the constraint that prevents the source absorbing arbitrary curve structure).
// R_tr[o] in (0,1) = subcritical per-outbreak transmission (propagated channel).
//   theta_o = 1 - R_tr[o]; logit(theta_o) partially pooled (non-centered).
//   informed=0 -> pooled to a grand mean (neutral); informed=1 -> mode-specific mean.
// =============================================================================
data {
  int<lower=1> O;                          // outbreaks
  int<lower=1> N;                          // total observations
  array[N] int<lower=0> I;                 // incidence (flattened, outbreak-major)
  array[O] int<lower=1> Tn;                // length per outbreak
  array[O] int<lower=1> start;             // 1-based start index of each outbreak in I
  vector[O] logA_center;                   // per-outbreak amplitude prior center = log(mean I)
  int<lower=1> K;                          // max GI lag (weeks)
  real<lower=0> gi_cv;                      // fixed GI coefficient of variation
  real gi_mean_prior_mu;                   // lognormal prior (log scale) for GI mean (weeks)
  real<lower=0> gi_mean_prior_sd;
  int<lower=0, upper=1> informed;          // 0 neutral, 1 article-informed (mode covariate)
  int<lower=1> n_modes;
  array[O] int<lower=1, upper=n_modes> mode;
}
parameters {
  real<lower=0> gi_mean;                   // generation-interval mean (weeks)
  real mu_theta;                           // grand mean (logit theta)
  vector[n_modes] mode_offset;             // mode-specific offsets (used iff informed)
  real<lower=0> tau_theta;                 // pooling sd (half-normal)
  vector[O] z_raw;                         // non-centered hierarchy
  vector[O] logA;                          // log source amplitude
  vector<lower=0>[O] pulse_mu;             // source pulse location (weeks)
  vector<lower=0>[O] pulse_sd;             // source pulse width (weeks)
  real<lower=0> phi;                       // NB dispersion
}
transformed parameters {
  vector[O] z;
  vector<lower=0, upper=1>[O] theta;
  vector<lower=0, upper=1>[O] R_tr;
  vector[K] w;
  vector[N] mu;
  vector[N] s;
  {
    real shape = 1 / square(gi_cv);
    real rate = shape / gi_mean;
    vector[K] wun;
    for (k in 1:K) wun[k] = gamma_cdf(k | shape, rate) - gamma_cdf(k - 1 | shape, rate);
    w = wun / sum(wun);
  }
  for (o in 1:O) {
    real m = mu_theta + (informed == 1 ? mode_offset[mode[o]] : 0);
    z[o] = m + tau_theta * z_raw[o];
    theta[o] = inv_logit(z[o]);
    R_tr[o] = 1 - theta[o];
  }
  for (o in 1:O) {
    int st = start[o];
    int T = Tn[o];
    real psh = square(pulse_mu[o] / pulse_sd[o]);
    real prt = pulse_mu[o] / square(pulse_sd[o]);
    vector[T] p;
    for (j in 1:T) p[j] = exp(gamma_lpdf(j | psh, prt));
    p /= sum(p);
    for (j in 1:T) {
      int g = st + j - 1;
      s[g] = exp(logA[o]) * p[j];
      real lam = 0;
      int kmax = min(K, j - 1);
      for (k in 1:kmax) lam += w[k] * I[g - k];
      mu[g] = s[g] + R_tr[o] * lam;
    }
  }
}
model {
  gi_mean ~ lognormal(gi_mean_prior_mu, gi_mean_prior_sd);
  // Weakly-informative source-leaning prior on logit(theta): typhoid ORI
  // transmission is the MINORITY channel (subcritical, source-dominated regime).
  // This is a structural/domain assumption, NOT the article-mode classification
  // (the neutral/informed contrast concerns the mode covariate only).
  mu_theta ~ normal(0.85, 0.8);            // theta centered ~0.70
  mode_offset ~ normal(0, 1);
  tau_theta ~ normal(0, 0.7);              // half-normal; moderate pooling
  z_raw ~ std_normal();
  logA ~ normal(logA_center, 1.5);         // anchor source amplitude
  pulse_mu ~ lognormal(log(4), 0.8);
  pulse_sd ~ lognormal(log(4), 0.9);       // allow a broad / continuous source to explain humps
  phi ~ normal(0, 5);                      // half-normal (lower=0)
  I ~ neg_binomial_2(mu, phi);
}
generated quantities {
  vector[O] theta_realized;
  array[N] int y_rep;
  vector[N] log_lik;
  for (o in 1:O) {
    real ssum = 0;
    real msum = 0;
    for (j in 1:Tn[o]) {
      int g = start[o] + j - 1;
      ssum += s[g];
      msum += mu[g];
    }
    theta_realized[o] = ssum / msum;
  }
  for (n in 1:N) {
    y_rep[n] = neg_binomial_2_rng(mu[n], phi);
    log_lik[n] = neg_binomial_2_lpmf(I[n] | mu[n], phi);
  }
}
