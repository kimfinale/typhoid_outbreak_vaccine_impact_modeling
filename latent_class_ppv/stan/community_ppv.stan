// Partially-pooled COMMUNITY-SURVEILLANCE PPV from suspected-case culture positivity.
//   k_o ~ Binomial(n_o, pi_o * Se_BC)        (Sp_BC ~ 1)
//   logit(pi_o) = mu + sigma * z_o           (partial pooling across settings)
// Se_BC is NOT identifiable from positivity alone (product pi*Se) -> pinned by an
// INFORMATIVE prior (community/mild value). Consequence: the LEVEL of pi tracks the
// Se_BC prior; the SPREAD sigma (setting-to-setting variation) is data-driven.
data {
  int<lower=1> O;
  array[O] int<lower=1> n;
  array[O] int<lower=0> k;
  real se_a; real se_b;            // Se_BC ~ Beta (community, informative)
}
parameters {
  real mu;                          // mean logit community PPV
  real<lower=0> sigma;              // between-setting SD (logit)
  vector[O] z;                      // non-centred
  real<lower=0, upper=1> Se_BC;     // community blood-culture sensitivity (prior-pinned)
}
transformed parameters {
  vector<lower=0, upper=1>[O] pi = inv_logit(mu + sigma * z);
}
model {
  mu ~ normal(0, 1.5);
  sigma ~ normal(0, 1);             // half-normal
  z ~ std_normal();
  Se_BC ~ beta(se_a, se_b);
  for (o in 1:O) k[o] ~ binomial(n[o], pi[o] * Se_BC);
}
generated quantities {
  real pi_mean = inv_logit(mu);                          // typical community PPV
  real pi_new  = inv_logit(mu + sigma * normal_rng(0, 1)); // predictive PPV, new outbreak (transferable)
}
