// Community/surveillance PPV across endemicity strata (extends community_ppv.stan).
//   k_o ~ Binomial(n_o, pi_o * Se_BC)                 (Sp ~ 1; positivity among tested)
//   logit(pi_o) = mu + beta * outbreak_o + sigma * z_o
// mu           = mean logit PPV in the ENDEMIC-SURVEILLANCE stratum (outbreak=0)
// mu + beta    = mean logit PPV in the OUTBREAK stratum (outbreak=1)
// beta         = the endemicity gradient (log-odds); exp(beta) = OR outbreak-vs-surveillance
// Se_BC is not identifiable from positivity alone (product pi*Se) -> pinned by prior;
// its LEVEL scales pi, its shared value cancels in the OR (gradient is Se-robust).
data {
  int<lower=1> O;
  array[O] int<lower=1> n;               // blood-cultured (tested)
  array[O] int<lower=0> k;               // culture-confirmed
  vector<lower=0, upper=1>[O] outbreak;  // 1 = outbreak stratum, 0 = endemic surveillance
  real se_a; real se_b;                  // Se_BC ~ Beta (informative)
}
parameters {
  real mu;
  real beta;
  real<lower=0> sigma;
  vector[O] z;
  real<lower=0, upper=1> Se_BC;
}
transformed parameters {
  vector<lower=0, upper=1>[O] pi = inv_logit(mu + beta * outbreak + sigma * z);
}
model {
  mu    ~ normal(0, 1.5);
  beta  ~ normal(0, 1.5);
  sigma ~ normal(0, 1);                   // half-normal
  z     ~ std_normal();
  Se_BC ~ beta(se_a, se_b);
  for (o in 1:O) k[o] ~ binomial(n[o], pi[o] * Se_BC);
}
generated quantities {
  real pi_surv = inv_logit(mu);
  real pi_out  = inv_logit(mu + beta);
  real OR_out_vs_surv = exp(beta);
  real pi_new_out  = inv_logit(mu + beta + sigma * normal_rng(0, 1));  // predictive, new outbreak
  real pi_new_surv = inv_logit(mu +        sigma * normal_rng(0, 1));  // predictive, new surveillance
}
