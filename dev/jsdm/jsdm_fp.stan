// Prototype (dev/jsdm/): multispecies false-positive occupancy with latent
// factors on occupancy, outside the package. Data are aggregated by unit
// (site x species): y detections in n automated occasions, v detections in
// nv verified occasions (no false positives).
//
//   logit psi_ij = beta0_j + beta1_j x_i + sigma_L * (Lambda z_i)_j
//   occupied:   y ~ Bin(n, p11_j [* site det effect]),  v ~ Bin(nv, pv_j)
//   unoccupied: y ~ Bin(n, p10_j),                       v = 0
//   logit p10_j = logit p11_j - gap_j, gap_j > 0 (labels cannot switch)
//
// Lambda is J x K with zeros above the diagonal. diag_pos = 1 constrains its
// diagonal to be positive; diag_pos = 0 leaves it free and fixes the sign of
// each factor after sampling, draw by draw, as jsdmstan does.
// det_re = 1 adds a site effect on the logit of p11 shared by all species
// (correlated detection), with standard deviation sigma_d.
data {
  int<lower=1> S;                 // sites
  int<lower=1> J;                 // species
  int<lower=0> K;                 // latent factors
  int<lower=0, upper=1> diag_pos;
  int<lower=0, upper=1> det_re;
  vector[S] x;
  array[S, J] int<lower=0> n;
  array[S, J] int<lower=0> y;
  array[S, J] int<lower=0> nv;
  array[S, J] int<lower=0> v;
}
transformed data {
  int M = K * J - (K * (K - 1)) %/% 2;   // free loadings (lower triangle)
}
parameters {
  real mu_b0;
  real<lower=0> sd_b0;
  vector[J] z_b0;
  real mu_b1;
  real<lower=0> sd_b1;
  vector[J] z_b1;
  real mu_a;
  real<lower=0> sd_a;
  vector[J] z_a;
  real mu_lg;
  real<lower=0> sd_lg;
  vector[J] z_lg;
  vector[J] av;                   // logit of verified detection
  vector[M] L;
  real<lower=0> sigma_L;
  matrix[K, S] Z;
  vector[det_re ? S : 0] z_d;
  array[det_re] real<lower=0> sigma_d;
}
transformed parameters {
  matrix[J, K] Lambda_raw = rep_matrix(0, J, K);
  {
    int idx = 0;
    for (k in 1:K) {
      for (j in k:J) {
        idx += 1;
        Lambda_raw[j, k] = (j == k && diag_pos == 1) ? exp(L[idx]) : L[idx];
      }
    }
  }
}
model {
  vector[J] b0 = mu_b0 + sd_b0 * z_b0;
  vector[J] b1 = mu_b1 + sd_b1 * z_b1;
  vector[J] a = mu_a + sd_a * z_a;
  vector[J] gap = exp(mu_lg + sd_lg * z_lg);
  matrix[S, J] lf = rep_matrix(0, S, J);
  if (K > 0) lf = (sigma_L * Lambda_raw * Z)';

  mu_b0 ~ normal(0, 1.5);  sd_b0 ~ normal(0, 1);  z_b0 ~ std_normal();
  mu_b1 ~ normal(0, 1);    sd_b1 ~ normal(0, 1);  z_b1 ~ std_normal();
  mu_a ~ normal(0, 1.5);   sd_a ~ normal(0, 1);   z_a ~ std_normal();
  mu_lg ~ normal(log(2), 0.5); sd_lg ~ normal(0, 0.5); z_lg ~ std_normal();
  av ~ normal(0, 1.5);
  // with diag_pos, L of a diagonal entry is the log of the loading
  L ~ std_normal();
  sigma_L ~ normal(0, 1);
  to_vector(Z) ~ std_normal();
  if (det_re) { z_d ~ std_normal(); sigma_d ~ normal(0, 1); }

  for (i in 1:S) {
    real d = det_re ? sigma_d[1] * z_d[i] : 0;
    for (j in 1:J) {
      real eta = b0[j] + b1[j] * x[i] + lf[i, j];
      real lp1 = log_inv_logit(eta)
                 + binomial_logit_lpmf(y[i, j] | n[i, j], a[j] + d)
                 + binomial_logit_lpmf(v[i, j] | nv[i, j], av[j]);
      real lp0 = log1m_inv_logit(eta)
                 + binomial_logit_lpmf(y[i, j] | n[i, j], a[j] + d - gap[j]);
      if (v[i, j] > 0) target += lp1;
      else target += log_sum_exp(lp1, lp0);
    }
  }
}
generated quantities {
  // sign of each factor fixed by its diagonal loading; the products
  // Lambda Lambda' do not depend on it
  matrix[J, K] Lambda = sigma_L * Lambda_raw;
  for (k in 1:K) {
    if (Lambda[k, k] < 0) Lambda[, k] = -Lambda[, k];
  }
  matrix[J, J] LL = Lambda * Lambda';
}
