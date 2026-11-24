// Royle-Nichols occupancy with false positives and a site effect on the
// per-individual detection r (prototype for step 1, dev/roadmap-v2.md).
// Occasions are exchangeable, so each site enters through its number of
// automated detections k (out of nocc) and verified detections kv (out of nv).
// N is marginalized over 0..Nmax, and the site effect e_i ~ N(0, sigma_r) by
// Gauss-Hermite quadrature (nodes gh_x, log weights gh_logw for a standard
// normal), so no site-level parameter is sampled.
functions {
  // binomial log-pmf from the log of the failure probability, stable when
  // the success probability is close to 0 or to 1 (the constant is dropped)
  real binom_logmiss(int k, int n, real log_miss) {
    real out = (n - k) * log_miss;
    if (k > 0) out += k * log1m_exp(log_miss);
    return out;
  }
}
data {
  int<lower=1> nsite;
  int<lower=1> nocc;
  vector[nsite] x;
  array[nsite] int<lower=0, upper=nocc> k;
  array[nsite] int<lower=0> nv;
  array[nsite] int<lower=0> kv;
  int<lower=1> Nmax;
  int<lower=0, upper=1> site_effect;   // 0 fixes sigma_r at 0
  int<lower=0, upper=1> lam_effect;    // 1 adds a site effect on log lambda
  int<lower=1> nq;
  vector[nq] gh_x;
  vector[nq] gh_logw;
}
transformed data {
  vector[Nmax + 1] Nvec = linspaced_vector(Nmax + 1, 0, Nmax);
}
parameters {
  real b0;
  real b1;
  real a0;
  array[site_effect] real<lower=0> sigma_r;
  array[lam_effect] real<lower=0> sigma_lam;
  vector[lam_effect ? nsite : 0] u_raw;
  real logit_p10;
  real logit_rv;
}
model {
  real sr = site_effect ? sigma_r[1] : 0;
  b0 ~ normal(0, 2);
  b1 ~ normal(0, 1);
  a0 ~ normal(0, 1.5);
  sigma_r ~ normal(0, 1);
  sigma_lam ~ normal(0, 1);
  u_raw ~ std_normal();
  logit_p10 ~ normal(-3, 1);
  logit_rv ~ normal(0, 1.5);
  {
    real log1m_p10 = log1m_inv_logit(logit_p10);
    real log1m_rv = log1m_inv_logit(logit_rv);
    for (i in 1:nsite) {
      real lam = exp(b0 + b1 * x[i] + (lam_effect ? sigma_lam[1] * u_raw[i] : 0));
      vector[Nmax + 1] lp_n;
      vector[nq] lq;
      // terms that do not depend on e
      for (n in 0:Nmax) {
        lp_n[n + 1] = poisson_lpmf(n | lam);
        if (nv[i] > 0) {
          if (n == 0) {
            lp_n[n + 1] += kv[i] == 0 ? 0 : negative_infinity();
          } else {
            lp_n[n + 1] += binom_logmiss(kv[i], nv[i], n * log1m_rv);
          }
        }
      }
      for (q in 1:nq) {
        real log1m_r = log1m_inv_logit(a0 + sr * gh_x[q]);
        vector[Nmax + 1] lp = lp_n;
        for (n in 0:Nmax) {
          // P(no automated detection) = (1 - r)^N (1 - p10)
          real log_miss = n * log1m_r + log1m_p10;
          lp[n + 1] += binom_logmiss(k[i], nocc, log_miss);
        }
        lq[q] = gh_logw[q] + log_sum_exp(lp);
      }
      target += log_sum_exp(lq);
    }
  }
}
generated quantities {
  real p10 = inv_logit(logit_p10);
  real rv = inv_logit(logit_rv);
  real sigma_r_out = site_effect ? sigma_r[1] : 0;
  real sigma_lam_out = lam_effect ? sigma_lam[1] : 0;
}
