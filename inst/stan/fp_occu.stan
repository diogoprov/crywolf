// ============================================================================
// Bayesian false-positive occupancy model, multiple-methods (site-confirmation)
// design of Miller et al. (2011): an automated method that can return false
// positives, plus verified detections that cannot.
//
//   z_i        ~ Bernoulli(psi_i),        logit(psi_i)  = X_i beta
//   y2_ik      ~ Bernoulli(z_i p11_ik + (1 - z_i) p10_i)
//                logit(p11_ik) = Z_ik alpha + u_rec[ik]
//                logit(p10_i)  = (alpha_1 - fp_gap) + Xfp_i gamma_fp
//   y1_im      ~ Bernoulli(z_i r11),      logit(r11) = delta0
//
// The latent z is marginalized. fp_gap > 0 keeps p11 > p10 at the reference
// (all covariates at 0), which resolves label switching. Column 1 of Z must be
// the intercept. Kfp = 0 gives the constant-rate model; estimate_fp = 0 gives
// the standard occupancy model (no false positives).
//
// Automated occasions enter as rows of identical covariates within a site
// (n2 trials, y2 detections), which is exact for a product of Bernoulli
// trials and much faster when detection covariates are constant within
// sites. Rows must be sorted by site (fp_occu() does both), so that each
// site's rows are contiguous and can be summed with segment().
// ============================================================================
functions {
  // Marginal log-likelihood of each site, with the latent state z summed out.
  vector site_loglik(vector lpsi, vector lp11_fixed, real gamma0, matrix Xfp,
                     vector gamma_fp, vector u, array[] int rec2,
                     int estimate_re, int estimate_fp, int Kfp,
                     vector y2v, vector f2v, array[] int site2, array[] int start2,
                     array[] int len2, array[] int any_det2,
                     array[] int y1, array[] int site1, real delta0,
                     array[] int known_present) {
    int nsite = rows(lpsi);
    int N2 = rows(y2v);
    vector[N2] lp11 = lp11_fixed;
    if (estimate_re == 1) lp11 += u[rec2];
    vector[nsite] lp10 = rep_vector(gamma0, nsite);
    if (Kfp > 0) lp10 += Xfp * gamma_fp;

    // per-occasion log-likelihoods, vectorized
    vector[N2] l11 = y2v .* log_inv_logit(lp11) + f2v .* log1m_inv_logit(lp11);
    vector[N2] l10;
    if (estimate_fp == 1) {
      vector[N2] lp10_n = lp10[site2];
      l10 = y2v .* log_inv_logit(lp10_n) + f2v .* log1m_inv_logit(lp10_n);
    }

    vector[nsite] occ = rep_vector(0, nsite);
    vector[nsite] unocc = rep_vector(0, nsite);
    for (i in 1:nsite) {
      if (len2[i] > 0) {
        occ[i] = sum(segment(l11, start2[i], len2[i]));
        if (estimate_fp == 1) unocc[i] = sum(segment(l10, start2[i], len2[i]));
      }
      if (estimate_fp == 0 && any_det2[i] == 1) unocc[i] = negative_infinity();
    }
    for (n in 1:size(y1))
      occ[site1[n]] += bernoulli_logit_lpmf(y1[n] | delta0);

    vector[nsite] ll;
    for (i in 1:nsite) {
      real occ_branch = log_inv_logit(lpsi[i]) + occ[i];
      ll[i] = known_present[i] == 1
              ? occ_branch
              : log_sum_exp(occ_branch, log1m_inv_logit(lpsi[i]) + unocc[i]);
    }
    return ll;
  }
}

data {
  int<lower=1> nsite;
  int<lower=1> Kbeta;
  matrix[nsite, Kbeta] X;                    // occupancy design (with intercept)

  int<lower=0> N2;                           // rows of automated occasions
  array[N2] int<lower=1> n2;                 // occasions in each row
  array[N2] int<lower=0> y2;                 // detections in each row
  array[N2] int<lower=1, upper=nsite> site2;  // sorted by site
  array[nsite] int<lower=1> start2;          // first row of each site in y2
  array[nsite] int<lower=0> len2;            // number of rows of each site
  int<lower=1> Kalpha;
  matrix[N2, Kalpha] Z;                      // detection design (column 1 = intercept)
  int<lower=1> nrec;
  array[N2] int<lower=1, upper=nrec> rec2;
  int<lower=0, upper=1> estimate_re;         // recorder random intercept on p11

  int<lower=0> Kfp;
  matrix[nsite, Kfp] Xfp;                    // false-positive design (no intercept)
  int<lower=0, upper=1> estimate_fp;

  int<lower=0> N1;                           // verified occasions (long format)
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=nsite> site1;
  array[nsite] int<lower=0, upper=1> known_present;

  real<lower=0> prior_sd_int;                // sd of Normal priors on intercepts
  real<lower=0> prior_sd_slope;              // sd of Normal priors on slopes
  real gap_mean;                             // Normal(gap_mean, gap_sd) prior on fp_gap, truncated at 0
  real<lower=0> gap_sd;
  real<lower=0> recorder_sd;                 // half-Normal(0, recorder_sd) on sigma_u
}

transformed data {
  vector[N2] y2v = to_vector(y2);
  vector[N2] f2v = to_vector(n2) - y2v;      // non-detections
  // sites with an automated detection cannot be unoccupied without false positives
  array[nsite] int any_det2 = rep_array(0, nsite);
  for (n in 1:N2) {
    if (y2[n] > n2[n]) reject("y2 cannot exceed n2");
    if (y2[n] > 0) any_det2[site2[n]] = 1;
  }
}

parameters {
  vector[Kbeta] beta;
  vector[Kalpha] alpha;
  vector[nrec] u_raw;
  real<lower=0> sigma_u;
  real<lower=0> fp_gap;
  vector[Kfp] gamma_fp;
  real delta0;
}

transformed parameters {
  vector[nrec] u = sigma_u * u_raw;
  real gamma0 = alpha[1] - fp_gap;
}

model {
  beta[1]        ~ normal(0, prior_sd_int);
  if (Kbeta > 1)  beta[2:Kbeta]   ~ normal(0, prior_sd_slope);
  alpha[1]       ~ normal(0, prior_sd_int);
  if (Kalpha > 1) alpha[2:Kalpha] ~ normal(0, prior_sd_slope);
  u_raw          ~ std_normal();
  sigma_u        ~ normal(0, recorder_sd);
  fp_gap         ~ normal(gap_mean, gap_sd);
  gamma_fp       ~ normal(0, prior_sd_slope);
  delta0         ~ normal(0, prior_sd_int);

  {
    vector[nsite] ll_site = site_loglik(X * beta, Z * alpha, gamma0, Xfp, gamma_fp,
                                        u, rec2, estimate_re, estimate_fp, Kfp,
                                        y2v, f2v, site2, start2, len2, any_det2,
                                        y1, site1, delta0, known_present);
    target += sum(ll_site);
  }
}

generated quantities {
  vector[nsite] log_lik;
  real p11_ref = inv_logit(alpha[1]);
  real p10_ref = estimate_fp == 1 ? inv_logit(gamma0) : 0.0;
  real p10_bar = 0;
  real r11 = inv_logit(delta0);
  {
    vector[nsite] lp10 = rep_vector(gamma0, nsite);
    if (Kfp > 0) lp10 += Xfp * gamma_fp;
    if (estimate_fp == 1) p10_bar = mean(inv_logit(lp10));
    log_lik = site_loglik(X * beta, Z * alpha, gamma0, Xfp, gamma_fp,
                          u, rec2, estimate_re, estimate_fp, Kfp,
                          y2v, f2v, site2, start2, len2, any_det2,
                          y1, site1, delta0, known_present);
  }
}
