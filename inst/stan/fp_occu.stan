// ============================================================================
// Bayesian false-positive occupancy model, multiple-methods (site-confirmation)
// design of Miller et al. (2011): an automated method that can return false
// positives, plus verified detections that cannot.
//
//   z_i        ~ Bernoulli(psi_i),        logit(psi_i)  = X_i beta
//   y2_ik      ~ Bernoulli(z_i p11_ik + (1 - z_i) p10_i)
//                logit(p11_ik) = Z_ik alpha + u_rec[ik]
//                logit(p10_i)  = (alpha_1 - fp_gap) + Xfp_i gamma_fp [+ e_i]
//   y1_im      ~ Bernoulli(z_i r11),      logit(r11) = delta0
//
// The latent z is marginalized. fp_gap > 0 keeps p11 > p10 at the reference
// (all covariates at 0), which resolves label switching. Column 1 of Z must be
// the intercept. Kfp = 0 gives the constant-rate model; estimate_fp = 0 gives
// the standard occupancy model (no false positives). estimate_fp_re = 1 adds a
// site random effect e_i ~ Normal(0, sigma_fp) to the logit of p10.
//
// Observation-confirmation design (Chambert et al. 2015), when obs_confirm = 1:
// the observations of some automated occasions are checked, giving a state
// v3 = 0 (nothing), 1 (false observations only), 2 (true only), 3 (both). The
// detection and false-positive predictors are then those of s1 and s0, the
// probabilities of at least one true and at least one false observation on an
// occasion, and p11 = 1 - (1 - s1)(1 - s0), p10 = s0 for unchecked occasions.
//
// Automated occasions enter as rows of identical covariates within a site
// (n2 trials, y2 detections), which is exact for a product of Bernoulli
// trials and much faster when detection covariates are constant within
// sites. Rows must be sorted by site (fp_occu() does both), so that each
// site's rows are contiguous and can be summed with segment().
//
// The sum over sites runs with reduce_sum(), which splits the sites among
// threads when the model is compiled with stan_threads (threads_per_chain in
// fp_occu()). Each slice of sites only touches its own rows of the automated
// and confirmed occasions. With one thread the result is the same sum.
// ============================================================================
functions {
  // Marginal log-likelihood of sites a to b, with the latent state z summed
  // out. The automated (and confirmed) rows of these sites are contiguous.
  vector site_loglik(int a, int b, matrix X, vector beta, matrix Z, vector alpha,
                     real gamma0, matrix Xfp, vector gamma_fp, vector e_fp,
                     vector u, array[] int rec2, int estimate_re,
                     int estimate_fp, int Kfp,
                     vector y2v, vector n2v, array[] int site2,
                     array[] int start2, array[] int len2, array[] int any_det2,
                     vector n1v, vector s1v, real delta0,
                     array[] int known_present, int obs_confirm,
                     array[] int v3, array[] int site3, matrix Z3,
                     array[] int rec3, array[] int start3, array[] int len3) {
    int ns = b - a + 1;
    vector[ns] lpsi = X[a:b] * beta;
    vector[ns] lp10 = rep_vector(gamma0, ns);
    if (Kfp > 0) lp10 += Xfp[a:b] * gamma_fp;
    lp10 += e_fp[a:b];
    vector[ns] occ = rep_vector(0, ns);
    vector[ns] unocc = rep_vector(0, ns);

    // automated occasions, as binomial rows; log_inv_logit(x) = x - log1p_exp(x)
    // and log1m_inv_logit(x) = -log1p_exp(x)
    int r0 = start2[a];
    int nr = start2[b] + len2[b] - r0;
    if (nr > 0) {
      int r1 = r0 + nr - 1;
      vector[nr] y = y2v[r0:r1];
      vector[nr] n = n2v[r0:r1];
      vector[nr] lp11 = Z[r0:r1] * alpha;
      if (estimate_re == 1) lp11 += u[rec2[r0:r1]];
      array[nr] int loc;
      for (k in 1:nr) loc[k] = site2[r0 + k - 1] - a + 1;
      vector[nr] l11;
      vector[nr] l10;
      if (estimate_fp == 1) {
        vector[nr] lp10_n = lp10[loc];
        l10 = y .* lp10_n - n .* log1p_exp(lp10_n);
        if (obs_confirm == 1) {
          // log(1 - p11) = log(1 - s1) + log(1 - s0)
          vector[nr] lq11 = -log1p_exp(lp11) - log1p_exp(lp10_n);
          l11 = y .* log1m_exp(lq11) + (n - y) .* lq11;
        }
      }
      if (obs_confirm == 0) l11 = y .* lp11 - n .* log1p_exp(lp11);
      for (i in 1:ns) {
        int s = a + i - 1;
        if (len2[s] > 0) {
          occ[i] = sum(segment(l11, start2[s] - r0 + 1, len2[s]));
          if (estimate_fp == 1) unocc[i] = sum(segment(l10, start2[s] - r0 + 1, len2[s]));
        }
      }
    }
    if (estimate_fp == 0) {
      for (i in 1:ns) if (any_det2[a + i - 1] == 1) unocc[i] = negative_infinity();
    }

    // verified occasions, as counts per site
    occ += s1v[a:b] * log_inv_logit(delta0) + (n1v[a:b] - s1v[a:b]) * log1m_inv_logit(delta0);

    // confirmed occasions (observation-confirmation design)
    if (obs_confirm == 1) {
      int q0 = start3[a];
      int nq = start3[b] + len3[b] - q0;
      if (nq > 0) {
        int q1 = q0 + nq - 1;
        vector[nq] lp11_3 = Z3[q0:q1] * alpha;
        if (estimate_re == 1) lp11_3 += u[rec3[q0:q1]];
        for (k in 1:nq) {
          int i = site3[q0 + k - 1] - a + 1;
          int v = v3[q0 + k - 1];
          real ls1 = log_inv_logit(lp11_3[k]);
          real lq1 = log1m_inv_logit(lp11_3[k]);
          real ls0 = log_inv_logit(lp10[i]);
          real lq0 = log1m_inv_logit(lp10[i]);
          occ[i] += (v >= 2 ? ls1 : lq1) + (v % 2 == 1 ? ls0 : lq0);
          unocc[i] += v >= 2 ? negative_infinity() : (v == 1 ? ls0 : lq0);
        }
      }
    }

    vector[ns] ll;
    for (i in 1:ns) {
      real occ_branch = log_inv_logit(lpsi[i]) + occ[i];
      ll[i] = known_present[a + i - 1] == 1
              ? occ_branch
              : log_sum_exp(occ_branch, log1m_inv_logit(lpsi[i]) + unocc[i]);
    }
    return ll;
  }

  // partial sum over a slice of sites, for reduce_sum()
  real partial_loglik(array[] int slice, int a, int b, matrix X, vector beta,
                      matrix Z, vector alpha, real gamma0, matrix Xfp,
                      vector gamma_fp, vector e_fp, vector u, array[] int rec2,
                      int estimate_re, int estimate_fp, int Kfp,
                      vector y2v, vector n2v, array[] int site2,
                      array[] int start2, array[] int len2, array[] int any_det2,
                      vector n1v, vector s1v, real delta0,
                      array[] int known_present, int obs_confirm,
                      array[] int v3, array[] int site3, matrix Z3,
                      array[] int rec3, array[] int start3, array[] int len3) {
    return sum(site_loglik(a, b, X, beta, Z, alpha, gamma0, Xfp, gamma_fp, e_fp, u,
                           rec2, estimate_re, estimate_fp, Kfp, y2v, n2v, site2,
                           start2, len2, any_det2, n1v, s1v, delta0,
                           known_present, obs_confirm, v3, site3, Z3, rec3,
                           start3, len3));
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

  int<lower=0, upper=1> obs_confirm;         // observation-confirmation design
  int<lower=0> N3;                           // confirmed occasions
  array[N3] int<lower=0, upper=3> v3;        // confirmation state
  array[N3] int<lower=1, upper=nsite> site3;
  matrix[N3, Kalpha] Z3;                     // detection design of confirmed occasions
  array[N3] int<lower=1, upper=nrec> rec3;

  real<lower=0> prior_sd_int;                // sd of Normal priors on intercepts
  real<lower=0> prior_sd_slope;              // sd of Normal priors on slopes
  real gap_mean;                             // Normal(gap_mean, gap_sd) prior on fp_gap, truncated at 0
  real<lower=0> gap_sd;
  real<lower=0> recorder_sd;                 // half-Normal(0, recorder_sd) on sigma_u

  int<lower=0, upper=1> estimate_fp_re;      // site random effect on logit(p10)
  real<lower=0> fp_site_sd;                  // half-Normal(0, fp_site_sd) on sigma_fp
}

transformed data {
  vector[N2] y2v = to_vector(y2);
  vector[N2] n2v = to_vector(n2);
  // sites with an automated detection cannot be unoccupied without false positives
  array[nsite] int any_det2 = rep_array(0, nsite);
  for (n in 1:N2) {
    if (y2[n] > n2[n]) reject("y2 cannot exceed n2");
    if (y2[n] > 0) any_det2[site2[n]] = 1;
  }
  // verified occasions and detections per site
  vector[nsite] n1v = rep_vector(0, nsite);
  vector[nsite] s1v = rep_vector(0, nsite);
  for (n in 1:N1) {
    n1v[site1[n]] += 1;
    s1v[site1[n]] += y1[n];
  }
  // first row and number of rows of each site among the confirmed occasions
  array[nsite] int len3 = rep_array(0, nsite);
  for (n in 1:N3) {
    if (n > 1 && site3[n] < site3[n - 1]) reject("confirmed occasions must be sorted by site");
    len3[site3[n]] += 1;
  }
  array[nsite] int start3;
  start3[1] = 1;
  for (i in 2:nsite) start3[i] = start3[i - 1] + len3[i - 1];
  // sites, sliced among threads by reduce_sum()
  array[nsite] int site_index = linspaced_int_array(nsite, 1, nsite);
}

parameters {
  vector[Kbeta] beta;
  vector[Kalpha] alpha;
  vector[nrec] u_raw;
  real<lower=0> sigma_u;
  real<lower=0> fp_gap;
  vector[Kfp] gamma_fp;
  real delta0;
  vector[estimate_fp_re == 1 ? nsite : 0] e_raw;
  array[estimate_fp_re] real<lower=0> sigma_fp;
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
  e_raw          ~ std_normal();
  sigma_fp       ~ normal(0, fp_site_sd);

  // site effects on the false-positive rate (zero when estimate_fp_re = 0)
  vector[nsite] e_fp = rep_vector(0, nsite);
  if (estimate_fp_re == 1) e_fp = sigma_fp[1] * e_raw;
  target += reduce_sum(partial_loglik, site_index, 1, X, beta, Z, alpha, gamma0,
                       Xfp, gamma_fp, e_fp, u, rec2, estimate_re, estimate_fp, Kfp,
                       y2v, n2v, site2, start2, len2, any_det2, n1v, s1v, delta0,
                       known_present, obs_confirm, v3, site3, Z3, rec3, start3,
                       len3);
}

generated quantities {
  vector[nsite] log_lik;
  real p11_ref = obs_confirm == 1
                 ? 1 - (1 - inv_logit(alpha[1])) * (1 - inv_logit(gamma0))
                 : inv_logit(alpha[1]);
  real p10_ref = estimate_fp == 1 ? inv_logit(gamma0) : 0.0;
  real p10_bar = 0;
  real r11 = inv_logit(delta0);
  {
    vector[nsite] e_fp = rep_vector(0, nsite);
    if (estimate_fp_re == 1) e_fp = sigma_fp[1] * e_raw;
    vector[nsite] lp10 = rep_vector(gamma0, nsite);
    if (Kfp > 0) lp10 += Xfp * gamma_fp;
    lp10 += e_fp;
    if (estimate_fp == 1) p10_bar = mean(inv_logit(lp10));
    log_lik = site_loglik(1, nsite, X, beta, Z, alpha, gamma0, Xfp, gamma_fp, e_fp, u,
                          rec2, estimate_re, estimate_fp, Kfp, y2v, n2v, site2,
                          start2, len2, any_det2, n1v, s1v, delta0,
                          known_present, obs_confirm, v3, site3, Z3, rec3,
                          start3, len3);
  }
}
