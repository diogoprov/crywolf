// ============================================================================
// Royle-Nichols occupancy model with false positives (fp_occu_rn()).
//
//   N_i        ~ Poisson(lambda_i),        log(lambda_i) = X_i beta [+ v_i]
//   y2_ik      ~ Bernoulli(1 - (1 - r_ik)^N_i (1 - p10_i))
//                logit(r_ik)  = Z_ik alpha + u_rec[ik] [+ d_i]
//                logit(p10_i) = gamma0 + Xfp_i gamma_fp [+ e_i]
//   y1_im      ~ Bernoulli(1 - (1 - rv_im)^N_i),
//                logit(rv_im) = delta0 + X1s_im delta_s
//
// r is the detection probability of one individual on one occasion, rv that
// of the verified method; the site is occupied when N_i > 0, so
// psi_i = 1 - exp(-lambda_i). False positives arise on every occasion, also
// at occupied sites (with N_i = 0 the automated detection probability is
// p10_i). d_i ~ Normal(0, sigma_det) is a site effect on the per-individual
// detection (estimate_det_re = 1), v_i ~ Normal(0, sigma_lam) one on log
// abundance (estimate_abund_re = 1), e_i ~ Normal(0, sigma_fp) one on the
// logit of p10 (estimate_fp_re = 1). estimate_fp = 0 gives the
// Royle-Nichols model without false positives. Unlike fp_occu.stan, the
// false-positive intercept is not tied to the detection intercept by a gap:
// detection at an occupied site, 1 - (1 - r)^N (1 - p10), always exceeds
// p10, so the labels cannot switch, and p10 may exceed the per-individual r.
// Its prior is the marginal prior the gap implies in fp_occu.stan,
// Normal(-gap_mean, sqrt(prior_sd_int^2 + gap_sd^2)), without truncation.
//
// Observation-confirmation design (Chambert et al. 2015), when obs_confirm =
// 1: the observations of confirmed occasions give a state v3 = 0 (nothing),
// 1 (false observations only), 2 (true only), 3 (both), with probabilities of
// at least one true observation s1 = 1 - (1 - r)^N and at least one false
// observation s0 = p10.
//
// N is summed out over 0..Nmax. The automated rows of a site are sorted by
// site (fp_occu_rn() builds them with the data of fp_occu()); rows with
// identical covariates within a site are collapsed into binomial rows (n2
// occasions, y2 detections), and confirmed occasions with identical
// covariates and state into rows of w3 occasions. Every term for N >= 1 is
// written from log(1 - p), so that detection probabilities close to 0 or 1
// keep finite gradients. The sum over sites runs with reduce_sum().
// ============================================================================
functions {
  // log-likelihood of y detections in n trials given log(1 - p) for each
  // value of N (the binomial constant is dropped)
  vector binom_logmiss(real y, real n, vector log_miss) {
    vector[rows(log_miss)] out = (n - y) * log_miss;
    if (y > 0) out += y * log1m_exp(log_miss);
    return out;
  }

  // marginal log-likelihood of sites a to b
  vector site_loglik(int a, int b, vector beta, vector alpha, real gamma0,
                     vector gamma_fp, vector e_fp, vector u, vector d_det,
                     vector v_ab, real delta0, vector delta_s,
                     matrix X, matrix Z, matrix Xfp, matrix X1s,
                     array[] int y2, array[] int n2, array[] int rec2,
                     array[] int start2, array[] int len2,
                     array[] int any_det2, vector Ysite, vector Nsite,
                     array[] int y1, array[] int start1, array[] int len1,
                     vector n1v, vector s1v, array[] int known_present,
                     array[] int v3, array[] int w3, array[] int rec3, matrix Z3,
                     array[] int start3, array[] int len3,
                     vector Nvec, vector lgamma_n1, array[] int options) {
    int estimate_re = options[1];
    int estimate_fp = options[2];
    int Kfp = options[3];
    int obs_confirm = options[4];
    int Kver_s = options[5];
    int ns = b - a + 1;
    int Nmax = rows(Nvec);
    vector[ns] ll;
    for (i in 1:ns) {
      int s = a + i - 1;
      real lp10 = gamma0 + e_fp[s];
      if (Kfp > 0) lp10 += Xfp[s] * gamma_fp;
      real ls0 = log_inv_logit(lp10);
      real lq0 = log1m_inv_logit(lp10);
      real l10 = estimate_fp == 1 ? lq0 : 0;
      // N = 0: false positives only
      real em0 = 0;
      if (estimate_fp == 1) {
        em0 = Ysite[s] * lp10 - Nsite[s] * log1p_exp(lp10);
      } else if (any_det2[s] == 1) {
        em0 = negative_infinity();
      }
      if (s1v[s] > 0 || known_present[s] == 1) em0 = negative_infinity();
      // N = 1..Nmax
      vector[Nmax] em = rep_vector(0, Nmax);
      for (q in start2[s]:(start2[s] + len2[s] - 1)) {
        real eta = Z[q] * alpha + d_det[s];
        if (estimate_re == 1) eta += u[rec2[q]];
        // log P(no detection on one occasion | N) = N log(1 - r) + log(1 - p10)
        em += binom_logmiss(y2[q], n2[q], Nvec * log1m_inv_logit(eta) + l10);
      }
      if (Kver_s == 0) {
        if (n1v[s] > 0) {
          em += binom_logmiss(s1v[s], n1v[s], Nvec * log1m_inv_logit(delta0));
        }
      } else {
        for (w in start1[s]:(start1[s] + len1[s] - 1)) {
          real eta1 = delta0 + X1s[w] * delta_s;
          em += binom_logmiss(y1[w], 1, Nvec * log1m_inv_logit(eta1));
        }
      }
      if (obs_confirm == 1) {
        for (q in start3[s]:(start3[s] + len3[s] - 1)) {
          real eta3 = Z3[q] * alpha + d_det[s];
          if (estimate_re == 1) eta3 += u[rec3[q]];
          vector[Nmax] lq1 = Nvec * log1m_inv_logit(eta3);
          int v = v3[q];
          em += w3[q] * ((v >= 2 ? log1m_exp(lq1) : lq1) + (v % 2 == 1 ? ls0 : lq0));
          em0 += v >= 2 ? negative_infinity() : w3[q] * (v == 1 ? ls0 : lq0);
        }
      }
      real llam = X[s] * beta + v_ab[s];
      real lp_0 = -exp(llam);
      vector[Nmax] lp_n = Nvec * llam - exp(llam) - lgamma_n1;
      ll[i] = log_sum_exp(lp_0 + em0, log_sum_exp(lp_n + em));
    }
    return ll;
  }

  real partial_loglik(array[] int slice, int a, int b, vector beta, vector alpha,
                      real gamma0, vector gamma_fp, vector e_fp, vector u,
                      vector d_det, vector v_ab, real delta0, vector delta_s,
                      matrix X, matrix Z, matrix Xfp, matrix X1s,
                      array[] int y2, array[] int n2, array[] int rec2,
                      array[] int start2, array[] int len2,
                      array[] int any_det2, vector Ysite, vector Nsite,
                      array[] int y1, array[] int start1, array[] int len1,
                      vector n1v, vector s1v, array[] int known_present,
                      array[] int v3, array[] int w3, array[] int rec3, matrix Z3,
                      array[] int start3, array[] int len3,
                      vector Nvec, vector lgamma_n1, array[] int options) {
    return sum(site_loglik(a, b, beta, alpha, gamma0, gamma_fp, e_fp, u, d_det,
                           v_ab, delta0, delta_s, X, Z, Xfp, X1s, y2, n2, rec2,
                           start2, len2, any_det2, Ysite, Nsite, y1, start1,
                           len1, n1v, s1v, known_present, v3, w3, rec3, Z3, start3,
                           len3, Nvec, lgamma_n1, options));
  }
}

data {
  int<lower=1> nsite;
  int<lower=1> Kbeta;
  matrix[nsite, Kbeta] X;                    // abundance design (log lambda, with intercept)

  int<lower=0> N2;                           // rows of automated occasions
  array[N2] int<lower=1> n2;
  array[N2] int<lower=0> y2;
  array[N2] int<lower=1, upper=nsite> site2;  // sorted by site
  array[nsite] int<lower=1> start2;
  array[nsite] int<lower=0> len2;
  int<lower=1> Kalpha;
  matrix[N2, Kalpha] Z;                      // per-individual detection design (column 1 = intercept)
  int<lower=1> nrec;
  array[N2] int<lower=1, upper=nrec> rec2;
  int<lower=0, upper=1> estimate_re;

  int<lower=0> Kfp;
  matrix[nsite, Kfp] Xfp;                    // false-positive design (site level, no intercept)
  int<lower=0, upper=1> estimate_fp;

  int<lower=0> N1;                           // verified occasions (long format, sorted by site)
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=nsite> site1;
  int<lower=0> Kver_s;
  matrix[N1, Kver_s] X1s;
  array[nsite] int<lower=0, upper=1> known_present;

  int<lower=0, upper=1> obs_confirm;
  int<lower=0> N3;
  array[N3] int<lower=0, upper=3> v3;
  array[N3] int<lower=1> w3;                 // confirmed occasions in each row
  array[N3] int<lower=1, upper=nsite> site3;
  matrix[N3, Kalpha] Z3;
  array[N3] int<lower=1, upper=nrec> rec3;

  int<lower=1> Nmax;                         // largest abundance summed over

  real<lower=0> prior_sd_int;
  real<lower=0> prior_sd_slope;
  real gap_mean;
  real<lower=0> gap_sd;
  real<lower=0> recorder_sd;
  int<lower=0, upper=1> estimate_fp_re;
  real<lower=0> fp_site_sd;
  int<lower=0, upper=1> estimate_det_re;
  real<lower=0> det_site_sd;
  int<lower=0, upper=1> estimate_abund_re;
  real<lower=0> abund_site_sd;
  int<lower=0, upper=1> prior_only;
}

transformed data {
  array[nsite] int any_det2 = rep_array(0, nsite);
  vector[nsite] Ysite = rep_vector(0, nsite);
  vector[nsite] Nsite = rep_vector(0, nsite);
  for (n in 1:N2) {
    if (y2[n] > n2[n]) reject("y2 cannot exceed n2");
    if (n > 1 && site2[n] < site2[n - 1]) reject("automated rows must be sorted by site");
    if (y2[n] > 0) any_det2[site2[n]] = 1;
    Ysite[site2[n]] += y2[n];
    Nsite[site2[n]] += n2[n];
  }
  if (obs_confirm == 1 && estimate_fp == 0) {
    reject("the observation-confirmation design needs a false-positive process (estimate_fp = 1)");
  }
  vector[nsite] n1v = rep_vector(0, nsite);
  vector[nsite] s1v = rep_vector(0, nsite);
  array[nsite] int len1 = rep_array(0, nsite);
  for (n in 1:N1) {
    if (n > 1 && site1[n] < site1[n - 1]) reject("verified occasions must be sorted by site");
    n1v[site1[n]] += 1;
    s1v[site1[n]] += y1[n];
    len1[site1[n]] += 1;
  }
  array[nsite] int start1;
  start1[1] = 1;
  for (i in 2:nsite) start1[i] = start1[i - 1] + len1[i - 1];
  array[nsite] int len3 = rep_array(0, nsite);
  for (n in 1:N3) {
    if (n > 1 && site3[n] < site3[n - 1]) reject("confirmed occasions must be sorted by site");
    len3[site3[n]] += 1;
  }
  array[nsite] int start3;
  start3[1] = 1;
  for (i in 2:nsite) start3[i] = start3[i - 1] + len3[i - 1];
  vector[Nmax] Nvec = linspaced_vector(Nmax, 1, Nmax);
  vector[Nmax] lgamma_n1;
  for (n in 1:Nmax) lgamma_n1[n] = lgamma(n + 1);
  array[nsite] int site_index = linspaced_int_array(nsite, 1, nsite);
  vector[nsite] zero_site = rep_vector(0, nsite);
  array[5] int options = {estimate_re, estimate_fp, Kfp, obs_confirm, Kver_s};
}

parameters {
  vector[Kbeta] beta;
  vector[Kalpha] alpha;
  vector[nrec] u_raw;
  real<lower=0> sigma_u;
  real gamma0;
  vector[Kfp] gamma_fp;
  real delta0;
  vector[Kver_s] delta_s;
  vector[estimate_fp_re == 1 ? nsite : 0] e_raw;
  array[estimate_fp_re] real<lower=0> sigma_fp;
  vector[estimate_det_re == 1 ? nsite : 0] d_raw;
  array[estimate_det_re] real<lower=0> sigma_det;
  vector[estimate_abund_re == 1 ? nsite : 0] v_raw;
  array[estimate_abund_re] real<lower=0> sigma_lam;
}

transformed parameters {
  vector[nrec] u = sigma_u * u_raw;
}

model {
  beta[1]        ~ normal(0, prior_sd_int);
  if (Kbeta > 1)  beta[2:Kbeta]   ~ normal(0, prior_sd_slope);
  alpha[1]       ~ normal(0, prior_sd_int);
  if (Kalpha > 1) alpha[2:Kalpha] ~ normal(0, prior_sd_slope);
  u_raw          ~ std_normal();
  sigma_u        ~ normal(0, recorder_sd);
  gamma0         ~ normal(-gap_mean, sqrt(square(prior_sd_int) + square(gap_sd)));
  gamma_fp       ~ normal(0, prior_sd_slope);
  delta0         ~ normal(0, prior_sd_int);
  delta_s        ~ normal(0, prior_sd_slope);
  e_raw          ~ std_normal();
  sigma_fp       ~ normal(0, fp_site_sd);
  d_raw          ~ std_normal();
  sigma_det      ~ normal(0, det_site_sd);
  v_raw          ~ std_normal();
  sigma_lam      ~ normal(0, abund_site_sd);

  if (prior_only == 0) {
    vector[nsite] e_fp = estimate_fp_re == 1 ? sigma_fp[1] * e_raw : zero_site;
    vector[nsite] d_det = estimate_det_re == 1 ? sigma_det[1] * d_raw : zero_site;
    vector[nsite] v_ab = estimate_abund_re == 1 ? sigma_lam[1] * v_raw : zero_site;
    target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                         e_fp, u, d_det, v_ab, delta0, delta_s, X, Z, Xfp, X1s, y2, n2,
                         rec2, start2, len2, any_det2, Ysite, Nsite, y1, start1, len1,
                         n1v, s1v, known_present, v3, w3, rec3, Z3, start3, len3, Nvec,
                         lgamma_n1, options);
  }
}

generated quantities {
  vector[nsite] log_lik;
  real lambda_ref = exp(beta[1]);
  real psi_ref = -expm1(-lambda_ref);
  real r_ref = inv_logit(alpha[1]);
  real p10_ref = estimate_fp == 1 ? inv_logit(gamma0) : 0.0;
  real r_verified = inv_logit(delta0);
  real lambda_bar;
  real psi_bar;
  // largest probability, over sites, that N exceeds Nmax (the mass the sum
  // over 0..Nmax leaves out)
  real p_above_nmax;
  {
    vector[nsite] e_fp = estimate_fp_re == 1 ? sigma_fp[1] * e_raw : zero_site;
    vector[nsite] d_det = estimate_det_re == 1 ? sigma_det[1] * d_raw : zero_site;
    vector[nsite] v_ab = estimate_abund_re == 1 ? sigma_lam[1] * v_raw : zero_site;
    vector[nsite] lam = exp(X * beta + v_ab);
    lambda_bar = mean(lam);
    psi_bar = mean(-expm1(-lam));
    p_above_nmax = 0;
    for (i in 1:nsite) p_above_nmax = fmax(p_above_nmax, exp(poisson_lccdf(Nmax | lam[i])));
    log_lik = site_loglik(1, nsite, beta, alpha, gamma0, gamma_fp, e_fp, u, d_det, v_ab,
                          delta0, delta_s, X, Z, Xfp, X1s, y2, n2, rec2, start2, len2,
                          any_det2, Ysite, Nsite, y1, start1, len1, n1v, s1v,
                          known_present, v3, w3, rec3, Z3, start3, len3, Nvec, lgamma_n1,
                          options);
  }
}
