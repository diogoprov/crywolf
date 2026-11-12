// ============================================================================
// Multispecies false-positive occupancy model (fp_occu_ms()): the model of
// fp_occu.stan for each species, with the coefficients of the species drawn
// from distributions of the community.
//
// Each site x species is a unit, numbered (site - 1) * nspecies + species,
// and the data of the units enter as those of the sites of fp_occu.stan.
// For unit i of species s:
//
//   z_i     ~ Bernoulli(psi_i),     logit(psi_i) = X_i beta_s
//   y2_ik   ~ Bernoulli(z_i p11_ik + (1 - z_i) p10_ik)
//             logit(p11_ik) = alpha_s + Z_ik[-1] alpha[-1] + u_rec[ik]
//             logit(p10_ik) = alpha_s - gap_s + Xfp_ik gamma_fp
//   y1_im   ~ Bernoulli(z_i r11_im), as in fp_occu.stan
//
//   beta_s  = mu_beta + sigma_beta .* zb_s,        zb_s ~ Normal(0, 1)
//   alpha_s = alpha[1] + sigma_alpha * za_s,       za_s ~ Normal(0, 1)
//   gap_s   = mu_gap * exp(sigma_gap * zg_s),      zg_s ~ Normal(0, 1)
//
// so that p11 > p10 at the reference for every species (gap_s > 0), and
// mu_gap is the median gap of the community. With by_det = 1, by_fp = 1 or
// by_ver = 1, the slopes of detection (alpha[-1]), the slopes of false
// positives (gamma_fp) or the intercept of the verified method (delta0) are
// also those of each species, drawn around the community values:
//
//   alpha_sk = alpha[k] + sigma_det_s[k] * zds_ks,  gamma_sk = gamma_fp[k] + sigma_fp_s[k] * zfs_ks,
//   delta_s  = delta0 + sigma_ver * zv_s,           zds, zfs, zv ~ Normal(0, 1)
//
// and the emissions add the deviations of each species to the shared terms
// (alpha_dev, gamma_dev, d_ver). Otherwise these terms and the recorder
// effects are shared by the species. The intercepts of the species reach the shared
// emissions (fp_emissions.stanfunctions) through the offsets of each unit
// on the logit of p11 (d_det) and of p10 (e_fp), the same offsets as the
// site effects of fp_occu.stan; gamma0 there is the community intercept
// alpha[1] - mu_gap, so e_fp = (alpha_s - alpha[1]) - gap_s + mu_gap.
//
// Random effects on logit(p10), added to e_fp: with estimate_fp_group = 1,
// v_c ~ Normal(0, sigma_v) for each group of sites x species (cell c; the
// groups are those of fp_occu_ms(fp_group = ), e.g. the physical site when
// the sites are site-months), and with estimate_fp_re = 1,
// t_i ~ Normal(0, sigma_t) for each unit:
//
//   logit(p10_ik) = alpha_s - gap_s + v_cell(i) + t_i + Xfp_ik gamma_fp
//
// With centered = 1 the same model has another parametrization, for species
// and recorders with many occasions, whose narrow posteriors make the
// non-centred one take very small steps: za holds alpha_s, zg holds
// log(gap_s), u_raw holds u, zds, zfs and zv hold the slopes and verified
// intercepts of the species, and v_raw and t_raw hold v and t.
// The sum over units runs with reduce_sum(), as over sites in fp_occu.stan.
// ============================================================================
functions {
#include fp_emissions.stanfunctions

  // Marginal log-likelihood of units a to b, with z summed out; lpsi holds
  // the logit of the occupancy probability of every unit.
  vector unit_loglik(int a, int b, vector lpsi, vector alpha, real gamma0,
                     vector gamma_fp, vector e_fp, vector u, vector d_det, real delta0,
                     vector delta_s, matrix alpha_dev, matrix gamma_dev, vector d_ver,
                     array[] int sp_unit,
                     tuple(matrix, matrix, matrix, matrix, matrix) design,
                     tuple(array[] int, array[] int, array[] int, array[] int,
                           vector, vector, array[] int, array[] int,
                           array[] int, array[] int, vector, vector, vector) automated,
                     tuple(vector, vector, array[] int, matrix, vector,
                           array[] int, array[] int) verified,
                     tuple(array[] int, array[] int, array[] int, matrix,
                           array[] int, array[] int, matrix) confirmed,
                     array[] int options) {
    int nu = b - a + 1;
    tuple(vector[nu], vector[nu]) em = unit_emissions(a, b, alpha, gamma0, gamma_fp, e_fp, u,
                                                      d_det, delta0, delta_s, alpha_dev,
                                                      gamma_dev, d_ver, sp_unit, design,
                                                      automated, verified, confirmed, options);
    vector[nu] ll;
    for (i in 1:nu) {
      real lp = lpsi[a + i - 1];
      real occ_branch = log_inv_logit(lp) + em.1[i];
      ll[i] = verified.3[a + i - 1] == 1
              ? occ_branch
              : log_sum_exp(occ_branch, log1m_inv_logit(lp) + em.2[i]);
    }
    return ll;
  }

  // partial sum over a slice of units, for reduce_sum()
  real partial_loglik(array[] int slice, int a, int b, vector lpsi, vector alpha,
                      real gamma0, vector gamma_fp, vector e_fp, vector u, vector d_det,
                      real delta0, vector delta_s, matrix alpha_dev, matrix gamma_dev,
                      vector d_ver, array[] int sp_unit,
                      tuple(matrix, matrix, matrix, matrix, matrix) design,
                      tuple(array[] int, array[] int, array[] int, array[] int,
                            vector, vector, array[] int, array[] int,
                            array[] int, array[] int, vector, vector, vector) automated,
                      tuple(vector, vector, array[] int, matrix, vector,
                            array[] int, array[] int) verified,
                      tuple(array[] int, array[] int, array[] int, matrix,
                            array[] int, array[] int, matrix) confirmed,
                      array[] int options) {
    return sum(unit_loglik(a, b, lpsi, alpha, gamma0, gamma_fp, e_fp, u, d_det, delta0,
                           delta_s, alpha_dev, gamma_dev, d_ver, sp_unit, design,
                           automated, verified, confirmed, options));
  }
}

data {
  int<lower=1> nsite;
  int<lower=1> nspecies;
  int<lower=1> nunit;                        // site x species: nsite * nspecies
  array[nunit] int<lower=1, upper=nspecies> sp_unit;  // species of each unit
  int<lower=1> Kbeta;
  matrix[nunit, Kbeta] X;                    // occupancy design (with intercept)

  int<lower=0> N2;                           // rows of automated occasions
  array[N2] int<lower=1> n2;                 // occasions in each row
  array[N2] int<lower=0> y2;                 // detections in each row
  array[N2] int<lower=1, upper=nunit> site2;  // sorted by site
  array[nunit] int<lower=1> start2;          // first row of each site in y2
  array[nunit] int<lower=0> len2;            // number of rows of each site
  int<lower=1> Kalpha;
  matrix[N2, Kalpha] Z;                      // detection design (column 1 = intercept)
  int<lower=1> nrec;
  array[N2] int<lower=1, upper=nrec> rec2;
  int<lower=0, upper=1> estimate_re;         // recorder random intercept on p11

  int<lower=0> Kfp;
  int<lower=0, upper=1> fp_by_row;           // false-positive covariates vary by occasion
  matrix[nunit, fp_by_row == 1 ? 0 : Kfp] Xfp;  // false-positive design (no intercept)
  matrix[N2, fp_by_row == 1 ? Kfp : 0] Xfp2;    // ... by automated row
  int<lower=0, upper=1> estimate_fp;

  int<lower=0> N1;                           // verified occasions (long format)
  array[N1] int<lower=0, upper=1> y1;
  array[N1] int<lower=1, upper=nunit> site1;
  int<lower=0> Kver_s;                       // slopes of the verified detection
  matrix[N1, Kver_s] X1s;                    // their design (verified rows sorted by site)
  array[nunit] int<lower=0, upper=1> known_present;

  int<lower=0, upper=1> obs_confirm;         // observation-confirmation design
  int<lower=0> N3;                           // confirmed occasions
  array[N3] int<lower=0, upper=3> v3;        // confirmation state
  array[N3] int<lower=1, upper=nunit> site3;
  matrix[N3, Kalpha] Z3;                     // detection design of confirmed occasions
  array[N3] int<lower=1, upper=nrec> rec3;
  matrix[N3, fp_by_row == 1 ? Kfp : 0] Xfp3; // false-positive design of confirmed occasions

  real<lower=0> prior_sd_int;                // sd of Normal priors on intercepts
  real<lower=0> prior_sd_slope;              // sd of Normal priors on slopes
  real gap_mean;                             // Normal(gap_mean, gap_sd) prior on fp_gap, truncated at 0
  real<lower=0> gap_sd;
  real<lower=0> recorder_sd;                 // half-Normal(0, recorder_sd) on sigma_u
  real<lower=0> species_sd;                  // half-Normal(0, species_sd) on the sds among species
  int<lower=0, upper=1> centered;            // parametrization of alpha_sp, gap_sp, u and the species slopes
  int<lower=0, upper=1> by_det;              // detection slopes by species
  int<lower=0, upper=1> by_fp;               // false-positive slopes by species
  int<lower=0, upper=1> by_ver;              // verified intercept by species
  int<lower=0, upper=1> estimate_fp_group;   // effect on logit(p10) per group of sites x species
  int<lower=1> ngroup;                       // groups of sites (1 without them)
  array[nsite] int<lower=1, upper=ngroup> group_site;  // group of each site
  int<lower=0, upper=1> estimate_fp_re;      // effect on logit(p10) per unit
  real<lower=0> fp_site_sd;                  // half-Normal(0, fp_site_sd) on sigma_v and sigma_t
}

transformed data {
  if (nunit != nsite * nspecies) reject("nunit must be nsite * nspecies");
  int Kds = by_det == 1 ? Kalpha - 1 : 0;    // detection slopes by species
  int Kfs = by_fp == 1 ? Kfp : 0;            // false-positive slopes by species
  int nsv = by_ver == 1 ? nspecies : 0;      // verified intercepts by species
  if ((estimate_fp_group == 1 || estimate_fp_re == 1) && estimate_fp == 0) {
    reject("random effects on p10 need a false-positive process (estimate_fp = 1)");
  }
  // cell (group of sites x species) of each unit
  int ncell = ngroup * nspecies;
  array[nunit] int cell_unit;
  for (i in 1:nunit) {
    cell_unit[i] = (group_site[(i - 1) %/% nspecies + 1] - 1) * nspecies + sp_unit[i];
  }
  vector[N2] y2v = to_vector(y2);
  vector[N2] n2v = to_vector(n2);
  // sites with an automated detection cannot be unoccupied without false positives
  array[nunit] int any_det2 = rep_array(0, nunit);
  for (n in 1:N2) {
    if (y2[n] > n2[n]) reject("y2 cannot exceed n2");
    if (n > 1 && site2[n] < site2[n - 1]) reject("automated rows must be sorted by site");
    if (y2[n] > 0) any_det2[site2[n]] = 1;
  }
  if (obs_confirm == 1 && estimate_fp == 0) {
    reject("the observation-confirmation design needs a false-positive process (estimate_fp = 1)");
  }
  // verified occasions and detections per site
  vector[nunit] n1v = rep_vector(0, nunit);
  vector[nunit] s1v = rep_vector(0, nunit);
  for (n in 1:N1) {
    n1v[site1[n]] += 1;
    s1v[site1[n]] += y1[n];
  }
  // with covariates on the verified detection, its rows are used one by one
  vector[N1] y1v = to_vector(y1);
  array[nunit] int len1 = rep_array(0, nunit);
  for (n in 1:N1) {
    if (Kver_s > 0 && n > 1 && site1[n] < site1[n - 1]) reject("verified occasions must be sorted by site");
    len1[site1[n]] += 1;
  }
  array[nunit] int start1;
  start1[1] = 1;
  for (i in 2:nunit) start1[i] = start1[i - 1] + len1[i - 1];
  // first row and number of rows of each site among the confirmed occasions
  array[nunit] int len3 = rep_array(0, nunit);
  for (n in 1:N3) {
    if (n > 1 && site3[n] < site3[n - 1]) reject("confirmed occasions must be sorted by site");
    len3[site3[n]] += 1;
  }
  array[nunit] int start3;
  start3[1] = 1;
  for (i in 2:nunit) start3[i] = start3[i - 1] + len3[i - 1];
  // units, sliced among threads by reduce_sum()
  array[nunit] int unit_index = linspaced_int_array(nunit, 1, nunit);
  // totals of automated occasions and detections per site (unoccupied branch)
  vector[nunit] Ysite = rep_vector(0, nunit);
  vector[nunit] Nsite = rep_vector(0, nunit);
  for (n in 1:N2) {
    Ysite[site2[n]] += y2[n];
    Nsite[site2[n]] += n2[n];
  }
  // the occupied branch uses one GLM call per site when detection has
  // covariates (otherwise the rows are already collapsed) and the design is
  // the multiple-methods one; Z_slopes is Z without its intercept column, and lchoose_site
  // the binomial constants that the GLM includes and the likelihood omits
  int use_glm = (Kalpha > 1 && obs_confirm == 0) ? 1 : 0;
  matrix[N2, Kalpha - 1] Z_slopes = rep_matrix(0, N2, Kalpha - 1);
  if (Kalpha > 1) Z_slopes = block(Z, 1, 2, N2, Kalpha - 1);
  vector[nunit] lchoose_site = rep_vector(0, nunit);
  for (n in 1:N2) lchoose_site[site2[n]] += lchoose(n2[n], y2[n]);
  // recorder of each site, when it does not change within sites
  int rec_by_site = 1;
  array[nunit] int rec_site = rep_array(1, nunit);
  for (n in 1:N2) {
    if (n == 1 || site2[n] != site2[n - 1]) rec_site[site2[n]] = rec2[n];
    else if (rec2[n] != rec_site[site2[n]]) rec_by_site = 0;
  }
  // the data of the likelihood, grouped for site_loglik() and reduce_sum();
  // the order of the fields is the order of the unpacking in site_loglik()
  tuple(matrix[nunit, Kbeta], matrix[N2, Kalpha], matrix[N2, Kalpha - 1],
        matrix[nunit, fp_by_row == 1 ? 0 : Kfp], matrix[N2, fp_by_row == 1 ? Kfp : 0])
    design = (X, Z, Z_slopes, Xfp, Xfp2);
  tuple(array[N2] int, array[N2] int, array[N2] int, array[N2] int, vector[N2],
        vector[N2], array[nunit] int, array[nunit] int, array[nunit] int,
        array[nunit] int, vector[nunit], vector[nunit], vector[nunit])
    automated = (y2, n2, site2, rec2, y2v, n2v, start2, len2, any_det2, rec_site,
                 Ysite, Nsite, lchoose_site);
  tuple(vector[nunit], vector[nunit], array[nunit] int, matrix[N1, Kver_s],
        vector[N1], array[nunit] int, array[nunit] int)
    verified = (n1v, s1v, known_present, X1s, y1v, start1, len1);
  tuple(array[N3] int, array[N3] int, array[N3] int, matrix[N3, Kalpha],
        array[nunit] int, array[nunit] int, matrix[N3, fp_by_row == 1 ? Kfp : 0])
    confirmed = (v3, site3, rec3, Z3, start3, len3, Xfp3);
  array[8] int options = {estimate_re, estimate_fp, Kfp, obs_confirm, use_glm,
                          rec_by_site, Kver_s, fp_by_row};
}

parameters {
  vector[Kbeta] mu_beta;                     // community means of the occupancy coefficients
  vector<lower=0>[Kbeta] sigma_beta;         // and their sds among species
  matrix[Kbeta, nspecies] zb;
  vector[Kalpha] alpha;                      // alpha[1]: community mean detection intercept
  real<lower=0> sigma_alpha;
  vector[nspecies] za;
  real<lower=0> mu_gap;                      // median gap of the community
  real<lower=0> sigma_gap;                   // sd of log(gap) among species
  vector[nspecies] zg;
  vector[nrec] u_raw;
  real<lower=0> sigma_u;
  vector[Kfp] gamma_fp;
  real delta0;
  vector[Kver_s] delta_s;
  vector<lower=0>[Kds] sigma_det_s;          // sds of the detection slopes among species
  matrix[Kds, nspecies] zds;
  vector<lower=0>[Kfs] sigma_fp_s;           // sds of the false-positive slopes among species
  matrix[Kfs, nspecies] zfs;
  vector<lower=0>[by_ver] sigma_ver;         // sd of the verified intercept among species
  vector[nsv] zv;
  vector<lower=0>[estimate_fp_group] sigma_v;   // sd of the effects of the cells on logit(p10)
  vector[estimate_fp_group == 1 ? ncell : 0] v_raw;
  vector<lower=0>[estimate_fp_re] sigma_t;      // sd of the effects of the units on logit(p10)
  vector[estimate_fp_re == 1 ? nunit : 0] t_raw;
}

transformed parameters {
  vector[nrec] u = centered == 1 && estimate_re == 1 ? u_raw : sigma_u * u_raw;
  real gamma0 = alpha[1] - mu_gap;
  matrix[Kbeta, nspecies] beta_sp = rep_matrix(mu_beta, nspecies) + diag_pre_multiply(sigma_beta, zb);
  vector[nspecies] alpha_sp = centered == 1 ? za : alpha[1] + sigma_alpha * za;
  vector[nspecies] gap_sp = centered == 1 ? exp(zg) : mu_gap * exp(sigma_gap * zg);
  // deviations of the species from the community slopes and verified
  // intercept, one row per species
  matrix[nspecies, Kds] alpha_dev;
  matrix[nspecies, Kfs] gamma_dev;
  vector[nsv] delta_dev;
  if (Kds > 0) {
    alpha_dev = centered == 1 ? (zds - rep_matrix(alpha[2:Kalpha], nspecies))'
                              : diag_pre_multiply(sigma_det_s, zds)';
  }
  if (Kfs > 0) {
    gamma_dev = centered == 1 ? (zfs - rep_matrix(gamma_fp, nspecies))'
                              : diag_pre_multiply(sigma_fp_s, zfs)';
  }
  if (nsv > 0) delta_dev = centered == 1 ? zv - delta0 : sigma_ver[1] * zv;
}

model {
  mu_beta[1]     ~ normal(0, prior_sd_int);
  if (Kbeta > 1)  mu_beta[2:Kbeta] ~ normal(0, prior_sd_slope);
  sigma_beta     ~ normal(0, species_sd);
  to_vector(zb)  ~ std_normal();
  alpha[1]       ~ normal(0, prior_sd_int);
  if (Kalpha > 1) alpha[2:Kalpha] ~ normal(0, prior_sd_slope);
  sigma_alpha    ~ normal(0, species_sd);
  if (centered == 1) za ~ normal(alpha[1], sigma_alpha);
  else za ~ std_normal();
  mu_gap         ~ normal(gap_mean, gap_sd);
  sigma_gap      ~ normal(0, species_sd);
  if (centered == 1 && estimate_fp == 1) zg ~ normal(log(mu_gap), sigma_gap);
  else zg ~ std_normal();
  if (centered == 1 && estimate_re == 1) u_raw ~ normal(0, sigma_u);
  else u_raw ~ std_normal();
  sigma_u        ~ normal(0, recorder_sd);
  gamma_fp       ~ normal(0, prior_sd_slope);
  delta0         ~ normal(0, prior_sd_int);
  delta_s        ~ normal(0, prior_sd_slope);
  sigma_det_s    ~ normal(0, species_sd);
  sigma_fp_s     ~ normal(0, species_sd);
  sigma_ver      ~ normal(0, species_sd);
  sigma_v        ~ normal(0, fp_site_sd);
  sigma_t        ~ normal(0, fp_site_sd);
  if (estimate_fp_group == 1) {
    if (centered == 1) v_raw ~ normal(0, sigma_v[1]);
    else v_raw ~ std_normal();
  }
  if (estimate_fp_re == 1) {
    if (centered == 1) t_raw ~ normal(0, sigma_t[1]);
    else t_raw ~ std_normal();
  }
  if (centered == 1) {
    if (Kds > 0) to_vector(zds) ~ normal(to_vector(rep_matrix(alpha[2:Kalpha], nspecies)),
                                         to_vector(rep_matrix(sigma_det_s, nspecies)));
    if (Kfs > 0) to_vector(zfs) ~ normal(to_vector(rep_matrix(gamma_fp, nspecies)),
                                         to_vector(rep_matrix(sigma_fp_s, nspecies)));
    if (nsv > 0) zv ~ normal(delta0, sigma_ver[1]);
  } else {
    to_vector(zds) ~ std_normal();
    to_vector(zfs) ~ std_normal();
    zv ~ std_normal();
  }
  {
    vector[nunit] lpsi = rows_dot_product(X, beta_sp[, sp_unit]');
    vector[nunit] d_det = (alpha_sp - alpha[1])[sp_unit];
    vector[nunit] e_fp = (alpha_sp - alpha[1] - gap_sp + mu_gap)[sp_unit];
    if (estimate_fp_group == 1) e_fp += (centered == 1 ? v_raw : sigma_v[1] * v_raw)[cell_unit];
    if (estimate_fp_re == 1) e_fp += centered == 1 ? t_raw : sigma_t[1] * t_raw;
    vector[nsv > 0 ? nunit : 0] d_ver;
    if (nsv > 0) d_ver = delta_dev[sp_unit];
    target += reduce_sum(partial_loglik, unit_index, 1, lpsi, alpha, gamma0, gamma_fp,
                         e_fp, u, d_det, delta0, delta_s, alpha_dev, gamma_dev, d_ver,
                         sp_unit, design, automated, verified, confirmed, options);
  }
}

generated quantities {
  vector[nunit] log_lik;
  // detection and false-positive rates of each species at the reference
  // (all covariates at 0), and of the community (species effects at 0)
  vector[nspecies] p11_sp = inv_logit(alpha_sp);
  vector[nspecies] p10_sp = estimate_fp == 1 ? inv_logit(alpha_sp - gap_sp) : rep_vector(0, nspecies);
  real p11_ref = inv_logit(alpha[1]);
  real p10_ref = estimate_fp == 1 ? inv_logit(gamma0) : 0.0;
  real r11 = inv_logit(delta0);
  // slopes and verified detection of each species
  matrix[Kds, nspecies] alpha_sp_s;
  matrix[Kfs, nspecies] gamma_sp;
  vector[nsv] r11_sp;
  vector[nsv] delta_sp;
  if (Kds > 0) alpha_sp_s = rep_matrix(alpha[2:Kalpha], nspecies) + alpha_dev';
  if (Kfs > 0) gamma_sp = rep_matrix(gamma_fp, nspecies) + gamma_dev';
  if (nsv > 0) {
    delta_sp = delta0 + delta_dev;
    r11_sp = inv_logit(delta_sp);
  }
  if (obs_confirm == 1) {
    p11_sp = 1 - (1 - p11_sp) .* (1 - p10_sp);
    p11_ref = 1 - (1 - p11_ref) * (1 - p10_ref);
  }
  {
    vector[nunit] lpsi = rows_dot_product(X, beta_sp[, sp_unit]');
    vector[nunit] d_det = (alpha_sp - alpha[1])[sp_unit];
    vector[nunit] e_fp = (alpha_sp - alpha[1] - gap_sp + mu_gap)[sp_unit];
    if (estimate_fp_group == 1) e_fp += (centered == 1 ? v_raw : sigma_v[1] * v_raw)[cell_unit];
    if (estimate_fp_re == 1) e_fp += centered == 1 ? t_raw : sigma_t[1] * t_raw;
    vector[nsv > 0 ? nunit : 0] d_ver;
    if (nsv > 0) d_ver = delta_dev[sp_unit];
    log_lik = unit_loglik(1, nunit, lpsi, alpha, gamma0, gamma_fp, e_fp, u, d_det, delta0,
                          delta_s, alpha_dev, gamma_dev, d_ver, sp_unit, design,
                          automated, verified, confirmed, options);
  }
}
