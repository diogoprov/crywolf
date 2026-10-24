# The Stan model

[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) fits
one Stan program for every model in the package: the multiple-methods
design of Miller et al. (2011), its observation-confirmation variant
(Chambert et al., 2015), and the standard occupancy model. This article
walks through it for readers who want to know exactly what is estimated,
or who want to change it. The model itself is described in
[`vignette("crywolf", package = "crywolf")`](https://provetelab.org/crywolf/articles/crywolf.md);
here the focus is on how it is coded. The program is a file installed
with the package:

``` r

fp_stan_file()
```

[`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md)
copies it to the user cache directory and compiles it the first time it
is needed, so later sessions reuse the executable.

## Data

[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) turns
the `fp_data` object and the formulas into the list that Stan reads,
which is kept in the fitted object as `m$stan_data`. The data block
declares it:

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
      int<lower=0, upper=1> fp_by_row;           // false-positive covariates vary by occasion
      matrix[nsite, fp_by_row == 1 ? 0 : Kfp] Xfp;  // false-positive design (no intercept)
      matrix[N2, fp_by_row == 1 ? Kfp : 0] Xfp2;    // ... by automated row
      int<lower=0, upper=1> estimate_fp;

      int<lower=0> N1;                           // verified occasions (long format)
      array[N1] int<lower=0, upper=1> y1;
      array[N1] int<lower=1, upper=nsite> site1;
      int<lower=0> Kver_s;                       // slopes of the verified detection
      matrix[N1, Kver_s] X1s;                    // their design (verified rows sorted by site)
      array[nsite] int<lower=0, upper=1> known_present;

      int<lower=0, upper=1> obs_confirm;         // observation-confirmation design
      int<lower=0> N3;                           // confirmed occasions
      array[N3] int<lower=0, upper=3> v3;        // confirmation state
      array[N3] int<lower=1, upper=nsite> site3;
      matrix[N3, Kalpha] Z3;                     // detection design of confirmed occasions
      array[N3] int<lower=1, upper=nrec> rec3;
      matrix[N3, fp_by_row == 1 ? Kfp : 0] Xfp3; // false-positive design of confirmed occasions

      real<lower=0> prior_sd_int;                // sd of Normal priors on intercepts
      real<lower=0> prior_sd_slope;              // sd of Normal priors on slopes
      real gap_mean;                             // Normal(gap_mean, gap_sd) prior on fp_gap, truncated at 0
      real<lower=0> gap_sd;
      real<lower=0> recorder_sd;                 // half-Normal(0, recorder_sd) on sigma_u

      int<lower=0, upper=1> estimate_fp_re;      // site random effect on logit(p10)
      real<lower=0> fp_site_sd;                  // half-Normal(0, fp_site_sd) on sigma_fp
      int<lower=0, upper=1> estimate_det_re;     // site random effect on logit(p11)
      real<lower=0> det_site_sd;                 // half-Normal(0, det_site_sd) on sigma_det
    }

The occupancy design `X`, the detection design `Z`, and the
false-positive design `Xfp` are model matrices built from the formulas.
`Xfp` has no intercept column, because the intercept of the
false-positive rate is not a free parameter (see below). Switches turn
parts of the model on and off: `estimate_fp = 0` gives the standard
occupancy model (`fp = NULL`), `Kfp = 0` the constant false-positive
rate, `estimate_re = 1` the recorder random intercept, and
`obs_confirm = 1` the observation-confirmation design. The prior scales
are data too, which is why
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)
can change them without recompiling.

### Binomial rows

Automated occasions do not enter one by one. Occasions of the same site
with identical detection covariates (and recorder) are collapsed into
one row with `n2` trials and `y2` detections. Given the occupancy state,
the occasions of a site are independent Bernoulli trials with the same
probability, so the product of their likelihoods equals the binomial
likelihood without its coefficient, which does not depend on the
parameters. The posterior is therefore exactly the same, and the model
is much faster when detection covariates are constant within sites. With
occasion-level covariates, such as the acoustic indices of the
[`vignette("curiango", package = "crywolf")`](https://provetelab.org/crywolf/articles/curiango.md)
data, each occasion is its own row. Rows are sorted by site, and
`start2` and `len2` give the first row and the number of rows of each
site.

## Parameters

    parameters {
      vector[Kbeta] beta;
      vector[Kalpha] alpha;
      vector[nrec] u_raw;
      real<lower=0> sigma_u;
      real<lower=0> fp_gap;
      vector[Kfp] gamma_fp;
      real delta0;
      vector[Kver_s] delta_s;
      vector[estimate_fp_re == 1 ? nsite : 0] e_raw;
      array[estimate_fp_re] real<lower=0> sigma_fp;
      vector[estimate_det_re == 1 ? nsite : 0] d_raw;
      array[estimate_det_re] real<lower=0> sigma_det;
    }
    transformed parameters {
      vector[nrec] u = sigma_u * u_raw;
      real gamma0 = alpha[1] - fp_gap;
    }

The recorder effects are written as `sigma_u * u_raw`, with `u_raw`
standard normal (the non-centred parametrization), which samples better
when there are few recorders or little information about their variance.

With `random_fp = TRUE` (`estimate_fp_re = 1`), the site effects on the
false-positive rate are written in the same way, `sigma_fp[1] * e_raw`,
and are added to its linear predictor; with the option off, `e_raw` and
`sigma_fp` have length zero and the model is the same as without them.
`random_det = TRUE` (`estimate_det_re = 1`) does the same for true
detection, with `sigma_det[1] * d_raw` added to the logit of p\_{11}.
The likelihood receives the site effects that are off as a data vector
of zeros (`zero_site`), so that `reduce_sum()` does not differentiate
them.

The false-positive intercept is not a parameter of its own:
`gamma0 = alpha[1] - fp_gap`, with `fp_gap` constrained to be positive.
Because column 1 of `Z` is the intercept, `alpha[1]` is the logit of
true detection at the reference values of the covariates (all zero), and
the gap keeps the false-positive rate below it there. Without this
constraint, the likelihood has a second mode in which occupied and
unoccupied sites swap roles, and the sampler could visit either.

## Likelihood

The occupancy state z_i is discrete, so Stan cannot sample it. The
function `loglik_sites()` sums it out for each site in a range of sites,
from `a` to `b`. It is called through `site_loglik()`, which receives
the data grouped in tuples (built at the end of `transformed data`) and
hands their fields on as arguments; Stan passes function arguments by
reference, so nothing is copied. The tuples need CmdStan 2.40 or later:

    functions {
      // Marginal log-likelihood of sites a to b, with the latent state z summed
      // out. The automated (and confirmed) rows of these sites are contiguous.
      // Called through site_loglik(), which takes the data grouped in tuples.
      vector loglik_sites(int a, int b, matrix X, vector beta, matrix Z, vector alpha,
                          real gamma0, matrix Xfp, vector gamma_fp, vector e_fp,
                          vector u, vector d_det, array[] int rec2, int estimate_re,
                          int estimate_fp, int Kfp,
                          vector y2v, vector n2v, array[] int site2,
                          array[] int start2, array[] int len2, array[] int any_det2,
                          vector n1v, vector s1v, real delta0,
                          array[] int known_present, int obs_confirm,
                          array[] int v3, array[] int site3, matrix Z3,
                          array[] int rec3, array[] int start3, array[] int len3,
                          vector Ysite, vector Nsite, int use_glm, matrix Z_slopes,
                          array[] int y2, array[] int n2, vector lchoose_site,
                          int rec_by_site, array[] int rec_site,
                          int Kver_s, matrix X1s, vector delta_s, vector y1v,
                          array[] int start1, array[] int len1,
                          int fp_by_row, matrix Xfp2, matrix Xfp3) {
        int ns = b - a + 1;
        vector[ns] lpsi = X[a:b] * beta;
        vector[ns] lp10 = rep_vector(gamma0, ns);
        if (Kfp > 0 && fp_by_row == 0) lp10 += Xfp[a:b] * gamma_fp;
        lp10 += e_fp[a:b];
        vector[ns] occ = rep_vector(0, ns);
        vector[ns] unocc = rep_vector(0, ns);

        // unoccupied branch: the false-positive rate is constant within a site,
        // so the automated rows of a site collapse into its totals of occasions
        // (Nsite) and detections (Ysite); log1m_inv_logit(x) = -log1p_exp(x)
        vector[ns] log1p_exp_lp10;
        int r0 = start2[a];
        int nr = start2[b] + len2[b] - r0;
        // with occasion-level false-positive covariates, the rate of each row
        vector[fp_by_row == 1 ? nr : 0] log1p_exp_lp10_row;
        if (estimate_fp == 1 && fp_by_row == 0) {
          log1p_exp_lp10 = log1p_exp(lp10);
          unocc = Ysite[a:b] .* lp10 - Nsite[a:b] .* log1p_exp_lp10;
        } else if (estimate_fp == 1 && nr > 0) {
          int r1 = r0 + nr - 1;
          array[nr] int loc;
          for (k in 1:nr) loc[k] = site2[r0 + k - 1] - a + 1;
          vector[nr] lp10_row = lp10[loc] + Xfp2[r0:r1] * gamma_fp;
          log1p_exp_lp10_row = log1p_exp(lp10_row);
          vector[nr] l10 = y2v[r0:r1] .* lp10_row - n2v[r0:r1] .* log1p_exp_lp10_row;
          for (i in 1:ns) {
            int s = a + i - 1;
            if (len2[s] > 0) unocc[i] = sum(segment(l10, start2[s] - r0 + 1, len2[s]));
          }
        }

        // occupied branch, automated occasions
        if (nr > 0 && use_glm == 1) {
          // one GLM call per site: the product Z * alpha, the binomial
          // log-likelihood and its gradient are computed in C++ without
          // intermediate autodiff variables; lchoose_site removes the binomial
          // constant. When the recorder is constant within sites the intercept
          // is a scalar, which avoids one autodiff variable per row (1.6 times
          // faster per gradient on the curiango data than the vector form)
          vector[rows(alpha) - 1] slopes = alpha[2:rows(alpha)];
          for (i in 1:ns) {
            int s = a + i - 1;
            if (len2[s] > 0) {
              int q0 = start2[s];
              int q1 = q0 + len2[s] - 1;
              if (estimate_re == 1 && rec_by_site == 0) {
                occ[i] = binomial_logit_glm_lpmf(y2[q0:q1] | n2[q0:q1], Z_slopes[q0:q1],
                                                 alpha[1] + d_det[s] + u[rec2[q0:q1]], slopes)
                         - lchoose_site[s];
              } else {
                real intercept = d_det[s] + (estimate_re == 1 ? alpha[1] + u[rec_site[s]] : alpha[1]);
                occ[i] = binomial_logit_glm_lpmf(y2[q0:q1] | n2[q0:q1], Z_slopes[q0:q1],
                                                 intercept, slopes) - lchoose_site[s];
              }
            }
          }
        } else if (nr > 0) {
          // vectorized rows; log_inv_logit(x) = x - log1p_exp(x)
          int r1 = r0 + nr - 1;
          vector[nr] y = y2v[r0:r1];
          vector[nr] n = n2v[r0:r1];
          vector[nr] lp11 = Z[r0:r1] * alpha;
          if (estimate_re == 1) lp11 += u[rec2[r0:r1]];
          lp11 += d_det[site2[r0:r1]];
          vector[nr] l11;
          if (obs_confirm == 1) {
            // log(1 - p11) = log(1 - s1) + log(1 - s0)
            vector[nr] lq11;
            if (fp_by_row == 1) {
              lq11 = -log1p_exp(lp11) - log1p_exp_lp10_row;
            } else {
              array[nr] int loc;
              for (k in 1:nr) loc[k] = site2[r0 + k - 1] - a + 1;
              lq11 = -log1p_exp(lp11) - log1p_exp_lp10[loc];
            }
            l11 = y .* log1m_exp(lq11) + (n - y) .* lq11;
          } else {
            l11 = y .* lp11 - n .* log1p_exp(lp11);
          }
          for (i in 1:ns) {
            int s = a + i - 1;
            if (len2[s] > 0) occ[i] = sum(segment(l11, start2[s] - r0 + 1, len2[s]));
          }
        }
        if (estimate_fp == 0) {
          for (i in 1:ns) if (any_det2[a + i - 1] == 1) unocc[i] = negative_infinity();
        }

        // verified occasions: counts per site, or rows when the verified method
        // has covariates (sorted by site, like the automated rows)
        if (Kver_s == 0) {
          occ += s1v[a:b] * log_inv_logit(delta0) + (n1v[a:b] - s1v[a:b]) * log1m_inv_logit(delta0);
        } else {
          int w0 = start1[a];
          int nw = start1[b] + len1[b] - w0;
          if (nw > 0) {
            int w1 = w0 + nw - 1;
            vector[nw] eta1 = delta0 + X1s[w0:w1] * delta_s;
            vector[nw] l1 = y1v[w0:w1] .* eta1 - log1p_exp(eta1);
            for (i in 1:ns) {
              int s = a + i - 1;
              if (len1[s] > 0) occ[i] += sum(segment(l1, start1[s] - w0 + 1, len1[s]));
            }
          }
        }

        // confirmed occasions (observation-confirmation design)
        if (obs_confirm == 1) {
          int q0 = start3[a];
          int nq = start3[b] + len3[b] - q0;
          if (nq > 0) {
            int q1 = q0 + nq - 1;
            vector[nq] lp11_3 = Z3[q0:q1] * alpha;
            if (estimate_re == 1) lp11_3 += u[rec3[q0:q1]];
            lp11_3 += d_det[site3[q0:q1]];
            for (k in 1:nq) {
              int i = site3[q0 + k - 1] - a + 1;
              int v = v3[q0 + k - 1];
              real ls1 = log_inv_logit(lp11_3[k]);
              real lq1 = log1m_inv_logit(lp11_3[k]);
              real lp10_k = fp_by_row == 1 ? lp10[i] + Xfp3[q0 + k - 1] * gamma_fp : lp10[i];
              real ls0 = log_inv_logit(lp10_k);
              real lq0 = log1m_inv_logit(lp10_k);
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

      // The data of the likelihood grouped in tuples (built at the end of
      // transformed data). Function arguments are passed by reference, so
      // handing the fields on to loglik_sites() copies nothing.
      vector site_loglik(int a, int b, vector beta, vector alpha, real gamma0,
                         vector gamma_fp, vector e_fp, vector u, vector d_det, real delta0,
                         vector delta_s,
                         tuple(matrix, matrix, matrix, matrix, matrix) design,
                         tuple(array[] int, array[] int, array[] int, array[] int,
                               vector, vector, array[] int, array[] int,
                               array[] int, array[] int, vector, vector, vector) automated,
                         tuple(vector, vector, array[] int, matrix, vector,
                               array[] int, array[] int) verified,
                         tuple(array[] int, array[] int, array[] int, matrix,
                               array[] int, array[] int, matrix) confirmed,
                         array[] int options) {
        // design: X, Z, Z_slopes, Xfp, Xfp2
        // automated: y2, n2, site2, rec2, y2v, n2v, start2, len2, any_det2,
        //            rec_site, Ysite, Nsite, lchoose_site
        // verified: n1v, s1v, known_present, X1s, y1v, start1, len1
        // confirmed: v3, site3, rec3, Z3, start3, len3, Xfp3
        // options: estimate_re, estimate_fp, Kfp, obs_confirm, use_glm,
        //          rec_by_site, Kver_s, fp_by_row
        return loglik_sites(a, b, design.1, beta, design.2, alpha, gamma0, design.4,
                            gamma_fp, e_fp, u, d_det, automated.4, options[1], options[2],
                            options[3], automated.5, automated.6, automated.3,
                            automated.7, automated.8, automated.9, verified.1,
                            verified.2, delta0, verified.3, options[4], confirmed.1,
                            confirmed.2, confirmed.4, confirmed.3, confirmed.5,
                            confirmed.6, automated.11, automated.12, options[5],
                            design.3, automated.1, automated.2, automated.13,
                            options[6], automated.10, options[7], verified.4,
                            delta_s, verified.5, verified.6, verified.7, options[8],
                            design.5, confirmed.7);
      }

      // partial sum over a slice of sites, for reduce_sum()
      real partial_loglik(array[] int slice, int a, int b, vector beta, vector alpha,
                          real gamma0, vector gamma_fp, vector e_fp, vector u, vector d_det,
                          real delta0, vector delta_s,
                          tuple(matrix, matrix, matrix, matrix, matrix) design,
                          tuple(array[] int, array[] int, array[] int, array[] int,
                                vector, vector, array[] int, array[] int,
                                array[] int, array[] int, vector, vector, vector) automated,
                          tuple(vector, vector, array[] int, matrix, vector,
                                array[] int, array[] int) verified,
                          tuple(array[] int, array[] int, array[] int, matrix,
                                array[] int, array[] int, matrix) confirmed,
                          array[] int options) {
        return sum(site_loglik(a, b, beta, alpha, gamma0, gamma_fp, e_fp, u, d_det, delta0,
                               delta_s, design, automated, verified, confirmed, options));
      }
    }

For site i, `occ[i]` accumulates the log-likelihood of all its data if
the site is occupied and `unocc[i]` if it is not; the site’s
contribution is then

\log\left\[\psi_i \exp(\text{occ}\_i) + (1 - \psi_i)
\exp(\text{unocc}\_i)\right\],

computed with `log_sum_exp()` for numerical stability. Three cases
shortcut it:

- A site with a verified detection (`known_present`) must be occupied,
  so only the first term is kept.
- Under the standard model (`estimate_fp = 0`), an automated detection
  also proves occupancy, so `unocc[i]` is set to -\infty.
- Under the observation-confirmation design, an occasion confirmed to
  hold a true observation (v_3 \geq 2) rules out the unoccupied state in
  the same way.

The two branches are computed differently. The false-positive rate is
constant within a site, so the unoccupied branch needs only the totals
of automated occasions and detections of each site (`Nsite`, `Ysite`,
built in `transformed data`): \text{unocc}\_i = Y_i \eta_i - N_i
\log(1 + e^{\eta_i}), with \eta_i the linear predictor of the rate. This
is the binomial log-likelihood without its constant, written as y \eta -
n \log(1 + e^{\eta}) rather than y \log p + (n - y) \log(1 - p); the
first form needs one transcendental function per term instead of two. In
the occupied branch the detection probability can change from one
occasion to the next, so when `det` has covariates each site’s rows go
through one call of Stan’s `binomial_logit_glm_lpmf()`, which evaluates
the linear predictor, the log-likelihood and its gradient in one pass
without creating an autodiff variable per row (`lchoose_site` removes
the binomial constant that this function includes). When the rows are
already collapsed (`det = ~ 1`), or under the observation-confirmation
design, the per-row terms are computed as vectors and summed by site
with `segment()`, which works because the rows of each site are
contiguous. The verified occasions enter as the number of occasions and
detections of each site, which gives the same sum as one Bernoulli term
per occasion.

On the
[`vignette("curiango", package = "crywolf")`](https://provetelab.org/crywolf/articles/curiango.md)
data, these two changes made each gradient evaluation three to four
times faster than computing every row as a vector, with one thread.
Threads, described below, shorten the fits further, though by less than
their number.

### Parallel likelihood

Most of the computing time goes to the automated occasions: in the
[`vignette("curiango", package = "crywolf")`](https://provetelab.org/crywolf/articles/curiango.md)
data, with 97,806 one-minute occasions and a recorder effect, they take
most of the time spent evaluating the log density. Because sites are
independent given the parameters, the model block adds them up with
`reduce_sum()`, which splits the sites into slices and evaluates the
slices on different threads. `partial_loglik()` is the sum of
`site_loglik()` over one slice; each slice reads only its own rows of
the automated and confirmed occasions, located with `start2` and `len2`
(and their counterparts for the confirmed occasions, built in
`transformed data`).

[`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md)
builds the model with threading enabled, and
`fp_occu(threads_per_chain = )` sets the number of threads of each
chain. With one thread the sum is computed serially, as before. With
more, only the order of the additions changes, so the log density
differs from the serial one by floating-point rounding (around 10^{-16}
in relative terms). Chains with the same seed then follow the same path
for a while and drift apart later, as they would on another computer,
and the posterior is the same. The gain is largest when there are many
occasions: keep `parallel_chains * threads_per_chain` at or below the
number of physical cores.

## Priors and model block

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
      delta_s        ~ normal(0, prior_sd_slope);
      e_raw          ~ std_normal();
      sigma_fp       ~ normal(0, fp_site_sd);
      d_raw          ~ std_normal();
      sigma_det      ~ normal(0, det_site_sd);

      // The four calls differ only in the site effects on the false-positive
      // and true-detection rates: parameters when estimated, otherwise the data
      // vector zero_site, which reduce_sum() neither copies nor differentiates
      // for each slice.
      if (estimate_fp_re == 1 && estimate_det_re == 1) {
        vector[nsite] e_fp = sigma_fp[1] * e_raw;
        vector[nsite] d_det = sigma_det[1] * d_raw;
        target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                             e_fp, u, d_det, delta0, delta_s, design, automated, verified,
                             confirmed, options);
      } else if (estimate_fp_re == 1) {
        vector[nsite] e_fp = sigma_fp[1] * e_raw;
        target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                             e_fp, u, zero_site, delta0, delta_s, design, automated, verified,
                             confirmed, options);
      } else if (estimate_det_re == 1) {
        vector[nsite] d_det = sigma_det[1] * d_raw;
        target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                             zero_site, u, d_det, delta0, delta_s, design, automated, verified,
                             confirmed, options);
      } else {
        target += reduce_sum(partial_loglik, site_index, 1, beta, alpha, gamma0, gamma_fp,
                             zero_site, u, zero_site, delta0, delta_s, design, automated,
                             verified, confirmed, options);
      }
    }

All priors are on the logit scale and are listed with their defaults in
[`?fp_priors`](https://provetelab.org/crywolf/reference/fp_priors.md).
The model block adds the sum over sites, computed by `reduce_sum()`, to
the log density.

## Generated quantities

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
        vector[nsite] d_det = rep_vector(0, nsite);
        if (estimate_det_re == 1) d_det = sigma_det[1] * d_raw;
        vector[nsite] lp10 = rep_vector(gamma0, nsite);
        if (Kfp > 0 && fp_by_row == 0) lp10 += Xfp * gamma_fp;
        lp10 += e_fp;
        if (estimate_fp == 1 && fp_by_row == 0) p10_bar = mean(inv_logit(lp10));
        // with occasion-level covariates: mean rate over the automated occasions
        if (estimate_fp == 1 && fp_by_row == 1 && N2 > 0) {
          p10_bar = dot_product(n2v, inv_logit(lp10[site2] + Xfp2 * gamma_fp)) / sum(n2v);
        }
        log_lik = site_loglik(1, nsite, beta, alpha, gamma0, gamma_fp, e_fp, u, d_det, delta0,
                              delta_s, design, automated, verified, confirmed, options);
      }
    }

`log_lik` holds one value per site, with z_i summed out, so
[`loo()`](https://mc-stan.org/loo/reference/loo.html) on a fitted model
gives leave-one-site-out cross-validation: a site with all its occasions
is the unit left out. `p11_ref` and `p10_ref` are the true- and
false-detection probabilities at the reference values of the covariates,
`p10_bar` is the false-positive rate averaged over the sites in the
data, and `r11` is the detection probability of the verified method
(`p_verified` in the output of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)).

## Modifying a copy

To change the model, for example its priors’ families or the structure
of a submodel, edit a copy of the file and fit it with ‘cmdstanr’
directly, using the data list that
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
prepared:

``` r

file.copy(fp_stan_file(), "my_fp_occu.stan")
# ... edit my_fp_occu.stan ...
sim <- simulate_fp(seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)  # for the data list
my_model <- cmdstanr::cmdstan_model("my_fp_occu.stan")
my_fit <- my_model$sample(data = m$stan_data, chains = 4, parallel_chains = 4,
                          adapt_delta = 0.95)
my_fit$summary(c("beta", "alpha", "fp_gap", "p10_ref"))
```

A copy compiled without `cpp_options = list(stan_threads = TRUE)` runs
the likelihood on one thread per chain, which gives the same results.
The result is a `CmdStanMCMC` object with the Stan names of the
parameters (`beta`, `alpha`, `gamma0`, …), not an `fp_occu` object, so
the methods of crywolf do not apply to it. If the edit adds data, add
the new elements to the list before sampling.
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md)
recomputes the likelihood in R from the posterior draws, so it describes
the original model, not an edited one.

## References

Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339. <https://doi.org/10.1890/14-1507.1>

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
<https://doi.org/10.1890/10-1396.1>
