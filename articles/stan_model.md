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
    }
    transformed parameters {
      vector[nrec] u = sigma_u * u_raw;
      real gamma0 = alpha[1] - fp_gap;
    }

The recorder effects are written as `sigma_u * u_raw`, with `u_raw`
standard normal (the non-centred parametrization), which samples better
when there are few recorders or little information about their variance.

The false-positive intercept is not a parameter of its own:
`gamma0 = alpha[1] - fp_gap`, with `fp_gap` constrained to be positive.
Because column 1 of `Z` is the intercept, `alpha[1]` is the logit of
true detection at the reference values of the covariates (all zero), and
the gap keeps the false-positive rate below it there. Without this
constraint, the likelihood has a second mode in which occupied and
unoccupied sites swap roles, and the sampler could visit either.

## Likelihood

The occupancy state z_i is discrete, so Stan cannot sample it. The
function `site_loglik()` sums it out for each site:

    functions {
      // linear predictor of true detection on the confirmed occasions
      vector confirmed_lp11(matrix Z3, vector alpha, vector u, array[] int rec3,
                    int estimate_re) {
        vector[rows(Z3)] lp = Z3 * alpha;
        if (estimate_re == 1) lp += u[rec3];
        return lp;
      }

      // Marginal log-likelihood of each site, with the latent state z summed out.
      vector site_loglik(vector lpsi, vector lp11_fixed, real gamma0, matrix Xfp,
                         vector gamma_fp, vector u, array[] int rec2,
                         int estimate_re, int estimate_fp, int Kfp,
                         vector y2v, vector f2v, array[] int site2, array[] int start2,
                         array[] int len2, array[] int any_det2,
                         array[] int y1, array[] int site1, real delta0,
                         array[] int known_present,
                         int obs_confirm, array[] int v3, array[] int site3,
                         vector lp11_3) {
        int nsite = rows(lpsi);
        int N2 = rows(y2v);
        vector[N2] lp11 = lp11_fixed;
        if (estimate_re == 1) lp11 += u[rec2];
        vector[nsite] lp10 = rep_vector(gamma0, nsite);
        if (Kfp > 0) lp10 += Xfp * gamma_fp;

        // per-occasion log-likelihoods, vectorized
        vector[N2] l11;
        vector[N2] l10;
        if (estimate_fp == 1) {
          vector[N2] lp10_n = lp10[site2];
          l10 = y2v .* log_inv_logit(lp10_n) + f2v .* log1m_inv_logit(lp10_n);
        }
        if (obs_confirm == 1) {
          // log(1 - p11) = log(1 - s1) + log(1 - s0)
          vector[N2] lq11 = log1m_inv_logit(lp11) + log1m_inv_logit(lp10[site2]);
          l11 = y2v .* log1m_exp(lq11) + f2v .* lq11;
        } else {
          l11 = y2v .* log_inv_logit(lp11) + f2v .* log1m_inv_logit(lp11);
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
        if (obs_confirm == 1) {
          for (n in 1:size(v3)) {
            int i = site3[n];
            real ls1 = log_inv_logit(lp11_3[n]);
            real lq1 = log1m_inv_logit(lp11_3[n]);
            real ls0 = log_inv_logit(lp10[i]);
            real lq0 = log1m_inv_logit(lp10[i]);
            occ[i] += (v3[n] >= 2 ? ls1 : lq1) + (v3[n] % 2 == 1 ? ls0 : lq0);
            unocc[i] += v3[n] >= 2 ? negative_infinity() : (v3[n] == 1 ? ls0 : lq0);
          }
        }

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

The per-occasion terms are computed as vectors and summed by site with
`segment()`, which works because the rows of each site are contiguous.

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

      {
        vector[nsite] ll_site = site_loglik(X * beta, Z * alpha, gamma0, Xfp, gamma_fp,
                                            u, rec2, estimate_re, estimate_fp, Kfp,
                                            y2v, f2v, site2, start2, len2, any_det2,
                                            y1, site1, delta0, known_present,
                                            obs_confirm, v3, site3, confirmed_lp11(Z3, alpha, u, rec3, estimate_re));
        target += sum(ll_site);
      }
    }

All priors are on the logit scale and are listed with their defaults in
[`?fp_priors`](https://provetelab.org/crywolf/reference/fp_priors.md).
The model block calls `site_loglik()` once and adds the sum over sites
to the log density.

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
        vector[nsite] lp10 = rep_vector(gamma0, nsite);
        if (Kfp > 0) lp10 += Xfp * gamma_fp;
        if (estimate_fp == 1) p10_bar = mean(inv_logit(lp10));
        log_lik = site_loglik(X * beta, Z * alpha, gamma0, Xfp, gamma_fp,
                              u, rec2, estimate_re, estimate_fp, Kfp,
                              y2v, f2v, site2, start2, len2, any_det2,
                              y1, site1, delta0, known_present,
                              obs_confirm, v3, site3, confirmed_lp11(Z3, alpha, u, rec3, estimate_re));
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
