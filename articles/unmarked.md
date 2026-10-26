# Coming from unmarked::occuFP()

The ‘unmarked’ package (Fiske & Chandler, 2011; Kellner et al., 2023)
fits the false-positive occupancy models of Royle and Link (2006) and
Miller et al. (2011) by maximum likelihood with
[`occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html).
crywolf fits the multiple-methods model of Miller et al. (2011) in a
Bayesian framework. This article shows how data and formulas prepared
for
[`occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html)
translate to
[`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md)
and [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
fits both to the same simulated data, and lists what differs between the
two.

## The data layout

[`occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html)
takes an `unmarkedFrameOccuFP`: one detection matrix `y` with sites in
rows and all occasions in columns, plus a vector `type` with the number
of occasions of each data type, in this order:

1.  occasions without false positives (the verified method);
2.  occasions with false positives (the automated method);
3.  occasions with certain and uncertain detections (detections coded 2
    and 1).

The multiple-methods design uses types 1 and 2 only, so
`type = c(n_verified, n_automated, 0)`. crywolf keeps the two methods in
separate matrices: `fp_data_wide(y_auto, y_verified)` takes the type-2
columns as `y_auto` and the type-1 columns as `y_verified`. Rows of `NA`
in `y_verified` mark sites that were not verified, as in ‘unmarked’.

| unmarked                        | crywolf                             |
|---------------------------------|-------------------------------------|
| `y[, 1:n1]` (type 1)            | `y_verified`                        |
| `y[, n1 + 1:n2]` (type 2)       | `y_auto`                            |
| `siteCovs`                      | `site_covs`                         |
| `obsCovs` (automated occasions) | `obs_covs` (named list of matrices) |
| `type = c(n1, n2, 0)`           | implied by the two matrices         |
| type-3 data                     | not supported                       |

[`fp_to_unmarked()`](https://provetelab.org/crywolf/reference/fp_to_unmarked.md)
goes the other way: it turns any `fp_data` object into the
`unmarkedFrameOccuFP` of
[`occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html),
with the verified columns first and a `METH` covariate that tells the
two methods apart, so that the same data can be fitted both ways.

## The same model fitted both ways

We simulate 300 sites with 10 automated occasions each, a third of them
verified on three occasions, and a constant false-positive rate of 0.05.

``` r

library(crywolf)
sim <- simulate_fp(nsite = 300, nocc = 10, frac_verified = 0.3, nver = 3,
                   seed = 2026)

# site-by-occasion matrices, the layout of 'unmarked'
to_wide <- function(df, nsite, nocc) {
  m <- matrix(NA_integer_, nsite, nocc)
  m[cbind(df$site, df$occasion)] <- df$y
  m
}
y_auto <- to_wide(sim$data$auto, 300, 10)
y_ver <- to_wide(sim$data$verified, 300, 3)
x <- sim$data$sites$x
```

In
[`occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html),
true detection on both methods is described by `detformula`, so the two
methods need an observation covariate that tells them apart. Here `METH`
is `"verified"` on the type-1 columns and `"auto"` on the type-2
columns; `"auto"` comes first alphabetically and becomes the reference
level.

``` r

library(unmarked)
meth <- matrix(rep(c("verified", "auto"), times = c(3, 10)),
               nrow = 300, ncol = 13, byrow = TRUE)
umf <- unmarkedFrameOccuFP(y = cbind(y_ver, y_auto),
                           siteCovs = data.frame(x = x),
                           obsCovs = list(METH = meth),
                           type = c(3, 10, 0))
#> Warning: obsCovs contains characters. Converting them to factors.
m_um <- occuFP(detformula = ~ METH, FPformula = ~ 1, stateformula = ~ x,
               data = umf)
m_um
#> 
#> Call:
#> occuFP(detformula = ~METH, FPformula = ~1, stateformula = ~x, 
#>     data = umf)
#> 
#> Occupancy (logit-scale):
#>             Estimate    SE     z  P(>|z|)
#> (Intercept)    0.156 0.125  1.25 2.12e-01
#> x             -0.766 0.139 -5.51 3.68e-08
#> 
#> Detection (logit-scale):
#>              Estimate     SE    z  P(>|z|)
#> (Intercept)     0.430 0.0526 8.17 3.12e-16
#> METHverified    0.234 0.1771 1.32 1.86e-01
#> 
#> false positive (logit-scale):
#>  Estimate    SE     z   P(>|z|)
#>     -2.99 0.131 -22.8 1.56e-115
#> 
#> AIC: 3266.285 
#> Number of sites: 300
```

In crywolf the verified method has its own detection probability, so no
method covariate is needed:

``` r

d <- fp_data_wide(y_auto, y_ver, site_covs = data.frame(x = x))
m_cw <- fp_occu(d, occ = ~ x, fp = ~ 1, threads_per_chain = 1, seed = 1)
#> Running MCMC with 4 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 3 finished in 45.7 seconds.
#> Chain 2 finished in 45.9 seconds.
#> Chain 4 finished in 47.2 seconds.
#> Chain 1 finished in 49.1 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 47.0 seconds.
#> Total execution time: 49.3 seconds.
coef(m_cw)
#>      occ_(Intercept)                occ_x      det_(Intercept) 
#>           0.15887315          -0.76977980           0.42722960 
#> verified_(Intercept)           p_verified              p11_ref 
#>           0.66219155           0.65878164           0.60514411 
#>       fp_(Intercept)               fp_gap              p10_ref 
#>          -2.98190153           3.40913113           0.04860622 
#>              p10_bar 
#>           0.04860622
```

Both fits recover the generating values (occupancy intercept 0, slope
−0.8, true detection logit(0.6) = 0.41 for both methods, false positives
logit(0.05) = −2.94), and with this much data the maximum-likelihood
estimates and the posterior means nearly coincide:

| parameter                         | truth | unmarked | crywolf |
|:----------------------------------|------:|---------:|--------:|
| occupancy intercept               |  0.00 |     0.16 |    0.16 |
| occupancy slope                   | -0.80 |    -0.77 |   -0.77 |
| true detection, automated (logit) |  0.41 |     0.43 |    0.43 |
| true detection, verified (logit)  |  0.41 |     0.66 |    0.66 |
| false positive (logit)            | -2.94 |    -2.99 |   -2.98 |

Maximum-likelihood estimates from occuFP() and posterior means from
fp_occu(), logit scale. {.table}

## Formulas

| [`occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html) | [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) | Notes |
|----|----|----|
| `stateformula` | `occ` | site covariates |
| `detformula` | `det` | true detection of the automated method; site or observation covariates |
| (method term in `detformula`) | `verified_(Intercept)` | crywolf gives the verified method a constant detection probability |
| `FPformula` | `fp` | crywolf: site covariates only |
| `Bformula` | none | type-3 data only |

Two differences in the parametrization matter when comparing output:

- In crywolf the intercept of the false-positive rate is tied to the
  intercept of true detection by a positive gap, \gamma_0 = \alpha_0 -
  \Delta with \Delta \> 0, so that at the reference values of the
  covariates a true detection is always more likely than a false
  positive. Without such a constraint the likelihood has a second mode
  in which the labels of occupied and unoccupied sites are swapped; an
  optimizer such as the one in
  [`occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html)
  returns whichever mode it reaches from its starting values, so check
  that the estimated true-detection probability exceeds the
  false-positive one. Because the reference is where all covariates are
  zero, centre or standardize them before fitting with crywolf.
- [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
  places priors on all parameters
  ([`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)),
  which matter little with as much data as above but keep the estimates
  finite in small or sparse datasets.

## What crywolf adds

- **Posterior distributions** of every parameter and of derived
  quantities such as the mean false-positive rate over sites
  (`p10_bar`), with convergence diagnostics.
- **A recorder random intercept** on true detection
  (`random_recorder = TRUE`).
- **The observation-confirmation design** of Chambert et al. (2015), in
  which detections on a subset of automated occasions are checked
  (`fp_data(confirmed = )`).
- **Cross-validation** by leave-one-site-out LOO
  ([`loo()`](https://mc-stan.org/loo/reference/loo.html)) in place of
  AIC, and **posterior predictive checks**
  ([`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md)).
- **Calibration checks by simulation**
  ([`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md),
  [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md)),
  which ask whether the model’s intervals are trustworthy for a given
  design and false-positive structure
  ([`vignette("calibration", package = "crywolf")`](https://provetelab.org/crywolf/articles/calibration.md)).

‘unmarked’ is the better choice when you need type-3 data,
false-positive covariates that vary by occasion, or a fast
maximum-likelihood fit inside a larger workflow built on its other
models.

## Predictions

Predictions follow the same logic in both packages:

``` r

nd <- data.frame(x = c(-1, 0, 1))
predict(m_um, type = "state", newdata = nd)
#>   Predicted         SE     lower     upper
#> 1 0.7154715 0.03888113 0.6336171 0.7852370
#> 2 0.5388560 0.03097511 0.4778700 0.5987006
#> 3 0.3519150 0.04165223 0.2751642 0.4371623
predict(m_cw, newdata = nd, process = "occ")
#>    x  estimate     lower     upper
#> 1 -1 0.7152729 0.6394499 0.7873743
#> 2  0 0.5394848 0.4822685 0.6003674
#> 3  1 0.3529699 0.2733537 0.4364194
```

## References

Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339. <https://doi.org/10.1890/14-1507.1>

Fiske, I., & Chandler, R. (2011). unmarked: An R package for fitting
hierarchical models of wildlife occurrence and abundance. *Journal of
Statistical Software*, 43(10), 1–23.
<https://doi.org/10.18637/jss.v043.i10>

Kellner, K. F., Smith, A. D., Royle, J. A., Kéry, M., Belant, J. L., &
Chandler, R. B. (2023). The unmarked R package: Twelve years of advances
in occurrence and abundance modelling in ecology. *Methods in Ecology
and Evolution*, 14, 1408–1415. <https://doi.org/10.1111/2041-210X.14123>

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
<https://doi.org/10.1890/10-1396.1>

Royle, J. A., & Link, W. A. (2006). Generalized site occupancy models
allowing for false positive and false negative errors. *Ecology*, 87,
835–841.
<https://doi.org/10.1890/0012-9658(2006)87%5B835:GSOMAF%5D2.0.CO;2>
