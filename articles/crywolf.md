# Getting started with crywolf

Automated detectors, such as the acoustic classifiers that scan
thousands of hours of recordings, make two kinds of mistakes: they miss
calls, and they label other sounds as the target species. Occupancy
models deal with the first kind by repeating surveys, but false
positives inflate occupancy and, worse, can distort the effects of
environmental covariates. **crywolf** fits the multiple-methods
false-positive occupancy model of Miller et al. (2011), in which a
subset of detections is checked by an expert, and helps you ask whether
the correction is good enough for your design.

## The design

At every site, the classifier scans several recordings (the *automated*
occasions). At a subset of sites, someone also listens to a few clips
and records whether the species was really there (the *verified*
occasions). Verified detections are certain, so they tell the model how
often the classifier is wrong.

## Installation

crywolf fits its models with [CmdStan](https://mc-stan.org/cmdstanr/)
through **cmdstanr**. Install both once:

``` r

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                       getOption("repos")))
cmdstanr::install_cmdstan()
remotes::install_github("diogoprov/crywolf")
```

``` r

library(crywolf)
```

## Your data

[`fp_data()`](https://diogoprov.github.io/crywolf/reference/fp_data.md)
takes three tables in long format: the automated detections (one row per
site and occasion), the sites with their covariates, and, optionally,
the verified detections. Here we simulate data instead, with a
false-positive rate that rises along a gradient correlated with the
occupancy covariate, as happens when a soundalike species becomes
commoner where the target is rare.

``` r

sim <- simulate_fp(nsite = 300, nocc = 20, psi = 0.5, b1 = -0.8,
                   p11 = 0.6, p10 = 0.05,
                   conf_effect = 1.5, conf_cor = 0.7,
                   frac_verified = 0.3, seed = 2026)
sim
#> <fp_sim> structured false positives
#>   truth: psi = 0.50, b1 = -0.80, p11 = 0.60, p10 = 0.050 (realized mean 0.106), r11 = 0.60
#>   confusion: effect 1.50, correlation with x 0.70 (realized 0.76)
#>   realized occupancy: 0.53
#> <fp_data>
#>   sites:                  300
#>   automated occasions:    6000 (naive occupancy 0.87)
#>   verified sites:         85 (255 occasions, 39 with a confirmed detection)
#>   recorders:              5
#>   site covariates:        x, w 
#>   observation covariates: bio, aci
head(sim$data$auto)
#>   site occasion y        bio         aci recorder
#> 1    1        1 0 -1.3818134 -0.36543963        5
#> 2    2        1 1  0.4355887  1.57117657        5
#> 3    3        1 0  0.3721785 -2.45643384        4
#> 4    4        1 0 -0.7030971  1.26171428        5
#> 5    5        1 1 -1.0114298 -1.02790936        4
#> 6    6        1 1 -1.6309700 -0.06734317        5
head(sim$data$sites)
#>   site          x           w
#> 1    1  0.5035356  0.40752869
#> 2    2 -1.1789366 -0.01491001
#> 3    3  0.1025980  0.66579653
#> 4    4 -0.1328932 -0.22992891
#> 5    5 -0.7446706 -1.44037323
#> 6    6 -2.6891101 -3.06909421
```

If your detections are in site-by-occasion matrices, as in **unmarked**,
[`fp_data_wide()`](https://diogoprov.github.io/crywolf/reference/fp_data_wide.md)
builds the same object.

## Fitting the model

[`fp_occu()`](https://diogoprov.github.io/crywolf/reference/fp_occu.md)
takes a formula for each process. With `fp = ~ 1` the false-positive
rate is constant, which is what most applications assume.

``` r

m_const <- fp_occu(sim$data, occ = ~ x, det = ~ 1, fp = ~ 1,
                   random_recorder = FALSE, seed = 1)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 2 finished in 3.6 seconds.
#> Chain 4 finished in 3.9 seconds.
#> Chain 3 finished in 4.2 seconds.
#> Chain 1 finished in 4.2 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 4.0 seconds.
#> Total execution time: 4.3 seconds.
m_const
#> Bayesian false-positive occupancy model (multiple-methods design)
#> 
#>   occupancy:       ~x 
#>   true detection:  ~1  
#>   false positive:  ~1 
#>   data:            300 sites, 6000 automated and 255 verified occasions
#> 
#>              variable    mean      sd   lower    upper rhat ess_bulk ess_tail
#>       occ_(Intercept)  0.4498 0.12468  0.2116  0.69666    1     4426     2919
#>                 occ_x -0.2553 0.12958 -0.5145 -0.00286    1     4307     2970
#>       det_(Intercept)  0.3506 0.03735  0.2783  0.42288    1     4655     3015
#>  verified_(Intercept) -0.0379 0.17200 -0.3740  0.28997    1     4521     3034
#>            p_verified  0.4906 0.04269  0.4076  0.57199    1     4521     3034
#>               p11_ref  0.5867 0.00905  0.5691  0.60417    1     4655     3015
#>        fp_(Intercept) -2.3967 0.09269 -2.5853 -2.22202    1     5456     3203
#>                fp_gap  2.7473 0.09158  2.5763  2.93488    1     4991     3308
#>               p10_ref  0.0837 0.00707  0.0701  0.09779    1     5456     3203
#>               p10_bar  0.0837 0.00707  0.0701  0.09779    1     5456     3203
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.003, min bulk/tail ESS 4307/2919
```

The true occupancy slope is −0.8. Because the false positives track a
gradient correlated with `x`, the constant-rate model attributes part of
the false positives to the occupancy process, and the slope is pulled
towards zero. The site where the classifier is fooled most looks like a
site where the species is present.

The practical remedy is to let the false-positive rate depend on the
occupancy covariate, which you always have:

``` r

m_x <- fp_occu(sim$data, occ = ~ x, det = ~ 1, fp = ~ x,
               random_recorder = FALSE, seed = 1)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 1 finished in 4.9 seconds.
#> Chain 2 finished in 5.1 seconds.
#> Chain 4 finished in 5.1 seconds.
#> Chain 3 finished in 5.3 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 5.1 seconds.
#> Total execution time: 5.4 seconds.
coef(m_const, "occ")
#> occ_(Intercept)           occ_x 
#>       0.4497632      -0.2553304
coef(m_x, "occ")
#> occ_(Intercept)           occ_x 
#>       0.2907979      -0.5595734
```

and, if you have measured a covariate of whatever fools the classifier,
adding it as well. In a simulation we know that variable, `w`, so we can
fit the generating model as a positive control:

``` r

m_w <- fp_occu(sim$data, occ = ~ x, det = ~ 1, fp = ~ w,
               random_recorder = FALSE, seed = 1)
#> Running MCMC with 4 parallel chains...
#> 
#> Chain 4 finished in 5.0 seconds.
#> Chain 2 finished in 5.2 seconds.
#> Chain 1 finished in 5.4 seconds.
#> Chain 3 finished in 5.7 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 5.3 seconds.
#> Total execution time: 5.8 seconds.
summary(m_w)
#> Posterior summary (95% equal-tailed intervals)
#> 
#>              variable    mean      sd    lower   upper rhat ess_bulk ess_tail
#>       occ_(Intercept)  0.1070 0.12816 -0.14845  0.3582    1     4902     2716
#>                 occ_x -0.7834 0.14643 -1.06963 -0.5012    1     4428     2897
#>       det_(Intercept)  0.4323 0.03794  0.35785  0.5080    1     4759     2744
#>  verified_(Intercept)  0.3571 0.17975  0.00826  0.7079    1     5243     3156
#>            p_verified  0.5877 0.04326  0.50207  0.6699    1     5243     3156
#>               p11_ref  0.6064 0.00905  0.58852  0.6243    1     4759     2744
#>        fp_(Intercept) -2.9041 0.10097 -3.10811 -2.7116    1     3614     3034
#>                  fp_w  1.4724 0.07767  1.32306  1.6257    1     3760     3075
#>                fp_gap  3.3364 0.10796  3.12574  3.5459    1     3410     2793
#>               p10_ref  0.0522 0.00499  0.04277  0.0623    1     3614     3034
#>               p10_bar  0.1072 0.00466  0.09826  0.1164    1     4598     3225
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.003, min bulk/tail ESS 3410/2716
```

## Comparing models

The pointwise log-likelihood is computed per site, with the latent state
marginalized, so [`loo()`](https://mc-stan.org/loo/reference/loo.html)
gives leave-one-site-out cross-validation.

``` r

loo::loo_compare(list(constant = loo(m_const), fp_x = loo(m_x), fp_w = loo(m_w)))
#>     model elpd_diff se_diff p_worse diag_diff diag_elpd
#>      fp_w       0.0     0.0      NA                    
#>      fp_x     -89.5    15.9    1.00                    
#>  constant    -135.9    18.2    1.00
```

## Predictions

[`predict()`](https://rdrr.io/r/stats/predict.html) returns the
occupancy probability, or the false-positive rate, with credible
intervals.

``` r

nd <- data.frame(x = seq(-2, 2, by = 1))
predict(m_x, newdata = nd, process = "occ")
#>    x  estimate     lower     upper
#> 1 -2 0.7994252 0.6964307 0.8827909
#> 2 -1 0.6992305 0.6237379 0.7708294
#> 3  0 0.5719128 0.5110966 0.6314804
#> 4  1 0.4338290 0.3438219 0.5285994
#> 5  2 0.3080601 0.1877566 0.4461158
predict(m_x, newdata = nd, process = "fp")
#>    x    estimate       lower      upper
#> 1 -2 0.009661468 0.005558678 0.01561041
#> 2 -1 0.025501627 0.017722220 0.03492059
#> 3  0 0.066168481 0.054623116 0.07795714
#> 4  1 0.162114250 0.139852866 0.18448956
#> 5  2 0.346338560 0.278565834 0.41050497
```

## Covariates on detection and recorders

The true-detection formula accepts observation covariates (columns of
`auto`) and site covariates, and `random_recorder = TRUE` adds a random
intercept per recorder.

``` r

fp_occu(sim$data, occ = ~ x, det = ~ bio + aci, fp = ~ x,
        random_recorder = TRUE)
```

## Is the correction calibrated for my design?

A single fit cannot tell you whether the credible intervals are honest.
[`fp_calibrate()`](https://diogoprov.github.io/crywolf/reference/fp_calibrate.md)
and
[`fp_verification_curve()`](https://diogoprov.github.io/crywolf/reference/fp_verification_curve.md)
answer that by simulation; see
[`vignette("calibration", package = "crywolf")`](https://diogoprov.github.io/crywolf/articles/calibration.md).

## References

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
<https://doi.org/10.1890/10-1396.1>
