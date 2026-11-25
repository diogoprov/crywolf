# Model weights by stacking

Weights of several fitted models for combining their predictions, by
stacking of their leave-one-site-out predictive densities (Yao et al.
2018), or by pseudo-BMA. A thin wrapper of
[`loo::loo_model_weights()`](https://mc-stan.org/loo/reference/loo_model_weights.html)
that accepts fits or
[`loo()`](https://mc-stan.org/loo/reference/loo.html) results.

## Usage

``` r
fp_loo_weights(..., method = c("stacking", "pseudobma"))
```

## Arguments

- ...:

  Named `fp_occu` or `fp_occu_rn` fits, or the results of
  [`loo()`](https://mc-stan.org/loo/reference/loo.html) on them.

- method:

  `"stacking"` (default) or `"pseudobma"`.

## Value

A vector of weights summing to one, named by the arguments in `...`.

## References

Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using stacking
to average Bayesian predictive distributions. *Bayesian Analysis*, 13,
917–1007. [doi:10.1214/17-BA1091](https://doi.org/10.1214/17-BA1091)

## Examples

``` r
# \donttest{
sim <- simulate_fp(nsite = 200, nocc = 10, conf_effect = 1.5,
                   conf_cor = 0.7, seed = 1)
m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
              chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 20.7 seconds.
#> Chain 2 finished in 26.5 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 23.6 seconds.
#> Total execution time: 26.6 seconds.
#> 
m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE,
              chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 23.2 seconds.
#> Chain 2 finished in 25.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 24.4 seconds.
#> Total execution time: 25.8 seconds.
#> 
fp_loo_weights(constant = m0, fp_on_x = m1)
#>  constant   fp_on_x 
#> 0.1786987 0.8213013 
# }
```
