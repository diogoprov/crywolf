# Posterior predictive check for a false-positive occupancy model

Replicates the detection data from the posterior of an
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) fit
and compares observed and replicated data with a discrepancy statistic,
giving a Bayesian p-value. Detections are summed per site, so the check
asks whether the model reproduces how many detections each site gets.

## Usage

``` r
fp_ppc(
  object,
  method = c("automated", "verified"),
  fit_stat = c("freeman-tukey", "chi-squared"),
  ndraws = 1000,
  seed = NULL
)
```

## Arguments

- object:

  An `fp_occu` object; an `fp_occu_dyn` object, for which detections are
  summed per site and season and the latent states are drawn jointly
  over seasons by forward filtering, backward sampling; or an
  `fp_occu_ms` object, for which detections are summed per site and
  species, the discrepancy is summed over species, and each species also
  gets its own p-value.

- method:

  Which data to check: `"automated"` (default; every site) or
  `"verified"` (verified sites only).

- fit_stat:

  Discrepancy statistic: `"freeman-tukey"` (default), \\\sum_i
  (\sqrt{y_i} - \sqrt{E_i})^2\\, or `"chi-squared"`, \\\sum_i (y_i -
  E_i)^2 / E_i\\, where \\y_i\\ is the number of detections at site
  \\i\\ and \\E_i\\ its expected value given the parameters and the
  latent state.

- ndraws:

  Number of posterior draws to use, sampled without replacement. `NULL`
  uses all of them.

- seed:

  Optional integer seed for the subsampling and replication.

## Value

An object of class `fp_ppc`, a list with `y` (observed detections per
site), `yrep` (matrix of replicated detections, draws in rows and sites
in columns), `fit_obs` and `fit_rep` (discrepancy of the observed and
replicated data, one value per draw), `p_value`, `site` (the site
identifiers of the columns), `method`, and `fit_stat`; for a dynamic
model, also `season`, the season of each column; for a multispecies
model, also `species`, the species of each column, and `p_species`, the
Bayesian p-value of each species.

## Details

The latent occupancy state is marginalized in the Stan model, so for
each posterior draw \\z_i\\ is first drawn from its conditional
posterior given the data and the parameters, \$\$\Pr(z_i = 1 \mid
\mathbf{y}\_i, \theta) = \frac{\psi_i\\ f(\mathbf{y}\_i \mid z_i = 1,
\theta)}{\psi_i\\ f(\mathbf{y}\_i \mid z_i = 1, \theta) + (1 - \psi_i)\\
f(\mathbf{y}\_i \mid z_i = 0, \theta)},\$\$ which is 1 at sites with a
verified detection. Given \\z_i\\, automated detections are replicated
with probability \\z_i p\_{11} + (1 - z_i) p\_{10}\\ and verified
detections with probability \\z_i r\_{11}\\. The discrepancy is computed
for the observed and the replicated data with the same \\z\\ and
parameters, and the Bayesian p-value is the proportion of draws in which
the replicated discrepancy exceeds the observed one. Values near 0 or 1
indicate lack of fit; values near 0.5 indicate that the observed data
are typical of the fitted model.

Under the observation-confirmation design, the conditional probability
of \\z_i\\ also uses the confirmed occasions, and `"automated"` checks
all automated occasions, confirmed or not, as detections or
non-detections. Confirmed occasions are usually chosen among the
detections, so the unconfirmed ones alone would hold too few detections
for the model to reproduce.

A p-value near 0.5 does not show that the false-positive process is
correctly specified: a check of detection counts can have little power
against a false-positive rate that tracks an unmeasured gradient. Use
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
to ask whether the intervals of a given model are calibrated for a
design.

The elements `y` and `yrep` follow the conventions of the 'bayesplot'
package, so its `ppc_*()` functions can be used directly, for example
`bayesplot::ppc_bars(ppc$y, ppc$yrep)`.

## See also

[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)

## Examples

``` r
# \donttest{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
             chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 21.7 seconds.
#> Chain 2 finished in 21.9 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 21.8 seconds.
#> Total execution time: 22.1 seconds.
#> 
ppc <- fp_ppc(m, ndraws = 500, seed = 1)
ppc
#> Posterior predictive check (automated detections per site)
#> 
#>   statistic:          freeman-tukey
#>   draws, sites:       500, 200
#>   mean discrepancy:   observed 55.8, replicated 49.1
#>   Bayesian p-value:   0.026
#> 
#>   The p-value is far from 0.5: the model may not fit these data.
bayesplot::ppc_bars(ppc$y, ppc$yrep)

# }
```
