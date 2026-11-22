# Site-level estimates

For every site: the occupancy probability given its covariates, the
probability that it is occupied given its covariates and its data
(detections and verification), and its false-positive rate. The second
is the quantity to map or rank sites by: a site with automated
detections but a high false-positive rate can still be unlikely to be
occupied.

## Usage

``` r
fp_sites(object, prob = 0.95, ndraws = NULL, seed = NULL)

# S3 method for class 'fp_sites'
plot(x, coords = NULL, ...)
```

## Arguments

- object:

  An `fp_occu` object, an `fp_occu_dyn` object or an `fp_occu_ms` object
  (see Details).

- prob:

  Width of the credible intervals.

- ndraws:

  Number of posterior draws to use (all by default).

- seed:

  Optional seed for the choice of draws.

- x:

  An `fp_sites` object.

- coords:

  Optional names of two columns of `x` with the coordinates of the sites
  (x then y, e.g. `c("lon", "lat")`). With them the plot is a map of the
  probability that each site is occupied; without them, the sites are
  ranked by it.

- ...:

  Unused.

## Value

A data frame of class `fp_sites`, one row per site (per site-season for
a dynamic model, per site and species for a multispecies model), with
the site columns of the data, `psi` (occupancy probability, with
`psi_lower` and `psi_upper`), `occupied` (posterior probability that the
site is occupied given its data), `occupied_lower` and `occupied_upper`
(quantiles over the posterior draws of that conditional probability;
they show how much it depends on the uncertain parameters, not the
uncertainty of the occupancy state itself, which `occupied` already
gives), `p10` (posterior mean false-positive rate, when the model has
false positives), `detections` (number of automated detections), and
`verified` (whether a verified detection proves occupancy).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it.

## Details

For a dynamic model
([`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md))
the result has one row per site and season, with the site and season
columns and the season covariates. `psi` is then the probability of
occupancy in that season given the covariates alone (initial occupancy
carried forward by colonization and extinction), and `occupied` the
probability that the site-season is occupied given all the data of the
site, in every season (computed by the forward-backward algorithm for
each posterior draw). `detections`, `verified` and `p10` refer to the
site-season.

For a multispecies model
([`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md))
the result has one row per site and species, with the species column of
the data after the site columns, and every column refers to that species
at that site.

## Examples

``` r
# \donttest{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE,
             chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 23.3 seconds.
#> Chain 2 finished in 23.5 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 23.4 seconds.
#> Total execution time: 23.6 seconds.
#> 
st <- fp_sites(m)
head(st)
#>   site          x w       psi psi_lower psi_upper     occupied occupied_lower
#> 1    1 -0.7125125 0 0.6068485 0.5183068 0.6960445 0.0003687517   1.353295e-04
#> 2    2  0.1594060 0 0.4958390 0.4200741 0.5672305 0.0003008647   6.733963e-05
#> 3    3 -0.9376502 0 0.6339738 0.5321936 0.7281057 0.0004153213   1.498936e-04
#> 4    4  1.6787706 0 0.3119244 0.1911776 0.4483337 0.9999993569   9.999953e-01
#> 5    5  0.3164019 0 0.4755301 0.3960121 0.5506958 0.9999999442   9.999997e-01
#> 6    6 -0.9213331 0 0.6320399 0.5314427 0.7259646 0.0004117265   1.484120e-04
#>   occupied_upper detections verified        p10
#> 1   0.0007946070          0    FALSE 0.05184956
#> 2   0.0008708825          1    FALSE 0.05422322
#> 3   0.0009027371          0    FALSE 0.05138950
#> 4   0.9999999978          8    FALSE 0.06091389
#> 5   0.9999999973          8    FALSE 0.05475647
#> 6   0.0008908884          0    FALSE 0.05142083
plot(st)

# }
```
