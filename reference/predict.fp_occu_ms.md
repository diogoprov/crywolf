# Predicted probabilities of a fitted multispecies model

The occupancy probability or the false-positive rate of each species of
an
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md)
fit, from the coefficients of that species, as
[`predict.fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
gives them for one species. Random effects on the false-positive rate
(`fp_group`, `random_fp`) are left out, as in
[`predict.fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu-methods.md).
**In development**, as
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md).

## Usage

``` r
# S3 method for class 'fp_occu_ms'
predict(
  object,
  newdata,
  process = c("occ", "fp"),
  species = NULL,
  prob = 0.95,
  summary = TRUE,
  ...
)
```

## Arguments

- object:

  An `fp_occu_ms` object.

- newdata:

  Data frame with the covariates of the requested process. Without a
  column of species (the species column of the data), every row is
  predicted for every species in `species`; with it, each row for its
  own species. If missing, the sites used in the fit, for every species;
  for a false-positive rate with covariates measured on each occasion,
  the automated occasions of each species, with the covariates of their
  sites.

- process:

  `"occ"` for the occupancy probability or `"fp"` for the false-positive
  rate, both at the site level.

- species:

  Names of the species to predict when `newdata` has no column of
  species. The default, `NULL`, is every species.

- prob:

  Width of the equal-tailed credible interval.

- summary:

  Logical. Return a summary per row (`TRUE`) or the matrix of posterior
  draws (rows are draws, columns are rows of `newdata`)?

- ...:

  Unused, or passed on (see each method).

## Value

A data frame with the rows of `newdata` (repeated for each species when
it has no column of species), the species column, and the posterior mean
and credible interval of the probability; or, with `summary = FALSE`,
the matrix of draws, with rows as draws and columns as the rows of that
data frame.

## See also

[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md),
[`predict.fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)

## Examples

``` r
# \donttest{
s <- simulate_fp_ms(nspecies = 4, nsite = 60, seed = 1)
m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE,
                chains = 2, threads_per_chain = 1)
#> Parametrization (4 species at 60 sites): centered = FALSE, centered_occ = FALSE; see ?fp_occu_ms.
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 100.3 seconds.
#> Chain 2 finished in 108.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 104.5 seconds.
#> Total execution time: 108.8 seconds.
#> 
#> Warning: Largest R-hat is 1.013 (> 1.01); run longer chains.
predict(m, newdata = data.frame(x = c(-1, 0, 1)))
#>     x species  estimate     lower     upper
#> 1  -1     sp1 0.6258650 0.4211492 0.8241446
#> 2   0     sp1 0.6047237 0.4508931 0.7697389
#> 3   1     sp1 0.5797430 0.3880339 0.7635818
#> 4  -1     sp2 0.3987673 0.2216843 0.5537199
#> 5   0     sp2 0.4374394 0.3033104 0.5580321
#> 6   1     sp2 0.4789149 0.3065087 0.6602128
#> 7  -1     sp3 0.4504723 0.2919431 0.5988758
#> 8   0     sp3 0.4330593 0.3086327 0.5486169
#> 9   1     sp3 0.4174396 0.2597201 0.5787053
#> 10 -1     sp4 0.4878284 0.3390190 0.6395466
#> 11  0     sp4 0.4898734 0.3762007 0.6003460
#> 12  1     sp4 0.4920795 0.3512278 0.6434019
# }
```
