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
if (FALSE) { # \dontrun{
s <- simulate_fp_ms(nspecies = 4, nsite = 60, seed = 1)
m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE)
predict(m, newdata = data.frame(x = c(-1, 0, 1)))
} # }
```
