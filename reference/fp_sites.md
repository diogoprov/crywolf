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

  An `fp_occu` object.

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

A data frame of class `fp_sites`, one row per site, with the site
columns of the data, `psi` (occupancy probability, with `psi_lower` and
`psi_upper`), `occupied` (posterior probability that the site is
occupied given its data), `occupied_lower` and `occupied_upper`
(quantiles over the posterior draws of that conditional probability;
they show how much it depends on the uncertain parameters, not the
uncertainty of the occupancy state itself, which `occupied` already
gives), `p10` (posterior mean false-positive rate, when the model has
false positives), `detections` (number of automated detections), and
`verified` (whether a verified detection proves occupancy).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
st <- fp_sites(m)
head(st)
plot(st)
} # }
```
