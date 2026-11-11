# Species richness of each site from a multispecies model

Counts, for each site, the species of an
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md)
model that are present, from the posterior of the occupancy state of
every species at the site. **In development**, as
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md).

## Usage

``` r
fp_richness(object, prob = 0.95, ndraws = NULL, seed = NULL)
```

## Arguments

- object:

  An `fp_occu_ms` object.

- prob:

  Width of the credible intervals.

- ndraws:

  Number of posterior draws to use (all by default).

- seed:

  Optional seed for the choice of draws.

## Value

A data frame, one row per site, with the site columns of the data,
`richness` (posterior mean number of species present, with
`richness_lower` and `richness_upper`, quantiles of its posterior
distribution), `richness_psi` (expected number of species from the
occupancy covariates alone, before the data of the site), `detected`
(number of species with automated detections) and `verified` (number of
species whose presence a verified detection proves).

## Details

For each posterior draw, the state of each species at each site is drawn
from its conditional posterior given the data of that species at that
site (as in
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)),
and the richness of the site is the number of species drawn as present.
The species are those of the model: a species never detected anywhere is
not in the list and is not counted.

## See also

[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md),
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)

## Examples

``` r
if (FALSE) { # \dontrun{
s <- simulate_fp_ms(nspecies = 5, nsite = 60, seed = 1)
m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE)
r <- fp_richness(m)
plot(rowSums(s$truth$z), r$richness)
} # }
```
