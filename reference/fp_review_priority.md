# Review priority of each site, for a directed review round of sheepdog

Gives every site (every site and species, for a multispecies model) a
size for
[`sheepdog::select_directed()`](https://provetelab.org/sheepdog/reference/select_directed.html):
the variance of its occupancy state given the data, \\q(1 - q)\\, where
\\q\\ is the posterior probability that the site is occupied (`occupied`
of
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)).
The size is largest, 0.25, where the model is least sure whether the
species is present, and 0 where the data settle it. sheepdog then draws
units with probability proportional to size and reviews one occasion of
each.

## Usage

``` r
fp_review_priority(
  object,
  species_code = NULL,
  validation = NULL,
  ndraws = NULL,
  seed = NULL
)
```

## Arguments

- object:

  An `fp_occu` or `fp_occu_ms` object.

- species_code:

  The label of the species in the predictions of sheepdog (its
  `species_code`, such as `"Nyctidromus albicollis_Common Pauraque"`).
  For an `fp_occu` model, a single label, required: without it, sheepdog
  would draw occasions of every species. For an `fp_occu_ms` model,
  optional: a character vector named by the species of the model, for
  when the model used other names (such as the short `species` of
  [`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md));
  by default the species names of the model are taken as the labels.

- validation:

  Optional data frame of reviewed clips, such as
  `fp_read_sheepdog()$validation`, with the columns `site`, `verdict`
  (`TRUE` for a confirmed detection) and, for a multispecies model,
  `species_code`. Units with a confirmed detection are left out.

- ndraws:

  Number of posterior draws to use (all by default).

- seed:

  Optional seed for the choice of draws.

## Value

A data frame with one row per unit, sorted by decreasing size: `site`,
`species_code` and `size`, the columns that
[`sheepdog::select_directed()`](https://provetelab.org/sheepdog/reference/select_directed.html)
expects in its `units`. Sites and labels are text.

## Details

**Which units.** The units are the sites, or the sites and species, the
level of the occupancy state; sheepdog draws one occasion among the
candidate occasions of each drawn unit. Units where the species is
already known to be present are left out: those with a verified
detection in the data of the fit, and those with a confirmed clip in
`validation` (which may hold rounds reviewed after the fit). Without
that, the floor of sheepdog (`floor`) would still draw them.

**Why the variance.** The size ranks units by how uncertain their state
is, not by the exact gain of one review, which would also depend on the
detection and false-positive rates of the occasion drawn. It is a simple
rule, and it needs nothing beyond
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md).

**Fitting a model to the directed rounds.** The sizes depend only on
data already in the fit, and sheepdog adds a floor, draws units at
random given the sizes, and draws the occasion among those of the unit
with a prediction at or above its `threshold`, detections the fit also
holds when it uses the same threshold. A model fit to the data of all
rounds (the earlier ones and the directed one) therefore conditions on
everything that chose the units, and the selection is ignorable for that
model:
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) can
use the clips without weighting them. This holds only when the size of
every directed round came from a model of data that the new fit also
holds. Design-based estimates from the clips alone, such as the
precision of the classifier, still need the weights `1 / inclusion_prob`
that sheepdog records.

## Experimental

This function is experimental (new after 2.3.0, part of the integration
with 'sheepdog'): its arguments and results may change in a minor
version, as
[`sheepdog::select_directed()`](https://provetelab.org/sheepdog/reference/select_directed.html)
settles. See
[crywolf-api](https://provetelab.org/crywolf/reference/crywolf-api.md).

## See also

[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md),
[`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md)

## Examples

``` r
# \donttest{
# the small export of sheepdog shipped with crywolf
dir <- system.file("extdata", "sheepdog-export", package = "crywolf")
x <- suppressWarnings(fp_read_sheepdog(dir, species = "phy_cuv"))
auto <- fp_aggregate(x$occasions, x$detections, threshold = 0.5,
                     by = c("site", "occasion_id"), validation = x$validation)
m <- fp_occu(fp_data(auto, x$sites), random_recorder = FALSE,
             chains = 2, threads_per_chain = 1)
#> Warning: No verified data: the false-positive rate is identified only by the priors. Treat the estimates with caution.
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.7 seconds.
#> Chain 2 finished in 0.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.7 seconds.
#> Total execution time: 0.9 seconds.
#> 
units <- fp_review_priority(m, species_code = unique(x$detections$species_code),
                            validation = x$validation)
head(units)
#>   site                species_code       size
#> 1  BAI Physalaemus cuvieri_phy-cuv 0.24625664
#> 2  ARA Physalaemus cuvieri_phy-cuv 0.23085412
#> 3  XAR Physalaemus cuvieri_phy-cuv 0.03659857
#> 4  NHE Physalaemus cuvieri_phy-cuv 0.01294631
# then, in sheepdog, with the predictions of the survey:
# sel <- sheepdog::select_directed(pred, units, n_units = 10, threshold = 0.5)
# }
```
