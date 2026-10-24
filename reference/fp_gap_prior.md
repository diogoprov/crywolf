# A prior for the gap from validated detections

Turns the precision of the classifier, measured on reviewed detections,
into a prior for the gap \\\Delta\\ between the true-detection and
false-positive intercepts, and returns it in
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md).

## Usage

``` r
fp_gap_prior(n_true, n_false, p11, psi, ndraws = 4000, ...)
```

## Arguments

- n_true, n_false:

  Numbers of reviewed detections confirmed as true and as false. They
  must be a random sample of the detections at the threshold of the
  analysis; detections picked for their high confidence overstate the
  precision.

- p11:

  Probability that an occasion at an occupied site is a detection, at
  the reference values of the covariates. A single value, or a vector of
  draws that expresses its uncertainty.

- psi:

  Occupancy probability at the reference values of the covariates. A
  single value, or a vector of draws.

- ndraws:

  Number of draws used to propagate the uncertainty.

- ...:

  Further arguments for
  [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md).

## Value

An `fp_priors` object; the derived values are in the attribute
`"gap_draws"`.

## Details

In the model, an occasion at an occupied site is a detection with
probability \\p\_{11}\\, which includes false alarms; an occasion at an
unoccupied site, with probability \\p\_{10}\\. Assuming that false
alarms occur at the same rate \\p\_{10}\\ at occupied sites,
independently of true detections, the probability that an occasion holds
a true detection is \\s = (p\_{11} - p\_{10}) / (1 - p\_{10})\\, and the
share of automated detections that are true is \$\$\pi = \frac{\psi
s}{\psi p\_{11} + (1 - \psi) p\_{10}}.\$\$ For given \\\pi\\,
\\p\_{11}\\ and \\\psi\\, this has a single solution \\p\_{10}\\ between
0 and \\p\_{11}\\, the smaller root of a quadratic. The function draws
\\\pi\\ from its Beta(1 + `n_true`, 1 + `n_false`) posterior, solves for
\\p\_{10}\\, computes \\\Delta = \mathrm{logit}(p\_{11}) -
\mathrm{logit}(p\_{10})\\, and uses the mean and standard deviation of
the result as the parameters of the Normal prior on \\\Delta\\. The
model truncates that prior at zero; when most of the derived values are
far from zero, as is usual, the truncation changes little. Guesses of
`p11` and `psi` usually come from a pilot fit or from earlier studies;
give them as vectors when they are uncertain, so that the prior widens
accordingly.

## See also

[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md),
[`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)

## Examples

``` r
# 46 of 50 randomly drawn detections were true; p11 and psi from a pilot
fp_gap_prior(n_true = 46, n_false = 4, p11 = 0.3, psi = 0.4)
#> <fp_priors> (logit scale)
#>   intercepts      ~ Normal(0, 1.5)
#>   slopes          ~ Normal(0, 2.5)
#>   fp_gap          ~ Normal(3.50353, 0.478711), truncated at 0
#>   sd_recorder     ~ half-Normal(0, 2)
#>   sd_fp_site      ~ half-Normal(0, 1)
#>   sd_det_site     ~ half-Normal(0, 1)
# with uncertainty in p11
fp_gap_prior(46, 4, p11 = rbeta(1000, 30, 70), psi = 0.4)
#> <fp_priors> (logit scale)
#>   intercepts      ~ Normal(0, 1.5)
#>   slopes          ~ Normal(0, 2.5)
#>   fp_gap          ~ Normal(3.50994, 0.488533), truncated at 0
#>   sd_recorder     ~ half-Normal(0, 2)
#>   sd_fp_site      ~ half-Normal(0, 1)
#>   sd_det_site     ~ half-Normal(0, 1)
```
