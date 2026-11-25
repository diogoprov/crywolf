# Detection that depends on abundance: the Royle-Nichols model

In [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
the probability of detecting the species at an occupied site is the same
function of covariates at every site. In reality, a site with many
calling individuals is detected more often than a site with one.
Variation in detection among sites then has two possible sources:
abundance, and everything else (distance to the recorder, background
noise, vegetation).
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) with
`random_det = TRUE` captures the sum of both, but cannot tell them
apart.

[`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md)
fits the model of Royle and Nichols (2003), in which detection depends
on abundance, with the false positives and verification of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md). A
site effect on the per-individual detection rate then measures the
variation in detection that abundance does not explain. This article
describes the model, shows how to fit and check it, and reports the
simulations that show what it can and cannot recover. The fitting code
is shown but not run when the article is built.

## The model

The number of individuals at site i is N_i \sim
\mathrm{Poisson}(\lambda_i), with \log \lambda_i = X_i \beta + v_i. The
site is occupied when N_i \> 0, so the occupancy probability is \psi_i =
1 - e^{-\lambda_i}. Each individual is detected on an occasion with
probability r, and the automated method detects at least one with
probability 1 - (1 - r\_{ik})^{N_i}, where \mathrm{logit}(r\_{ik}) =
Z\_{ik} \alpha + d_i. False positives occur on every automated occasion
with probability p\_{10}, also at occupied sites. An occasion is a
detection when either a true or a false observation occurs: \Pr(y\_{ik}
= 1 \mid N_i) = 1 - (1 - r\_{ik})^{N_i} (1 - p\_{10,i}). A verified
occasion has no false positives and is a detection with probability 1 -
(1 - r^{v})^{N_i}. The abundance N_i is summed out of the likelihood
over 0 to `n_max`.

The two site effects are the point of the model:

- v_i, with standard deviation `sd_abund_site` (`random_abund = TRUE`),
  is variation in abundance beyond the covariates of `occ`.
- d_i, with standard deviation `sd_det_site` (`random_det = TRUE`), is
  variation in the per-individual detection rate.

Both raise the variance of detections among occupied sites. They differ
in what else they predict. Abundance also changes the share of occupied
sites and the detection by the verified method. A detection effect
changes only the automated detections. Verified or confirmed occasions
are therefore what separates the two.

There are two ways to provide that information:

- **Verified occasions**: detections by an independent method without
  false positives, at some sites (the `verified` table of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)).
- **Confirmed occasions**: automated occasions whose detections were
  reviewed and split into true and false observations (the `confirmed`
  argument of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
  the observation-confirmation design of Chambert et al. 2015). A fully
  annotated survey confirms every occasion.

Unlike
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md), the
false-positive intercept is not tied to the detection intercept by a
gap. Detection at an occupied site always exceeds p\_{10}, so the labels
of true and false detections cannot switch.

## Fitting and checking

[`simulate_fp_rn()`](https://provetelab.org/crywolf/reference/simulate_fp_rn.md)
generates data from the model. Here 100 sites with 15 automated
occasions each, verified occasions at 30% of the sites, and a site
effect on detection:

``` r

s <- simulate_fp_rn(nsite = 100, nocc = 15, psi = 0.6, b1 = 0.5, r = 0.3,
                    p10 = 0.05, sigma_r = 1, sigma_lam = 0.5,
                    frac_verified = 0.3, nver = 3, seed = 1)
```

The fit with both site effects:

``` r

m <- fp_occu_rn(s$data, occ = ~ x, random_det = TRUE, random_abund = TRUE,
                n_max = 25, iter_warmup = 1500,
                iter_sampling = 1500, seed = 1)
m
```

Three settings matter more than in
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md):

- `n_max` must be large enough that the Poisson mass above it is
  negligible at every site. The fit warns when it is not, from the
  derived quantity `p_above_nmax`. The cost of the fit grows with
  `n_max`.
- The default `adapt_delta` is 0.99, above the 0.95 of
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  The site effects create funnels, and at 0.95 divergences appeared in
  most designs of the simulations below.
- With many sites, the standard deviations of the site effects mix
  slowly. At 100 sites, 1000 + 1000 iterations left a largest R-hat of
  1.01-1.03 in half or more of the simulated fits. Longer chains are the
  remedy.

Before the estimates, check convergence (R-hat, bulk ESS and divergences
in `m$diagnostics`) and the fit to the data:

``` r

fp_ppc(m)                        # automated detections per site
fp_ppc(m, method = "verified")   # verified detections per site
```

[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)
gives, for each site, the posterior probability that it is occupied and
its expected abundance given its data:

``` r

head(fp_sites(m))   # occupied, abundance, abundance_lower, abundance_upper, ...
```

[`loo()`](https://mc-stan.org/loo/reference/loo.html) estimates the
leave-one-site-out predictive density, with the site effects integrated
out of each site’s likelihood. Without that integration most sites get
high Pareto k. The comparison of interest is usually a model with and
without the detection effect:

``` r

m0 <- update(m, random_det = FALSE)
fp_loo_weights(with_det = m, without_det = m0)
```

## Calibration

Simulation-based calibration (Talts et al. 2018) checks that the program
fits the model the simulator generates. Parameters are drawn from the
priors, data are simulated from them, and the rank of each true value
among the posterior draws should be uniform. Two designs, 500 datasets
each, had 50 sites and 10 occasions: design A confirmed every occasion
and had no verified data, and design B confirmed 20% of the occasions
and verified 3 occasions at 30% of the sites.

| design | parameter      | datasets | p (uniformity) | mean rank (of 0-100) |
|:-------|:---------------|---------:|---------------:|---------------------:|
| A      | alpha\[1\]     |      500 |          0.280 |                 49.7 |
| A      | beta\[1\]      |      500 |          0.260 |                 50.6 |
| A      | beta\[2\]      |      500 |          0.028 |                 45.8 |
| A      | gamma0         |      500 |          0.370 |                 51.7 |
| A      | sigma_det\[1\] |      500 |          0.900 |                 49.7 |
| A      | sigma_lam\[1\] |      500 |          0.260 |                 49.7 |
| B      | alpha\[1\]     |      500 |          0.650 |                 50.3 |
| B      | beta\[1\]      |      500 |          0.830 |                 50.3 |
| B      | beta\[2\]      |      500 |          0.070 |                 48.6 |
| B      | delta0         |      500 |          0.810 |                 50.4 |
| B      | gamma0         |      500 |          0.850 |                 48.7 |
| B      | sigma_det\[1\] |      500 |          0.480 |                 48.0 |
| B      | sigma_lam\[1\] |      500 |          0.210 |                 51.1 |

The ranks are close to uniform. The exception is the slope of log
abundance (`beta[2]`) in design A. Its mean rank is 45.8, about 3
standard errors below 50, and its posterior median lies on average 0.02
above the true value. With 13 tests, a single p-value of 0.03 can arise
by chance, and the departure is small. Divergences occurred in 0.2%-0.4%
of the fits.

## Recovery across designs

The general grid crossed the number of sites (25, 50, 100), occasions
(10, 20), the annotation scheme, \sigma_r (0 or 1) and \sigma\_\lambda
(0 or 0.7), with 100 datasets per cell. The annotation schemes were:

- **verified**: 3 verified occasions at 30% of the sites;
- **partial**: the same, plus 20% of automated occasions confirmed;
- **full**: every occasion confirmed, no verified data.

Other values were fixed (\psi = 0.6, slope 0.5, r = 0.3, p\_{10} =
0.05), and the fits used the default priors. A signal is a lower 5%
bound of `sd_det_site` above 0.1. Where \sigma_r = 0, a signal is false.
Where \sigma_r = 1, it is power.

| sites x occasions | power: verified | partial | full | false: verified | partial | full |
|:------------------|:----------------|:--------|:-----|:----------------|:--------|:-----|
| 25 x 10           | 51%             | 46%     | 55%  | 11%             | 7%      | 12%  |
| 25 x 20           | 51%             | 59%     | 70%  | 6%              | 8%      | 8%   |
| 50 x 10           | 68%             | 69%     | 70%  | 7%              | 5%      | 12%  |
| 50 x 20           | 79%             | 84%     | 80%  | 7%              | 7%      | 9%   |
| 100 x 10          | 83%             | 86%     | 84%  | 7%              | 11%     | 11%  |
| 100 x 20          | 93%             | 93%     | 94%  | 5%              | 6%      | 1%   |

Signal for a site effect on detection, with sigma_lambda = 0.7. {.table
style="width:100%;"}

Three results stand out:

- Power grows with the number of sites and occasions, from about half at
  25 sites and 10 occasions to above 90% at 100 sites and 20 occasions.
- Beyond the smallest design, the annotation scheme matters little.
  Verifying 30% of the sites gave about the same power as confirming
  every occasion. At 25 sites and 10 occasions without variation in
  abundance, verified data alone gave less power (41%) than partial
  confirmation (54%).
- Variation in abundance raises the false signal only moderately.
  Without it (\sigma\_\lambda = 0) the false signal was 0%-6%; with
  \sigma\_\lambda = 0.7 it was 1%-12%.

The other quantities were well recovered. The 90% intervals covered the
true value in 85%-97% of the fits for p\_{10}, 81%-99% for r and 85%-99%
for \psi. The standard deviation of abundance was weakly identified.
When it was 0, its posterior median was 0.24-0.63, the scale of its
prior.

## Few sites and the prior

A second grid had 21 sites with 11-19 occasions, every occasion
confirmed and no verified data, at three occupancy levels. Each dataset
was fitted under four priors on `sd_det_site`: half-Normal with scale
0.5, 1 (the default) and 2, and half-Student-t with 3 degrees of
freedom.

| prior on sd_det_site | false signal | posterior median |
|:---------------------|:-------------|:-----------------|
| half-Normal(0, 0.5)  | 0%-1%        | 0.3-0.39         |
| half-Normal(0, 1)    | 11%-31%      | 0.46-0.59        |
| half-Normal(0, 2)    | 20%-43%      | 0.55-0.7         |
| half-t(3, 0, 1)      | 9%-30%       | 0.44-0.58        |

21 sites, true sd_det_site = 0, over occupancy and p10. {.table}

With 21 sites the posterior of `sd_det_site` near 0 follows its prior.
The false signal ranges from almost none under the narrowest prior to
above 40% under the widest. With few sites, report `sd_det_site` under
more than one prior, through `fp_priors(det_site_sd = )`. A signal that
holds only under wide priors is not evidence of variation in detection.

## What the model can and cannot say

- A large `sd_det_site` with many sites is good evidence of variation in
  detection beyond abundance.
- Evidence for the absence of variation in detection needs a large
  design. Where \sigma_r = 0, the Savage-Dickey Bayes factor for
  \sigma_r = 0 exceeded 3 in 0%-3% of the fits at 25 sites and 10
  occasions, and in 24%-77% at 100 sites and 20 occasions. Variation in
  abundance lowers these shares. A small `sd_det_site` from a small
  design is not evidence that detection does not vary.
- `sd_abund_site` is weakly identified and leans on its prior when the
  variation in abundance is small.
- Abundance in this model is a latent quantity that links detection to
  occupancy. Its scale depends on the assumption that individuals are
  detected independently with the same r. It should not be read as a
  count of animals without independent support.

## References

Chambert, T., Miller, D. A. W. and Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332-339. <https://doi.org/10.1890/14-1507.1>

Royle, J. A. and Nichols, J. D. (2003). Estimating abundance from
repeated presence-absence data or point counts. *Ecology*, 84, 777-790.
<https://doi.org/10.1890/0012-9658(2003)084%5B0777:EAFRPA%5D2.0.CO;2>

Talts, S., Betancourt, M., Simpson, D., Vehtari, A. and Gelman, A.
(2018). Validating Bayesian inference algorithms with simulation-based
calibration. arXiv:1804.06788.
