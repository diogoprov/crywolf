# Step 4 of the Royle-Nichols model: SBC and simulation grid

Validation of `fp_occu_rn()` by simulation, in two parts.

## 1. Simulation-based calibration (`run_sbc.R`)

Checks that the Stan program and `simulate_fp_rn()` describe the same model.
If they do, the rank of each true value among the posterior draws is uniform
(Talts et al. 2018, arXiv:1804.06788).

- Parameters are drawn from the priors of the fit:
  `fp_priors(intercept_sd = 0.75, slope_sd = 0.3, gap_mean = 2, gap_sd = 0.5,
  det_site_sd = 1, abund_site_sd = 0.3)`.
- These priors are narrower than the package defaults so that abundance
  rarely reaches `n_max = 50`. In a prior check of 4000 datasets, 0.1% had a
  site with P(N > 50) > 0.001. With the first choice (sd 1, 0.5 and 0.5, and
  `n_max = 40`) the share was 7.7%, enough to bias the ranks.
- Design A: 50 sites, 10 occasions, every occasion confirmed, no verified
  data (the ACAD scheme).
- Design B: 50 sites, 10 occasions, 20% of occasions confirmed, 3 verified
  occasions at 30% of the sites.
- Model: `occ = ~ x`, `random_det`, `random_abund`.
- 500 replicates per design. Ranks are computed over 100 draws thinned from
  4000. The ESS is kept so the thinning can be checked.
- Parameters checked: `beta[1]`, `beta[2]`, `alpha[1]`, `gamma0`, `delta0`
  (design B only), `sigma_det`, `sigma_lam`.

## 2. General grid (`run_grid.R`)

Measures bias, interval coverage and convergence over a wider space than the
ACAD grid (`../acad-real/`). The main question is how much annotation is
needed to tell a site effect on detection (sigma_r) from a site effect on
abundance (sigma_lam).

- 72 cells, 100 replicates each:
  - nsite: 25, 50, 100
  - nocc: 10, 20
  - annotation: verified, partial, full (defined below)
  - sigma_r: 0, 1
  - sigma_lam: 0, 0.7
- The annotation schemes:
  - verified: 3 verified occasions at 30% of the sites, nothing confirmed.
  - partial: the same, plus 20% of automated occasions confirmed.
  - full: every occasion confirmed and no verified data.
- Fixed values: psi 0.6, b1 0.5, r 0.3, p10 0.05, rv 0.5.
- Fits use the package default priors and `n_max = 25`. Both site effects
  are always estimated.
- Each replicate keeps:
  - summaries of the main quantities
  - diagnostics
  - 1000 draws of sigma_r and sigma_lam, for the Bayes factor and for
    false-signal rates as in the ACAD grid

## Running

Both scripts follow the AWS recipe in `dev/dyn-sim/aws`. One job file is
saved per replicate, and finished replicates are skipped on a rerun. Local
tests use `LOAD_ALL=1`, `REPS`, `GRID_CELLS` and `MAX_TASKS`.

## Pilot and cost

Local pilot on 24/11/2026, run on an Apple Silicon Mac with 10 fits at a
time. The times are minutes per fit:

| design | min/fit | divergences | max R-hat |
|---|---|---|---|
| grid 25 x 10, verified | 0.6-0.8 | 0-3 | 1.003 |
| grid 50 x 10, full | 1.3-1.4 | 0-3 | 1.009-1.014 |
| grid 100 x 20, partial | 2.8-3.6 | 1-2 | 1.007-1.010 |
| grid 100 x 20, full | 2.7-3.0 | 6-20 | 1.021-1.024 |
| SBC A and B (50 x 10) | 1.5-2.0 | 0-18 | 1.004-1.005 |

Sum over the grid: about 125 min per replicate of all 72 cells, so about
210 core-hours for 100 replicates. The SBC takes about 30 core-hours for
2 x 500 replicates.

Divergences show up in most designs, and R-hat reaches 1.02 with full
annotation at 100 sites. The full run therefore uses adapt_delta 0.99, which
is assumed to cost about 1.5 times more. That gives about 360 core-hours.

## Run

AWS spot, three c7a.16xlarge with 64 workers each, started 25/11/2026 at
02:58 UTC with no interruptions:

| run | end (UTC) |
|---|---|
| SBC | 04:21 |
| grid cells 1-48 | 06:37 |
| grid cells 49-72 | 08:44 |

The three runs took about 10.9 instance-hours, about US$10 at the spot price
of US$0.91-0.95 per hour. Summaries are made with `summarise_sbc.R`
(`sbc.rds`, `sbc-ranks.png`) and `summarise_grid.R` (`grid.rds`).

## Results

### SBC

The ranks are close to uniform: every histogram in `sbc-ranks.png` stays
inside the 99% binomial band. The chi-squared p-values range from 0.07 to 0.90,
except for `beta[2]` in design A. Its p-value is 0.03 and its mean rank is
45.8 against 50, which is 3.2 standard errors low.

- Design B shows the same direction for `beta[2]`, but weaker (1.1 standard
  errors).
- On average the posterior median of `beta[2]` lies 0.02 above the true value.
- One cause was checked and ruled out: the site covariate matches the sites
  in the Stan data.
- With 13 tests, a single p-value of 0.03 can arise by chance.
- Fits with divergences (0.2-0.4%) do not drive the result. R-hat was 1.01 or
  more in 1-2% of fits, and no fit had truncated abundance.

### Grid

**Sigma_r.** The lower 5% bound of sigma_r clears 0.1 at the following
rates. Where sigma_r = 0 this is a false signal; where sigma_r = 1 it is
power.

| sigma_lam | false signal | power |
|---|---|---|
| 0 | 0-6% | 41-54% at 25 x 10; 90-95% at 100 x 20 |
| 0.7 | 5-12% | 46-55% at 25 x 10; 93-94% at 100 x 20 |

- Power grows with sites and occasions.
- Power depends little on the annotation scheme. Verified data at 30% of the
  sites do about as well as confirming every occasion.
- Variation in abundance raises the false signal only moderately.
- When sigma_r = 1, the 90% coverage of sigma_r is 0.84-1.00.

**Bayes factor.** BF01 < 1/3 in about half the fits where sigma_r = 1
(median over designs 0.42-0.55). Where sigma_r = 0, the share of fits with
BF01 > 3 grows with the design: 0-3% at 25 x 10 and 24-77% at 100 x 20. It
is lower with sigma_lam = 0.7 (24-30% at 100 x 20) than without (66-77%).
Small designs therefore cannot confirm that sigma_r is absent.

**Sigma_lam.** It is overestimated when its true value is 0 (median 0.43-0.47,
the half-Normal(0, 1) prior speaking) and slightly underestimated when it is
0.7 (median bias -0.11 to -0.13).

**Other parameters.**
- Coverage of the 90% interval across cells is 0.85-0.97 for p10, 0.81-0.99
  for r and 0.85-0.99 for psi at x = 0.
- Median biases are within 0.06 for r and psi.

**Convergence.** The share of fits with R-hat < 1.01 and no divergences falls
with design size, from 0.77-0.92 at 25 sites to 0.19-0.53 at 100 sites, and is
lowest with full annotation.

- Divergences are rare: at most 17% of fits in a cell have any.
- The failures are mostly a largest R-hat of 1.01-1.03, with a bulk ESS of
  170-1100, in the site-effect standard deviations.
- 1000 + 1000 iterations are short for these designs. Analyses of real data
  with many sites should use longer chains.
