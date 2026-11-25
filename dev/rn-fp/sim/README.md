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
