# Royle–Nichols with false positives: step 1 (dev/roadmap-v2.md)

Question: is a site effect on the per-individual detection r recoverable on
top of site variation in abundance (lambda), without full annotation?

- `simulate_rn_fp.R`: N_i ~ Poisson(lambda_i), log lambda_i = b0 + b1 x_i + u_i,
  u_i ~ N(0, sigma_lam); logit r_i = a0 + e_i, e_i ~ N(0, sigma_r);
  p11_i = 1 - (1 - r_i)^N_i; false positives p10 on every occasion; verified
  occasions (by ear, no false positives) at a share of the sites, with
  s_i = 1 - (1 - rv)^N_i. The verified occasions depend on N but not on r:
  they are what can separate sigma_r from sigma_lam.
- `rn_fp.stan`: N marginalized over 0..Nmax; e_i integrated by Gauss–Hermite
  quadrature; u_i sampled (non-centered). Binomials written from log(1 - p),
  because `binomial_lpmf` lost log(1 - p) when p was close to 1 (non-finite
  gradients and most transitions divergent in the first version).
- `pilot.R`: three data sets, one fit each (200 sites, 10 occasions, 30% of
  sites verified with 3 occasions).

## Pilot, 23/11/2026 (one replicate per scenario; not a result)

| scenario | true sigma_r | sigma_r (90%) | true sigma_lam | sigma_lam (90%) | divergences | time |
|---|---|---|---|---|---|---|
| A | 1 | 0.88 (0.56–1.20) | 0 | 0.39 (0.04–0.84) | 4% | 142 s |
| B | 0 | 0.34 (0.05–0.66) | 0.7 | 0.35 (0.04–0.75), R-hat 1.10 | 6% | 142 s |
| C | 1 | 0.65 (0.14–1.13) | 0.7 | 0.91 (0.48–1.25) | 0.5% | 224 s |

Reading: sigma_r was recovered when it was 1 and did not absorb the variation
in lambda when it was 0. sigma_lam was poorly identified on its own (B). The
divergences come from the funnel of u_i; the grid should raise adapt_delta or
integrate u_i too. One replicate per cell says nothing about coverage or
power: that needs the grid.
