# Exporting occupancy draws to external JSDMs

Feasibility check, 25/11/2026. Instead of a JSDM inside crywolf (archived;
see dev/roadmap-v2.md), the posterior draws of the occupancy states of an
`fp_occu_ms()` fit could be passed to an external JSDM (gllvm, Hmsc) by
multiple imputation. Each draw is a binary site × species matrix, so
both packages take it as it is.

## Files

- `z_draws.R`: draws the states from an `fp_occu_ms()` fit, as
  `fp_richness()` does internally. The result is an array of draws × sites ×
  species. Presence is fixed where a verified detection proves it. It uses
  internal helpers, so it needs `devtools::load_all()`.
- `simulate_jsdm_fp.R`: copy of the simulator of the JSDM prototype
  (branch `jsdm-proto`, `dev/jsdm/`).
- `pilot.R`: the check below. `pilot.rds`: its results. `z_draws.rds`: 200
  state draws from the fit, kept so that the fit is not needed for reruns.

## Design

There was one dataset from the JSDM simulator, so the numbers below show
feasibility, not performance. It had 100 sites, 10 species and 10
automated occasions, with two latent factors (mean |ΛΛ′| 0.69). Verified
data were available at 30% of the sites.

`fp_occu_ms(occ = ~ x, det = ~ 1, fp = ~ 1)` was fitted with 4 chains, and
20 state draws were taken from it. `gllvm()` used two latent variables and a
logit link, with `x` as covariate. It was fitted to:

- the true states (oracle);
- the naive states (at least one automated detection);
- each of the 20 draws.

ΛΛ′ was then averaged over the draws. `Hmsc` used probit, a site random
level and 2 chains of 1000 samples. It was fitted to the true states and to
10 draws.

## Results

| Input to gllvm | Correlation with true ΛΛ′ | Mean \|ΛΛ′\| |
|---|---|---|
| True states | 0.96 | 0.53 |
| Naive states | 0.82 | 0.34 |
| Mean over 20 draws | 0.91 | 0.28 |
| Single draw (median) | 0.82 | 0.32 |

- **It works and is cheap.** The `fp_occu_ms()` fit took 13 min. Each gllvm
  fit took under a second, and each Hmsc fit about 6 s.
- **The pattern improves on naive detections.** Averaged over draws, the
  correlation with the true ΛΛ′ was 0.91, against 0.82 for naive states.
- **The strength is attenuated.** Mean |ΛΛ′| was 0.28, below both the oracle
  (0.53) and the naive states (0.34). Our interpretation: `fp_occu_ms()`
  draws the states of each species independently. Where the data leave a
  state uncertain, the draws carry no correlation, and they dilute the
  correlation of the analysis model. The imputation model is simpler than
  the analysis model, a known source of bias towards zero in multiple
  imputation.
- **Hmsc was inconclusive.** Even on the true states, its residual
  correlation matched the truth at only 0.55, so the comparison used here is
  not adequate for it. The 10 draws gave 0.13. The workflow runs, but Hmsc
  needs its own setup (number of factors, priors) before any comparison.

## Decision

Exporting state draws is technically simple and gives a better pattern of
residual co-occurrence than naive detections. It underestimates the
strength of the associations. It also inherits the confusion bias of
`fp_occu_ms()` (dev/jsdm/ on branch `jsdm-proto`, pilot E0–E2). A function
in the package would document both limits.

Before that:

- a grid of datasets (sizes of design, strength of ΛΛ′), with gllvm only,
  locally;
- an adequate Hmsc setup.

## Run

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  NDRAWS=20 Rscript dev/jsdm-export/pilot.R
```
