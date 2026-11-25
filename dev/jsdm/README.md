# JSDM prototype: latent factors on false-positive occupancy

Feasibility prototype for a joint species distribution model (JSDM) in
crywolf (dev/roadmap-v2.md, JSDM assessment of 25/11/2026). It lives outside
the package and touches no exported code.

## Files

- `jsdm_fp.stan`: multispecies false-positive occupancy with K latent factors
  on the occupancy logit, `sigma_L * Lambda z_i`. Lambda is J x K with zeros
  above the diagonal (GLLVM parameterisation, as in jsdmstan). Data are
  aggregated by site x species: automated detections (`y` of `n`) and verified
  detections (`v` of `nv`, no false positives). `logit p10 = logit p11 - gap`
  with `gap > 0`, so the labels cannot switch. Options: `diag_pos` (positive
  diagonal of Lambda, or free diagonal with the factor signs fixed after
  sampling, draw by draw) and `det_re` (site effect on the detection logit,
  shared by all species).
- `simulate_jsdm_fp.R`: matching simulator, with an optional shared site
  effect on detection (`sigma_d`).
- `pilot.R`: pilot below. `pilot.rds`: its results.

Priors: community means and SDs for intercepts, slopes, detection and log gap;
`L ~ N(0, 1)`, `sigma_L ~ N(0, 1)`, `Z ~ N(0, 1)`, `sigma_d ~ N(0, 1)`.

## Pilot (25/11/2026)

Design: 100 sites, 10 species, 10 automated occasions. Verified data (3
occasions) were available at about 30% of the sites. Each fit used 4 chains
of 1000 warm-up and 1000 sampling iterations, with `adapt_delta = 0.95`.
There were 3 datasets per scenario, so the results below are a pilot, not a
calibration.

| Scenario | Truth | Fit |
|---|---|---|
| A1 | K = 2 | K = 2, free diagonal, sign fixed after sampling |
| A2 | same data as A1 | K = 2, positive diagonal |
| B1 | K = 0, detection correlated among species (shared site effect, sd 1) | K = 1, no detection effect |
| B2 | same data as B1 | K = 1, with the shared detection effect |
| C | K = 0, independent detection | K = 1 (null) |

Each fit took 1.7–2.6 min with sequential chains. No fit had divergent
transitions or hit the maximum tree depth.

### Identification and mixing (A1 vs A2)

- **A1**: ΛΛ′ mixed in all three datasets (largest R-hat 1.00–1.02, smallest
  bulk ESS 308–1009). The sign-fixed Λ also mixed (R-hat ≤ 1.02).
- The largest R-hat over all parameters was 1.17–1.53. It presumably comes
  from quantities that change sign with the factors (raw loadings and Z). This
  was not checked parameter by parameter.
- **A2**: in one dataset ΛΛ′ did not mix (R-hat 1.52, ESS 7), and in another
  its ESS was 39. The positive diagonal splits the posterior into modes that
  the chains do not cross.
- The free diagonal with signs fixed after sampling is the parameterisation
  to keep.

### Recovery of ΛΛ′ (A1)

- The posterior mean of the off-diagonal entries correlated with the truth at
  0.94, 0.76 and 0.77.
- Magnitude was not recovered. Mean |ΛΛ′| (estimated vs true) was 1.06 vs
  0.69, 1.36 vs 1.04, and 0.64 vs 1.20. The errors went in both directions,
  and the RMSE was 0.67–1.14, the same order as the entries themselves.
- With 100 sites and 10 species, the model recovers the pattern of residual
  co-occurrence but not its strength.

### Correlated detection (B1, B2, C)

- **C (null)**: mean |ΛΛ′| was 0.03–0.06.
- **B1**: in two of three datasets the factor absorbed the shared detection
  effect as occupancy correlation, with mean |ΛΛ′| of 1.05 and 0.47. This
  confirms risk 2 in the roadmap.
- **B2**: with the shared detection effect in the model, mean |ΛΛ′| fell to
  0.03–0.15, and `sigma_d` was estimated at 0.86–0.92 (true 1).
- A factor model without a site effect on detection should not be offered.
  With it, the spurious signal disappeared in this design.

## Decision

Feasible, with two conditions for a package version:

1. Free diagonal and signs fixed after sampling (A1), not a positive
   diagonal.
2. A site effect on detection, shared by species, included by default
   whenever `n_factors > 0`.

Before an implementation in the package:

- A grid over sites (100, 200, 400) and species (10, 20) to establish when the
  strength of ΛΛ′, and not only its pattern, is recovered.
- An SBC run for the factor model.

Both need AWS (an estimate goes to the user before launch).

## Run

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  REPS=3 WORKERS=10 Rscript dev/jsdm/pilot.R
```
