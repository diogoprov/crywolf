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
- `pilot_confusion.R`: classifier-confusion pilot below.
  `pilot_confusion.rds`: its results.
- `pilot_confusion_pkg.R`: the same confusion in `fp_occu()` and
  `fp_occu_ms()`. `pilot_confusion_pkg.rds`: its results.

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

### Classifier confusion (D1–D3)

Truth: no factors and no shared detection effect. Species 1, 3 and 5
received automated detections from the presence of species 2, 4 and 6, with
per-occasion probability `conf`. Verified occasions had no confusion. Every
fit used K = 1 with the shared detection effect, the configuration proposed
for the package.

| Scenario | `conf` | Verified sites |
|---|---|---|
| D1 | 0.15 | 30% |
| D2 | 0.30 | 30% |
| D3 | 0.30 | 10% |

Each fit took 2.6–3.3 min. No fit had divergent transitions, and ΛΛ′ mixed
in all of them (R-hat ≤ 1.02, bulk ESS ≥ 199).

- **Spurious association.** The posterior mean of ΛΛ′ for the confused pairs
  was positive in all nine fits. It was 0.09–1.19 in D1, 0.97–4.53 in D2 and
  1.25–3.51 in D3. The null scenario C had mean |ΛΛ′| of 0.03–0.06.
- **Spillover.** The factor also distorted the other pairs, with mean |ΛΛ′|
  of 0.12–1.12.
- **The shared detection effect does not help.** `sigma_d` was estimated at
  0.05–0.19 (true 0), so it did not absorb the confusion.
- **Verified data at 30% of the sites did not prevent the bias.**
- **Occupancy bias.** The intercepts of the confused species were biased
  upwards, by 0.18–0.61 logit in D1, 0.70–1.49 in D2, and 0.00–1.47 in D3.
  For the other species the bias was −0.16 to 0.30.

Confusion makes the false-positive rate of one species depend on the presence
of another. The model assumes a false-positive rate that does not vary with
site, so the factor absorbs the confusion as co-occurrence. The upward
occupancy bias suggests the same assumption matters outside the JSDM, which
the next pilot tests.

### Classifier confusion in the package models (E0–E2)

Data from `simulate_fp_ms()` with 6 species, 100 sites and its defaults
(verified data at 30% of the sites). Species 1, 3 and 5 received automated
detections from the presence of species 2, 4 and 6. At `conf` 0.30, the
detection rate of species 1 at sites where it was absent and species 2 was
present rose from 0.09 to 0.38 per occasion, against 0.61 where species 1
was present. Each dataset was fitted with `fp_occu_ms()` (all species) and
with `fp_occu()` (species 1 alone), using `occ = ~ x, det = ~ bio`.

| Scenario | `conf` | `fp_occu_ms()`, confused species | `fp_occu_ms()`, other species | `fp_occu()`, species 1 |
|---|---|---|---|---|
| E0 | 0 | −0.03 to 0.01 | 0.03 to 0.22 | −0.19 to 0.18 |
| E1 | 0.15 | −0.09 to 0.13 | −0.22 to 0.23 | −0.43 to 0.33 |
| E2 | 0.30 | 0.56 to 0.96 | −0.04 to 0.09 | −0.06 to 1.72 |

The values are the occupancy intercept bias in logit (posterior mean minus
truth), over three datasets per scenario.

No fit had divergent transitions, and the largest R-hat was 1.01.
`fp_occu()` took about 1 min per fit; `fp_occu_ms()` took 33–42 min with
sequential chains.

- **Moderate confusion was absorbed.** At `conf` 0.15, the bias was within
  the range of the control in both models.
- **Strong confusion biased occupancy upwards.** At `conf` 0.30, the bias was
  0.56–0.96 logit in all three multispecies fits. In the single-species
  fits it was about 1.7 logit in two of three datasets.
- **The confusion and the JSDM interact.** The JSDM prototype was already
  biased at `conf` 0.15 (0.18–0.61 logit). Our interpretation, not tested:
  the factor turns the confused detections into a shared signal and
  amplifies the bias.

This is a limitation of the package models, not only of the JSDM. They
assume that the false-positive rate of a species does not depend on which
other species are present. With verified data at 30% of the sites, that
assumption held up under moderate confusion but not when the confused
detections were as frequent as true ones.

### What ΛΛ′ is not

ΛΛ′ is residual co-occurrence, which is not evidence of species interactions
(Blanchet et al. 2020). In a competitive metacommunity model, Barbier et al.
(2025) found the opposite split to this pilot. Individual interactions could
not be inferred from cross-species correlations, but the mean and variance of
the interaction distribution could. The two results are not comparable. The
failure here is statistical (a fitted model equal to the generating one, with
few sites), whereas theirs comes from mapping interactions onto correlations.
Any package documentation of ΛΛ′ should carry this caveat.

## Decision

Feasible in the narrow sense: the model fits, mixes and, without
confusion, separates correlated detection from co-occurrence.

Not ready for the package. Classifier confusion produced a spurious
association that the shared detection effect did not remove, and confusion
between similar calls is common in automated acoustic data. A package
version needs three conditions:

1. Free diagonal and signs fixed after sampling (A1), not a positive
   diagonal.
2. A site effect on detection, shared by species, included by default
   whenever `n_factors > 0`.
3. A confusion term: the false-positive rate of species a increases with the
   presence of species b, for pairs the user declares. The verified data
   would identify it.

Next steps, in order:

- Document in `fp_occu()` and `fp_occu_ms()` the assumption that the
  false-positive rate does not depend on the presence of other species, and
  its consequence under strong confusion (pilot E0–E2).
- Prototype the confusion term and repeat D1–D3.
- Only then, a grid over sites (50, 100, 200, 400) and species (10, 20) and an
  SBC run. Both need AWS (an estimate goes to the user before launch).

## Run

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  REPS=3 WORKERS=10 Rscript dev/jsdm/pilot.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  REPS=3 WORKERS=9 Rscript dev/jsdm/pilot_confusion.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  REPS=3 WORKERS=9 Rscript dev/jsdm/pilot_confusion_pkg.R
```

## References

- Barbier, M., Bunin, G. & Leibold, M. A. (2025). Getting more by asking for
  less: Linking species interactions to species co-distributions in
  metacommunities. *Peer Community Journal*, 5, e1.
  https://doi.org/10.24072/pcjournal.483
- Blanchet, F. G., Cazelles, K. & Gravel, D. (2020). Co-occurrence is not
  evidence of ecological interactions. *Ecology Letters*, 23, 1050–1063.
  https://doi.org/10.1111/ele.13525
