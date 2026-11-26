# Exporting occupancy draws to external JSDMs

Instead of a JSDM inside crywolf (archived; see dev/roadmap-v2.md), the
posterior draws of the occupancy states of an `fp_occu_ms()` fit can be
passed to an external JSDM (gllvm, Hmsc) by multiple imputation. Each draw
is a binary site × species matrix, so both packages take it as it is. The
package function is `fp_state_draws()` (experimental, from 26/11/2026). It
draws the states as `z_draws.R` does.

## Files

- `z_draws.R`: the state draws, as an array of draws × sites × species.
  Presence is fixed where a verified detection proves it. Same draws as
  `fp_state_draws()`; it uses internal helpers, so it needs
  `devtools::load_all()`.
- `simulate_jsdm_fp.R`: copy of the simulator of the JSDM prototype
  (branch `jsdm-proto`, `dev/jsdm/`).
- `pilot.R`, `pilot.rds`: feasibility check on one dataset (25/11/2026).
  `z_draws.rds`: 200 state draws from its fit, kept for reruns.
- `grid.R`, `run-grid.sh`, `grid/`: the grid of 30 datasets with gllvm
  (26/11/2026), one result file per dataset. `grid_summary.R`,
  `grid_summary.rds`: its summary.
- `hmsc_setup.R`, `hmsc_setup.rds`: the Hmsc setup, on the pilot dataset.

**Species order.** `fp_occu_ms()` sorts the species as text (sp1, sp10,
sp2, …). Until 26/11/2026, the Hmsc part of `pilot.R` compared the draws
with the truth in that order. Its Hmsc result for the draws (0.13) was
therefore wrong. The scripts now reorder the species by name; the gllvm
results of the pilot did not change, because gllvm placed the species by
name.

## Pilot (one dataset)

100 sites, 10 species, 10 automated occasions, two latent factors (mean
|ΛΛ′| 0.69), verified data at 30% of the sites.
`fp_occu_ms(occ = ~ x, det = ~ 1, fp = ~ 1)` with 4 chains (13 min), 20
state draws. `gllvm()` with two latent variables, logit link and `x` as
covariate, fitted to the true states (oracle), the naive states (at least
one automated detection) and each draw; ΛΛ′ averaged over the draws.

| Input to gllvm | Correlation with true ΛΛ′ | Mean \|ΛΛ′\| |
|---|---|---|
| True states | 0.96 | 0.53 |
| Naive states | 0.82 | 0.34 |
| Mean over 20 draws | 0.91 | 0.28 |
| Single draw (median) | 0.82 | 0.31 |

## Grid with gllvm

`grid.R` crossed 50, 100 and 200 sites with two strengths of the latent
factors (loadings with SD 0.5 or 1; mean true |ΛΛ′| about 0.23 and 0.93),
with 5 replicates each: 30 datasets of 10 species and 10 occasions, verified
data at 30% of the sites. The fits and gllvm were as in the pilot, with 20
state draws per dataset. It ran on the Mac in about 3 h (2 fits at a time,
4 chains each), with another heavy session sharing the machine.

All 30 `fp_occu_ms()` fits converged: no divergences, largest R-hat 1.009,
smallest bulk ESS 676. gllvm fitted every draw. Means over the 5 replicates:

| Sites | Loadings SD | Correlation with true ΛΛ′: oracle / naive / draws | Mean \|ΛΛ′\|: oracle / naive / draws |
|---|---|---|---|
| 50 | 0.5 | 0.31 / 0.07 / 0.29 | 0.35 / 0.34 / 0.26 |
| 50 | 1 | 0.68 / 0.45 / 0.66 | 0.66 / 0.53 / 0.45 |
| 100 | 0.5 | 0.62 / 0.32 / 0.58 | 0.29 / 0.18 / 0.21 |
| 100 | 1 | 0.86 / 0.75 / 0.83 | 0.70 / 0.35 / 0.51 |
| 200 | 0.5 | 0.58 / 0.34 / 0.51 | 0.15 / 0.14 / 0.11 |
| 200 | 1 | 0.89 / 0.71 / 0.88 | 0.64 / 0.29 / 0.40 |

- **Pattern.** The mean over draws came within 0.01 to 0.07 of the oracle
  in every cell, and above the naive states by 0.08 to 0.26. It beat the naive
  states in 28 of the 30 datasets. Averaging over draws was better than a
  single draw in every cell.
- **Strength.** The associations were weaker than the oracle's in every
  cell. With strong factors, the draws gave 62 to 73% of the oracle's mean
  |ΛΛ′|. Against the naive states the comparison went both ways: the draws
  gave more with 100 sites (both strengths) and with 200 sites and strong
  factors, less in the other three cells.
- **Weak factors and few sites.** With 50 sites and weak factors, even the
  oracle recovered little (0.31): the data hold little about the
  associations, whatever the input.
- Our interpretation of the attenuation, as in the pilot: `fp_occu_ms()`
  draws the states of each species independently. Where the data leave a
  state uncertain, the draws carry no co-occurrence, and the analysis model
  sees weaker associations. It is the known bias of multiple imputation with
  an imputation model simpler than the analysis model.

## Hmsc setup

On the true states of the pilot dataset, `hmsc_setup.R` compared the scale of
the comparison and two settings of Hmsc (probit, a site-level random effect,
2 chains). With the residual covariance `Omega` against the true ΛΛ′, Hmsc
recovered the pattern at 0.85 with its default prior on the number of
factors, 0.86 with two factors fixed, and 0.85 with chains twice as long.
The pilot compared the residual correlation with the correlation of ΛΛ′,
which is rank-deficient; that comparison is less adequate. The default setup
is enough.

| Input to Hmsc (default setup) | Correlation with true ΛΛ′ | Mean \|Ω\| |
|---|---|---|
| True states | 0.86 | 0.14 |
| Naive states | 0.77 | 0.04 |
| Mean over 20 draws | 0.83 | 0.04 |

The result matches gllvm: the draws keep the pattern close to the oracle,
better than naive detections, and the associations are weaker. Hmsc was not
run on the grid.

## Decision

`fp_state_draws()` enters the package as experimental. Its documentation
shows the loops for gllvm and Hmsc as code that is not run, so neither
package is a dependency. It states the two limits: the strength of the
associations is underestimated, and confusion between species is not
corrected (pilot E0–E2 of the JSDM prototype, `dev/jsdm/` on branch
`jsdm-proto`).

## Run

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  NDRAWS=20 Rscript dev/jsdm-export/pilot.R
sh dev/jsdm-export/run-grid.sh          # about 3 h; resumes where it stopped
Rscript dev/jsdm-export/grid_summary.R
Rscript dev/jsdm-export/hmsc_setup.R
```
