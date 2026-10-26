# Stability of the crywolf interface

What crywolf promises not to change between versions, and which
functions are still experimental.

## What is frozen

From crywolf 2.0.0, for the stable functions listed below:

- **Arguments**: their names, their order and their default values, and
  those of the methods of crywolf's classes. A minor version may add an
  argument after the existing named ones (and before a trailing `...`),
  with a default that keeps the earlier behaviour.

- **Classes**: the crywolf classes their results inherit from
  (`fp_data`, `fp_occu`, `fp_priors`, `fp_sim`, `fp_calibration`,
  `fp_sbc`, `fp_ppc`, `fp_prior_predictive`, `fp_effects`, `fp_sites`,
  `fp_threshold_sweep`, `fp_sheepdog`, `summary.fp_occu`,
  `summary.fp_calibration`, `summary.fp_sbc`), and the class of the
  design warnings of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)
  (`crywolf_design_warning`), which code may catch. Classes of other
  packages (the `CmdStanModel` of
  [`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md),
  the `psis_loo` of
  [`loo()`](https://mc-stan.org/loo/reference/loo.html)) are theirs to
  keep.

- **Results**: the names documented in
  [crywolf-output](https://provetelab.org/crywolf/reference/crywolf-output.md)
  (parameters, derived quantities, columns of summaries and of the other
  results described there) and the components and columns described in
  the Value section of each stable function, for example the tables of
  [`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md),
  [`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md)
  and
  [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md).
  New columns or components may be added. Columns that a reader passes
  on from its input unchanged (the selection-table columns of
  [`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md);
  the `scores` and `sites` tables and the `schema` of
  [`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md))
  follow BirdNET and 'sheepdog' and are not covered.

These change only in a major version, announced in the NEWS of the
previous minor version. Not covered: arguments passed through `...` to
'cmdstanr', the fields of `stan_data` and `design`, raw Stan names, the
text of messages and warnings, and the numerical results themselves,
which may move within Monte Carlo error when the sampler or its settings
change (for example `threads_per_chain`, whose default became `"auto"`
in 2.0.0). A test of the package (run on its continuous integration, not
on CRAN) records the arguments of every stable function and method, so
that a change cannot pass unnoticed.

## Stable functions

Data and reading:
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md),
[`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md),
[`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md),
[`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md),
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md),
[`fp_threshold_sweep()`](https://provetelab.org/crywolf/reference/fp_threshold_sweep.md).

Fitting and checking:
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md),
[`fp_gap_prior()`](https://provetelab.org/crywolf/reference/fp_gap_prior.md),
[`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md),
[`fp_stan_file()`](https://provetelab.org/crywolf/reference/fp_stan_file.md),
[`fp_prior_predictive()`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md),
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md),
[`fp_loo_weights()`](https://provetelab.org/crywolf/reference/fp_loo_weights.md),
the methods in
[fp_occu-methods](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
([`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`loo()`](https://mc-stan.org/loo/reference/loo.html),
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html))
and
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md).

Results and simulation:
[`fp_effects()`](https://provetelab.org/crywolf/reference/fp_effects.md),
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md),
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md),
[`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md),
[`fp_calibrate_fit()`](https://provetelab.org/crywolf/reference/fp_calibrate_fit.md),
[`fp_sbc()`](https://provetelab.org/crywolf/reference/fp_sbc.md).

The [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods of the
classes above are stable too.

## Experimental functions

[`fp_occu_mc()`](https://provetelab.org/crywolf/reference/fp_occu_mc.md),
[`fp_draws_data()`](https://provetelab.org/crywolf/reference/fp_draws_data.md)
and
[`fp_to_unmarked()`](https://provetelab.org/crywolf/reference/fp_to_unmarked.md),
with the `fp_occu_mc` class and its methods, are new in 2.0.0. They
propagate a classifier calibration from 'sheepdog' into the model or
pass data to 'unmarked', and are still being evaluated. Their arguments,
classes and results may change in a minor version; the NEWS will say so.
Each says so in its help.
[`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md),
by contrast, reads a versioned export format (`sheepdog-occupancy`,
schema 1) and is stable.
