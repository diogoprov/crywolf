# Performance

This page collects what is known about how long crywolf takes, what
makes it faster, and what was tried without success. All times were
measured on one 18-core laptop (Apple M5 Pro) and will differ on other
computers; the ratios are the useful part. The scripts and raw results
are in `dev/profile/` of the source repository, and each change is
described in `NEWS.md`.

## A large dataset over the versions

The `curiango` data have 97,806 one-minute automated occasions in 1,753
site-nights. The model used for timing has three occupancy covariates,
two detection covariates, a constant false-positive rate, and a recorder
random effect, fitted with four chains of 1000 warm-up and 1000 sampling
iterations and the default sampler settings.

| Version | Change | Wall time |
|----|----|----|
| 1.0.0 | one thread per chain | 70 min |
| 1.1.0 | likelihood summed over sites in parallel (`reduce_sum()`), four threads per chain | 16.5 min |
| 1.3.0 | per-site totals for unoccupied sites, one GLM call per site for occupied ones |  |
| 1.5.0 | (same likelihood as 1.3.0), four threads per chain | 6.6 min |

The 1.3.0 change was measured per gradient evaluation rather than per
fit: on these data a gradient took 8.5 ms in 1.2.0 and 2.6 ms in 1.3.0
with one thread. In every version the posterior is the same; only the
arithmetic is arranged differently, and the log density and its gradient
agree up to floating-point rounding.

## What you can do

**Use more threads per chain.** The likelihood is a sum over sites,
which Stan can split among cores. With `threads_per_chain = "auto"`,
crywolf divides the physical cores among the chains that run in
parallel. On the `curiango` data, 150 iterations of one chain took:

| Threads | 1      | 2      | 4      | 8      |
|---------|--------|--------|--------|--------|
| Time    | 55.1 s | 28.8 s | 16.1 s | 10.3 s |

The default is one thread, because with more than one the partial sums
are added in an order that changes from run to run. A fit with several
threads is therefore not reproduced exactly by its seed, even on the
same computer, although its posterior is the same. Use threads for
exploratory fits and large datasets, and one thread when exact
reproduction matters. A fit whose chains took more than a minute on one
core prints a reminder.

**Run chains in parallel.** `parallel_chains` defaults to the number of
chains, so the four chains of a fit run at the same time when there are
four free cores.

**In calibrations, use fewer chains per fit.**
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
and its relatives run many small fits in parallel workers, with the
chains of each fit one after another. There the cost is proportional to
the number of chains: two chains take half the time of four, and still
gave an effective sample size of about 2000 per parameter in the default
design. Keep one thread per chain in parallel calibrations, since the
workers already use the cores.

**Simplify the detection model when you can.** Occasions of a site with
the same detection covariates are collapsed into one binomial row, so a
model without occasion-level detection covariates evaluates far fewer
terms. With covariates on every occasion, the per-site GLM call of 1.3.0
is what keeps the cost down.

## Small fits

For simulated datasets of the sizes used in calibrations, the fixed cost
of starting CmdStan and reading its output was small next to sampling
(medians of three fits, chains in parallel):

| Design                   | Chains | Wall time | Sampler time | Fixed cost |
|--------------------------|--------|-----------|--------------|------------|
| 100 sites, 6 occasions   | 2      | 6.6 s     | 6.2 s        | 0.4 s      |
| 300 sites, 20 occasions  | 2      | 21.2 s    | 20.8 s       | 0.5 s      |
| 1000 sites, 30 occasions | 2      | 66.1 s    | 65.4 s       | 0.7 s      |
| 100 sites, 6 occasions   | 4      | 7.7 s     | 7.1 s        | 0.6 s      |
| 300 sites, 20 occasions  | 4      | 24.4 s    | 23.4 s       | 1.0 s      |
| 1000 sites, 30 occasions | 4      | 71.3 s    | 70.2 s       | 1.1 s      |

## Tried without success

On the `curiango` data each sampling iteration takes about 100 leapfrog
steps, and these were the attempts to reduce that number or the cost of
each step. None gave a consistent gain, and run-to-run differences of up
to 30% in effective sample size per second made smaller effects
impossible to detect.

- **A dense mass matrix** (`metric = "dense_e"`, which can be passed to
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
  through `...`) cut the wall time from 393 s to 161 s but produced
  divergent transitions and halved the effective sample size of the
  recorder standard deviation. It is not a safe default with a
  hierarchical recorder effect.
- **Sampling the false-positive intercept with the gap.** With a
  recorder effect, the mean detection intercept is the uncertain
  quantity and its posterior correlation with the gap reached 0.99.
  Sampling the false-positive intercept instead was meant to avoid that
  correlation, but the number of leapfrog steps did not fall.
- **A centred recorder effect**, which suits recorders with many
  occasions each, alone or with the previous change: no combination was
  best in both the `curiango` data and simulations with 5 to 30
  recorders.
- **`reduce_sum_static()`**, in the hope that a fixed partition of the
  sites would make draws with several threads reproducible. It did not.
- **Initialization by Pathfinder** was not pursued, because warm-up is
  only about 40% of the time, and **a QR decomposition** of the design
  matrices was not, because no posterior correlation among covariate
  effects exceeded 0.5.
- **Longer warm-up or a higher `adapt_delta`** with `random_fp = TRUE`
  took longer without a clear gain: 2 and 3 of 20 fits had an R-hat
  above 1.01, against 3 of 20 with the defaults.
- Earlier: a larger `grainsize` for `reduce_sum()` (no gain), and moving
  the computation to another language (the time is spent inside
  CmdStan).
