# Contributing to crywolf

Thank you for considering a contribution. crywolf fits Bayesian
false-positive occupancy models, so changes to the model or its defaults
are reviewed for their statistical consequences as well as for their
code.

## Reporting a problem

Open an issue at <https://github.com/diogoprov/crywolf/issues>. A useful
report has:

- a minimal reproducible example, ideally with simulated data from
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  (the [reprex](https://reprex.tidyverse.org) package helps);
- the output of
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) and of
  [`cmdstanr::cmdstan_version()`](https://mc-stan.org/cmdstanr/reference/set_cmdstan_path.html);
- for a fit that behaves oddly, the sampler diagnostics printed by
  `print(fit)` and the call that produced it.

Questions about how to use the models for a given design are welcome as
issues too.

## Proposing a change

1.  Open an issue first for anything beyond a typo or a small bug fix,
    so that we can agree on the approach before you write code.
2.  Fork the repository and create a branch from `main`.
3.  Follow the existing style: tidyverse style, the native pipe `|>`,
    and roxygen2 documentation with Markdown. Internal functions get a
    short roxygen comment with `@noRd`.
4.  Add tests with testthat for every change in behaviour. Tests that
    fit a model need CmdStan and are skipped on CRAN; the slow
    statistical tests run only with `CRYWOLF_EXTENDED_TESTS=true` (see
    `tests/README.md`).
5.  Run `devtools::document()`, `devtools::test()` and
    `devtools::check()` before opening the pull request, and add a
    bullet to the development section of `NEWS.md`.

## Changes to the Stan programs

The Stan programs live in `src/stan/`. A change to a likelihood or a
prior should come with a test that recovers known parameters from
simulated data and, when it changes results, a note in `NEWS.md`. The
interface promises in
[`help("crywolf-api")`](https://provetelab.org/crywolf/reference/crywolf-api.md)
apply: stable arguments, classes and result names change only in a major
version.

## Code of conduct

Please be respectful and constructive in issues and pull requests.
Contributors are expected to follow the [Contributor
Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
