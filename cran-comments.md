## Resubmission

This is a resubmission. The version moves from 2.1.0 to 2.3.0 because two
releases (2.2.0 and 2.3.0) were made on GitHub in between; NEWS.md lists
their changes, which add a dynamic and a multispecies model. The
Description now names them.

In response to the review of 19 November 2026:

* `\dontrun{}` is no longer used. Examples that fit models are in
  `\donttest{}` and run only when CmdStan is installed
  (`@examplesIf instantiate::stan_cmdstan_exists()`), since fitting needs
  this external software. Examples that read files use a small export
  shipped in `inst/extdata` and run unwrapped.
* Nothing is written to the user's filespace in examples, vignettes or
  tests: they write only to `tempdir()`. When the model has to be compiled
  after installation (binary packages), the executable is kept in
  `tools::R_user_dir("crywolf", "cache")`, as the CRAN policy allows; it
  holds only the Stan program and its executable, which a new version
  replaces. Under `R CMD check` this cache is in `tempdir()`.
* Examples and vignettes use at most two cores: fits in examples set
  `chains = 2, threads_per_chain = 1`, fits in vignettes set
  `parallel_chains = 2, threads_per_chain = 1`, and under `R CMD check`
  every fit is limited to two cores in total.
* Tessa Rhinehart is no longer listed as an author.

crywolf fits Bayesian false-positive occupancy models written in Stan. The
model is compiled through 'instantiate' when the package is installed from
source and CmdStan is available; without CmdStan (as on CRAN) the package
installs, the Stan program is copied, and the functions that fit models
stop with a message on how to install CmdStan. 'cmdstanr' is in Suggests
and is available from the repository in `Additional_repositories`
(https://stan-dev.r-universe.dev), as for other packages built on
'instantiate'. 'sheepdog' (Suggests) is available from
https://diogoprov.r-universe.dev.

Tests that need CmdStan are skipped, and vignettes evaluate model fits only when CmdStan is found (the
long ones show precomputed results). Under `R CMD check` a fit uses at
most two cores.

## Test environments

* GitHub Actions: macOS, Ubuntu and Windows (R release), without CmdStan
  and without 'cmdstanr', with `--as-cran`; macOS and Ubuntu with CmdStan.
* Local: macOS (R 4.6), with CmdStan 2.40.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission; 'cmdstanr' and 'sheepdog' in Suggests are not on CRAN
  (available from the repositories in `Additional_repositories`).
