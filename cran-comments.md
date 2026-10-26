## Submission

This is a new submission.

crywolf fits Bayesian false-positive occupancy models written in Stan. The
model is compiled through 'instantiate' when the package is installed from
source and CmdStan is available; without CmdStan (as on CRAN) the package
installs, the Stan program is copied, and the functions that fit models
stop with a message on how to install CmdStan. 'cmdstanr' is in Suggests
and is available from the repository in `Additional_repositories`
(https://stan-dev.r-universe.dev), as for other packages built on
'instantiate'. 'sheepdog' (Suggests) is available from
https://diogoprov.r-universe.dev.

Examples that fit models are in `\dontrun{}`, tests that need CmdStan are
skipped, and vignettes evaluate model fits only when CmdStan is found (the
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
