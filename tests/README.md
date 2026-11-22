# Tests of crywolf

Most tests run in seconds and need nothing beyond the package. Tests that fit
models need CmdStan 2.40 or later and are skipped without it and on CRAN
(`skip_if_no_cmdstan()` in `tests/testthat/helper.R`).

Longer statistical tests run only when the environment variable
`CRYWOLF_EXTENDED_TESTS` is `"true"` (`skip_if_not_extended()`):

```r
Sys.setenv(CRYWOLF_EXTENDED_TESTS = "true")
devtools::test()
```

They check parameter recovery over ten simulated datasets, the stability of
the estimates to small noise in the covariates and to the random seed, the
agreement with `unmarked::occuFP()`, that precision improves and run time
grows about linearly with the number of sites, that the scale of a covariate
changes only the scale of its effect, and that replicated detections match
the observed scale (`tests/testthat/test-extended.R`). They take about
half an hour, and run weekly on GitHub Actions
(`.github/workflows/extended-tests.yaml`). The tests that always run with CmdStan include a
recovery test on one dataset and a check that sampling with
`fp_occu(prior_only = TRUE)` recovers the priors.

Model checks that take hours (simulation-based calibration and the coverage of
the calibration functions) are not tests: their results are precomputed by the
scripts in `dev/` and shown in the vignettes.
