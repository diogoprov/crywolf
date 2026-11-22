skip_if_no_cmdstan <- function() {
  testthat::skip_on_cran()
  if (!crywolf:::.has_cmdstan()) testthat::skip("CmdStan 2.40 or later not installed")
}

# the long statistical tests (recovery over several seeds, noise, seeds) run
# only when CRYWOLF_EXTENDED_TESTS is "true"; see tests/README.md
#' @srrstats {G5.10} Extended tests run only with CRYWOLF_EXTENDED_TESTS=true (weekly in .github/workflows/extended-tests.yaml); tests/README.md describes them (G5.12).
#' @srrstats {G5.12} tests/README.md gives the requirements, run times and how to run the extended tests.
skip_if_not_extended <- function() {
  skip_if_no_cmdstan()
  if (!identical(tolower(Sys.getenv("CRYWOLF_EXTENDED_TESTS")), "true")) {
    testthat::skip("extended tests: set CRYWOLF_EXTENDED_TESTS=true")
  }
}
