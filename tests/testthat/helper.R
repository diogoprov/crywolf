skip_if_no_cmdstan <- function() {
  testthat::skip_on_cran()
  if (!crywolf:::.has_cmdstan()) testthat::skip("CmdStan 2.40 or later not installed")
}

# the long statistical tests (recovery over several seeds, noise, seeds) run
# only when CRYWOLF_EXTENDED_TESTS is "true"; see tests/README.md
skip_if_not_extended <- function() {
  skip_if_no_cmdstan()
  if (!identical(tolower(Sys.getenv("CRYWOLF_EXTENDED_TESTS")), "true")) {
    testthat::skip("extended tests: set CRYWOLF_EXTENDED_TESTS=true")
  }
}
