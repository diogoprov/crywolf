skip_if_no_cmdstan <- function() {
  testthat::skip_on_cran()
  if (!crywolf:::.has_cmdstan()) testthat::skip("CmdStan 2.40 or later not installed")
}
