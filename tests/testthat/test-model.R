
test_that("the compile cache is in tempdir() under R CMD check", {
  old <- Sys.getenv("_R_CHECK_PACKAGE_NAME_", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("_R_CHECK_PACKAGE_NAME_") else
    Sys.setenv("_R_CHECK_PACKAGE_NAME_" = old))
  Sys.setenv("_R_CHECK_PACKAGE_NAME_" = "crywolf")
  # normalised: dirname() and tempdir() use different separators on Windows
  norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  expect_identical(norm(dirname(crywolf:::.cache_dir())), norm(tempdir()))
  Sys.unsetenv("_R_CHECK_PACKAGE_NAME_")
  expect_identical(crywolf:::.cache_dir(), tools::R_user_dir("crywolf", which = "cache"))
})
