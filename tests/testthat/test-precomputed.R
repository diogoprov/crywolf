test_that("precomputed results warn only when made by another Stan program", {
  x <- list(a = 1)
  expect_match(crywolf:::.stale_note(x, "x.rds"), "may be out of date")
  stamped <- crywolf:::.stamp(x)
  expect_identical(crywolf:::.stale_note(stamped, "x.rds"), "")
  expect_identical(attr(stamped, "crywolf_stan_md5"), crywolf:::.stan_md5())
  attr(stamped, "crywolf_stan_md5") <- "0"
  expect_match(crywolf:::.stale_note(stamped, "x.rds"), "x.rds")
})
