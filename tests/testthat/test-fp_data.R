test_that("fp_data validates its inputs", {
  sites <- data.frame(site = 1:3, x = c(-1, 0, 1))
  auto <- data.frame(site = rep(1:3, each = 2), y = c(0, 1, 0, 0, 1, 1))
  expect_s3_class(fp_data(auto, sites), "fp_data")
  expect_error(fp_data(transform(auto, y = y * 2), sites), "0/1")
  expect_error(fp_data(auto, rbind(sites, sites)), "one row per site")
  expect_error(fp_data(transform(auto, site = site + 1), sites), "missing from")
  expect_error(fp_data(auto[, "site", drop = FALSE], sites), "lacks column")
})

test_that("missing detections are dropped and known_present is set", {
  sites <- data.frame(site = 1:3)
  auto <- data.frame(site = 1:3, y = c(1, NA, 0))
  ver <- data.frame(site = c(1, 1, 3), y = c(0, 1, 0))
  d <- fp_data(auto, sites, verified = ver)
  expect_equal(nrow(d$auto), 2)
  expect_equal(d$known_present, c(1L, 0L, 0L))
})

test_that("fp_data_wide matches the long layout", {
  ya <- matrix(c(1, 0, NA, 1, 0, 0), 3, 2)
  yv <- matrix(c(1, NA, NA, 0, NA, NA), 3, 2)
  d <- fp_data_wide(ya, yv, site_covs = data.frame(x = 1:3),
                    obs_covs = list(bio = matrix(1:6, 3, 2)),
                    recorder = c("a", "b", "a"))
  expect_equal(nrow(d$auto), 5)
  expect_equal(nrow(d$verified), 2)
  expect_true(all(c("bio", "recorder") %in% names(d$auto)))
  expect_equal(d$known_present, c(1L, 0L, 0L))
  expect_error(fp_data_wide(ya, obs_covs = list(bio = matrix(1, 2, 2))),
               "same dimensions")
})

test_that("confirmed occasions are validated and mark sites as occupied", {
  auto <- data.frame(site = c(1, 1, 2, 2, 3), y = c(1, 0, 1, 0, 1),
                     n_true = c(1, 0, 0, NA, NA), n_false = c(2, 0, 1, NA, NA))
  sites <- data.frame(site = 1:3)
  d <- fp_data(auto, sites, confirmed = c("n_true", "n_false"))
  expect_equal(d$known_present, c(1L, 0L, 0L))
  expect_output(print(d), "confirmed occasions:    3 \\(1 with a true")

  bad <- auto; bad$n_false[4] <- 0
  expect_error(fp_data(bad, sites, confirmed = c("n_true", "n_false")), "same occasions")
  bad <- auto; bad$n_true[2] <- 1
  expect_error(fp_data(bad, sites, confirmed = c("n_true", "n_false")), "detection must")
  bad <- auto; bad$n_true[1] <- -1
  expect_error(fp_data(bad, sites, confirmed = c("n_true", "n_false")), "non-negative")
  expect_error(fp_data(auto, sites, confirmed = "n_true"), "two columns")
})
