test_that("fp_to_unmarked() undoes fp_data_wide()", {
  skip_if_not_installed("unmarked")
  set.seed(1)
  ya <- matrix(stats::rbinom(40, 1, 0.3), 10, 4)
  yv <- matrix(NA, 10, 2)
  yv[1:3, ] <- stats::rbinom(6, 1, 0.5)
  noise <- matrix(stats::rnorm(40), 10, 4)
  d <- fp_data_wide(ya, yv, site_covs = data.frame(x = stats::rnorm(10)),
                    obs_covs = list(noise = noise))
  umf <- fp_to_unmarked(d, obs_covs = "noise", occasion = "occasion")
  expect_s4_class(umf, "unmarkedFrameOccuFP")
  expect_equal(unname(umf@y), unname(cbind(yv, ya)))
  expect_equal(umf@type, c(2, 4, 0))
  expect_equal(umf@siteCovs$x, d$sites$x)
  # obsCovs run site by site, the verified columns first
  oc <- umf@obsCovs
  expect_equal(as.character(oc$METH[1:6]), c("verified", "verified", rep("auto", 4)))
  expect_equal(oc$noise[3:6], noise[1, ])
  expect_equal(oc$noise[1:2], c(0, 0))
  # an occasion missing in the middle keeps its column
  ya[2, 2] <- NA
  d2 <- fp_data_wide(ya, yv)
  expect_equal(unname(fp_to_unmarked(d2, occasion = "occasion")@y[2, 3:6]), ya[2, ])
  expect_error(fp_to_unmarked(d, obs_covs = "METH"), "data\\$auto|METH")
})

test_that("sites with fewer occasions are padded at the end", {
  skip_if_not_installed("unmarked")
  auto <- data.frame(site = c(1, 1, 1, 2), y = c(1, 0, 1, 0))
  ver <- data.frame(site = 1, y = 1)
  d <- suppressWarnings(fp_data(auto, data.frame(site = 1:2), ver))
  umf <- fp_to_unmarked(d)
  expect_equal(unname(umf@y), rbind(c(1, 1, 0, 1), c(NA, 0, NA, NA)))
})

test_that("unsupported data are refused or flagged", {
  skip_if_not_installed("unmarked")
  s <- suppressWarnings(simulate_fp(nsite = 30, nocc = 3, frac_confirmed = 0.3,
                                    frac_verified = 0, seed = 1))
  expect_error(fp_to_unmarked(s$data), "observation-confirmation")
  s2 <- simulate_fp(nsite = 30, nocc = 3, frac_verified = 0, seed = 1)
  expect_warning(fp_to_unmarked(suppressWarnings(s2$data)), "No verified occasions")
})

test_that("occuFP() fits the converted data as the frame built by hand", {
  skip_if_not_installed("unmarked")
  sim <- simulate_fp(nsite = 300, nocc = 10, frac_verified = 0.3, nver = 3, seed = 2026)
  # the frame of vignette("unmarked")
  to_wide <- function(df, nocc) {
    m <- matrix(NA_integer_, 300, nocc)
    m[cbind(df$site, df$occasion)] <- df$y
    m
  }
  meth <- matrix(rep(c("verified", "auto"), times = c(3, 10)), nrow = 300, ncol = 13,
                 byrow = TRUE)
  by_hand <- suppressWarnings(unmarked::unmarkedFrameOccuFP(
    y = cbind(to_wide(sim$data$verified, 3), to_wide(sim$data$auto, 10)),
    siteCovs = data.frame(x = sim$data$sites$x), obsCovs = list(METH = meth),
    type = c(3, 10, 0)))
  umf <- fp_to_unmarked(sim$data, occasion = "occasion")
  expect_equal(unname(umf@y), unname(by_hand@y))
  fit <- function(u) {
    unmarked::coef(unmarked::occuFP(detformula = ~ METH, FPformula = ~ 1,
                                    stateformula = ~ x, data = u))
  }
  expect_equal(fit(umf), fit(by_hand))
  # a covariate of the automated occasions keeps the verified data in the fit
  sim$data$auto$noise <- stats::rnorm(nrow(sim$data$auto))
  umf2 <- fp_to_unmarked(sim$data, obs_covs = "noise", occasion = "occasion")
  m2 <- unmarked::occuFP(detformula = ~ METH + noise, FPformula = ~ 1,
                         stateformula = ~ x, data = umf2)
  expect_equal(unname(unmarked::coef(m2)["p(METHverified)"]),
               unname(fit(by_hand)["p(METHverified)"]), tolerance = 0.1)
})
