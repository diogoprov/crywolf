# Which changes does the posterior layer of the regression test flag? Refits
# two designs of tests/testthat/test-regression.R with altered priors and
# applies the test's tolerance. Result of 22/10/2026: a slope prior sd of 0.3
# (instead of 2.5) moved occ_x by 0.20 (tolerance 0.04) and was flagged; a
# recorder prior sd of 0.3 (instead of 2) moved sd_recorder by -0.43
# (tolerance 0.19) and was flagged; a gap prior mean of 3 (instead of 2)
# moved fp_gap by 0.02 (tolerance 0.05) and was not, because the data
# settle the gap in that design.
# Run from the package root: Rscript dev/check-regression-sensitivity.R
devtools::load_all(quiet = TRUE)
source(file.path("tests", "testthat", "helper-regression.R"))
refs <- readRDS(file.path("tests", "testthat", "fixtures", "regression-refs.rds"))
designs <- .regression_designs()
check <- function(name, priors, label) {
  design <- designs[[name]]
  design$priors <- priors
  fit <- .regression_fit(design)
  now <- as.data.frame(posterior::summarise_draws(fit$draws, "mean", "mcse_mean"))
  both <- merge(refs[[name]]$summary, now, by = "variable", suffixes = c("_ref", "_now"))
  tol <- 6 * sqrt(2) * pmax(both$mcse_mean_ref, both$mcse_mean_now) + 1e-8
  shift <- both$mean_now - both$mean_ref
  i <- which.max(abs(shift) / tol)
  cat(sprintf("%-16s %-24s flagged: %-12s largest shift %s = %.3f (tolerance %.3f)\n",
              name, label, paste(both$variable[abs(shift) > tol], collapse = ","),
              both$variable[i], shift[i], tol[i]))
}
check("chambert", fp_priors(gap_mean = 3), "gap_mean 2 -> 3")
check("chambert", fp_priors(slope_sd = 0.3), "slope_sd 2.5 -> 0.3")
check("miller_recorder", fp_priors(recorder_sd = 0.3), "recorder_sd 2 -> 0.3")
