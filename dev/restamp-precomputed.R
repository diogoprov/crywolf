# Marks the precomputed results of the vignettes as current for the Stan
# program in use, WITHOUT recomputing them. Only for changes that leave the
# posterior unchanged: the script first runs the posterior regression test
# (tests/testthat/test-regression.R) and stops if it fails. Otherwise use
# dev/precompute-all.R.
# Run from the package root: Rscript dev/restamp-precomputed.R
devtools::load_all(quiet = TRUE)
res <- as.data.frame(devtools::test(filter = "regression", reporter = "summary"))
if (any(res$failed > 0) || any(res$error)) {
  stop("The regression test failed: the posterior changed. Run dev/precompute-all.R.",
       call. = FALSE)
}
if (!nrow(res) || any(res$skipped)) {
  stop("The regression test did not run in full (no CmdStan?).", call. = FALSE)
}
files <- list.files("vignettes", pattern = "\\.rds$", full.names = TRUE)
purrr::walk(files, function(f) {
  saveRDS(crywolf:::.stamp(readRDS(f)), f, compress = "xz")
  message("restamped ", f)
})
