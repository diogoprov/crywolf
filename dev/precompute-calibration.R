# Pre-computes the two calibration runs shown in vignettes/calibration.Rmd
# (800 Stan fits, about 20 minutes on 8 workers), so that the vignette can
# show their results without running them.
#
# Run from the package root after installing the package, and rerun after
# any change to the Stan model or to fp_occu():
#   Rscript dev/precompute-calibration.R

library(crywolf)
fp_compile()
future::plan(future::multisession, workers = 8)

models <- list(constant = ~ 1, fp_on_x = ~ x)

design <- expand.grid(
  nsite = 300, nocc = 20, p10 = 0.05, frac_verified = 0.3,
  conf_effect = c(0, 0.5, 1.5), conf_cor = 0.7
)
cal_confusion <- fp_calibrate(design, models = models, occ = ~ x, reps = 50,
                              parallel = TRUE)

cal_verification <- fp_verification_curve(
  nsite = 300, nocc = 20, p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7,
  frac_verified = c(0.1, 0.3, 0.5, 0.7, 0.9),
  models = models, reps = 50, parallel = TRUE
)

saveRDS(list(confusion = cal_confusion, verification = cal_verification,
             info = list(date = Sys.Date(),
                         crywolf = as.character(utils::packageVersion("crywolf")),
                         cmdstan = as.character(cmdstanr::cmdstan_version()))),
        file.path("vignettes", "calibration-results.rds"), compress = "xz")
