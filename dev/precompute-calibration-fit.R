# Pre-computes the example of fp_calibrate_fit() in vignettes/calibration.Rmd
# ("Calibrating a fitted model"): a model whose false-positive rate depends
# on a gradient w, fitted to one simulated dataset with structured confusion,
# and data simulated from that fit refitted with a constant rate and with
# fp = ~ w (2 models x 50 replicates, about 15 minutes on 12 workers).
# Run from the package root after installing the package, and rerun after
# any change to the Stan model or to fp_occu():
#   WORKERS=12 Rscript dev/precompute-calibration-fit.R
library(crywolf)
fp_compile()
t0 <- Sys.time()
sim <- simulate_fp(nsite = 300, nocc = 20, p10 = 0.05, frac_verified = 0.3,
                   conf_effect = 1.5, conf_cor = 0.7, seed = 1)
m_w <- fp_occu(sim$data, occ = ~ x, fp = ~ w, random_recorder = FALSE,
               threads_per_chain = 1, seed = 1, refresh = 0, show_messages = FALSE)
future::plan(future::multisession, workers = as.integer(Sys.getenv("WORKERS", "8")))
cal_fit <- fp_calibrate_fit(m_w, models = list(constant = ~ 1, fp_w = ~ w),
                            reps = 50, parallel = TRUE)
out <- list(fit = summary(m_w)$table, calibration = cal_fit,
            minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
saveRDS(crywolf:::.stamp(out), file.path("vignettes", "calibration-fit-results.rds"),
        compress = "xz")
print(summary(cal_fit))
