# Follow-up of vignettes/threshold.Rmd: at a threshold of 0.2 the occupancy
# effect of temperature changes sign with a constant false-positive rate.
# Does a false-positive rate that depends on the conditions of the night
# (fp ~ temp + rh) remove the change? Fits thresholds 0.2 and 0.4 with
# fp ~ temp + rh; compare with threshold-results.rds (fp ~ 1).
# About 30 minutes with 4 chains x 2 threads (8 cores). When other heavy jobs
# share the machine, keep to their core limit and set
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1.
# Install the package first, then run from the package root:
#   Rscript dev/precompute-threshold-fp.R
library(crywolf)
source(file.path("dev", "curiango-common.R"))
t0 <- Sys.time()
sweep <- fp_threshold_sweep(cur_occasions, cur_detections, thresholds = c(0.2, 0.4),
                            sites = cur_sites, verified = cur_verified,
                            site = "site_night", recorder = "recorder",
                            occ = ~ temp + rh + moon, det = ~ bio + aci,
                            fp = ~ temp + rh, chains = 4, parallel_chains = 4,
                            threads_per_chain = 2, seed = 1, refresh = 0,
                            show_messages = FALSE)
attr(sweep, "minutes") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(crywolf:::.stamp(sweep), file.path("vignettes", "threshold-fp-results.rds"),
        compress = "xz")
