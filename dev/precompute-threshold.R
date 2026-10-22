# Pre-computes the results of vignettes/threshold.Rmd: the same model fitted
# to the curiango data at four thresholds of the custom classifier. About 30
# minutes on 16 cores. Run from the package root after installing the
# package: Rscript dev/precompute-threshold.R
library(crywolf)
source(file.path("dev", "curiango-common.R"))
t0 <- Sys.time()
sweep <- fp_threshold_sweep(cur_occasions, cur_detections,
                            thresholds = c(0.2, 0.4, 0.6, 0.8),
                            sites = cur_sites, verified = cur_verified,
                            site = "site_night", recorder = "recorder",
                            occ = ~ temp + rh + moon, det = ~ bio + aci, fp = ~ 1,
                            chains = 4, parallel_chains = 4, threads_per_chain = 4,
                            seed = 1, refresh = 0, show_messages = FALSE)
attr(sweep, "minutes") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(sweep, file.path("vignettes", "threshold-results.rds"), compress = "xz")
