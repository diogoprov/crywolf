# Simulation-based calibration of the site random effect on true detection
# (fp_occu(random_det = TRUE)), with the design of dev/precompute-sbc.R
# (80 sites, 6 occasions, 200 replicates). Fits are about three times
# slower than those of the constant-rate model.
# Run from the package root after installing the package:
#   WORKERS=14 Rscript dev/precompute-sbc-random-det.R
library(crywolf)
fp_compile()
future::plan(future::multisession, workers = as.integer(Sys.getenv("WORKERS", "8")))
t0 <- Sys.time()
s <- simulate_fp(nsite = 80, nocc = 6, seed = 1)
sbc <- fp_sbc(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE, random_det = TRUE,
              reps = 200, parallel = TRUE,
              fit_args = list(chains = 2, parallel_chains = 1, iter_warmup = 400,
                              iter_sampling = 500, threads_per_chain = 1))
attr(sbc, "minutes") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(crywolf:::.stamp(sbc), file.path("vignettes", "sbc-random-det-results.rds"),
        compress = "xz")
