# Pre-computes the simulation-based calibration shown in
# vignettes/calibration.Rmd (200 replicates, about 2 minutes on 8 workers).
# Run from the package root after installing the package:
#   Rscript dev/precompute-sbc.R
library(crywolf)
fp_compile()
future::plan(future::multisession, workers = 8)
s <- simulate_fp(nsite = 80, nocc = 6, seed = 1)
sbc <- fp_sbc(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE, reps = 200,
              parallel = TRUE,
              fit_args = list(chains = 2, parallel_chains = 1, iter_warmup = 400,
                              iter_sampling = 500))
saveRDS(crywolf:::.stamp(sbc), file.path("vignettes", "sbc-results.rds"), compress = "xz")
