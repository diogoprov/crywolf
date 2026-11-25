# Classifier confusion in the package models (README.md). Data from
# simulate_fp_ms() with 6 species and 100 sites; species 1, 3 and 5 receive
# automated detections from the presence of species 2, 4 and 6, with
# per-occasion probability `conf`. Verified occasions have no confusion.
#   E0  conf 0 (control)
#   E1  conf 0.15
#   E2  conf 0.30
# Each dataset is fitted with fp_occu_ms() (all species) and with fp_occu()
# (species 1 alone). Occupancy bias: posterior mean of the species occupancy
# intercept minus the truth.
#   REPS=3 WORKERS=9 Rscript dev/jsdm/pilot_confusion_pkg.R
reps <- as.integer(Sys.getenv("REPS", "3"))
workers <- as.integer(Sys.getenv("WORKERS", "9"))
out_file <- Sys.getenv("OUT", "dev/jsdm/pilot_confusion_pkg.rds")

scen <- data.frame(id = c("E0", "E1", "E2"), conf = c(0, 0.15, 0.3))
tasks <- expand.grid(model = c("ms", "single"), rep = seq_len(reps),
                     s = seq_len(nrow(scen)), stringsAsFactors = FALSE)
confused <- c("sp1", "sp3", "sp5")
source_sp <- c("sp2", "sp4", "sp6")

run <- function(model, rep, s) {
  suppressMessages(devtools::load_all(".", compile = FALSE, quiet = TRUE))
  sc <- scen[s, ]
  seed <- 5000 + 100 * s + rep
  sim <- simulate_fp_ms(nspecies = 6, nsite = 100, seed = seed)
  set.seed(seed + 1)
  dat <- sim$data
  z <- sim$truth$z
  for (k in seq_along(confused)) {
    rows <- dat$auto$species == confused[k]
    hit <- z[dat$auto$site[rows], source_sp[k]] == 1
    dat$auto$y[rows] <- pmax(dat$auto$y[rows],
                             stats::rbinom(sum(rows), 1, sc$conf * hit))
  }
  t0 <- Sys.time()
  if (model == "ms") {
    fit <- fp_occu_ms(dat, occ = ~ x, det = ~ bio, fp = ~ 1,
                      random_recorder = FALSE, chains = 4, parallel_chains = 1,
                      threads_per_chain = 1, seed = seed, refresh = 0)
    pars <- paste0("occ_(Intercept)[", paste0("sp", 1:6), "]")
    est <- colMeans(as.data.frame(as_draws_df(fit))[, pars])
    bias <- est - sim$truth$pars[pars]
    b_conf <- mean(bias[1:6 %in% c(1, 3, 5)])
    b_other <- mean(bias[1:6 %in% c(2, 4, 6)])
  } else {
    a <- dat$auto[dat$auto$species == "sp1", setdiff(names(dat$auto), "species")]
    v <- dat$verified[dat$verified$species == "sp1",
                      setdiff(names(dat$verified), "species")]
    d1 <- fp_data(auto = a, sites = dat$sites, verified = v)
    fit <- fp_occu(d1, occ = ~ x, det = ~ bio, fp = ~ 1,
                   random_recorder = FALSE, chains = 4, parallel_chains = 1,
                   threads_per_chain = 1, seed = seed, refresh = 0)
    est <- mean(as.data.frame(as_draws_df(fit))[["occ_(Intercept)"]])
    b_conf <- est - sim$truth$pars[["occ_(Intercept)[sp1]"]]
    b_other <- NA
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  dg <- fit$fit$diagnostic_summary(quiet = TRUE)
  data.frame(scenario = sc$id, model = model, rep = rep, seed = seed,
             minutes = secs / 60, divergent = sum(dg$num_divergent),
             rhat_max = max(fit$fit$summary()$rhat, na.rm = TRUE),
             occ_bias_conf = b_conf, occ_bias_other = b_other)
}

# compile once here: concurrent load_all() calls in the workers can block
suppressMessages(devtools::load_all(".", quiet = TRUE))
fp_compile()
.compile_model("fp_occu_ms")
future::plan(future::multisession, workers = workers)
res <- do.call(rbind, furrr::future_pmap(tasks, run,
  .options = furrr::furrr_options(seed = NULL, scheduling = Inf)))
saveRDS(res, out_file)
print(res, digits = 3, row.names = FALSE)
