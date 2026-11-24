# =============================================================================
# dev/rn-fp/acad-real/run_grid.R
# FrogNet request: realistic ACAD design fitted with fp_occu_rn(), and the
# sensitivity of sigma_r to its prior.
#
#   design 21 sites, 11-19 automated occasions per site (one fixed draw, seed
#          500000), every occasion confirmed (frac_confirmed = 1), no verified
#          data; r = 0.27, b1 = 0.5 on a standardized covariate x,
#          sigma_lam = 0.7.
#   cells  psi {0.97, 0.8, 0.3} x p10 {0.05, 0.15} x sigma_r {0, 1}: 12 cells,
#          REPS replicates each (default 100). Seed 600000 + 1000 x cell + rep.
#   fits   each dataset under four priors on sigma_r: half-Normal(0, 0.5),
#          half-Normal(0, 1), half-Normal(0, 2) and half-Student-t(3, 0, 1)
#          (fp_occu_rn_t3.stan, a copy of the package program with only that
#          line changed). occ = ~ x, random_det, random_abund, n_max 40 for
#          psi 0.97 and 30 otherwise, 4 chains of 1000 + 1000 run one after
#          the other, adapt_delta 0.95 (the default).
#
# Rerun: RERUN = an .rds of (cell, rep, prior) refits only those fits, with
# the same data and seed, at ADAPT_DELTA (e.g. 0.99), one job file per fit
# (cellXX-repYYYY-<prior>.rds) in JOBS.
#
# Each task is saved in JOBS (default dev/rn-fp/acad-real/jobs) as soon as it
# ends and skipped on a rerun. Run from the package root with crywolf installed:
#   WORKERS=64 Rscript dev/rn-fp/acad-real/run_grid.R
# =============================================================================
# LOAD_ALL = 1 (local test, with WORKERS = 1): the checkout instead of the installed package
if (Sys.getenv("LOAD_ALL") == "1") devtools::load_all(quiet = TRUE) else library(crywolf)
reps <- as.integer(Sys.getenv("REPS", "100"))
workers <- as.integer(Sys.getenv("WORKERS", "8"))
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.95"))
rerun <- Sys.getenv("RERUN", "")
job_dir <- Sys.getenv("JOBS", "dev/rn-fp/acad-real/jobs")
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(500000)
nocc <- sample(11:19, 21, replace = TRUE)
cells <- expand.grid(sigma_r = c(0, 1), p10 = c(0.05, 0.15), psi = c(0.97, 0.8, 0.3))
cells$n_max <- ifelse(cells$psi > 0.9, 40L, 30L)
priors <- c("hn0.5", "hn1", "hn2", "t3")
tasks <- expand.grid(rep = seq_len(reps), cell = seq_len(nrow(cells)))
tasks$seed <- 600000 + 1000 * tasks$cell + tasks$rep
tasks$prior <- "all"
if (nzchar(rerun)) {
  tasks <- readRDS(rerun)[, c("rep", "cell", "prior")]
  tasks$seed <- 600000 + 1000 * tasks$cell + tasks$rep
}
tasks$file <- file.path(job_dir, ifelse(tasks$prior == "all",
  sprintf("cell%02d-rep%04d.rds", tasks$cell, tasks$rep),
  sprintf("cell%02d-rep%04d-%s.rds", tasks$cell, tasks$rep, tasks$prior)))
tasks <- tasks[!file.exists(tasks$file), ]
if (nzchar(Sys.getenv("MAX_TASKS"))) tasks <- head(tasks, as.integer(Sys.getenv("MAX_TASKS")))
message(nrow(tasks), " tasks to run on ", workers, " workers")

t3_file <- "dev/rn-fp/acad-real/fp_occu_rn_t3.stan"
# compile both programs once, before the workers
invisible(crywolf:::.compile_model("fp_occu_rn"))
invisible(cmdstanr::cmdstan_model(t3_file, cpp_options = list(stan_threads = TRUE)))
pars <- c("sigma_det[1]", "sigma_lam[1]", "alpha[1]", "gamma0", "beta[1]",
          "beta[2]", "p10_ref", "r_ref", "psi_bar", "p_above_nmax")

fit_one <- function(s, cl, prior, seed) {
  sd <- c(hn0.5 = 0.5, hn1 = 1, hn2 = 2, t3 = 1)[[prior]]
  if (prior == "t3") {
    t3 <- cmdstanr::cmdstan_model(t3_file, cpp_options = list(stan_threads = TRUE))
    orig <- crywolf:::.compile_model
    assignInNamespace(".compile_model", function(...) t3, "crywolf")
    on.exit(assignInNamespace(".compile_model", orig, "crywolf"))
  }
  out_dir <- tempfile("rn"); dir.create(out_dir)
  t0 <- Sys.time()
  m <- suppressWarnings(fp_occu_rn(
    s$data, occ = ~ x, random_det = TRUE, random_abund = TRUE,
    n_max = cl$n_max, priors = fp_priors(det_site_sd = sd), chains = 4,
    parallel_chains = 1, threads_per_chain = 1, seed = seed, refresh = 0,
    show_messages = FALSE, output_dir = out_dir, adapt_delta = adapt_delta))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  d <- m$fit$diagnostic_summary(quiet = TRUE)
  sm <- as.data.frame(m$fit$summary(pars, "mean", "median",
                                    ~quantile(.x, c(0.05, 0.95)), "rhat",
                                    "ess_bulk", "ess_tail"))
  rhat_max <- max(m$fit$summary()$rhat, na.rm = TRUE)
  sig <- as.vector(m$fit$draws("sigma_det[1]"))
  unlink(out_dir, recursive = TRUE)
  list(prior = prior, summary = sm, secs = secs, rhat_max = rhat_max,
       divergent = sum(d$num_divergent), treedepth = sum(d$num_max_treedepth),
       ebfmi = d$ebfmi, sigma_r_draws = sig)
}

run_task <- function(cell, rep, seed, file, prior) {
  cl <- cells[cell, ]
  s <- simulate_fp_rn(nsite = 21, nocc = nocc, psi = cl$psi, b1 = 0.5, r = 0.27,
                      p10 = cl$p10, sigma_r = cl$sigma_r, sigma_lam = 0.7,
                      frac_verified = 0, frac_confirmed = 1, seed = seed)
  use <- if (prior == "all") priors else prior
  fits <- lapply(use, function(p) fit_one(s, cl, p, seed))
  names(fits) <- use
  res <- list(cell = cell, rep = rep, seed = seed, design = cl,
              truth = s$truth[c("b0", "b1", "a0", "r", "p10", "sigma_r",
                                "sigma_lam", "psi")],
              n_occupied = sum(s$truth$z), adapt_delta = adapt_delta, fits = fits)
  tmp <- paste0(file, ".tmp")
  saveRDS(res, tmp)
  file.rename(tmp, file)
  invisible(NULL)
}

if (nrow(tasks) > 0) {
  if (workers > 1) future::plan(future::multisession, workers = workers) else future::plan(future::sequential)
  furrr::future_pwalk(tasks, run_task,
                      .options = furrr::furrr_options(seed = NULL, scheduling = Inf))
}
message("done: ", length(list.files(job_dir, "\\.rds$")), " job files")
