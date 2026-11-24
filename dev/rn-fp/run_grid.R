# =============================================================================
# dev/rn-fp/run_grid.R
# Step 1 grid of the Royle-Nichols model with false positives (README.md):
# is the site effect on the per-individual detection r recoverable on top of
# site variation in abundance, and how does it depend on the share of
# verified sites?
#
#   cells  sigma_r {0, 0.5, 1} x sigma_lam {0, 0.7} x frac_verified
#          {0.1, 0.3, 0.6}: 18 cells, 20 replicates each (360 fits).
#          200 sites, 10 automated occasions, 3 verified occasions at the
#          verified sites; the other parameters at the defaults of
#          simulate_rn_fp() (b0 = 0, b1 = 0.5, a0 = -1, p10 = 0.05, rv = 0.5).
#   fit    rn_fp.stan with site effects on r and on log lambda, Nmax = 20,
#          9 Gauss-Hermite nodes, 4 chains of 1000 + 1000 iterations run one
#          after the other on one core, adapt_delta = ADAPT_DELTA (default
#          0.99). Seed 300000 + 1000 x cell + replicate, for data and fit.
#
# Each task is saved in dev/rn-fp/jobs/ as soon as it ends and skipped on a
# rerun. Run from the package root:
#   WORKERS=16 Rscript dev/rn-fp/run_grid.R
# REPS sets the replicates per cell (REPS = 1 for a test); ITER the
# iterations; GRID_CELLS (e.g. "1,7") keeps some cells. Summaries: summarise_grid.R.
# =============================================================================

source("dev/rn-fp/simulate_rn_fp.R")

reps <- as.integer(Sys.getenv("REPS", "20"))
workers <- as.integer(Sys.getenv("WORKERS", "16"))
iter <- as.integer(Sys.getenv("ITER", "1000"))
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.99"))
cells_env <- Sys.getenv("GRID_CELLS", "")
# DESIGN "acad" (set by run_grid_acad.R): the design of the FrogNet ACAD data
design <- Sys.getenv("DESIGN", "step1")
job_dir <- Sys.getenv("JOB_DIR", file.path("dev", "rn-fp",
                                           if (design == "acad") "jobs-acad" else "jobs"))
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)

if (design == "acad") {
  # psi at the mean of x through b0 = log(-log(1 - psi)); nver = 17 verifies
  # every occasion of every site, nver = 3 three occasions at every site
  cells <- expand.grid(sigma_r = c(0, 0.5, 1), sigma_lam = c(0, 0.7),
                       psi = c(0.8, 0.3), nver = c(17L, 3L))
  cells$nsite <- 21L; cells$nocc <- 17L; cells$frac_verified <- 1
  cells$b0 <- log(-log(1 - cells$psi)); cells$Nmax <- 30L; seed_base <- 400000
} else {
  cells <- expand.grid(sigma_r = c(0, 0.5, 1), sigma_lam = c(0, 0.7),
                       frac_verified = c(0.1, 0.3, 0.6))
  cells$nsite <- 200L; cells$nocc <- 10L; cells$nver <- 3L; cells$b0 <- 0
  cells$Nmax <- 20L; seed_base <- 300000
}
keep <- if (nzchar(cells_env)) as.integer(strsplit(cells_env, ",")[[1]]) else
  seq_len(nrow(cells))
tasks <- expand.grid(rep = seq_len(reps), cell = keep)
tasks$seed <- seed_base + 1000 * tasks$cell + tasks$rep
tasks$file <- file.path(job_dir, sprintf("cell%02d-rep%04d.rds", tasks$cell, tasks$rep))
tasks <- tasks[!file.exists(tasks$file), ]
message(nrow(tasks), " tasks to run on ", workers, " workers")

# Gauss-Hermite nodes and log weights for a standard normal (Golub-Welsch)
gh <- function(nq) {
  J <- matrix(0, nq, nq)
  b <- sqrt(seq_len(nq - 1) / 2)
  J[cbind(1:(nq - 1), 2:nq)] <- b
  J[cbind(2:nq, 1:(nq - 1))] <- b
  ev <- eigen(J, symmetric = TRUE)
  list(nq = nq, gh_x = sqrt(2) * ev$values, gh_logw = log(ev$vectors[1, ]^2))
}
GH <- gh(9)
stan_file <- "dev/rn-fp/rn_fp.stan"
invisible(cmdstanr::cmdstan_model(stan_file))   # compile once, before the workers
pars <- c("b0", "b1", "a0", "sigma_r_out", "sigma_lam_out", "p10", "rv")

run_task <- function(cell, rep, seed, file) {
  cl <- cells[cell, ]
  s <- simulate_rn_fp(nsite = cl$nsite, nocc = cl$nocc, b0 = cl$b0,
                      sigma_r = cl$sigma_r, sigma_lam = cl$sigma_lam,
                      frac_verified = cl$frac_verified, nver = cl$nver,
                      seed = seed)
  mod <- cmdstanr::cmdstan_model(stan_file)
  out_dir <- tempfile("rn")
  dir.create(out_dir)
  t0 <- Sys.time()
  f <- mod$sample(data = c(s$data, list(Nmax = cl$Nmax, site_effect = 1L,
                                        lam_effect = 1L), GH),
                  chains = 4, parallel_chains = 1, iter_warmup = iter,
                  iter_sampling = iter, adapt_delta = adapt_delta, seed = seed,
                  refresh = 0, show_messages = FALSE, output_dir = out_dir)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  d <- f$diagnostic_summary(quiet = TRUE)
  sm <- as.data.frame(f$summary(pars, "mean", "median",
                                ~quantile(.x, c(0.05, 0.95)), "rhat",
                                "ess_bulk", "ess_tail"))
  sm$truth <- with(s$truth, c(b0, b1, a0, sigma_r, sigma_lam, p10, rv))
  res <- list(cell = cell, rep = rep, seed = seed, design = cl, summary = sm,
              secs = secs, divergent = sum(d$num_divergent),
              treedepth = sum(d$num_max_treedepth), ebfmi = d$ebfmi,
              sigma_r_draws = as.vector(f$draws("sigma_r_out")),
              sigma_lam_draws = as.vector(f$draws("sigma_lam_out")))
  unlink(out_dir, recursive = TRUE)
  tmp <- paste0(file, ".tmp")
  saveRDS(res, tmp)
  file.rename(tmp, file)
  invisible(NULL)
}

if (nrow(tasks) > 0) {
  future::plan(future::multisession, workers = workers)
  furrr::future_pwalk(tasks, run_task,
                      .options = furrr::furrr_options(seed = NULL, scheduling = Inf))
}
message("done: ", length(list.files(job_dir, "\\.rds$")), " job files")
