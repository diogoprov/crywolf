# Validation of the tuple refactor of the Stan program (crywolf 1.7.0): the
# likelihood data travel in tuples instead of ~45 positional arguments.
# Needs CmdStan >= 2.40 (tuples as reduce_sum() arguments do not compile in
# 2.39). 1. Draws identical to the previous program (same seed, 1 thread) in
# six designs. 2. Time on the curiango data along a fixed trajectory (same
# start, step size and metric, no adaptation), 1 and 4 threads.
# Run from the package root: Rscript dev/profile/validate-tuples.R
# Output: dev/profile/validate-tuples.txt
sink(file.path("dev", "profile", "validate-tuples.txt"), split = TRUE)
devtools::load_all(quiet = TRUE)
cat("CmdStan", cmdstanr::cmdstan_version(), "\n")
old_stan <- tempfile(fileext = ".stan")
writeLines(system2("git", c("show", "v1.6.0:inst/stan/fp_occu.stan"), stdout = TRUE), old_stan)
opts <- list(stan_threads = TRUE)
old <- cmdstanr::cmdstan_model(old_stan, cpp_options = opts, quiet = TRUE)
new <- fp_compile()

sampler <- list(chains = 2, parallel_chains = 2, threads_per_chain = 1, iter_warmup = 300,
                iter_sampling = 300, adapt_delta = 0.95, refresh = 0, show_messages = FALSE)
same <- function(label, data, ...) {
  sd <- .build_stan_data(data, ...)$stan_data
  f1 <- do.call(new$sample, c(list(data = sd, seed = 7), sampler))
  f0 <- do.call(old$sample, c(list(data = sd, seed = 7), sampler))
  v <- intersect(c("beta", "alpha", "fp_gap", "delta0", "delta_s", "gamma_fp",
                   "sigma_u", "sigma_fp", "log_lik"), f0$metadata()$stan_variables)
  d <- max(abs(unclass(posterior::as_draws_matrix(f1$draws(v))) -
                 unclass(posterior::as_draws_matrix(f0$draws(v)))))
  cat(sprintf("1. %-44s max |draw difference| = %g\n", label, d))
}
s1 <- simulate_fp(nsite = 150, nocc = 8, a_bio = 0.5, sd_recorder = 0.5, conf_effect = 1,
                  conf_cor = 0.5, seed = 3)
s2 <- simulate_fp(nsite = 120, nocc = 6, frac_confirmed = 0.4, frac_verified = 0.2, seed = 6)
s3 <- simulate_fp(nsite = 120, nocc = 6, seed = 8)
set.seed(9)
d4 <- s2$data
d4$auto$noise <- stats::rnorm(nrow(d4$auto))
d4$verified$listener <- stats::rnorm(nrow(d4$verified))
same("Miller, recorder (GLM path), fp ~ w", s1$data, ~ x, ~ bio, ~ w, TRUE, fp_priors())
same("Chambert, constant rate", s2$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
same("no false positives", s3$data, ~ x, ~ 1, NULL, FALSE, fp_priors())
same("random_fp", s3$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors(), random_fp = TRUE)
same("Chambert, fp ~ noise, ver ~ listener, markov", d4, ~ x, ~ 1, ~ noise, FALSE, fp_priors(),
     ver = ~ listener, markov = "occasion")
same("Miller, recorder, markov, det ~ bio", s1$data, ~ x, ~ bio, ~ 1, TRUE, fp_priors(),
     markov = "occasion")

# 2. curiango, fixed trajectory
source(file.path("dev", "curiango-common.R"))
auto <- cur_occasions
auto$y <- as.integer(cur_detections$confidence >= 0.40)
d <- suppressWarnings(fp_data(auto, cur_sites, cur_verified, site = "site_night",
                              recorder = "recorder"))
sd <- .build_stan_data(d, ~ temp + rh + moon, ~ bio + aci, ~ 1, TRUE, fp_priors())$stan_data
warm <- new$sample(data = sd, seed = 1, chains = 1, threads_per_chain = 4, iter_warmup = 300,
                   iter_sampling = 10, refresh = 0, show_messages = FALSE)
run <- function(mod, threads) {
  f <- mod$sample(data = sd, seed = 1, chains = 1, threads_per_chain = threads,
                  iter_warmup = 0, iter_sampling = 100, adapt_engaged = FALSE,
                  step_size = warm$metadata()$step_size_adaptation[1],
                  inv_metric = unname(warm$inv_metric(matrix = FALSE)[1]),
                  init = warm, refresh = 0, show_messages = FALSE)
  list(sec = f$time()$chains$sampling,
       draws = unclass(posterior::as_draws_matrix(f$draws(c("beta", "alpha", "lp__")))))
}
for (th in c(1, 4)) {
  r0 <- run(old, th); r1 <- run(new, th); r0b <- run(old, th); r1b <- run(new, th)
  cat(sprintf("2. curiango, %d thread(s): previous %.1f / %.1f s, tuples %.1f / %.1f s; same draws (1 thread only): %s\n",
              th, r0$sec, r0b$sec, r1$sec, r1b$sec,
              if (th == 1) identical(r0$draws, r1$draws) else "n/a"))
}
sink()
