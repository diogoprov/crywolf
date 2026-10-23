# Validation of the tuple refactor, part 2: the log-likelihood of the
# previous and the new Stan program at the same parameter values (40 draws
# of a fit of the new program), written with 15 significant digits. Draws of
# a whole run may diverge after many iterations when the compiler arranges
# the arithmetic differently; equal log-likelihoods show the posterior is the
# same. Run from the package root: Rscript dev/profile/validate-tuples-loglik.R
# Output: dev/profile/validate-tuples-loglik.txt
sink(file.path("dev", "profile", "validate-tuples-loglik.txt"), split = TRUE)
devtools::load_all(quiet = TRUE)
cat("CmdStan", cmdstanr::cmdstan_version(), "\n")
old_stan <- tempfile(fileext = ".stan")
writeLines(system2("git", c("show", "v1.6.0:inst/stan/fp_occu.stan"), stdout = TRUE), old_stan)
old <- cmdstanr::cmdstan_model(old_stan, cpp_options = list(stan_threads = TRUE), quiet = TRUE)
new <- fp_compile()
compare <- function(label, data, ...) {
  sd <- .build_stan_data(data, ...)$stan_data
  fit <- new$sample(data = sd, seed = 1, chains = 2, parallel_chains = 2, threads_per_chain = 1,
                    iter_warmup = 200, iter_sampling = 20, refresh = 0, show_messages = FALSE)
  pars <- intersect(c("beta", "alpha", "u_raw", "sigma_u", "fp_gap", "gamma_fp", "delta0",
                      "delta_s", "e_raw", "sigma_fp"), fit$metadata()$stan_variables)
  params <- fit$draws(variables = pars)
  ll <- function(mod) {
    gq <- mod$generate_quantities(fitted_params = params, data = sd, seed = 1,
                                  threads_per_chain = 1, sig_figs = 15)
    unclass(posterior::as_draws_matrix(gq$draws("log_lik")))
  }
  a <- ll(old); b <- ll(new)
  cat(sprintf("%-46s max |difference| %.2e, relative %.2e\n", label, max(abs(a - b)),
              max(abs(a - b) / pmax(1, abs(a)))))
}
s1 <- simulate_fp(nsite = 150, nocc = 8, a_bio = 0.5, sd_recorder = 0.5, conf_effect = 1,
                  conf_cor = 0.5, seed = 3)
s2 <- simulate_fp(nsite = 120, nocc = 6, frac_confirmed = 0.4, frac_verified = 0.2, seed = 6)
s3 <- simulate_fp(nsite = 120, nocc = 6, seed = 8)
set.seed(9)
d4 <- s2$data
d4$auto$noise <- stats::rnorm(nrow(d4$auto))
d4$verified$listener <- stats::rnorm(nrow(d4$verified))
compare("Miller, recorder (GLM path), fp ~ w", s1$data, ~ x, ~ bio, ~ w, TRUE, fp_priors())
compare("Chambert, constant rate", s2$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
compare("no false positives", s3$data, ~ x, ~ 1, NULL, FALSE, fp_priors())
compare("random_fp", s3$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors(), random_fp = TRUE)
compare("Chambert, fp ~ noise, ver ~ listener, markov", d4, ~ x, ~ 1, ~ noise, FALSE,
        fp_priors(), ver = ~ listener, markov = "occasion")
compare("Miller, recorder, markov, det ~ bio", s1$data, ~ x, ~ bio, ~ 1, TRUE, fp_priors(),
        markov = "occasion")
source(file.path("dev", "curiango-common.R"))
auto <- cur_occasions
auto$y <- as.integer(cur_detections$confidence >= 0.40)
d <- suppressWarnings(fp_data(auto, cur_sites, cur_verified, site = "site_night",
                              recorder = "recorder"))
compare("curiango, recorder, det ~ bio + aci", d, ~ temp + rh + moon, ~ bio + aci, ~ 1,
        TRUE, fp_priors())
sink()
