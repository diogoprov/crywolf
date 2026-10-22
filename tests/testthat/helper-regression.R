# Designs of the posterior regression test (test-regression.R). The
# reference values are made by dev/make-regression-refs.R, which sources this
# file; the data are simulated with fixed seeds, so both build the same data.

.regression_designs <- function() {
  s_miller <- simulate_fp(nsite = 150, nocc = 8, a_bio = 0.5, sd_recorder = 0.5,
                          conf_effect = 1, conf_cor = 0.5, seed = 3)
  s_chambert <- simulate_fp(nsite = 120, nocc = 6, frac_confirmed = 0.4,
                            frac_verified = 0.2, seed = 6)
  s_plain <- simulate_fp(nsite = 120, nocc = 6, seed = 8)
  s_options <- simulate_fp(nsite = 120, nocc = 6, frac_verified = 0.3, seed = 9)
  d_options <- s_options$data
  set.seed(9)
  d_options$auto$noise <- stats::rnorm(nrow(d_options$auto))
  d_options$verified$listener <- stats::rnorm(nrow(d_options$verified))
  list(
    miller_recorder = list(data = s_miller$data, occ = ~ x, det = ~ bio, fp = ~ w,
                           random_recorder = TRUE),
    chambert = list(data = s_chambert$data, occ = ~ x, fp = ~ 1,
                    random_recorder = FALSE),
    no_false_positives = list(data = s_plain$data, occ = ~ x, fp = NULL,
                              random_recorder = FALSE),
    all_options = list(data = d_options, occ = ~ x, fp = ~ noise, ver = ~ listener,
                       markov = "occasion", random_fp = TRUE,
                       random_recorder = FALSE)
  )
}

# one fit with fixed seed and one thread, so that a platform reproduces it
.regression_fit <- function(design) {
  suppressWarnings(do.call(fp_occu, c(design, list(
    seed = 1, chains = 4, parallel_chains = 4, threads_per_chain = 1,
    iter_warmup = 500, iter_sampling = 500, refresh = 0, show_messages = FALSE))))
}

# the data of the Stan program for a design, built as fp_occu() builds it
.regression_stan_data <- function(design) {
  crywolf:::.build_stan_data(
    design$data, design$occ, if (is.null(design$det)) ~ 1 else design$det,
    design$fp, design$random_recorder, fp_priors(),
    random_fp = isTRUE(design$random_fp),
    ver = if (is.null(design$ver)) ~ 1 else design$ver,
    markov = design$markov)$stan_data
}

# the parameters block of the Stan program, as written in the output files
.regression_pars <- function(fit) {
  intersect(c("beta", "alpha", "u_raw", "sigma_u", "fp_gap", "gamma_fp", "delta0",
              "delta_s", "e_raw", "sigma_fp"),
            fit$fit$metadata()$stan_variables)
}

# log_lik recomputed from given parameter draws, without sampling
.regression_log_lik <- function(stan_data, params) {
  gq <- fp_compile()$generate_quantities(fitted_params = params, data = stan_data,
                                         seed = 1, threads_per_chain = 1, sig_figs = 15,
                                         output_basename = crywolf:::.output_basename())
  unclass(posterior::as_draws_matrix(gq$draws("log_lik")))
}
