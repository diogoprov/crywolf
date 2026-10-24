# The R likelihood behind fp_ppc(), fp_sites(), simulate() and loo(reloo =
# TRUE) must equal the log_lik of the Stan program in every combination of
# the model options. Fits write 15 significant digits (sig_figs), so the two
# agree to rounding. Moved here from dev/validate-model-options.R; the
# identity of draws across versions is covered by test-regression.R.

ll_gap <- function(data, ...) {
  m <- suppressWarnings(fp_occu(data, occ = ~ x, ..., chains = 1, iter_warmup = 150,
                                iter_sampling = 40,
                                threads_per_chain = 1, seed = 3, refresh = 0,
                                show_messages = FALSE, sig_figs = 15))
  pars <- crywolf:::.ppc_pars(m$fit, m$stan_data, NULL)
  r <- crywolf:::.ppc_core(m$stan_data, pars, "automated", "freeman-tukey")
  ll <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
  max(abs(r$ll - ll) / pmax(1, abs(ll)))
}

test_that("the R likelihood matches the Stan log_lik with every option", {
  skip_if_no_cmdstan()
  s_miller <- simulate_fp(nsite = 80, nocc = 6, a_bio = 0.5, sd_recorder = 0.5,
                          conf_effect = 1, conf_cor = 0.5, seed = 3)
  s_chambert <- simulate_fp(nsite = 80, nocc = 6, frac_confirmed = 0.4,
                            frac_verified = 0.2, seed = 6)
  s_plain <- simulate_fp(nsite = 80, nocc = 6, seed = 8)
  set.seed(1)
  d_ver <- s_miller$data
  d_ver$verified$listener <- stats::rnorm(nrow(d_ver$verified))
  d_noise <- s_miller$data
  d_noise$auto$noise <- stats::rnorm(nrow(d_noise$auto))
  d_chambert <- s_chambert$data
  d_chambert$auto$noise <- stats::rnorm(nrow(d_chambert$auto))
  d_chambert$verified$listener <- stats::rnorm(nrow(d_chambert$verified))

  gaps <- c(
    verified_covariate = ll_gap(d_ver, det = ~ bio, fp = ~ 1, ver = ~ listener,
                                random_recorder = TRUE),
    fp_occasion_and_site = ll_gap(d_noise, det = ~ bio, fp = ~ noise + x,
                                  random_recorder = TRUE),
    fp_occasion_chambert = ll_gap(d_chambert, fp = ~ noise, random_recorder = FALSE),
    markov_glm = ll_gap(s_miller$data, det = ~ bio, fp = ~ 1, markov = "occasion",
                        random_recorder = TRUE),
    markov_random_fp = ll_gap(s_plain$data, fp = ~ 1, markov = "occasion",
                              random_fp = TRUE, random_recorder = FALSE),
    markov_noise_ver_chambert = ll_gap(d_chambert, fp = ~ noise, ver = ~ listener,
                                       markov = "occasion", random_recorder = FALSE),
    # site effect on true detection in the GLM, vectorized and confirmation
    # paths, alone and with the site effect on false positives
    random_det_glm = ll_gap(s_miller$data, det = ~ bio, fp = ~ 1, random_det = TRUE,
                            random_recorder = TRUE),
    random_det_rows = ll_gap(s_plain$data, fp = ~ 1, random_det = TRUE,
                             random_recorder = FALSE),
    random_det_chambert = ll_gap(d_chambert, fp = ~ noise, random_det = TRUE,
                                 random_recorder = FALSE),
    random_det_and_fp = ll_gap(s_plain$data, fp = ~ 1, random_det = TRUE,
                               random_fp = TRUE, random_recorder = FALSE)
  )
  purrr::iwalk(gaps, function(g, name) {
    expect_lt(g, 1e-8, label = paste(name, "relative difference R vs Stan"))
  })
})
