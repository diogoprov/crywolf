# Checks that the likelihood with per-site totals and binomial_logit_glm_lpmf()
# (crywolf >= 1.2.0.9000) gives the same log density, gradient, log_lik and
# p10_bar as the previous model, at the same parameter values, in every
# design the Stan program supports. Usage, from the package root:
#   Rscript dev/profile/validate-glm.R <previous fp_occu.stan>
# Result of 19/10/2026 in validate-glm.rds: all differences at or below 3e-15
# (relative) for the log density and its gradient, 0 for log_lik and p10_bar.
# Threads only change the order of the additions in reduce_sum() (see
# benchmark-curiango.R); the comparison here runs with one thread.
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
opts <- list(cpp_options = list(stan_threads = TRUE), quiet = TRUE,
             compile_model_methods = TRUE, force_recompile = TRUE)
previous <- do.call(cmdstanr::cmdstan_model,
                    c(list(commandArgs(trailingOnly = TRUE)[1]), opts))
current <- do.call(cmdstanr::cmdstan_model, c(list(fp_stan_file()), opts))

check <- function(label, sd) {
  set.seed(1)
  init <- list(beta = rnorm(sd$Kbeta, 0, 0.5), alpha = rnorm(sd$Kalpha, 0, 0.5),
               u_raw = rnorm(sd$nrec), sigma_u = 0.7, fp_gap = 2.5,
               gamma_fp = rnorm(sd$Kfp, 0, 0.3), delta0 = 0.3)
  if (isTRUE(sd$estimate_fp_re == 1)) {
    init$e_raw <- rnorm(sd$nsite)
    init$sigma_fp <- array(0.8)
  }
  fit <- function(mod) {
    mod$sample(data = sd, chains = 1, iter_sampling = 1, iter_warmup = 0,
               fixed_param = TRUE, init = list(init), seed = 1, refresh = 0,
               show_messages = FALSE, threads_per_chain = 1)
  }
  fp <- fit(previous)
  fc <- fit(current)
  fp$init_model_methods()
  fc$init_model_methods()
  up <- as.numeric(fp$unconstrain_draws(format = "draws_matrix")[1, ])
  gp <- fp$grad_log_prob(up)
  gc <- fc$grad_log_prob(up)
  data.frame(
    design = label, nsite = sd$nsite, N2 = sd$N2,
    use_glm = sd$Kalpha > 1 && sd$obs_confirm == 0,
    max_loglik_diff = max(abs(as.numeric(fp$draws("log_lik")) -
                              as.numeric(fc$draws("log_lik")))),
    lp_rel_diff = abs(fp$log_prob(up) - fc$log_prob(up)) / abs(fp$log_prob(up)),
    grad_rel_diff = max(abs(gp - gc)) / max(abs(gp)),
    p10bar_diff = abs(as.numeric(fp$draws("p10_bar")) -
                        as.numeric(fc$draws("p10_bar")))
  )
}

stan_data <- function(s, occ, det, fp, random_recorder, random_fp = FALSE) {
  crywolf:::.build_stan_data(s$data, occ, det, fp, random_recorder, fp_priors(),
                             random_fp = random_fp)$stan_data
}
# a recorder that changes within sites, to exercise the per-row intercepts
mixed_recorders <- function(s) {
  s$data$auto$recorder <- sample(1:5, nrow(s$data$auto), replace = TRUE)
  s
}

designs <- list(
  "curiango: Miller, det~bio+aci, RE (GLM, recorder constant within site)" =
    curiango_bench_data(),
  "Miller, det~bio, RE varying within site (GLM, vector intercept)" =
    stan_data(mixed_recorders(simulate_fp(nsite = 150, nocc = 8, a_bio = 0.5,
                                          sd_recorder = 0.5, seed = 3)),
              ~ x, ~ bio, ~ 1, TRUE),
  "Miller, det~1 collapsed, RE varying within site (vectorized, u[rec2])" =
    stan_data(mixed_recorders(simulate_fp(nsite = 150, nocc = 8,
                                          sd_recorder = 0.5, seed = 4)),
              ~ x, ~ 1, ~ 1, TRUE),
  "Chambert, det~bio, fp~w (vectorized, log1p_exp_lp10[loc])" =
    stan_data(simulate_fp(nsite = 120, nocc = 6, a_bio = 0.5, conf_effect = 1,
                          conf_cor = 0.5, frac_confirmed = 0.4,
                          frac_verified = 0.2, seed = 6),
              ~ x, ~ bio, ~ w, FALSE),
  "Chambert, det~bio, RE and random_fp (vectorized)" =
    stan_data(simulate_fp(nsite = 120, nocc = 6, a_bio = 0.5, sd_recorder = 0.5,
                          frac_confirmed = 0.4, frac_verified = 0.2, seed = 11),
              ~ x, ~ bio, ~ 1, TRUE, random_fp = TRUE),
  "Miller, det~1 collapsed, fp~x (vectorized, totals)" =
    stan_data(simulate_fp(nsite = 200, nocc = 10, conf_effect = 1,
                          conf_cor = 0.5, seed = 8),
              ~ x, ~ 1, ~ x, FALSE),
  "Miller, det~bio+aci, RE, random_fp (GLM)" =
    stan_data(simulate_fp(nsite = 100, nocc = 8, a_bio = 0.3, a_aci = -0.3,
                          sd_recorder = 0.4, sigma_p10 = 1, seed = 9),
              ~ x, ~ bio + aci, ~ 1, TRUE, random_fp = TRUE),
  "standard occupancy fp = NULL, det~bio (GLM)" =
    stan_data(simulate_fp(nsite = 100, nocc = 8, a_bio = 0.3, seed = 10),
              ~ x, ~ bio, NULL, FALSE),
  "standard occupancy fp = NULL, det~1 (vectorized)" =
    stan_data(simulate_fp(nsite = 100, nocc = 8, seed = 12),
              ~ x, ~ 1, NULL, FALSE)
)
out <- purrr::imap(designs, function(sd, label) check(label, sd)) |>
  purrr::list_rbind()
print(out, digits = 3)
saveRDS(out, "dev/profile/validate-glm.rds")
