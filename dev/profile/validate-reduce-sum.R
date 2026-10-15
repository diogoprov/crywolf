# Checks that the reduce_sum() version of the Stan model gives the same
# log density, gradient and per-site log_lik as the previous version (main).
# Usage, from the package root:
#   git show <ref>:inst/stan/fp_occu.stan > <dir>/old/fp_occu.stan
#   Rscript dev/profile/validate-reduce-sum.R <dir>/old/fp_occu.stan
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
args <- commandArgs(trailingOnly = TRUE)
old <- cmdstanr::cmdstan_model(args[1], compile_model_methods = TRUE,
                               force_recompile = TRUE, quiet = TRUE)
new_dir <- file.path(dirname(dirname(args[1])), "new")
dir.create(new_dir, showWarnings = FALSE)
file.copy(fp_stan_file(), new_dir, overwrite = TRUE)
new <- cmdstanr::cmdstan_model(file.path(new_dir, "fp_occu.stan"),
                               compile_model_methods = TRUE, force_recompile = TRUE,
                               cpp_options = list(stan_threads = TRUE), quiet = TRUE)
stopifnot(crywolf:::.has_threads(new))

sim <- simulate_fp(nsite = 120, nocc = 8, conf_effect = 1.5, conf_cor = 0.7,
                   sd_recorder = 0.5, seed = 3)
sim_chambert <- simulate_fp(nsite = 120, nocc = 8, frac_confirmed = 0.3, seed = 4)
b <- function(d, occ, det, fp, re) {
  crywolf:::.build_stan_data(d, occ, det, fp, re, fp_priors())$stan_data
}
cases <- list(
  miller_const = b(sim$data, ~ x, ~ 1, ~ 1, FALSE),
  miller_fp_x_re = b(sim$data, ~ x, ~ 1, ~ x, TRUE),
  miller_det_cov = b(sim$data, ~ x, ~ bio, ~ w, TRUE),
  standard = b(sim$data, ~ x, ~ 1, NULL, FALSE),
  chambert = b(sim_chambert$data, ~ x, ~ 1, ~ 1, FALSE),
  chambert_re = b(sim_chambert$data, ~ x, ~ 1, ~ x, TRUE),
  curiango = curiango_bench_data()
)

out_old <- tempfile("old"); out_new <- tempfile("new")
dir.create(out_old); dir.create(out_new)
compare <- function(sd, name, n_points = 20) {
  # a minimal run, only to get objects with the model methods
  fo <- old$sample(data = sd, chains = 1, iter_warmup = 1, iter_sampling = 1,
                   refresh = 0, show_messages = FALSE, seed = 1, init = 0, output_dir = out_old)
  fn <- new$sample(data = sd, chains = 1, iter_warmup = 1, iter_sampling = 1,
                   refresh = 0, show_messages = FALSE, seed = 1, init = 0,
                   threads_per_chain = 4, output_dir = out_new)
  fo$init_model_methods(); fn$init_model_methods()
  set.seed(1)
  res <- purrr::map(seq_len(n_points), function(k) {
    th <- stats::rnorm(with(sd, Kbeta + Kalpha + nrec + 3 + Kfp), 0, 0.7)
    go <- fo$grad_log_prob(th); gn <- fn$grad_log_prob(th)
    c(lp = abs(attr(go, "log_prob") - attr(gn, "log_prob")) /
        max(1, abs(attr(go, "log_prob"))),
      grad = max(abs(go - gn) / pmax(1, abs(go))))
  })
  r <- do.call(rbind, res)
  data.frame(case = name, max_rel_lp = max(r[, "lp"]),
             max_rel_grad = max(r[, "grad"]), finite = all(is.finite(r)))
}
tab <- purrr::imap(cases, compare) |> purrr::list_rbind()
print(tab, digits = 3)

# per-site log_lik on the same draws: fit the old model briefly and recompute
# its generated quantities with the new one
gq_check <- function(sd, name) {
  fo <- old$sample(data = sd, chains = 2, iter_warmup = 200, iter_sampling = 100,
                   refresh = 0, show_messages = FALSE, seed = 2,
                   parallel_chains = 2)
  gq <- new$generate_quantities(fo, data = sd, seed = 2, threads_per_chain = 2,
                                parallel_chains = 2)
  a <- posterior::as_draws_matrix(fo$draws("log_lik"))
  bb <- posterior::as_draws_matrix(gq$draws("log_lik"))
  data.frame(case = name, max_abs_diff_log_lik = max(abs(a - bb)))
}
gq <- purrr::imap(cases[names(cases) != "curiango"], gq_check) |> purrr::list_rbind()
print(gq, digits = 3)
saveRDS(list(grad = tab, log_lik = gq), "dev/profile/validate-reduce-sum.rds")
