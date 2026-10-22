# Pre-computes the results of vignettes/fp_variation.Rmd. Two simulated
# datasets: (1) false-positive rates that vary at random among sites,
# fitted with a constant rate and with random_fp = TRUE; (2) a rate that
# depends on background noise measured on each occasion, fitted with a
# constant rate and with fp = ~ noise. Run from the package root after
# installing the package: Rscript dev/precompute-fp-variation.R
library(crywolf)
fit <- function(d, ...) {
  t0 <- Sys.time()
  f <- suppressWarnings(fp_occu(d, occ = ~ x, random_recorder = FALSE, seed = 1,
                                refresh = 0, show_messages = FALSE, ...))
  list(table = summary(f)$table, diagnostics = f$diagnostics,
       sites = as.data.frame(fp_sites(f, ndraws = 1000, seed = 1)),
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
       loo = loo(f))
}

# 1. among-site variation with no covariate behind it
s1 <- simulate_fp(nsite = 200, nocc = 10, sigma_p10 = 1, frac_verified = 0.3, seed = 1)
among <- list(truth = s1$truth[c("p10", "p10_marginal", "sigma_p10", "p10_site")],
              constant = fit(s1$data, fp = ~ 1),
              random = fit(s1$data, fp = ~ 1, random_fp = TRUE))

# 2. occasion-level noise: simulated with the generator of simulate.fp_occu()
s2 <- simulate_fp(nsite = 200, nocc = 10, frac_verified = 0.3, seed = 2)
d2 <- s2$data
set.seed(2)
d2$auto$noise <- stats::rnorm(nrow(d2$auto))
sl <- crywolf:::.build_stan_data(d2, ~ x, ~ 1, ~ noise, FALSE, fp_priors())
gen <- list(data = d2, formulas = list(occ = ~ x, det = ~ 1, fp = ~ noise, ver = ~ 1),
            stan_data = sl$stan_data, design = sl$design, occasions = sl$occasions)
truth2 <- list(beta = c(0, -0.8), alpha = stats::qlogis(0.6), delta0 = stats::qlogis(0.6),
               delta_s = numeric(0), gamma0 = stats::qlogis(0.05), gamma_fp = 1)
set.seed(3)
sim2 <- crywolf:::.simulate_from(gen, truth2)
within <- list(truth = c(occ_x = -0.8, `fp_(Intercept)` = stats::qlogis(0.05), fp_noise = 1,
                         p10_bar = sim2$truth$p10_bar),
               constant = fit(sim2$data, fp = ~ 1),
               noise = fit(sim2$data, fp = ~ noise))
saveRDS(list(among = among, within = within),
        file.path("vignettes", "fp-variation-results.rds"), compress = "xz")
