# two species on the same sites and occasions: the data of one simulation,
# with the detections of the second species drawn again
ms_sim <- function(seed, ...) {
  sim <- simulate_fp(nsite = 40, nocc = 6, a_bio = 0.5, sd_recorder = 0.5, seed = seed, ...)
  d <- sim$data
  set.seed(seed)
  auto <- transform(d$auto, y = stats::rbinom(nrow(d$auto), 1, 0.3))
  verified <- transform(d$verified, y = stats::rbinom(nrow(d$verified), 1, 0.4))
  b <- suppressWarnings(fp_data(auto, d$sites, verified, recorder = d$recorder))
  tag <- function(x, sp) transform(x, sp = sp)
  ms <- suppressWarnings(fp_data(rbind(tag(d$auto, "a"), tag(auto, "b")), d$sites,
                                 rbind(tag(d$verified, "a"), tag(verified, "b")),
                                 recorder = d$recorder, species = "sp"))
  list(a = d, b = b, ms = ms)
}

test_that("fp_occu_ms() checks its data and arguments", {
  s <- ms_sim(1)
  expect_error(fp_occu_ms(s$a), "species column")
  expect_error(fp_occu_ms(s$ms, markov = "occasion"), "not available")
  expect_error(fp_occu_ms(s$ms, random_det = TRUE), "not available")
  sl <- crywolf:::.build_ms_stan_data(s$ms, ~ x, ~ bio, ~ 1, TRUE, fp_priors())
  sd <- sl$stan_data
  expect_equal(c(sd$nsite, sd$nspecies, sd$nunit), c(40, 2, 80))
  expect_equal(sd$sp_unit, rep(1:2, 40))
  # unit (site - 1) * 2 + species; the verified detections of site 1 by species
  kp <- s$ms$known_present
  expect_equal(sd$known_present, as.vector(t(kp)))
  expect_equal(sl$design$species, c("a", "b"))
  expect_equal(sd$X[1:2, 2], rep(s$ms$sites$x[1], 2))
})

test_that("the multispecies likelihood equals that of fp_occu() species by species", {
  skip_if_no_cmdstan()
  gq <- function(model, pars, stan_data) {
    # a folder per call, as cmdstanr names files by the minute
    out <- tempfile("gq")
    dir.create(out)
    g <- crywolf:::.compile_model(model)$generate_quantities(
      fitted_params = pars, data = stan_data, seed = 1, sig_figs = 18,
      threads_per_chain = 1, output_dir = out)
    unclass(posterior::as_draws_matrix(g$draws("log_lik")))
  }
  gap <- function(s, det, fp) {
    m <- suppressWarnings(fp_occu(s$a, occ = ~ x, det = det, fp = fp, chains = 1,
                                  iter_warmup = 100, iter_sampling = 20,
                                  threads_per_chain = 1, seed = 3, refresh = 0,
                                  show_messages = FALSE, sig_figs = 18))
    pa <- posterior::as_draws_df(m$fit$draws())
    pa <- pa[, !grepl("^(lp__|u|gamma0|log_lik|p11_ref|p10_ref|p10_bar|r11)(\\[|$)",
                      names(pa))]
    # species b: its own occupancy coefficients, detection intercept and gap
    set.seed(4)
    pb <- pa
    pb[["beta[1]"]] <- stats::rnorm(20)
    pb[["beta[2]"]] <- stats::rnorm(20)
    pb[["alpha[1]"]] <- stats::rnorm(20)
    pb$fp_gap <- stats::runif(20, 0.5, 3)
    sl_a <- crywolf:::.build_stan_data(s$a, ~ x, det, fp, TRUE, fp_priors())
    sl_b <- crywolf:::.build_stan_data(s$b, ~ x, det, fp, TRUE, fp_priors())
    ll_a <- gq("fp_occu", posterior::as_draws_array(pa), sl_a$stan_data)
    ll_b <- gq("fp_occu", posterior::as_draws_array(pb), sl_b$stan_data)
    # the multispecies parameters that give each species its own values:
    # community means 0 (gap 1), standard deviations 1
    pm <- pa[, !grepl("^(beta|fp_gap)", names(pa))]
    pm[["alpha[1]"]] <- 0
    pm[c("mu_beta[1]", "mu_beta[2]", "sigma_beta[1]", "sigma_beta[2]")] <-
      list(0, 0, 1, 1)
    pm[c("zb[1,1]", "zb[2,1]", "zb[1,2]", "zb[2,2]")] <-
      list(pa[["beta[1]"]], pa[["beta[2]"]], pb[["beta[1]"]], pb[["beta[2]"]])
    pm$sigma_alpha <- 1
    pm[c("za[1]", "za[2]")] <- list(pa[["alpha[1]"]], pb[["alpha[1]"]])
    pm[c("mu_gap", "sigma_gap")] <- list(1, 1)
    pm[c("zg[1]", "zg[2]")] <- list(log(pa$fp_gap), log(pb$fp_gap))
    # CmdStan reads the parameters in the order of their declaration
    first <- c("^mu_beta", "^sigma_beta", "^zb", "^alpha", "^sigma_alpha$", "^za",
               "^mu_gap$", "^sigma_gap$", "^zg", "^u_raw", "^sigma_u$", "^gamma_fp",
               "^delta0$", "^delta_s")
    pm <- pm[, c(unlist(lapply(first, grep, names(pm), value = TRUE)),
                 ".chain", ".iteration", ".draw")]
    sl_m <- crywolf:::.build_ms_stan_data(s$ms, ~ x, det, fp, TRUE, fp_priors())
    ll_m <- gq("fp_occu_ms", posterior::as_draws_array(pm), sl_m$stan_data)
    # the same values in the centred parametrization: with community means 0
    # (gap 1) and standard deviations 1, only the recorder effects change
    pc <- pm
    for (u in grep("^u_raw", names(pc), value = TRUE)) pc[[u]] <- pc[[u]] * pc$sigma_u
    sl_c <- crywolf:::.build_ms_stan_data(s$ms, ~ x, det, fp, TRUE, fp_priors(), centered = TRUE)
    ll_c <- gq("fp_occu_ms", posterior::as_draws_array(pc), sl_c$stan_data)
    a <- cbind(ll_a, ll_b)
    by_species <- c(seq(1, by = 2, length.out = 40), seq(2, by = 2, length.out = 40))
    max(abs(a - ll_m[, by_species]) / pmax(1, abs(a)),
        abs(a - ll_c[, by_species]) / pmax(1, abs(a)))
  }
  s <- ms_sim(3)
  gaps <- c(site_fp = gap(s, ~ bio, ~ x), occasion_fp = gap(s, ~ bio, ~ aci))
  purrr::iwalk(gaps, function(g, name) {
    expect_lt(g, 1e-12, label = paste(name, "relative difference multispecies vs fp_occu"))
  })
})

# two species with the sites, occasions and recorders of one simulation: the
# rows of species b are those of species a moved to other sites, so that
# each species is also a valid data set of fp_occu(), in either design
ms_sim_moved <- function(seed, ...) {
  sim <- simulate_fp(nsite = 40, nocc = 6, a_bio = 0.5, sd_recorder = 0.5, seed = seed, ...)
  d <- sim$data
  set.seed(seed)
  d$verified$listener <- stats::rnorm(nrow(d$verified))
  move <- sample(d$sites$site)
  moved <- function(x) transform(x, site = move[match(x$site, d$sites$site)])
  conf <- d$confirmed
  b <- suppressWarnings(fp_data(moved(d$auto), d$sites, moved(d$verified),
                                recorder = d$recorder, confirmed = conf))
  tag <- function(x, sp) transform(x, sp = sp)
  ms <- suppressWarnings(fp_data(rbind(tag(d$auto, "a"), tag(moved(d$auto), "b")), d$sites,
                                 rbind(tag(d$verified, "a"), tag(moved(d$verified), "b")),
                                 recorder = d$recorder, confirmed = conf, species = "sp"))
  list(a = d, b = b, ms = ms)
}

test_that("with by_species the likelihood equals that of fp_occu() species by species", {
  skip_if_no_cmdstan()
  gq <- function(model, pars, stan_data) {
    out <- tempfile("gq")
    dir.create(out)
    g <- crywolf:::.compile_model(model)$generate_quantities(
      fitted_params = pars, data = stan_data, seed = 1, sig_figs = 18,
      threads_per_chain = 1, output_dir = out)
    unclass(posterior::as_draws_matrix(g$draws("log_lik")))
  }
  first <- c("^mu_beta", "^sigma_beta", "^zb", "^alpha", "^sigma_alpha$", "^za",
             "^mu_gap$", "^sigma_gap$", "^zg", "^u_raw", "^sigma_u$", "^gamma_fp",
             "^delta0$", "^delta_s", "^sigma_det_s", "^zds", "^sigma_fp_s", "^zfs",
             "^sigma_ver", "^zv")
  ordered <- function(p) {
    p[, c(unlist(lapply(first, grep, names(p), value = TRUE)), ".chain", ".iteration", ".draw")]
  }
  gap <- function(s, det, fp, ver = ~ 1) {
    by <- c("det", "fp", "ver")
    m <- suppressWarnings(fp_occu(s$a, occ = ~ x, det = det, fp = fp, ver = ver, chains = 1,
                                  iter_warmup = 100, iter_sampling = 20,
                                  threads_per_chain = 1, seed = 3, refresh = 0,
                                  show_messages = FALSE, sig_figs = 18))
    pa <- posterior::as_draws_df(m$fit$draws())
    pa <- pa[, !grepl("^(lp__|u|gamma0|log_lik|p11_ref|p10_ref|p10_bar|r11)(\\[|$)",
                      names(pa))]
    # species b: its own values of every term of by_species too
    set.seed(4)
    pb <- pa
    for (v in c("beta[1]", "beta[2]", "alpha[1]", "alpha[2]", "gamma_fp[1]", "delta0")) {
      pb[[v]] <- stats::rnorm(20)
    }
    pb$fp_gap <- stats::runif(20, 0.5, 3)
    sl_a <- crywolf:::.build_stan_data(s$a, ~ x, det, fp, TRUE, fp_priors(), ver = ver)
    sl_b <- crywolf:::.build_stan_data(s$b, ~ x, det, fp, TRUE, fp_priors(), ver = ver)
    ll_a <- gq("fp_occu", posterior::as_draws_array(pa), sl_a$stan_data)
    ll_b <- gq("fp_occu", posterior::as_draws_array(pb), sl_b$stan_data)
    # community values 0 (gap 1) and standard deviations 1, so that the
    # standard normal draws of each species are its values
    pm <- pa[, !grepl("^(beta|fp_gap)", names(pa))]
    pm[c("mu_beta[1]", "mu_beta[2]", "sigma_beta[1]", "sigma_beta[2]")] <- list(0, 0, 1, 1)
    pm[c("zb[1,1]", "zb[2,1]", "zb[1,2]", "zb[2,2]")] <-
      list(pa[["beta[1]"]], pa[["beta[2]"]], pb[["beta[1]"]], pb[["beta[2]"]])
    pm[c("alpha[1]", "alpha[2]", "sigma_alpha")] <- list(0, 0, 1)
    pm[c("za[1]", "za[2]")] <- list(pa[["alpha[1]"]], pb[["alpha[1]"]])
    pm[c("mu_gap", "sigma_gap")] <- list(1, 1)
    pm[c("zg[1]", "zg[2]")] <- list(log(pa$fp_gap), log(pb$fp_gap))
    pm[c("sigma_det_s[1]", "zds[1,1]", "zds[1,2]")] <- list(1, pa[["alpha[2]"]], pb[["alpha[2]"]])
    pm[c("gamma_fp[1]", "sigma_fp_s[1]", "zfs[1,1]", "zfs[1,2]")] <-
      list(0, 1, pa[["gamma_fp[1]"]], pb[["gamma_fp[1]"]])
    pm[c("delta0", "sigma_ver[1]", "zv[1]", "zv[2]")] <- list(0, 1, pa$delta0, pb$delta0)
    pm <- ordered(pm)
    sl_m <- crywolf:::.build_ms_stan_data(s$ms, ~ x, det, fp, TRUE, fp_priors(), ver = ver,
                                          by_species = by)
    ll_m <- gq("fp_occu_ms", posterior::as_draws_array(pm), sl_m$stan_data)
    # centred: with these community values only the recorder effects change
    pc <- pm
    for (u in grep("^u_raw", names(pc), value = TRUE)) pc[[u]] <- pc[[u]] * pc$sigma_u
    sl_c <- crywolf:::.build_ms_stan_data(s$ms, ~ x, det, fp, TRUE, fp_priors(), ver = ver,
                                          centered = TRUE, by_species = by)
    ll_c <- gq("fp_occu_ms", posterior::as_draws_array(pc), sl_c$stan_data)
    # deviations 0 give the model with shared terms
    p0 <- pm
    p0[grep("^(zds|zfs|zv)", names(p0))] <- 0
    p0[c("alpha[2]", "gamma_fp[1]")] <- list(pa[["alpha[2]"]], pa[["gamma_fp[1]"]])
    ll_0 <- gq("fp_occu_ms", posterior::as_draws_array(p0), sl_m$stan_data)
    p_shared <- p0[, !grepl("^(sigma_det_s|zds|sigma_fp_s|zfs|sigma_ver|zv)", names(p0))]
    sl_s <- crywolf:::.build_ms_stan_data(s$ms, ~ x, det, fp, TRUE, fp_priors(), ver = ver)
    ll_s <- gq("fp_occu_ms", posterior::as_draws_array(p_shared), sl_s$stan_data)
    a <- cbind(ll_a, ll_b)
    by_species <- c(seq(1, by = 2, length.out = 40), seq(2, by = 2, length.out = 40))
    rel <- function(x, y) max(abs(x - y) / pmax(1, abs(x)))
    c(species = max(rel(a, ll_m[, by_species]), rel(a, ll_c[, by_species])),
      shared = rel(ll_s, ll_0))
  }
  designs <- list(
    multiple_methods = list(s = ms_sim_moved(3), det = ~ bio, fp = ~ aci),
    unit_fp = list(s = ms_sim_moved(3), det = ~ bio, fp = ~ x),
    verified_covariate = list(s = ms_sim_moved(3), det = ~ bio, fp = ~ aci, ver = ~ listener),
    confirmation = list(s = ms_sim_moved(6, frac_confirmed = 0.4, frac_verified = 0.2),
                        det = ~ bio, fp = ~ aci),
    confirmation_unit_fp = list(s = ms_sim_moved(6, frac_confirmed = 0.4, frac_verified = 0.2),
                                det = ~ bio, fp = ~ x)
  )
  purrr::iwalk(designs, function(d, name) {
    g <- gap(d$s, d$det, d$fp, if (is.null(d$ver)) ~ 1 else d$ver)
    expect_lt(g[["species"]], 1e-12, label = paste(name, ": by species vs fp_occu"))
    expect_lt(g[["shared"]], 1e-12, label = paste(name, ": zero deviations vs shared terms"))
  })
})

test_that("by_species is checked and named", {
  s <- ms_sim(1)
  expect_error(fp_occu_ms(s$ms, by_species = "occ"), "must name terms")
  expect_error(fp_occu_ms(s$ms, by_species = NA_character_), "must name terms")
  expect_identical(crywolf:::.check_by_species(NULL, s$ms), character())
  sl <- crywolf:::.build_ms_stan_data(s$ms, ~ x, ~ bio, ~ aci, TRUE, fp_priors(),
                                      by_species = c("det", "ver"))
  expect_equal(unlist(sl$stan_data[c("by_det", "by_fp", "by_ver")]),
               c(by_det = 1, by_fp = 0, by_ver = 1))
  no_ver <- s$ms
  no_ver$verified <- no_ver$verified[0, ]
  expect_error(crywolf:::.check_by_species("ver", no_ver), "needs verified")
})

test_that("a short multispecies fit has readable names and methods", {
  skip_if_no_cmdstan()
  s <- ms_sim(5)
  m <- suppressWarnings(fp_occu_ms(s$ms, occ = ~ x, det = ~ bio, fp = ~ aci,
                                   by_species = c("det", "fp", "ver"), chains = 1,
                                   iter_warmup = 150, iter_sampling = 100,
                                   threads_per_chain = 1, seed = 1, refresh = 0,
                                   show_messages = FALSE))
  expect_s3_class(m, "fp_occu_ms")
  v <- posterior::variables(m$draws)
  expect_true(all(c("occ_(Intercept)", "occ_x", "sd_occ_(Intercept)", "sd_occ_x",
                    "det_(Intercept)", "det_bio", "sd_det_(Intercept)",
                    "fp_(Intercept)", "fp_gap", "sd_log_fp_gap", "p11_ref", "p10_ref",
                    "sd_recorder", "verified_(Intercept)", "p_verified",
                    "occ_x[a]", "occ_(Intercept)[b]", "det_(Intercept)[a]",
                    "fp_gap[b]", "p10_ref[a]", "p11_ref[b]",
                    "sd_det_bio", "det_bio[a]", "sd_fp_aci", "fp_aci[b]",
                    "sd_verified_(Intercept)", "verified_(Intercept)[a]",
                    "p_verified[b]") %in% v))
  expect_false(any(grepl("^(mu_|z|beta|alpha|gap_sp)", v)))
  # every species keeps p11 above p10
  expect_true(all(m$draws[["p11_ref[a]"]] > m$draws[["p10_ref[a]"]]))
  expect_output(print(m), "2 species at 40 sites")
  expect_output(print(m), "intercept and slopes by species")
  expect_identical(m$by_species, c("det", "fp", "ver"))
  expect_s3_class(summary(m), "summary.fp_occu")
  expect_identical(as_draws_df(m), m$draws)
})

test_that("the centred parametrization gives the same posterior", {
  skip_if_no_cmdstan()
  s <- ms_sim(5)
  fit <- function(centered) {
    suppressWarnings(fp_occu_ms(s$ms, occ = ~ 1, det = ~ 1, random_recorder = TRUE,
                                centered = centered, chains = 2, iter_warmup = 500,
                                iter_sampling = 1000, threads_per_chain = 1, seed = 2,
                                refresh = 0, show_messages = FALSE))
  }
  a <- fit(FALSE)
  b <- fit(TRUE)
  expect_identical(b$stan_data$centered, 1L)
  v <- c("p11_ref[a]", "p10_ref[b]", "occ_(Intercept)[a]", "fp_gap[b]", "sd_recorder")
  da <- as.data.frame(a$draws)[, v]
  db <- as.data.frame(b$draws)[, v]
  gap <- (colMeans(da) - colMeans(db)) / vapply(da, stats::sd, numeric(1))
  expect_lt(max(abs(gap)), 0.3)
  expect_error(fp_occu_ms(s$ms, centered = 1), "must be TRUE or FALSE")
})

test_that("simulate_fp_ms() shares the design among species", {
  s <- simulate_fp_ms(nspecies = 3, nsite = 30, nocc = 5, frac_confirmed = 0.3,
                      sd_recorder = 0.5, seed = 1)
  expect_s3_class(s, "fp_sim_ms")
  expect_identical(s$data$species_names, c("sp1", "sp2", "sp3"))
  expect_equal(dim(s$truth$z), c(30, 3))
  # every occasion for every species, with the same covariates and checks
  by_sp <- split(s$data$auto, s$data$auto$species)
  expect_true(all(vapply(by_sp, nrow, 1L) == 150))
  cols <- c("site", "occasion", "bio", "recorder")
  expect_equal(by_sp$sp1[cols], by_sp$sp2[cols], ignore_attr = TRUE)
  expect_identical(is.na(by_sp$sp1$n_true), is.na(by_sp$sp3$n_true))
  expect_true(all(c("occ_x", "sd_occ_x", "fp_gap", "sd_log_fp_gap", "sd_recorder",
                    "occ_x[sp2]", "fp_gap[sp3]", "p10_ref[sp1]") %in% names(s$truth$pars)))
  # a verified detection only where the species is present
  v <- s$data$verified
  zv <- s$truth$z[cbind(v$site, match(v$species, colnames(s$truth$z)))]
  expect_true(all(zv[v$y == 1] == 1))
  expect_identical(simulate_fp_ms(nspecies = 3, nsite = 30, nocc = 5, seed = 1)$data,
                   simulate_fp_ms(nspecies = 3, nsite = 30, nocc = 5, seed = 1)$data)
  expect_output(print(s), "3 species at 30 sites")
  # p11_ref as fp_occu_ms() defines it under observation confirmation
  p11 <- stats::plogis(s$truth$pars[["det_(Intercept)[sp1]"]])
  expect_equal(s$truth$pars[["p11_ref[sp1]"]],
               1 - (1 - p11) * (1 - s$truth$pars[["p10_ref[sp1]"]]))
})

test_that("the outputs of fp_occu_ms() use the likelihood of each species", {
  skip_if_no_cmdstan()
  fit <- function(s, ...) {
    suppressWarnings(fp_occu_ms(s$data, occ = ~ x, det = ~ bio, ..., chains = 1,
                                iter_warmup = 150, iter_sampling = 30,
                                threads_per_chain = 1, seed = 1, refresh = 0,
                                show_messages = FALSE, sig_figs = 15))
  }
  # the log-likelihood of each species in R, from the code of fp_occu(),
  # against the log_lik of each unit of the Stan program
  gap <- function(m) {
    sp <- crywolf:::.ms_species(m, NULL)
    ll <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
    max(vapply(seq_along(sp), function(s) {
      r <- crywolf:::.ppc_core(sp[[s]]$stan_data, sp[[s]]$pars, "automated",
                               "freeman-tukey", ll_only = TRUE)$ll
      cols <- (seq_len(ncol(r)) - 1) * length(sp) + s
      max(abs(r - ll[, cols]) / pmax(1, abs(ll[, cols])))
    }, 1))
  }
  s1 <- simulate_fp_ms(nspecies = 3, nsite = 40, nocc = 6, sd_recorder = 0.5,
                       a_bio = 0.5, ngroup = 8, sd_fp_group = 0.5, seed = 1)
  s1$data$verified$listener <- stats::rnorm(nrow(s1$data$verified))
  m1 <- fit(s1, fp = ~ aci, ver = ~ listener, by_species = c("det", "fp", "ver"),
            fp_group = "group", random_fp = TRUE)
  s2 <- simulate_fp_ms(nspecies = 3, nsite = 40, nocc = 6, frac_confirmed = 0.4, seed = 2)
  m2 <- fit(s2, fp = ~ x, random_recorder = FALSE, by_species = c("det", "fp"))
  expect_lt(gap(m1), 1e-10, label = "multiple methods, R vs Stan")
  expect_lt(gap(m2), 1e-10, label = "observation confirmation, R vs Stan")

  st <- fp_sites(m1, ndraws = 20, seed = 1)
  expect_s3_class(st, "fp_sites")
  expect_equal(nrow(st), 40 * 3)
  expect_identical(st$species, rep(c("sp1", "sp2", "sp3"), 40))
  expect_true(all(st$occupied[st$verified] > 1 - 1e-12))
  r <- fp_richness(m1, ndraws = 20, seed = 1)
  expect_equal(nrow(r), 40)
  # the same draws: expected richness is the sum of the occupancy of the species
  expect_equal(r$richness_psi, as.vector(tapply(st$psi, st$site, sum)))
  expect_true(all(r$richness_lower >= r$verified))
  expect_true(all(r$richness <= 3 & r$detected <= 3))
  expect_true(all(r$richness_upper == round(r$richness_upper)))
  expect_error(fp_richness(list()), "fp_occu_ms")

  pp <- fp_ppc(m1, ndraws = 20, seed = 1)
  expect_s3_class(pp, "fp_ppc")
  expect_equal(ncol(pp$yrep), length(pp$site))
  expect_length(pp$species, length(pp$site))
  expect_named(pp$p_species, c("sp1", "sp2", "sp3"))
  expect_output(print(pp), "By species")
  expect_s3_class(fp_ppc(m1, "verified", ndraws = 20, seed = 1), "fp_ppc")

  new <- simulate(m2, nsim = 2, truth = "posterior", seed = 1)
  expect_length(new, 2)
  expect_s3_class(new[[1]]$data, "fp_data")
  expect_identical(new[[1]]$data$species_names, s2$data$species_names)
  expect_equal(nrow(new[[1]]$data$auto), nrow(s2$data$auto))
  expect_equal(dim(new[[1]]$truth$z), c(40, 3))
  expect_true(all(c("occ_x[sp1]", "det_bio[sp3]", "fp_x[sp2]") %in%
                    names(new[[1]]$truth$pars)))
})

test_that("the random effects on p10 enter the likelihood as in fp_occu()", {
  skip_if_no_cmdstan()
  gq <- function(model, pars, stan_data) {
    out <- tempfile("gq")
    dir.create(out)
    g <- crywolf:::.compile_model(model)$generate_quantities(
      fitted_params = posterior::as_draws_array(pars), data = stan_data, seed = 1,
      sig_figs = 18, threads_per_chain = 1, output_dir = out)
    unclass(posterior::as_draws_matrix(g$draws("log_lik")))
  }
  first <- c("^mu_beta", "^sigma_beta", "^zb", "^alpha", "^sigma_alpha$", "^za",
             "^mu_gap$", "^sigma_gap$", "^zg", "^u_raw", "^sigma_u$", "^gamma_fp",
             "^delta0$", "^delta_s", "^sigma_v", "^v_raw", "^sigma_t", "^t_raw")
  ordered <- function(p) {
    p[, c(unlist(lapply(first, grep, names(p), value = TRUE)), ".chain", ".iteration", ".draw")]
  }
  s <- ms_sim_moved(3)
  m <- suppressWarnings(fp_occu(s$a, occ = ~ x, det = ~ bio, fp = ~ aci, random_fp = TRUE,
                                chains = 1, iter_warmup = 100, iter_sampling = 20,
                                threads_per_chain = 1, seed = 3, refresh = 0,
                                show_messages = FALSE, sig_figs = 18))
  pa <- posterior::as_draws_df(m$fit$draws())
  pa <- pa[, !grepl("^(lp__|u|gamma0|log_lik|p11_ref|p10_ref|p10_bar|r11)(\\[|$)", names(pa))]
  set.seed(4)
  pb <- pa
  for (v in c("beta[1]", "beta[2]", "alpha[1]")) pb[[v]] <- stats::rnorm(20)
  for (v in grep("^e_raw", names(pb), value = TRUE)) pb[[v]] <- stats::rnorm(20)
  pb$fp_gap <- stats::runif(20, 0.5, 3)
  sl <- function(d) crywolf:::.build_stan_data(d, ~ x, ~ bio, ~ aci, TRUE, fp_priors(),
                                               random_fp = TRUE)$stan_data
  ll <- cbind(gq("fp_occu", pa, sl(s$a)), gq("fp_occu", pb, sl(s$b)))
  # community values 0 (gap 1), standard deviations 1; the effect of unit
  # (site, species) is that of the site in the fit of the species
  pm <- pa[, !grepl("^(beta|fp_gap|e_raw|sigma_fp)", names(pa))]
  pm[c("mu_beta[1]", "mu_beta[2]", "sigma_beta[1]", "sigma_beta[2]")] <- list(0, 0, 1, 1)
  pm[c("zb[1,1]", "zb[2,1]", "zb[1,2]", "zb[2,2]")] <-
    list(pa[["beta[1]"]], pa[["beta[2]"]], pb[["beta[1]"]], pb[["beta[2]"]])
  pm[c("alpha[1]", "sigma_alpha")] <- list(0, 1)
  pm[c("za[1]", "za[2]")] <- list(pa[["alpha[1]"]], pb[["alpha[1]"]])
  pm[c("mu_gap", "sigma_gap")] <- list(1, 1)
  pm[c("zg[1]", "zg[2]")] <- list(log(pa$fp_gap), log(pb$fp_gap))
  pm[["sigma_t[1]"]] <- pa[["sigma_fp[1]"]]
  for (i in 1:40) {
    pm[[sprintf("t_raw[%d]", 2 * i - 1)]] <- pa[[sprintf("e_raw[%d]", i)]]
    # species b has its own effects (pb keeps the sd of pa)
    pm[[sprintf("t_raw[%d]", 2 * i)]] <- pb[[sprintf("e_raw[%d]", i)]]
  }
  sl_m <- crywolf:::.build_ms_stan_data(s$ms, ~ x, ~ bio, ~ aci, TRUE, fp_priors(),
                                        random_fp = TRUE)$stan_data
  ll_m <- gq("fp_occu_ms", ordered(pm), sl_m)
  by_species <- c(seq(1, by = 2, length.out = 40), seq(2, by = 2, length.out = 40))
  rel <- function(x, y) max(abs(x - y) / pmax(1, abs(x)))
  expect_lt(rel(ll, ll_m[, by_species]), 1e-12, label = "random_fp vs fp_occu")

  # an effect per group of sites equals the effect per unit of the same value
  ms <- s$ms
  ms$sites$grp <- rep(1:8, length.out = 40)
  pg <- pm[, !grepl("^(sigma_t|t_raw)", names(pm))]
  pg[["sigma_v[1]"]] <- 1
  vg <- matrix(stats::rnorm(20 * 16), 20)        # draws x (group - 1) * 2 + species
  for (c in 1:16) pg[[sprintf("v_raw[%d]", c)]] <- vg[, c]
  pt <- pg[, !grepl("^(sigma_v|v_raw)", names(pg))]
  pt[["sigma_t[1]"]] <- 1
  for (i in 1:40) for (sp in 1:2) {
    pt[[sprintf("t_raw[%d]", (i - 1) * 2 + sp)]] <- vg[, (ms$sites$grp[i] - 1) * 2 + sp]
  }
  build <- function(...) crywolf:::.build_ms_stan_data(ms, ~ x, ~ bio, ~ aci, TRUE,
                                                       fp_priors(), ...)$stan_data
  ll_g <- gq("fp_occu_ms", ordered(pg), build(fp_group = "grp"))
  ll_t <- gq("fp_occu_ms", ordered(pt), build(random_fp = TRUE))
  expect_lt(rel(ll_t, ll_g), 1e-12, label = "fp_group vs the same effects by unit")
  # centred: v_raw holds v; with standard deviations 1 and community values 0
  # (gap 1) only the recorder effects change
  pc <- pg
  for (u in grep("^u_raw", names(pc), value = TRUE)) pc[[u]] <- pc[[u]] * pc$sigma_u
  pc[["sigma_v[1]"]] <- 2
  for (c in 1:16) pc[[sprintf("v_raw[%d]", c)]] <- 2 * vg[, c]
  pg2 <- pg
  pg2[["sigma_v[1]"]] <- 2
  ll_c <- gq("fp_occu_ms", ordered(pc), build(fp_group = "grp", centered = TRUE))
  ll_n <- gq("fp_occu_ms", ordered(pg2), build(fp_group = "grp"))
  expect_lt(rel(ll_n, ll_c), 1e-12, label = "centred fp_group")
})

test_that("random_fp and fp_group are checked, named and simulated", {
  s <- ms_sim(1)
  expect_error(fp_occu_ms(s$ms, random_fp = 1), "TRUE or FALSE")
  expect_error(fp_occu_ms(s$ms, fp_group = "nope"), "column of the sites")
  expect_error(fp_occu_ms(s$ms, fp = NULL, random_fp = TRUE), "need a false-positive")
  one <- s$ms
  one$sites$g1 <- 1
  expect_error(fp_occu_ms(one, fp_group = "g1"), "at least two groups")
  expect_error(fp_occu_ms(s$ms, fp_group = "site", random_fp = TRUE), "are the same")
  sd <- crywolf:::.build_ms_stan_data(s$ms, ~ x, ~ 1, ~ 1, TRUE, fp_priors(),
                                      fp_group = "x", random_fp = TRUE)$stan_data
  expect_equal(c(sd$estimate_fp_group, sd$estimate_fp_re), c(1, 1))
  expect_equal(sd$ngroup, length(unique(s$ms$sites$x)))
  sim <- simulate_fp_ms(nspecies = 3, nsite = 40, nocc = 5, ngroup = 8, sd_fp_group = 1,
                        sd_fp_site = 0.5, seed = 1)
  expect_true(all(table(sim$data$sites$group) == 5))
  expect_true(all(c("sd_fp_group", "sd_fp_site") %in% names(sim$truth$pars)))
  expect_error(simulate_fp_ms(sd_fp_group = 1), "needs `ngroup`")
})

test_that("simulate_fp_ms(psi_rare) sets the occupancy of the first species", {
  a <- simulate_fp_ms(nspecies = 3, nsite = 30, nocc = 5, seed = 4)
  b <- simulate_fp_ms(nspecies = 3, nsite = 30, nocc = 5, psi_rare = 0.05, seed = 4)
  expect_equal(unname(b$truth$pars["occ_(Intercept)[sp1]"]), qlogis(0.05))
  expect_equal(b$truth$pars["occ_x[sp1]"], a$truth$pars["occ_x[sp1]"])
  expect_equal(b$truth$pars["occ_(Intercept)[sp2]"], a$truth$pars["occ_(Intercept)[sp2]"])
})
