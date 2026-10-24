# Checks that fp_occu.stan defines the same posterior as the Stan programs of
# the fpcalib compendium (GitHub diogoprov/fpcalib, Zenodo
# doi:10.5281/zenodo.23379100): fp_occupancy_re_toggle.stan (main run) and
# fp_occupancy_fpcov.stan (covariate-dependent false-positive arms).
# Both simulators return the same datasets for the same seed
# (dev/compare-simulate-fp.R); the crywolf data go through .build_stan_data(),
# as in fp_occu(). The log densities (with Jacobian) are compared at random
# parameter values; the same posterior means a difference that does not vary
# with the parameters.
# Usage, from the package root: Rscript dev/compare-logdens-fpcalib.R <path to fpcalib>
# Writes dev/compare-logdens-fpcalib.txt.
fpc <- commandArgs(trailingOnly = TRUE)[1]
suppressMessages({
  library(cmdstanr)
  devtools::load_all(quiet = TRUE)
})

# simulador do fpcalib, isolado num ambiente
env <- new.env()
exprs <- parse(file.path(fpc, "R", "fp_sim_pilot.R"))
for (e in as.list(exprs)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("simulate_fp", "clampl")) eval(e, env)
}

# copias dos .stan, compiladas com os metodos do modelo (log_prob), para nao
# tocar nos executaveis do pacote e do compendio
sdir <- file.path(tempdir(), "compare-logdens")
dir.create(sdir, showWarnings = FALSE)
file.copy(c(file.path(fpc, "inst", "stan", c("fp_occupancy_re_toggle.stan",
                                              "fp_occupancy_fpcov.stan")),
            file.path("inst", "stan", "fp_occu.stan")), sdir, overwrite = TRUE)
compile <- function(f) cmdstan_model(file.path(sdir, f), compile_model_methods = TRUE,
                                     quiet = TRUE)
mod_const <- compile("fp_occupancy_re_toggle.stan")
mod_cov <- compile("fp_occupancy_fpcov.stan")
mod_cw <- compile("fp_occu.stan")

# um ajuste curtissimo so para expor log_prob() e unconstrain_variables()
methods_fit <- function(mod, data) {
  f <- suppressWarnings(suppressMessages(
    mod$sample(data = data, chains = 1, iter_warmup = 20, iter_sampling = 1,
               refresh = 0, show_messages = FALSE, show_exceptions = FALSE, seed = 1)))
  f$init_model_methods()
  f
}

map_args <- c(nocc = "K", frac_verified = "verified_coverage", nver = "n_verified",
              sd_recorder = "sig_u", nrec = "nobser", p10_conc = "p10_beta_conc")
to_cw <- function(a) {
  for (new in names(map_args)) {
    old <- map_args[[new]]
    if (!is.null(a[[old]])) { a[[new]] <- a[[old]]; a[[old]] <- NULL }
  }
  a
}

compare <- function(label, sim_args, fp_cov = NULL, seed = 2026, n_theta = 25) {
  A <- do.call(env$simulate_fp, c(sim_args, list(seed = seed)))
  B <- do.call(simulate_fp, c(to_cw(sim_args), list(seed = seed)))
  re <- A$stan_data$estimate_re
  sdA <- A$stan_data
  if (is.null(fp_cov)) {
    modA <- mod_const
    fp <- ~ 1
  } else {
    modA <- mod_cov
    W <- if (fp_cov == "x") A$umf$siteCovs$xs else A$w
    sdA$Kfp <- 1L
    sdA$Xfp <- matrix(W, ncol = 1)
    fp <- stats::as.formula(paste("~", fp_cov))
  }
  sdB <- .build_stan_data(B$data, occ = ~ x, det = ~ bio + aci, fp = fp,
                          random_recorder = re == 1, priors = fp_priors())$stan_data
  # sem efeito de gravador o crywolf usa 1 gravador e o fpcalib nobser, todos
  # so com a priori; os excedentes do fpcalib ficam em 0
  stopifnot(sdA$nobser == sdB$nrec || (re == 0 && sdB$nrec == 1))
  fA <- methods_fit(modA, sdA)
  fB <- methods_fit(mod_cw, sdB)
  set.seed(1)
  d <- replicate(n_theta, {
    beta <- rnorm(2)
    alpha <- rnorm(3, 0, 0.7)
    u_raw <- rnorm(sdB$nrec)
    sigma_u <- runif(1, 0.1, 2)
    fp_gap <- runif(1, 0.3, 4)
    delta0 <- rnorm(1)
    g <- array(rnorm(1, 0, 0.8), 1)
    thA <- list(beta = beta, alpha0 = alpha[1], a_bio = alpha[2], a_aci = alpha[3],
                u_raw = c(u_raw, rep(0, sdA$nobser - sdB$nrec)), sigma_u = sigma_u,
                fp_gap = fp_gap, delta0 = delta0)
    if (!is.null(fp_cov)) thA$gamma_fp <- g
    thB <- list(beta = beta, alpha = alpha, u_raw = u_raw, sigma_u = sigma_u,
                fp_gap = fp_gap, gamma_fp = if (is.null(fp_cov)) numeric(0) else g,
                delta0 = delta0, delta_s = numeric(0), e_raw = numeric(0),
                sigma_fp = numeric(0), d_raw = numeric(0), sigma_det = numeric(0))
    c(fA$log_prob(fA$unconstrain_variables(thA), jacobian = TRUE),
      fB$log_prob(fB$unconstrain_variables(thB), jacobian = TRUE))
  })
  diff <- d[1, ] - d[2, ]
  data.frame(design = label, estimate_re = re, Kfp = sdB$Kfp,
             logdens_min = min(d[1, ]), logdens_max = max(d[1, ]),
             mean_diff = mean(diff), max_dev_diff = max(abs(diff - mean(diff))))
}

base <- list(nsite = 300, K = 20, verified_coverage = 0.3, p10 = 0.05)
conf <- modifyList(base, list(conf_effect = 1.5, conf_cor = 0.7))
res <- rbind(
  compare("homogeneous", base),
  compare("homogeneous, p10 = 0.20, 10% verified",
          modifyList(base, list(p10 = 0.20, verified_coverage = 0.10))),
  compare("unstructured (Beta)", modifyList(base, list(p10_beta_conc = 2))),
  compare("structured, constant p10", conf),
  compare("structured, fp ~ x (practical fix)", conf, fp_cov = "x"),
  compare("structured, fp ~ w (positive control)", conf, fp_cov = "w"),
  compare("recorder effect (sig_u = 0.8)", modifyList(base, list(sig_u = 0.8)))
)

out <- c(
  paste("crywolf", as.character(utils::packageVersion("crywolf")),
        "commit", system("git rev-parse --short HEAD", intern = TRUE)),
  paste("fpcalib commit", system(paste("git -C", shQuote(fpc), "rev-parse --short HEAD"),
                                 intern = TRUE)),
  paste("CmdStan", cmdstan_version(), "| cmdstanr", as.character(utils::packageVersion("cmdstanr")),
        "|", R.version.string),
  paste("Run on", format(Sys.Date())),
  "",
  utils::capture.output(print(res, digits = 4, row.names = FALSE))
)
writeLines(out, file.path("dev", "compare-logdens-fpcalib.txt"))
writeLines(out)
