# Checks that crywolf::simulate_fp(seed = s) reproduces the datasets of the
# fpcalib compendium (R/fp_sim_pilot.R, v0.1.0, doi:10.5281/zenodo.23379100).
# Usage: Rscript dev/compare-simulate-fp.R <path to fpcalib>
devtools::load_all(quiet = TRUE)
src <- file.path(commandArgs(trailingOnly = TRUE)[1], "R", "fp_sim_pilot.R")
env <- new.env()
exprs <- parse(src)
is_def <- purrr::map_lgl(as.list(exprs), function(e) {
  is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("simulate_fp"))
})
eval(exprs[[which(is_def)]], env)
comp_sim <- env$simulate_fp

compare <- function(args, seed) {
  a <- do.call(comp_sim, c(args, list(seed = seed)))
  map_args <- c(nocc = "K", frac_verified = "verified_coverage", nver = "n_verified",
                sd_recorder = "sig_u", nrec = "nobser", p10_conc = "p10_beta_conc")
  pargs <- args
  for (new in names(map_args)) {
    old <- map_args[[new]]
    if (!is.null(pargs[[old]])) { pargs[[new]] <- pargs[[old]]; pargs[[old]] <- NULL }
  }
  b <- do.call(simulate_fp, c(pargs, list(seed = seed)))
  sd <- a$stan_data
  ver <- b$data$verified
  c(
    y_auto = identical(as.integer(sd$y2), as.integer(b$data$auto$y)),
    site_auto = identical(as.integer(sd$site2), as.integer(b$data$auto$site)),
    bio = isTRUE(all.equal(sd$bio2, b$data$auto$bio)),
    y_ver = identical(sd$y1, as.integer(ver$y)),
    site_ver = identical(sd$site1, as.integer(ver$site)),
    x = isTRUE(all.equal(a$umf$siteCovs$xs, b$data$sites$x)),
    w = isTRUE(all.equal(a$w, b$data$sites$w)),
    p10_marginal = isTRUE(all.equal(a$truth$p10_marginal, b$truth$p10_marginal))
  )
}

regimes <- list(
  homogeneous = list(nsite = 300, K = 20),
  unstructured_beta = list(nsite = 300, K = 20, p10_beta_conc = 2),
  unstructured_logit = list(nsite = 300, K = 20, sigma_p10 = 1),
  structured = list(nsite = 300, K = 20, p10 = 0.2, conf_effect = 1.5, conf_cor = 0.7,
                    verified_coverage = 0.9),
  covariates_re = list(nsite = 200, K = 10, a_bio = 0.5, a_aci = -0.3, sig_u = 0.5)
)
res <- purrr::imap(regimes, function(args, nm) {
  ok <- purrr::map(1:20, function(s) compare(args, s)) |> purrr::reduce(`&`)
  c(regime = nm, ok)
})
out <- do.call(rbind, res)
print(out, quote = FALSE)
cat("\nAll identical:", all(out[, -1] == "TRUE"), "\n")
