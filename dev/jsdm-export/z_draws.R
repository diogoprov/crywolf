# Posterior draws of the occupancy states of an fp_occu_ms() fit, as
# fp_richness() draws them internally: for each posterior draw, z[i, s] ~
# Bernoulli(P(z = 1 | data, draw)), with z = 1 where a verified detection
# proves presence. Returns an array [draw, site, species]. Needs the package
# loaded with devtools::load_all() (it uses internal helpers).
z_draws <- function(object, ndraws = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- .ms_species(object, ndraws)
  out <- NULL
  for (s in seq_along(sp)) {
    x <- sp[[s]]
    em <- .ppc_core(x$stan_data, x$pars, "automated", "freeman-tukey",
                    emissions_only = TRUE)
    lpsi <- x$pars$beta %*% t(x$stan_data$X)
    pz <- stats::plogis(lpsi + em$occ - em$unocc)
    pz[, x$stan_data$known_present == 1] <- 1
    if (is.null(out)) {
      out <- array(NA_integer_, c(nrow(pz), ncol(pz), length(sp)),
                   dimnames = list(NULL, object$data$sites[[object$data$site]],
                                   object$data$species_names))
    }
    out[, , s] <- matrix(stats::rbinom(length(pz), 1, pz), nrow(pz))
  }
  out
}
