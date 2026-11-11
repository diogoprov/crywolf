# Outputs of the multispecies model ------------------------------------------
#
# Given its parameters, each species of fp_occu_ms() is an fp_occu() model of
# its own rows (the likelihood test of tests/testthat/test-fp_occu_ms.R), so
# fp_sites(), fp_ppc(), simulate() and fp_richness() run the code of
# fp_occu() species by species, with the values of that species in each
# posterior draw. The draws are the same for every species, so that
# quantities summed over species (richness) keep their joint posterior.

#' Species richness of each site from a multispecies model
#'
#' Counts, for each site, the species of an [fp_occu_ms()] model that are
#' present, from the posterior of the occupancy state of every species at
#' the site. **In development**, as [fp_occu_ms()].
#'
#' For each posterior draw, the state of each species at each site is drawn
#' from its conditional posterior given the data of that species at that site
#' (as in [fp_sites()]), and the richness of the site is the number of
#' species drawn as present. The species are those of the model: a species
#' never detected anywhere is not in the list and is not counted.
#'
#' @param object An `fp_occu_ms` object.
#' @inheritParams fp_sites
#' @return A data frame, one row per site, with the site columns of the data,
#'   `richness` (posterior mean number of species present, with
#'   `richness_lower` and `richness_upper`, quantiles of its posterior
#'   distribution), `richness_psi` (expected number of species from the
#'   occupancy covariates alone, before the data of the site), `detected`
#'   (number of species with automated detections) and `verified` (number of
#'   species whose presence a verified detection proves).
#' @seealso [fp_occu_ms()], [fp_sites()]
#' @export
#' @examples
#' \dontrun{
#' s <- simulate_fp_ms(nspecies = 5, nsite = 60, seed = 1)
#' m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE)
#' r <- fp_richness(m)
#' plot(rowSums(s$truth$z), r$richness)
#' }
fp_richness <- function(object, prob = 0.95, ndraws = NULL, seed = NULL) {
  if (!inherits(object, "fp_occu_ms")) {
    stop("`object` must be an `fp_occu_ms` fit.", call. = FALSE)
  }
  if (!is.null(seed)) .set_seed(seed)
  sp <- .ms_species(object, ndraws)
  z <- psi <- 0
  for (x in sp) {
    r <- .ppc_core(x$stan_data, x$pars, "automated", "freeman-tukey")
    z <- z + r$z
    psi <- psi + stats::plogis(x$pars$beta %*% t(x$stan_data$X))
  }
  sites <- object$data$sites
  new <- c("richness", "richness_lower", "richness_upper", "richness_psi",
           "detected", "verified")
  clash <- intersect(new, names(sites))
  if (length(clash)) {
    stop("The sites table already has column(s) ", paste(clash, collapse = ", "),
         "; rename them before calling fp_richness().", call. = FALSE)
  }
  detected <- purrr::map(sp, function(x) {
    tabulate(rep(x$stan_data$site2, x$stan_data$y2), nbins = x$stan_data$nsite) > 0
  })
  a <- (1 - prob) / 2
  # quantiles that are counts of species, not interpolated between them
  count_q <- function(p) apply(z, 2, stats::quantile, probs = p, type = 1, names = FALSE)
  data.frame(sites,
             richness = colMeans(z), richness_lower = count_q(a),
             richness_upper = count_q(1 - a), richness_psi = colMeans(psi),
             detected = Reduce(`+`, detected),
             verified = rowSums(as.matrix(object$data$known_present) == 1),
             row.names = NULL)
}

# fp_sites() for the multispecies model: one row per site and species
.fp_sites_ms <- function(object, prob, ndraws, seed) {
  if (!is.null(seed)) .set_seed(seed)
  data <- object$data
  sp <- .ms_species(object, ndraws)
  if (data$species %in% names(data$sites)) {
    stop("The sites table already has a column `", data$species, "`, the species ",
         "column; rename it before calling fp_sites().", call. = FALSE)
  }
  out <- purrr::imap(sp, function(x, s) {
    r <- .ppc_core(x$stan_data, x$pars, "automated", "freeman-tukey")
    sites <- data$sites
    sites[[data$species]] <- data$species_names[s]
    .sites_table(sites, x$stan_data, x$pars, r$pz, prob)
  })
  out <- do.call(rbind, out)
  # site by site, the species in the order of the model
  out <- out[order(match(out[[data$site]], data$sites[[data$site]])), , drop = FALSE]
  rownames(out) <- NULL
  structure(out, class = c("fp_sites", "data.frame"), prob = prob)
}

# fp_ppc() for the multispecies model: the columns are site x species, the
# discrepancy is summed over species, and p_species holds the p-value of each
# species
.fp_ppc_ms <- function(object, method, fit_stat, ndraws, seed) {
  data <- object$data
  if (method == "verified" && object$stan_data$N1 == 0) {
    stop("The data have no verified occasions.", call. = FALSE)
  }
  if (!is.null(seed)) .set_seed(seed)
  sp <- .ms_species(object, ndraws)
  site_ids <- data$sites[[data$site]]
  res <- purrr::imap(sp, function(x, s) {
    if (method == "verified" && x$stan_data$N1 == 0) return(NULL)
    r <- .ppc_core(x$stan_data, x$pars, method, fit_stat)
    r$site <- site_ids[r$keep]
    r$species <- rep(data$species_names[s], length(r$keep))
    r$p <- mean(r$fit_rep > r$fit_obs)
    r
  })
  names(res) <- data$species_names
  res <- purrr::compact(res)
  fit_obs <- Reduce(`+`, purrr::map(res, "fit_obs"))
  fit_rep <- Reduce(`+`, purrr::map(res, "fit_rep"))
  structure(
    list(y = unlist(purrr::map(res, "y"), use.names = FALSE),
         yrep = do.call(cbind, purrr::map(res, "yrep")),
         fit_obs = fit_obs, fit_rep = fit_rep, p_value = mean(fit_rep > fit_obs),
         site = unlist(purrr::map(res, "site"), use.names = FALSE),
         species = unlist(purrr::map(res, "species"), use.names = FALSE),
         p_species = purrr::map_dbl(res, "p"),
         method = method, fit_stat = fit_stat),
    class = "fp_ppc"
  )
}

#' Simulate data from a fitted multispecies model
#'
#' Generates new detections from an [fp_occu_ms()] fit, keeping the
#' structure of the original data (sites, occasions and their covariates,
#' recorders, which occasions were verified or checked), as
#' [simulate.fp_occu()] does for one species. **In development**, as
#' [fp_occu_ms()].
#'
#' @param object An `fp_occu_ms` object.
#' @inheritParams simulate.fp_occu
#' @return A list of `nsim` elements, each a list with `data` (an `fp_data`
#'   object with the structure and the species column of the original data)
#'   and `truth`: `pars`, the generating values of each species named as in
#'   [simulate.fp_occu()] with the species in brackets (`occ_x[species]`,
#'   `det_(Intercept)[species]`, `fp_(Intercept)[species]`, ...); `z`, the
#'   latent states, a matrix of sites by species; `p10_site`, the
#'   false-positive rate of each site and species, a matrix of the same
#'   shape; and `p10_bar`, its mean over sites for each species. All the
#'   species of a dataset take their values from the same posterior draw
#'   (or from the posterior means).
#' @seealso [fp_occu_ms()], [simulate_fp_ms()]
#' @export
#' @examples
#' \dontrun{
#' s <- simulate_fp_ms(nspecies = 4, nsite = 60, seed = 1)
#' m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE)
#' new <- simulate(m, nsim = 2, seed = 1)
#' new[[1]]$truth$pars
#' }
simulate.fp_occu_ms <- function(object, nsim = 1, seed = NULL,
                                truth = c("mean", "posterior"), ...) {
  truth <- match.arg(truth)
  data <- object$data
  gens <- .ms_species_gens(object)
  dm <- .ms_draws(object, NULL)
  if (!is.null(seed)) .set_seed(seed)
  purrr::map(seq_len(nsim), function(i) {
    # one posterior draw (or the means) for all the species
    row <- if (truth == "mean") matrix(colMeans(dm), 1, dimnames = list(NULL, colnames(dm)))
           else dm[sample.int(nrow(dm), 1), , drop = FALSE]
    res <- purrr::imap(gens, function(g, s) {
      .simulate_from(g, .pick_pars(.ms_pars(row, object, g, s), "mean"))
    })
    names(res) <- data$species_names
    .ms_combine_sims(res, data)
  })
}

# the simulated data of each species (from .simulate_from()) as one
# multispecies data set and its truth
.ms_combine_sims <- function(res, data) {
  sp <- names(res)
  rows <- function(what) {
    tabs <- purrr::compact(purrr::map(res, function(r) r$data[[what]]))
    if (!length(tabs)) return(NULL)
    out <- do.call(rbind, tabs)
    rownames(out) <- NULL
    out
  }
  new <- withCallingHandlers(
    fp_data(auto = rows("auto"), sites = data$sites, verified = rows("verified"),
            site = data$site, y = data$y, recorder = data$recorder,
            confirmed = data$confirmed, species = data$species),
    crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  by_sp <- function(what) {
    m <- do.call(cbind, purrr::map(res, function(r) r$truth[[what]]))
    if (!is.null(m)) colnames(m) <- sp
    m
  }
  pars <- unlist(unname(purrr::imap(res, function(r, s) {
    stats::setNames(r$truth$pars, sprintf("%s[%s]", names(r$truth$pars), s))
  })))
  list(data = new,
       truth = list(pars = pars, z = by_sp("z"), p10_site = by_sp("p10_site"),
                    p10_bar = purrr::map_dbl(res, function(r) r$truth$p10_bar)))
}

# internal ---------------------------------------------------------------------

# the data of species s alone, as an fp_data of fp_occu(): its rows of the
# automated and verified tables (with the species column, kept for
# simulate()) and its column of known presences
.ms_species_data <- function(data, s) {
  sp <- data$species_names[s]
  keep <- function(df) {
    if (is.null(df)) return(NULL)
    df <- df[as.character(df[[data$species]]) == sp, , drop = FALSE]
    rownames(df) <- NULL
    df
  }
  d <- data
  d$auto <- keep(data$auto)
  d$verified <- keep(data$verified)
  d$known_present <- unname(as.matrix(data$known_present)[, s])
  d["species"] <- list(NULL)
  d["species_names"] <- list(NULL)
  d
}

# for each species, what .simulate_from() and .ppc_core() need: its data,
# the stan_data of fp_occu() for it, its design and occasions
.ms_species_gens <- function(object) {
  f <- object$formulas
  purrr::map(seq_along(object$data$species_names), function(s) {
    d <- .ms_species_data(object$data, s)
    no_ver <- is.null(d$verified) || nrow(d$verified) == 0
    sl <- .build_stan_data(d, f$occ, f$det, f$fp, object$random_recorder, object$priors,
                           ver = if (no_ver) ~ 1 else f$ver)
    list(data = d, formulas = f, stan_data = sl$stan_data, design = sl$design,
         occasions = sl$occasions)
  })
}

# the posterior draws of fp_occu_ms() that make the parameters of each
# species, a matrix with one row per draw (ndraws of them, if given)
.ms_draws <- function(object, ndraws) {
  sd <- object$stan_data
  fl <- .ms_flags(sd)
  vars <- c("beta_sp", "alpha", "alpha_sp", "delta0",
            if (fl$fp) "gap_sp", if (sd$estimate_re == 1) "u",
            if (fl$fp && sd$Kfp > 0) "gamma_fp", if (fl$det) "alpha_sp_s",
            if (fl$fp_s) "gamma_sp", if (fl$ver) "delta_sp",
            if (isTRUE(sd$Kver_s > 0)) "delta_s")
  dm <- unclass(posterior::as_draws_matrix(object$fit$draws(variables = vars)))
  if (!is.null(ndraws) && ndraws < nrow(dm)) {
    dm <- dm[sort(sample.int(nrow(dm), ndraws)), , drop = FALSE]
  }
  dm
}

# which terms vary by species and have values to vary
.ms_flags <- function(sd) {
  fp <- sd$estimate_fp == 1
  list(fp = fp,
       det = isTRUE(sd$by_det == 1) && sd$Kalpha > 1,
       fp_s = isTRUE(sd$by_fp == 1) && fp && sd$Kfp > 0,
       ver = isTRUE(sd$by_ver == 1))
}

# the parameters of species s in the format of .ppc_pars(), from the draws
# dm of fp_occu_ms(); gen is the generator of the species (.ms_species_gens()),
# whose recorders and verified terms may be a subset of those of the model
.ms_pars <- function(dm, object, gen, s) {
  sd <- object$stan_data
  des <- object$design
  fl <- .ms_flags(sd)
  S <- nrow(dm)
  cols <- function(v) dm[, grep(paste0("^", v, "(\\[|$)"), colnames(dm)), drop = FALSE]
  by_sp <- function(v, k) dm[, sprintf("%s[%d,%d]", v, k, s), drop = FALSE]
  alpha <- cols("alpha")
  alpha[, 1] <- dm[, sprintf("alpha_sp[%d]", s)]
  if (fl$det) alpha[, -1] <- by_sp("alpha_sp_s", seq_len(sd$Kalpha - 1))
  matched <- function(v, own, all) {
    j <- match(own, all)
    if (anyNA(j)) {
      stop("Species ", des$species[s], " has levels of `", v, "` the model does not have.",
           call. = FALSE)
    }
    cols(v)[, j, drop = FALSE]
  }
  gs <- gen$stan_data
  if (gs$Kalpha != sd$Kalpha || gs$Kfp != sd$Kfp) {
    stop("The rows of species ", des$species[s], " do not have every level of the ",
         "covariates of `det` or `fp`; give every species every occasion.", call. = FALSE)
  }
  list(
    beta = by_sp("beta_sp", seq_len(sd$Kbeta)),
    alpha = alpha,
    delta0 = if (fl$ver) dm[, sprintf("delta_sp[%d]", s)] else dm[, "delta0"],
    delta_s = if (isTRUE(gs$Kver_s > 0)) {
      matched("delta_s", gen$design$ver_names[-1], des$ver_names[-1])
    } else matrix(0, S, 0),
    u = if (sd$estimate_re == 1) matched("u", gen$design$recorders, des$recorders),
    gamma0 = if (fl$fp) dm[, sprintf("alpha_sp[%d]", s)] - dm[, sprintf("gap_sp[%d]", s)],
    gamma_fp = if (fl$fp_s) by_sp("gamma_sp", seq_len(sd$Kfp)) else
      if (fl$fp && sd$Kfp > 0) cols("gamma_fp") else matrix(0, S, 0)
  )
}

# each species with its generator and its parameters in the same draws
.ms_species <- function(object, ndraws) {
  dm <- .ms_draws(object, ndraws)
  purrr::imap(.ms_species_gens(object), function(g, s) {
    g$pars <- .ms_pars(dm, object, g, s)
    g
  })
}
