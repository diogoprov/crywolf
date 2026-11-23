#' Review priority of each site, for a directed review round of sheepdog
#'
#' Gives every site (every site and species, for a multispecies model) a
#' size for `sheepdog::select_directed()`: the variance of its occupancy
#' state given the data, \eqn{q(1 - q)}, where \eqn{q} is the posterior
#' probability that the site is occupied (`occupied` of [fp_sites()]). The
#' size is largest, 0.25, where the model is least sure whether the species
#' is present, and 0 where the data settle it. sheepdog then draws units
#' with probability proportional to size and reviews one occasion of each.
#'
#' @param object An `fp_occu` or `fp_occu_ms` object.
#' @param species_code The label of the species in the predictions of
#'   sheepdog (its `species_code`, such as `"Nyctidromus albicollis_Common
#'   Pauraque"`). For an `fp_occu` model, a single label, required: without
#'   it, sheepdog would draw occasions of every species. For an `fp_occu_ms`
#'   model, optional: a character vector named by the species of the model,
#'   for when the model used other names (such as the short `species` of
#'   [fp_read_sheepdog()]); by default the species names of the model are
#'   taken as the labels.
#' @param validation Optional data frame of reviewed clips, such as
#'   `fp_read_sheepdog()$validation`, with the columns `site`, `verdict`
#'   (`TRUE` for a confirmed detection) and, for a multispecies model,
#'   `species_code`. Units with a confirmed detection are left out.
#' @inheritParams fp_sites
#'
#' @details
#' **Which units.** The units are the sites, or the sites and species, the
#' level of the occupancy state; sheepdog draws one occasion among the
#' candidate occasions of each drawn unit. Units where the species is
#' already known to be present are left out: those with a verified detection
#' in the data of the fit, and those with a confirmed clip in `validation`
#' (which may hold rounds reviewed after the fit). Without that, the floor of
#' sheepdog (`floor`) would still draw them.
#'
#' **Why the variance.** The size ranks units by how uncertain their state
#' is, not by the exact gain of one review, which would also depend on the
#' detection and false-positive rates of the occasion drawn. It is a simple
#' rule, and it needs nothing beyond [fp_sites()].
#'
#' **Fitting a model to the directed rounds.** The sizes depend only on data
#' already in the fit, and sheepdog adds a floor, draws units at random given
#' the sizes, and draws the occasion among those of the unit with a
#' prediction at or above its `threshold`, detections the fit also holds when
#' it uses the same threshold. A model fit to the data of all rounds (the earlier ones and the
#' directed one) therefore conditions on everything that chose the units,
#' and the selection is ignorable for that model: [fp_occu()] can use the
#' clips without weighting them. This holds only when the size of every
#' directed round came from a model of data that the new fit also holds.
#' Design-based estimates from the clips alone, such as the precision of the
#' classifier, still need the weights `1 / inclusion_prob` that sheepdog
#' records.
#'
#' @return A data frame with one row per unit, sorted by decreasing size:
#'   `site`, `species_code` and `size`, the columns that
#'   `sheepdog::select_directed()` expects in its `units`. Sites and labels
#'   are text.
#' @section Experimental:
#' This function is experimental (new after 2.3.0, part of the integration
#' with 'sheepdog'): its arguments and results may change in a minor version,
#' as `sheepdog::select_directed()` settles.
#' See [crywolf-api].
#' @seealso [fp_sites()], [fp_read_sheepdog()]
#' @export
#' @examplesIf instantiate::stan_cmdstan_exists() && requireNamespace("jsonlite", quietly = TRUE)
#' \donttest{
#' # the small export of sheepdog shipped with crywolf
#' dir <- system.file("extdata", "sheepdog-export", package = "crywolf")
#' x <- suppressWarnings(fp_read_sheepdog(dir, species = "phy_cuv"))
#' auto <- fp_aggregate(x$occasions, x$detections, threshold = 0.5,
#'                      by = c("site", "occasion_id"), validation = x$validation)
#' m <- fp_occu(fp_data(auto, x$sites), random_recorder = FALSE,
#'              chains = 2, threads_per_chain = 1)
#' units <- fp_review_priority(m, species_code = unique(x$detections$species_code),
#'                             validation = x$validation)
#' head(units)
#' # then, in sheepdog, with the predictions of the survey:
#' # sel <- sheepdog::select_directed(pred, units, n_units = 10, threshold = 0.5)
#' }
fp_review_priority <- function(object, species_code = NULL, validation = NULL,
                               ndraws = NULL, seed = NULL) {
  if (!inherits(object, c("fp_occu", "fp_occu_ms"))) {
    stop("`object` must be an fp_occu or fp_occu_ms model.", call. = FALSE)
  }
  .check_review_args(object$data, species_code, validation)
  st <- fp_sites(object, ndraws = ndraws, seed = seed)
  .review_units(st, object$data, species_code, validation)
}

#' checks of fp_review_priority(), before the costly fp_sites()
#' @noRd
.check_review_args <- function(data, species_code, validation) {
  ms <- !is.null(data$species)
  if (!ms && (!is.character(species_code) || length(species_code) != 1 ||
              is.na(species_code))) {
    stop("`species_code` must be the label of the species in the predictions ",
         "of sheepdog, one string.", call. = FALSE)
  }
  if (ms && !is.null(species_code) &&
      (!is.character(species_code) || is.null(names(species_code)) ||
       !all(data$species_names %in% names(species_code)))) {
    stop("`species_code` must be a character vector named by the species of ",
         "the model: ", paste(data$species_names, collapse = ", "), ".",
         call. = FALSE)
  }
  if (is.null(validation)) return(invisible())
  miss <- setdiff(c("site", "verdict", if (ms) "species_code"), names(validation))
  if (length(miss)) {
    stop("`validation` lacks the column(s) ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  if (!ms && length(unique(validation$species_code)) > 1) {
    stop("`validation` holds clips of several species; keep those of the ",
         "species of the model (fp_read_sheepdog(species = )).", call. = FALSE)
  }
  labels <- if (is.null(species_code)) data$species_names else species_code
  if (ms && nrow(validation) && !any(validation$species_code %in% labels)) {
    warning("No `species_code` of `validation` is a species label of the units ",
            "(", paste(utils::head(labels, 3), collapse = ", "), if (length(labels) > 3)
              ", ...", "), so no unit is left out as confirmed; map the species ",
            "of the model to their labels with `species_code`.", call. = FALSE)
  }
  invisible()
}

#' the units of fp_review_priority() from the table of fp_sites()
#' @noRd
.review_units <- function(st, data, species_code, validation) {
  ms <- !is.null(data$species)
  label <- if (!ms) {
    rep(species_code, nrow(st))
  } else if (is.null(species_code)) {
    as.character(st[[data$species]])
  } else {
    unname(species_code[as.character(st[[data$species]])])
  }
  q <- st$occupied
  units <- data.frame(site = as.character(st[[data$site]]), species_code = label,
                      size = q * (1 - q), stringsAsFactors = FALSE)
  known <- st$verified
  if (!is.null(validation)) {
    hit <- validation[validation$verdict %in% TRUE, , drop = FALSE]
    key <- function(site, sp) paste(site, if (ms) sp, sep = "\r")
    known <- known | key(units$site, units$species_code) %in%
      key(as.character(hit$site), as.character(hit$species_code))
  }
  units <- units[!known, , drop = FALSE]
  units <- units[order(-units$size), , drop = FALSE]
  rownames(units) <- NULL
  units
}
