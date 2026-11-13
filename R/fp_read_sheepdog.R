#' Read an export of the sheepdog package
#'
#' Reads the folder written by `sheepdog::export_occupancy_data()`, which
#' follows the `sheepdog-occupancy` schema (version 1.x): every occasion of a
#' passive acoustic survey, the classifier predictions, the reviewed clips
#' with their sampling design, and the recordings annotated by ear. The
#' tables come back in the forms that [fp_aggregate()] and [fp_data()]
#' expect.
#'
#' @param path Folder of the export, holding `schema.json`.
#' @param species Optional species to keep, matched against `species_code`
#'   (the classifier label, such as `"Nyctidromus albicollis_Common
#'   Pauraque"`) and `species` (the short code that sheepdog derives from
#'   it). Without it, an export of several species is kept whole, and the
#'   verdicts and annotated recordings of all species are mixed: name the
#'   target species. For the multispecies model, name the species to keep
#'   (all of them: `species = schema$species`, from `schema.json`); every
#'   table then keeps its species columns, and [fp_draws_data()] builds the
#'   data of all of them.
#' @param keep_top Keep reviewed clips drawn with `sampling_mode = "top"`?
#'   See Details. Clips of `"score_bands"` are dropped either way.
#'
#' @details
#' **Which reviewed clips.** The false-positive models of crywolf treat the
#' checked detections as chosen independently of the classifier's score.
#' Clips of the `"random_detections"` design are; clips of the `"top"`
#' design, the highest-scoring detections, are not: they hold fewer false
#' positives than a random sample and would make the false-positive rate
#' look lower than it is. By default they are dropped with a warning
#' (`keep_top = FALSE`); they remain useful to calibrate the classifier
#' within sheepdog. Clips of the `"directed_pps"` design (sheepdog's
#' `select_directed()`) are kept with a warning: within each review unit
#' (for example site x species x month) the occasion is drawn independently
#' of the score, but the units are drawn with unequal probabilities, whose
#' size may come from the scores, and crywolf reads `inclusion_prob` without
#' weighting by it. Estimates from them are conditional on the units drawn;
#' weighting comes with a model of the scores.
#'
#' Clips of the `"score_bands"` design are always dropped, whatever
#' `keep_top`, with a warning. sheepdog draws them at random within bands of
#' score, so that thresholds and calibration curves see the whole score
#' range; being chosen by their score, they over-represent some bands and
#' would bias the false-positive rate in either direction. They serve
#' `sheepdog::threshold_bootstrap()`, `sheepdog::calibrate()` and
#' `sheepdog::fit_calibration()`, not these models. When every reviewed clip
#' is of this design and the export has no annotated recordings, there is
#' nothing left to estimate false positives from, and `fp_read_sheepdog()`
#' stops: review clips drawn with `"random_detections"` or annotate
#' recordings drawn with `"random_recordings"`.
#'
#' **Clips and recordings.** A positive clip proves that the species called
#' in its occasion; a negative clip says only that this prediction was
#' false. These go to `validation`, for the observation-confirmation design
#' through `fp_aggregate(validation = )`. A recording annotated by ear
#' (`"random_recordings"`) gives the state of the whole occasion; these go
#' to `verified`, the verified method of [fp_data()].
#'
#' **Censored scores.** Only predictions with a confidence of at least
#' `schema$min_confidence` were exported, so thresholds below it count
#' every exported prediction and no more.
#'
#' Occasions are identified by `site` and `occasion_id`; pass
#' `by = c("site", "occasion_id")` to [fp_aggregate()]. Identifiers, dates
#' and the keys of the calibration (such as `month`) are read as text. The
#' occasion covariates of the export (`schema$occasion_covariates`, from
#' schema version 1.1) are extra columns of `occasions`. The calibration draws
#' (`calibration_draws.parquet`) and the detection histories sheepdog
#' drew from them (`occasion_draws.parquet`, schema 1.2), when present, are
#' not read here; their
#' path is in `files`. When the calibration has continuous covariates
#' (`schema$calibration$covariates`, schema 1.3), the export holds only the
#' histories sheepdog drew, and no calibration draws.
#'
#' @return A list of class `fp_sheepdog` with the data frames `occasions`
#'   (`site`, `occasion_id`, `source_audio`, `recording_date`,
#'   `recording_hour`, `night`, `start`, `duration`, then any occasion
#'   covariates), `detections` (one row
#'   per prediction: `site`, `occasion_id`, `species_code`, `species`,
#'   `start`, `confidence`), `scores` (highest score and number of
#'   predictions per occasion and species), `validation` (one row per
#'   reviewed clip, with `start`, `end`, `confidence`, `sampling_mode`,
#'   `inclusion_prob`, `verdict` as `TRUE`, `FALSE` or `NA`), `verified`
#'   (one row per annotated recording, with `y`, 1 or 0, and
#'   `inclusion_prob`), `sites`, the parsed `schema`, and `files` (paths of
#'   the tables found). Tables absent from the export are `NULL`.
#' @seealso [fp_aggregate()], [fp_data()], [fp_read_validation()]
#' @export
#' @examples
#' \dontrun{
#' # the tables hold only the target species, so fp_aggregate() needs no
#' # `species`
#' x <- fp_read_sheepdog("occupancy_export", species = "compau")
#' auto <- fp_aggregate(x$occasions, x$detections, threshold = 0.5,
#'                      by = c("site", "occasion_id"), validation = x$validation)
#' d <- fp_data(auto, x$sites, x$verified, confirmed = c("n_true", "n_false"))
#' }
fp_read_sheepdog <- function(path, species = NULL, keep_top = FALSE) {
  .require_pkg("jsonlite")
  schema_file <- file.path(path, "schema.json")
  if (!file.exists(schema_file)) {
    stop("No schema.json in `path`; is it a folder written by ",
         "sheepdog::export_occupancy_data()?", call. = FALSE)
  }
  schema <- jsonlite::read_json(schema_file, simplifyVector = TRUE)
  if (!identical(schema$schema, "sheepdog-occupancy")) {
    stop("`path` holds a '", schema$schema, "' export; only ",
         "'sheepdog-occupancy' is supported.", call. = FALSE)
  }
  if (!identical(as.integer(schema$schema_version), 1L)) {
    stop("Version ", schema$schema_version, " of the sheepdog-occupancy ",
         "schema is not supported (only versions 1.x); update crywolf.",
         call. = FALSE)
  }
  tables <- c("occasions", "detections", "scores", "validation",
              "recordings_verified", "sites")
  files <- stats::setNames(file.path(path, paste0(tables, ".csv")), tables)
  files <- files[file.exists(files)]
  if (!"occasions" %in% names(files)) {
    stop("The export has no occasions.csv.", call. = FALSE)
  }
  read <- function(name, text = character()) {
    if (!name %in% names(files)) return(NULL)
    header <- names(utils::read.csv(files[[name]], nrows = 1, check.names = FALSE,
                                    fileEncoding = "UTF-8"))
    text <- intersect(c("site", "occasion_id", "species_code", "species", text), header)
    utils::read.csv(files[[name]], stringsAsFactors = FALSE, na.strings = "",
                    check.names = FALSE, fileEncoding = "UTF-8",
                    colClasses = stats::setNames(rep("character", length(text)), text))
  }
  # keys of the calibration beyond site, species and the grid are columns
  # of the occasions (schema 1.1), text like the identifiers
  keys <- setdiff(unlist(schema$calibration$keys), c("site", "species", "grid_index"))
  occasions <- read("occasions", c("source_audio", "recording_date",
                                   "recording_hour", "night", keys))
  detections <- read("detections")
  scores <- read("scores")
  validation <- read("validation", c("sampling_mode", "verdict", "reviewer"))
  verified <- read("recordings_verified", c("verdict", "reviewer"))
  sites <- read("sites")

  pick <- function(d) {
    if (is.null(species) || is.null(d)) return(d)
    cols <- intersect(c("species_code", "species"), names(d))
    hit <- purrr::reduce(purrr::map(cols, function(cl) d[[cl]] %in% species),
                         `|`, .init = rep(FALSE, nrow(d)))
    d[hit, , drop = FALSE]
  }
  detections <- pick(detections)
  scores <- pick(scores)
  validation <- pick(validation)
  verified <- pick(verified)
  several <- purrr::keep(list(validation = validation, verified = verified),
                         function(d) !is.null(d) && length(unique(d$species_code)) > 1)
  if (is.null(species) && length(several)) {
    warning("The export holds several species, so `",
            paste(names(several), collapse = "` and `"), "` mix the verdicts ",
            "of all of them. Name the target species in `species`, or all ",
            "the species to keep them for the multispecies model.", call. = FALSE)
  }

  if (!is.null(detections)) {
    names(detections)[names(detections) == "file_offset"] <- "start"
  }
  if (!is.null(validation)) {
    bands <- validation$sampling_mode %in% "score_bands"
    if (any(bands)) {
      if (all(bands) && (is.null(verified) || nrow(verified) == 0)) {
        stop("Every reviewed clip of the export was drawn with sampling_mode = ",
             "\"score_bands\", and there are no annotated recordings. These clips ",
             "were chosen by their score: they serve thresholds and calibration ",
             "in sheepdog, not the false-positive rate of crywolf's models. Review ",
             "clips drawn with \"random_detections\", or annotate recordings ",
             "drawn with \"random_recordings\". See ?fp_read_sheepdog.", call. = FALSE)
      }
      warning(sum(bands), " reviewed clip(s) drawn with sampling_mode = ",
              "\"score_bands\" were dropped: they were chosen by their score, for ",
              "thresholds and calibration in sheepdog, and would bias the ",
              "false-positive rate. See ?fp_read_sheepdog.", call. = FALSE)
      validation <- validation[!bands, , drop = FALSE]
    }
    top <- validation$sampling_mode %in% "top"
    if (any(top) && !keep_top) {
      warning(sum(top), " reviewed clip(s) drawn with sampling_mode = \"top\" ",
              "were dropped: they were chosen for their high scores and would ",
              "make the false-positive rate look lower than it is. ",
              "Use keep_top = TRUE to keep them.", call. = FALSE)
      validation <- validation[!top, , drop = FALSE]
    } else if (any(top)) {
      warning("Keeping ", sum(top), " clip(s) drawn with sampling_mode = ",
              "\"top\": the models of crywolf assume clips chosen ",
              "independently of the score.", call. = FALSE)
      # a clip drawn by both designs is one clip
      validation <- validation[!duplicated(validation[c("site", "occasion_id",
                                                        "species_code", "clip_start")]), ,
                               drop = FALSE]
    }
    pps <- validation$sampling_mode %in% "directed_pps"
    if (any(pps)) {
      warning(sum(pps), " reviewed clip(s) drawn with sampling_mode = \"directed_pps\": ",
              "their units were drawn with unequal probabilities (inclusion_prob), ",
              "possibly from the scores, and the models of crywolf do not weight ",
              "them. See ?fp_read_sheepdog.", call. = FALSE)
    }
    if (nrow(validation) == 0) validation <- NULL
  }
  if (!is.null(validation)) {
    names(validation)[names(validation) == "clip_start"] <- "start"
    names(validation)[names(validation) == "clip_end"] <- "end"
    names(validation)[names(validation) == "verdict"] <- "verdict_label"
    validation$verdict <- as.logical(validation$tp)
    validation$tp <- NULL
  }
  if (!is.null(verified)) {
    names(verified)[names(verified) == "verdict"] <- "verdict_label"
    names(verified)[names(verified) == "present"] <- "y"
    unknown <- is.na(verified$y)
    if (any(unknown)) {
      message(sum(unknown), " annotated recording(s) without a state were dropped.")
      verified <- verified[!unknown, , drop = FALSE]
    }
    verified$y <- as.integer(verified$y)
  }
  for (draws in c("calibration_draws", "occasion_draws")) {
    parquet <- file.path(path, paste0(draws, ".parquet"))
    if (file.exists(parquet)) files[[draws]] <- parquet
  }

  structure(list(occasions = occasions, detections = detections, scores = scores,
                 validation = validation, verified = verified, sites = sites,
                 schema = schema, files = files),
            class = "fp_sheepdog")
}

#' @export
print.fp_sheepdog <- function(x, ...) {
  s <- x$schema
  cat("<fp_sheepdog> export of sheepdog ", s$sheepdog_version,
      " (schema ", s$schema, " v", s$schema_version, ")\n", sep = "")
  if (!is.null(s$classifier)) cat("  classifier:          ", s$classifier, "\n")
  cat("  occasion unit:       ", s$occasion_unit, "; min_confidence ",
      s$min_confidence, "\n", sep = "")
  n <- function(d) if (is.null(d)) 0L else nrow(d)
  cat(sprintf("  occasions:            %d at %d site(s)\n", n(x$occasions),
              length(unique(x$occasions$site))))
  cat(sprintf("  predictions:          %d\n", n(x$detections)))
  if (!is.null(x$validation)) {
    cat(sprintf("  reviewed clips:       %d (%s)\n", n(x$validation),
                paste(names(table(x$validation$sampling_mode)), collapse = ", ")))
  }
  cat(sprintf("  annotated recordings: %d\n", n(x$verified)))
  if (!is.null(s$calibration)) {
    cat("  calibration:          ", s$calibration$formula, "\n", sep = "")
  }
  invisible(x)
}
