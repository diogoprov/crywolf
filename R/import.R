# Reading classifier output and building occasions -----------------------------

#' Read BirdNET selection tables
#'
#' Reads the selection tables that BirdNET writes in the Raven format (one
#' tab-separated `.txt` file per recording or per site), with the column
#' names standardized as in the birdnet-validation-sampler pipeline:
#' lower case, without units, and with spaces replaced by underscores
#' (`Begin Time (s)` becomes `begin_time`).
#'
#' @param path Paths of selection tables, or of directories that hold them.
#' @param pattern Regular expression for the table files inside directories.
#' @param recursive Search directories recursively?
#' @param tz Time zone of the recording times in the file names.
#'
#' @return A data frame with one row per detection, with the standardized
#'   columns of the tables plus `table` (the selection table it came from),
#'   `file` (the name of the audio file), `start` (start of the detection in
#'   the audio file, in seconds: `file_offset` when the table has it,
#'   otherwise `begin_time`), and `time` (recording start parsed from a
#'   `YYYYMMDD_HHMMSS` stamp in the file name, plus `start`; `NA` when the
#'   name has none).
#' @seealso [fp_aggregate()], [fp_read_validation()]
#' @export
#' @examples
#' tab <- tempfile(fileext = ".txt")
#' writeLines(c(
#'   paste("Selection", "View", "Channel", "Begin Time (s)", "End Time (s)",
#'         "Begin Path", "File Offset (s)", "Species Code", "Common Name",
#'         "Confidence", sep = "\t"),
#'   paste(1, "Spectrogram 1", 1, 0, 3, "/rec/SITE1_20240101_200000.wav", 0,
#'         "compau", "Common Pauraque", 0.81, sep = "\t")
#' ), tab)
#' fp_read_birdnet(tab)
fp_read_birdnet <- function(path, pattern = "\\.txt$", recursive = TRUE,
                            tz = "UTC") {
  files <- unlist(purrr::map(path, function(p) {
    if (dir.exists(p)) list.files(p, pattern, recursive = recursive, full.names = TRUE) else p
  }))
  if (!length(files)) stop("No selection table found.", call. = FALSE)
  out <- purrr::map(files, function(f) {
    # read as text, so that tables agree on column types before binding
    tb <- utils::read.delim(f, check.names = FALSE, colClasses = "character",
                            na.strings = c("", "NA"))
    names(tb) <- .std_names(names(tb))
    file_col <- intersect(c("begin_path", "begin_file"), names(tb))
    need <- c("begin_time", "end_time", "confidence")
    miss <- setdiff(need, names(tb))
    if (length(miss) || !length(file_col)) {
      stop("Table ", f, " lacks column(s): ",
           paste(c(miss, if (!length(file_col)) "begin_path or begin_file"),
                 collapse = ", "), ".", call. = FALSE)
    }
    # BirdNET writes a table with only the header for a recording without
    # detections
    if (nrow(tb) == 0) return(NULL)
    tb$table <- f
    tb$file <- basename(gsub("\\\\", "/", tb[[file_col[1]]]))
    tb$start <- if ("file_offset" %in% names(tb)) tb$file_offset else tb$begin_time
    tb
  }) |> purrr::compact()
  if (!length(out)) stop("The selection tables hold no detection.", call. = FALSE)
  cols <- purrr::reduce(purrr::map(out, names), union)
  out <- purrr::map(out, function(tb) {
    tb[setdiff(cols, names(tb))] <- NA_character_
    tb[cols]
  }) |> purrr::list_rbind()
  out[] <- purrr::map(out, utils::type.convert, as.is = TRUE)
  stamp <- regmatches(out$file, regexpr("\\d{8}_\\d{6}", out$file))
  out$time <- as.POSIXct(NA, tz = tz)
  has <- grepl("\\d{8}_\\d{6}", out$file)
  out$time[has] <- as.POSIXct(stamp, format = "%Y%m%d_%H%M%S", tz = tz) +
    out$start[has]
  out
}

#' Read the verdicts of the birdnet-validation-sampler
#'
#' Reads `validacao_por_clipe.csv`, the clip-by-clip result of the review
#' step of the birdnet-validation-sampler pipeline, and returns the verdicts
#' in a form that [fp_aggregate()] can attach to the occasions.
#'
#' @param path Path of `validacao_por_clipe.csv`.
#'
#' @details
#' The sampler picks clips among the detections of a species, by default
#' the ones with the highest confidence. The models in crywolf treat the
#' checked occasions as chosen independently of the classifier's score: clips
#' picked for their high scores hold fewer false positives than a random
#' sample, and would make the false-positive rate look lower than it is. Use
#' clips drawn at random among the detections for this purpose.
#'
#' @return A data frame with one row per clip: `site` (`local` in the
#'   sampler), `species_code`, `file` (the audio file the clip was cut from),
#'   `start` and `end` (in seconds within it), `confidence`, and `verdict`
#'   (`TRUE` for a true detection, `FALSE` for a false one, `NA` for clips
#'   not reviewed).
#' @seealso [fp_aggregate()]
#' @export
fp_read_validation <- function(path) {
  v <- utils::read.csv(path, stringsAsFactors = FALSE)
  need <- c("local", "species_code", "audio_origem", "inicio_clipe",
            "fim_clipe", "confidence", "veredito")
  miss <- setdiff(need, names(v))
  if (length(miss)) {
    stop("`path` lacks column(s): ", paste(miss, collapse = ", "),
         ". Is it the validacao_por_clipe.csv of the sampler?", call. = FALSE)
  }
  verdict <- c(positivo = TRUE, negativo = FALSE)[v$veredito]
  other <- setdiff(unique(v$veredito[is.na(verdict)]),
                   c("nao_validado", "clipe_ausente"))
  if (length(other)) {
    warning("Unknown verdict(s) treated as not reviewed: ",
            paste(other, collapse = ", "), ".", call. = FALSE)
  }
  data.frame(
    site = v$local, species_code = v$species_code,
    file = basename(gsub("\\\\", "/", v$audio_origem)),
    start = v$inicio_clipe, end = v$fim_clipe, confidence = v$confidence,
    verdict = unname(verdict), stringsAsFactors = FALSE
  )
}

#' Detections and confirmations per occasion
#'
#' Turns classifier detections into the automated occasions of [fp_data()]:
#' an occasion is a detection when the classifier found the species in it
#' with a confidence at or above `threshold`. Occasions come from a table you
#' provide, one row per occasion (for instance, one per recording), because
#' the classifier lists only the detections, not the recordings where it
#' found nothing.
#'
#' @param occasions Data frame with one row per occasion, with the columns
#'   named in `by`, the site, and any covariates.
#' @param detections Data frame of detections, such as the output of
#'   [fp_read_birdnet()], with the columns in `by` and `confidence`.
#' @param threshold Confidence at or above which a detection counts.
#' @param by Columns that identify an occasion in both tables.
#' @param species Optional species to keep, matched against
#'   `species_code` and `common_name` of `detections` (and `species_code` of
#'   `validation`).
#' @param validation Optional verdicts on detections, such as the output of
#'   [fp_read_validation()], with the columns in `by`, `confidence`, and
#'   `verdict`. They become counts of confirmed true and false observations
#'   for the observation-confirmation design (columns `n_true` and `n_false`,
#'   `NA` on occasions without a reviewed clip).
#' @param standardize Optional names of numeric columns of `occasions` to
#'   standardize (mean 0, standard deviation 1), as the model expects for
#'   detection covariates. Their means and standard deviations are kept in
#'   the attribute `"scaling"` of the result. Site covariates live in the
#'   sites table of [fp_data()]; standardize them there.
#'
#' @details
#' The observation-confirmation design assumes that every detection of a
#' checked occasion was checked. The function warns when an occasion has
#' more detections at or above the threshold than reviewed clips, and drops
#' verdicts on clips below the threshold, which are not detections at this
#' threshold.
#'
#' @return `occasions` with three or five more columns: `y` (1 for a
#'   detection), `score` (highest confidence in the occasion, 0 without any),
#'   `n_detections` (detections at or above the threshold), and, with
#'   `validation`, `n_true` and `n_false`. Pass it as `auto` to
#'   [fp_data()], with `confirmed = c("n_true", "n_false")` when it holds
#'   confirmations.
#' @seealso [fp_read_birdnet()], [fp_read_validation()], [fp_data()]
#' @export
#' @examples
#' occasions <- data.frame(site = c("A", "A", "B"),
#'                         file = c("a1.wav", "a2.wav", "b1.wav"))
#' detections <- data.frame(file = c("a1.wav", "a1.wav", "b1.wav"),
#'                          species_code = "compau",
#'                          confidence = c(0.9, 0.3, 0.2))
#' fp_aggregate(occasions, detections, threshold = 0.5)
fp_aggregate <- function(occasions, detections, threshold, by = "file",
                         species = NULL, validation = NULL, standardize = NULL) {
  .check_df(occasions, "occasions", by)
  .check_df(detections, "detections", c(by, "confidence"))
  if (anyDuplicated(occasions[by])) {
    stop("`occasions` must have one row per occasion; duplicated values of `",
         paste(by, collapse = "`, `"), "`.", call. = FALSE)
  }
  pick <- function(d) {
    if (is.null(species)) return(d)
    cols <- intersect(c("species_code", "common_name"), names(d))
    hit <- purrr::reduce(purrr::map(cols, function(cl) d[[cl]] %in% species),
                         `|`, .init = rep(FALSE, nrow(d)))
    d[hit, , drop = FALSE]
  }
  sp_cols <- intersect(c("species_code", "common_name"), names(detections))
  if (is.null(species) && length(sp_cols) &&
      length(unique(detections[[sp_cols[1]]])) > 1) {
    warning("`detections` holds several species and `species` is NULL, so a ",
            "detection of any of them counts. Name the target species.",
            call. = FALSE)
  }
  detections <- pick(detections)
  if (anyNA(detections$confidence)) {
    warning(sum(is.na(detections$confidence)), " detection(s) without a ",
            "confidence were dropped.", call. = FALSE)
    detections <- detections[!is.na(detections$confidence), , drop = FALSE]
  }
  key <- function(d) do.call(paste, c(unname(as.list(d[by])), sep = "\r"))
  k_occ <- key(occasions)
  k_det <- key(detections)
  unknown <- setdiff(unique(k_det), k_occ)
  if (length(unknown)) {
    warning(length(unknown), " occasion(s) with detections are not in ",
            "`occasions` and were ignored.", call. = FALSE)
  }
  above <- detections$confidence >= threshold
  occasions$n_detections <- as.integer(table(factor(k_det[above], levels = k_occ)))
  occasions$y <- as.integer(occasions$n_detections > 0)
  best <- tapply(detections$confidence, factor(k_det, levels = k_occ), max)
  best <- as.vector(best)
  occasions$score <- ifelse(is.na(best), 0, best)

  if (!is.null(validation)) {
    .check_df(validation, "validation", c(by, "confidence", "verdict"))
    validation <- pick(validation)
    validation <- validation[!is.na(validation$verdict), , drop = FALSE]
    below <- validation$confidence < threshold
    if (any(below)) {
      message(sum(below), " reviewed clip(s) below the threshold were dropped: ",
              "they are not detections at this threshold.")
      validation <- validation[!below, , drop = FALSE]
    }
    out_val <- setdiff(unique(key(validation)), k_occ)
    if (length(out_val)) {
      warning(length(out_val), " occasion(s) with reviewed clips are not in ",
              "`occasions` and were ignored.", call. = FALSE)
    }
    k_val <- factor(key(validation), levels = k_occ)
    n_true <- as.integer(table(k_val[validation$verdict]))
    n_false <- as.integer(table(k_val[!validation$verdict]))
    checked <- (n_true + n_false) > 0
    partial <- checked & (n_true + n_false) < occasions$n_detections
    if (any(partial)) {
      warning(sum(partial), " checked occasion(s) have more detections than ",
              "reviewed clips. The observation-confirmation design assumes ",
              "that every detection of a checked occasion was checked.",
              call. = FALSE)
    }
    occasions$n_true <- ifelse(checked, n_true, NA_integer_)
    occasions$n_false <- ifelse(checked, n_false, NA_integer_)
  }
  if (length(standardize)) {
    .check_df(occasions, "occasions", standardize)
    scaling <- data.frame(
      covariate = standardize,
      center = purrr::map_dbl(standardize, function(v) mean(occasions[[v]], na.rm = TRUE)),
      scale = purrr::map_dbl(standardize, function(v) stats::sd(occasions[[v]], na.rm = TRUE))
    )
    occasions[standardize] <- purrr::map(standardize, function(v) {
      (occasions[[v]] - scaling$center[scaling$covariate == v]) /
        scaling$scale[scaling$covariate == v]
    })
    attr(occasions, "scaling") <- scaling
  }
  occasions
}

.std_names <- function(x) {
  x <- tolower(x)
  x <- gsub("\\((s|hz)\\)", "", x)
  x <- trimws(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  sub("_$", "", x)
}
