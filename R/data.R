#' Passive acoustic monitoring of the Common Pauraque in the Pantanal
#'
#' Detections of the Common Pauraque (*Nyctidromus albicollis*) by two
#' automated classifiers in one-minute recordings made at seven sites in the
#' Pantanal wetlands and adjacent Cerrado, Mato Grosso do Sul, Brazil, between
#' November 2023 and February 2025, with an independent set of recordings
#' annotated by ear. The occupancy unit is the site-night, the occasion is the
#' one-minute recording, and the annotated recordings are the verified method
#' of the multiple-methods design.
#'
#' @format A list with three data frames:
#' \describe{
#'   \item{`minutes`}{One row per one-minute recording (97,806 rows):
#'     `site_night` (site-night identifier, the occupancy unit), `site`
#'     (recorder site, one of seven), `time` (start of the recording),
#'     `conf_custom` and `conf_birdnet` (confidence scores of a custom-trained
#'     classifier and of the BirdNET global model v2.4, 0 when the classifier
#'     returned nothing), and the acoustic indices `bio` (Bioacoustic Index) and
#'     `aci` (Acoustic Complexity Index).}
#'   \item{`nights`}{One row per site-night (1,753 rows): `site_night`, `site`,
#'     `night` (date of the evening), and the nightly means of temperature
#'     (`temp`, °C) and relative humidity (`rh`, %), and the fraction of the
#'     lunar disc illuminated (`moon`).}
#'   \item{`verified`}{One row per annotated recording (1,199 rows):
#'     `site_night`, `site`, `time`, and `y`, 1 if any 3-s segment of the
#'     recording held the species.}
#' }
#'
#' @details
#' The annotated recordings are an independent ground-truth test set of about
#' 1% of the recordings, chosen in proportion to the number of recordings at
#' each site, and annotated without reference to the classifiers (Piatti et
#' al., 2026). They can therefore serve as the verified method of
#' [fp_data()]. A detection by a classifier is a confidence score at or above a
#' threshold; Piatti et al. (2026) used 0.40 (custom) and 0.10 (BirdNET) to
#' maximize F1, and 0.95 for both to maximize precision.
#'
#' The source table has one row per detection window, with exact duplicates
#' within each recording; this dataset keeps one row per recording. The script
#' that builds it from the public deposit is `data-raw/curiango.R` in the
#' package repository.
#'
#' @source Piatti, L., de Sousa, D. L. H., Oliveira, B. S., de Melo, A. V.,
#'   Sugai, L. S. M., & Provete, D. B. (2026). Data and code for:
#'   Methodological choices influence ecological inference in passive acoustic
#'   monitoring of a Neotropical nightjar. Zenodo.
#'   \doi{10.5281/zenodo.21795898}. Licence CC BY 4.0.
#'
#' @references
#' Piatti, L., de Sousa, D. L. H., Oliveira, B. S., de Melo, A. V., Sugai,
#' L. S. M., & Provete, D. B. (2026). Methodological choices influence
#' ecological inference in passive acoustic monitoring of a Neotropical
#' nightjar. *Journal of Avian Biology*, in press.
#' @examples
#' str(curiango, max.level = 1)
#'
#' # detections by the custom classifier at the MaxF1 threshold
#' mean(curiango$minutes$conf_custom >= 0.40)
"curiango"
