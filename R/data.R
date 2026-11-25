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
#'     `site_night`, `site`, `time`, `y`, 1 if any 3-s segment of the
#'     recording held the species, and `test_set`, the test set the recording
#'     was drawn in (1, 2 or 3; see Details).}
#' }
#'
#' @details
#' The annotated recordings are a test set of about 1% of the recordings
#' (Piatti et al., 2026), assembled in two stages:
#'
#' * Test sets 1 and 2 (960 recordings, 43 with the species) are a weighted
#'   random draw of 1% of the night-time recordings of each site (16:00 to
#'   06:00). The weight of a recording was the product of a month weight (5
#'   from August to October, 3 from November to January, 1 otherwise) and an
#'   hour weight (4 from 18:00 to 22:00, 2 from 00:00 to 05:00, 1 otherwise).
#'   Drawn recordings without audio were dropped. Given year-month, hour and
#'   site, the classifier scores therefore do not predict inclusion.
#' * Test set 3 (239 recordings, 214 with the species) completes them with
#'   recordings chosen because a BirdNET run detected the species, plus some
#'   chosen by the custom classifier at one site (BEP), with a quota of about
#'   40 per site and none from MIM. That BirdNET run is not the one in
#'   `conf_birdnet`.
#'
#' The account of test sets 1 and 2 comes from the selection script of
#' L. Piatti. That of test set 3 comes from her selection folders and notes,
#' compared with the data.
#'
#' The verified method of [fp_data()] assumes that the verified occasions
#' were not chosen by what the automated method recorded. Test sets 1 and 2
#' meet this assumption once year-month, hour and site are accounted for.
#' Test set 3 does not. Its recordings hold the species more often than those
#' of sets 1 and 2 with the same classifier scores: among recordings with a
#' custom score below 0.3 and a BirdNET score above 0.1, 49 of 57 in set 3
#' and 8 of 39 in sets 1 and 2. The choice therefore used more than the
#' scores, and inclusion probabilities estimated from them would not remove
#' the bias. Leave set 3 out of the verified method, as in the example below.
#' Metrics of a classifier computed on all 1,199 recordings describe this
#' enriched set, not all recordings.
#'
#' A detection by a classifier is a confidence score at or above a threshold;
#' Piatti et al. (2026) used 0.40 (custom) and 0.10 (BirdNET) to maximize F1,
#' and 0.95 for both to maximize precision.
#'
#' The source table has one to 20 rows (detection windows) per recording,
#' and a classifier's score is non-zero only in the windows where it detected
#' the species. This dataset keeps one row per recording, with the maximum
#' score over its windows, so a recording is a detection when any of its
#' windows reaches the threshold. The script
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
#'
#' # share of annotated recordings with the species, by test set
#' tapply(curiango$verified$y, curiango$verified$test_set, mean)
#'
#' # the random draws only
#' ver <- curiango$verified[curiango$verified$test_set < 3, ]
"curiango"

#' Surveys of the American pika in five park units of the western United States
#'
#' Sign of the American pika (*Ochotona princeps*) recorded on repeated surveys
#' of 771 plots in five park units of Oregon, Idaho and California, from 2010
#' to 2021, with no surveys in 2015 and 2016. Plots are circles of 12 m radius,
#' searched for up to 30 minutes. Fresh scat, fresh haypiles and calls can be
#' false, whereas sightings and some fresh scat are unambiguous, so the data fit
#' a false-positive occupancy model with a verified method.
#'
#' @format A data frame with 5,126 rows, one per survey, and 29 columns:
#' \describe{
#'   \item{`plot`, `Plot_ID`}{Plot, with and without the park code.}
#'   \item{`Park`}{Park unit: Crater Lake (`CRLA`), Craters of the Moon
#'     (`CRMO`), Lava Beds (`LABE`), Lassen Volcanic (`LAVO`) or Newberry
#'     (`NEWB`, from 2017).}
#'   \item{`Survey_ID`, `Location_ID`}{Survey, and plot-year (4,012 values).}
#'   \item{`Field_Date`, `Year`, `Month`, `Yday`}{Date of the survey
#'     (day/month/year, as text), and its year, month and day of the year.}
#'   \item{`Season`, `ssn`}{Season of the survey, as text (`summer`, or `autumm`
#'     as spelled in the source) and as 0 or 1. The source describes `ssn` as
#'     the original survey (summer) or the resurvey (autumn), but 2,489 surveys
#'     are in autumn and 1,107 are resurveys by `resurv`.}
#'   \item{`FS`, `OS`, `FH`, `OH`}{Number of fresh and old scat and fresh and old
#'     haypile observations.}
#'   \item{`PC`, `PS`}{Number of calls and sightings.}
#'   \item{`FS_U`}{Number of fresh scat observations classed as unambiguous. See
#'     Details.}
#'   \item{`hassign`}{1 if the survey recorded any of the six signs above.}
#'   \item{`resurv`}{1 for a resurvey, 0 for the original survey. These 1,107 are
#'     the second surveys of a plot-year; the seven third surveys have 0.}
#'   \item{`PerClass_Rock`, `PerClass_Bare`, `PerClass_Grass`, `PerClass_Shrub`,
#'     `PerClass_Tree`, `PerClass_Forb`}{Cover class of rock, bare ground,
#'     grasses, shrubs, trees and forbs, as text. Goldman et al. (2023) give
#'     the method, but not the limits of each class or the meaning of `T`.}
#'   \item{`lava_complx`}{Rock complexity, 1 to 3. Goldman et al. (2023) give
#'     the levels as low, intermediate and high. It is missing on all surveys
#'     at Lassen Volcanic and on 294 others, although a note to their Table 1
#'     names Lava Beds.}
#'   \item{`tmaxmean`}{Mean daily maximum summer temperature (°C), from gridMET.}
#'   \item{`dailysumppt`}{Sum of daily precipitation over the previous water
#'     year (mm), from gridMET.}
#' }
#'
#' @details
#' Rows and columns are those of `post_process.csv` in the public deposit, with
#' empty cells, which the source defines as no data available, as `NA`. The script is `data-raw/pika.R` in the package
#' repository.
#'
#' Fresh scat is unambiguous when the surveyor's notes describe fresh mucous
#' holding it together, or when a fresh haypile was seen on the same survey
#' (Goldman et al., 2023). Those notes were kept consistently only in 2019 to
#' 2021. `FS_U` is `NA` on 682 surveys, all from 2010 to 2018, all with fresh
#' scat and none with a fresh haypile. `vignette("dynamic")` treats these `NA`
#' as 0, so the scat of those surveys is ambiguous.
#'
#' A verified method built from `FS_U` is not independent of the automated
#' occasions of the same survey. Every survey with both fresh scat and a fresh
#' haypile has `FS_U` above 0, so a haypile can turn scat into a verified
#' detection. The rule also differs between years, because of the notes.
#'
#' @source Goldman, M. R., Shinderman, M., Jeffress, M. R., Rodhouse, T. J., &
#'   Shoemaker, K. T. (2023). Data from: Integrating multiple sign types to
#'   improve occupancy estimation for inconspicuous species. Dryad.
#'   \doi{10.5061/dryad.ttdz08m2n}. Licence CC0. Data collected by the OSU
#'   Cascades Human and Ecosystem Resilience and Sustainability Lab with support
#'   from the National Park Service; uses of them should attribute that source.
#'
#' @references
#' Goldman, M. R., Shinderman, M., Jeffress, M. R., Rodhouse, T. J., &
#' Shoemaker, K. T. (2023). Integrating multiple sign types to improve occupancy
#' estimation for inconspicuous species. *Ecology and Evolution*, 13, e10019.
#' \doi{10.1002/ece3.10019}
#' @seealso `vignette("dynamic")` fits [fp_occu_dyn()] to these data.
#' @examples
#' table(pika$Park, pika$Year)
#'
#' # one automated occasion per ambiguous sign, and sightings and unambiguous
#' # fresh scat as the verified occasion, with years as seasons
#' u <- ifelse(is.na(pika$FS_U), 0, pika$FS_U)
#' year <- factor(pika$Year, levels = 2010:2021)
#' sign <- function(name, y) {
#'   data.frame(site = pika$plot, sign = name, y = as.integer(y), year = year)
#' }
#' auto <- rbind(sign("scat", pika$FS > 0 & u == 0), sign("haypile", pika$FH > 0),
#'               sign("call", pika$PC > 0))
#' verified <- data.frame(site = pika$plot, y = as.integer(pika$PS > 0 | u > 0),
#'                        year = year)
#' plots <- unique(pika[c("plot", "Park")])
#' d <- fp_data(auto, data.frame(site = plots$plot, park = plots$Park), verified,
#'              season = "year")
"pika"
