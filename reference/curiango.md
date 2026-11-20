# Passive acoustic monitoring of the Common Pauraque in the Pantanal

Detections of the Common Pauraque (*Nyctidromus albicollis*) by two
automated classifiers in one-minute recordings made at seven sites in
the Pantanal wetlands and adjacent Cerrado, Mato Grosso do Sul, Brazil,
between November 2023 and February 2025, with an independent set of
recordings annotated by ear. The occupancy unit is the site-night, the
occasion is the one-minute recording, and the annotated recordings are
the verified method of the multiple-methods design.

## Usage

``` r
curiango
```

## Format

A list with three data frames:

- `minutes`:

  One row per one-minute recording (97,806 rows): `site_night`
  (site-night identifier, the occupancy unit), `site` (recorder site,
  one of seven), `time` (start of the recording), `conf_custom` and
  `conf_birdnet` (confidence scores of a custom-trained classifier and
  of the BirdNET global model v2.4, 0 when the classifier returned
  nothing), and the acoustic indices `bio` (Bioacoustic Index) and `aci`
  (Acoustic Complexity Index).

- `nights`:

  One row per site-night (1,753 rows): `site_night`, `site`, `night`
  (date of the evening), and the nightly means of temperature (`temp`,
  °C) and relative humidity (`rh`, %), and the fraction of the lunar
  disc illuminated (`moon`).

- `verified`:

  One row per annotated recording (1,199 rows): `site_night`, `site`,
  `time`, `y`, 1 if any 3-s segment of the recording held the species,
  and `test_set`, the test set the recording was drawn in (1, 2 or 3;
  see Details).

## Source

Piatti, L., de Sousa, D. L. H., Oliveira, B. S., de Melo, A. V., Sugai,
L. S. M., & Provete, D. B. (2026). Data and code for: Methodological
choices influence ecological inference in passive acoustic monitoring of
a Neotropical nightjar. Zenodo.
[doi:10.5281/zenodo.21795898](https://doi.org/10.5281/zenodo.21795898) .
Licence CC BY 4.0.

## Details

The annotated recordings are a test set of about 1% of the recordings
(Piatti et al., 2026), assembled in two stages:

- Test sets 1 and 2 (960 recordings, 43 with the species) are random
  draws. About 1% of the recordings of each site were drawn. The
  inclusion probability also varies with the month and the hour: it is
  higher from August to January and from 18:00 to 22:00. Given
  year-month, hour and site, the classifier scores hardly predict
  inclusion.

- Test set 3 (239 recordings, 214 with the species) completes them with
  recordings chosen because a BirdNET run detected the species, plus
  some chosen by the custom classifier at one site (BEP), with a quota
  of about 40 per site and none from MIM. That BirdNET run is not the
  one in `conf_birdnet`.

The original selection script was lost. This account comes from the
selection folders and notes kept by L. Piatti, compared with the data.

The verified method of
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)
assumes that the verified occasions were not chosen by what the
automated method recorded. Test sets 1 and 2 meet this assumption once
year-month, hour and site are accounted for. Test set 3 does not. Its
recordings hold the species more often than those of sets 1 and 2 with
the same classifier scores: among recordings with a custom score below
0.3 and a BirdNET score above 0.1, 49 of 57 in set 3 and 8 of 39 in sets
1 and 2. The choice therefore used more than the scores, and inclusion
probabilities estimated from them would not remove the bias. Leave set 3
out of the verified method, as in the example below. Metrics of a
classifier computed on all 1,199 recordings describe this enriched set,
not all recordings.

A detection by a classifier is a confidence score at or above a
threshold; Piatti et al. (2026) used 0.40 (custom) and 0.10 (BirdNET) to
maximize F1, and 0.95 for both to maximize precision.

The source table has one to 20 rows (detection windows) per recording,
and a classifier's score is non-zero only in the windows where it
detected the species. This dataset keeps one row per recording, with the
maximum score over its windows, so a recording is a detection when any
of its windows reaches the threshold. The script that builds it from the
public deposit is `data-raw/curiango.R` in the package repository.

## References

Piatti, L., de Sousa, D. L. H., Oliveira, B. S., de Melo, A. V., Sugai,
L. S. M., & Provete, D. B. (2026). Methodological choices influence
ecological inference in passive acoustic monitoring of a Neotropical
nightjar. *Journal of Avian Biology*, in press.

## Examples

``` r
str(curiango, max.level = 1)
#> List of 3
#>  $ minutes :'data.frame':    97806 obs. of  7 variables:
#>  $ nights  :'data.frame':    1753 obs. of  6 variables:
#>  $ verified:'data.frame':    1199 obs. of  5 variables:

# detections by the custom classifier at the MaxF1 threshold
mean(curiango$minutes$conf_custom >= 0.40)
#> [1] 0.01795391

# share of annotated recordings with the species, by test set
tapply(curiango$verified$y, curiango$verified$test_set, mean)
#>          1          2          3 
#> 0.03773585 0.05175983 0.89539749 

# the random draws only
ver <- curiango$verified[curiango$verified$test_set < 3, ]
```
