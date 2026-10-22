# From BirdNET output to an occupancy model

A classifier such as BirdNET (Kahl et al., 2021) lists detections: for
each audio file, the time windows where it found a species, with a
confidence score. An occupancy model needs something else: one row per
sampling occasion, including the occasions where nothing was found, and
a separate record of which detections a person checked. This article
goes from the first to the second with
[`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md),
[`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md),
and
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md),
and fits the model to the result.

The example data are made up so that the article runs anywhere, but they
have the layout of real files: BirdNET selection tables in the Raven
format, and the clip-by-clip verdicts written by the ‘sheepdog’ package
(<https://provetelab.org/sheepdog/>), which cuts clips from BirdNET
detections for review and is installed from the same R-universe as
crywolf (<https://diogoprov.r-universe.dev>).

## Two kinds of checking

Two designs use checked data, and they need different files.

- In the **observation-confirmation design** (Chambert et al., 2015), a
  person listens to clips that the classifier flagged and says whether
  each holds the species. The verdicts refer to detections, so they
  enter as counts of confirmed true and false detections per occasion.
  This is what
  [`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)
  and the `validation` argument of
  [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
  produce.
- In the **multiple-methods design** (Miller et al., 2011), a person
  annotates whole recordings, whether or not the classifier flagged
  them. These are a second detection method, entered as the `verified`
  table of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md). In
  ‘sheepdog’, clips of whole recordings
  (`sampling_mode = "random_recordings"`) are of this kind and must not
  be passed to
  [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md).

Either way, the clips must be chosen without regard to the classifier’s
score. Clips picked because they have the highest scores contain fewer
false positives than the average detection and make the false-positive
rate look lower than it is. In ‘sheepdog’, draw them at random among the
detections (`sampling_mode = "random_detections"`).

## Example files

We start from a simulated dataset with known structure: 60 sites, each
recorded on 8 nights, with observation confirmation on about 40% of the
occasions.

``` r

library(crywolf)
sim <- simulate_fp(nsite = 60, nocc = 8, frac_confirmed = 0.4,
                   frac_verified = 0, seed = 1)
truth <- sim$data$auto
truth$file <- sprintf("S%02d_202401%02d_200000.wav", truth$site, truth$occasion)
head(truth[, c("site", "occasion", "file", "y", "n_true", "n_false")])
#>   site occasion                    file y n_true n_false
#> 1    1        1 S01_20240101_200000.wav 0     NA      NA
#> 2    2        1 S02_20240101_200000.wav 0      0       0
#> 3    3        1 S03_20240101_200000.wav 0     NA      NA
#> 4    4        1 S04_20240101_200000.wav 1     NA      NA
#> 5    5        1 S05_20240101_200000.wav 0     NA      NA
#> 6    6        1 S06_20240101_200000.wav 0     NA      NA
```

From it we write what the field data would look like. Each occasion with
a detection gets a BirdNET row with a high confidence, and some
occasions without one get a row below the threshold of 0.5 that we will
use. One selection table per site goes into a folder.

``` r

set.seed(2)
dir_tables <- file.path(tempdir(), "birdnet")
dir.create(dir_tables, showWarnings = FALSE)
low <- truth$y == 0 & stats::runif(nrow(truth)) < 0.2
det <- truth[truth$y == 1 | low, ]
det$confidence <- ifelse(det$y == 1, stats::runif(nrow(det), 0.5, 1),
                         stats::runif(nrow(det), 0.1, 0.5))
det$start <- 3 * sample(0:199, nrow(det), replace = TRUE)
purrr::walk(split(det, det$site), function(d) {
  tab <- data.frame(
    Selection = seq_len(nrow(d)), View = "Spectrogram 1", Channel = 1,
    `Begin Time (s)` = d$start, `End Time (s)` = d$start + 3,
    `Begin Path` = file.path("/audio", d$file), `File Offset (s)` = d$start,
    `Species Code` = "compau", `Common Name` = "Common Pauraque",
    Confidence = round(d$confidence, 4), check.names = FALSE)
  utils::write.table(tab, file.path(dir_tables, sprintf("S%02d.BirdNET.selection.table.txt",
                                                       d$site[1])),
                     sep = "\t", row.names = FALSE, quote = FALSE)
})
list.files(dir_tables)[1:3]
#> [1] "S01.BirdNET.selection.table.txt" "S02.BirdNET.selection.table.txt"
#> [3] "S03.BirdNET.selection.table.txt"
```

The verdicts follow the layout of `validation_by_clip.csv`: one row per
reviewed clip. A confirmed occasion contributes one positive clip if the
species was there and one negative clip if the classifier also fired on
something else.

``` r

conf <- truth[!is.na(truth$n_true) & truth$y == 1, ]
clips <- rbind(
  data.frame(conf[conf$n_true > 0, c("site", "file")], verdict = "positive"),
  data.frame(conf[conf$n_false > 0, c("site", "file")], verdict = "negative"))
clips <- merge(clips, det[, c("file", "start", "confidence")], by = "file")
validation <- data.frame(
  site = sprintf("S%02d", clips$site), species_code = "compau",
  source_audio = file.path("/audio", clips$file), clip_start = clips$start,
  clip_end = clips$start + 3, confidence = clips$confidence,
  sampling_mode = "random_detections", occasion_id = sub("\\.wav$", "", clips$file),
  verdict = clips$verdict)
f_validation <- file.path(tempdir(), "validation_by_clip.csv")
utils::write.csv(validation, f_validation, row.names = FALSE)
table(validation$verdict)
#> 
#> negative positive 
#>       12       56
```

## Reading the files

[`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md)
reads every selection table in the folder, standardizes the column names
(`Begin Time (s)` becomes `begin_time`), and adds the name of the audio
file and, when the file name carries a `YYYYMMDD_HHMMSS` stamp, the time
of each detection.

``` r

detections <- fp_read_birdnet(dir_tables)
head(detections[, c("file", "start", "species_code", "confidence", "time")])
#>                      file start species_code confidence                time
#> 1 S01_20240101_200000.wav   561       compau     0.3142 2024-01-01 20:09:21
#> 2 S01_20240103_200000.wav   330       compau     0.1231 2024-01-03 20:05:30
#> 3 S01_20240107_200000.wav   519       compau     0.1950 2024-01-07 20:08:39
#> 4 S02_20240104_200000.wav   189       compau     0.5073 2024-01-04 20:03:09
#> 5 S02_20240108_200000.wav   522       compau     0.3466 2024-01-08 20:08:42
#> 6 S03_20240103_200000.wav   237       compau     0.8389 2024-01-03 20:03:57
```

[`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)
reads the verdicts, in the format of ‘sheepdog’ 2.0.0 or in the
Portuguese format of the sampler up to version 1.2.0
(`validacao_por_clipe.csv`), and returns `TRUE` for a positive clip,
`FALSE` for a negative one, and `NA` for a clip not reviewed.

``` r

verdicts <- fp_read_validation(f_validation)
head(verdicts)
#>   site species_code                    file start end confidence verdict
#> 1  S02       compau S02_20240104_200000.wav   189 192  0.5073348   FALSE
#> 2  S03       compau S03_20240103_200000.wav   237 240  0.8389030    TRUE
#> 3  S03       compau S03_20240103_200000.wav   237 240  0.8389030   FALSE
#> 4  S03       compau S03_20240107_200000.wav   246 249  0.9667769    TRUE
#> 5  S05       compau S05_20240102_200000.wav   144 147  0.8415215    TRUE
#> 6  S05       compau S05_20240104_200000.wav   162 165  0.7387293    TRUE
#>           occasion_id
#> 1 S02_20240104_200000
#> 2 S03_20240103_200000
#> 3 S03_20240103_200000
#> 4 S03_20240107_200000
#> 5 S05_20240102_200000
#> 6 S05_20240104_200000
```

## Building the occasions

The classifier lists detections only, so the occasions come from a table
of all recordings: here, one per site and night, with the site
covariate.

``` r

occasions <- data.frame(site = truth$site, file = truth$file)
sites <- sim$data$sites
```

[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
marks an occasion as a detection when it has a detection at or above the
threshold, and turns the verdicts into the counts `n_true` and `n_false`
of the occasions that were checked. Verdicts on clips below the
threshold are dropped, since those clips are not detections at this
threshold.

``` r

auto <- fp_aggregate(occasions, detections, threshold = 0.5,
                     validation = verdicts)
head(auto)
#>   site                    file n_detections y  score n_true n_false
#> 1    1 S01_20240101_200000.wav            0 0 0.3142     NA      NA
#> 2    2 S02_20240101_200000.wav            0 0 0.0000     NA      NA
#> 3    3 S03_20240101_200000.wav            0 0 0.0000     NA      NA
#> 4    4 S04_20240101_200000.wav            1 1 0.8801     NA      NA
#> 5    5 S05_20240101_200000.wav            0 0 0.0000     NA      NA
#> 6    6 S06_20240101_200000.wav            0 0 0.0000     NA      NA
c(detections = sum(auto$y), checked = sum(!is.na(auto$n_true)))
#> detections    checked 
#>        160         60
all(auto$y == truth$y)
#> [1] TRUE
```

The occasions are the same as in the simulated truth.
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
warns when an occasion has more detections at or above the threshold
than reviewed clips, because the confirmation design assumes that every
detection of a checked occasion was checked.

## Fitting the model

The result goes to
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) as
the automated table, with the names of the confirmation counts.

``` r

d <- fp_data(auto, sites, confirmed = c("n_true", "n_false"))
d
#> <fp_data>
#>   sites:                  60
#>   automated occasions:    480 (naive occupancy 0.65)
#>   verified sites:         0 (0 occasions; verified detections at 0 sites)
#>   confirmed occasions:    60 (56 with a true observation, at 24 sites)
#>   sites known occupied:   24
#>   site covariates:        x, w 
#>   observation covariates: file, n_detections, score
```

``` r

m <- fp_occu(d, occ = ~ x, fp = ~ 1, random_recorder = FALSE, seed = 1,
             refresh = 0, show_messages = FALSE)
summary(m)
#> Posterior summary (95% equal-tailed intervals)
#> 
#>         variable    mean     sd   lower   upper rhat ess_bulk ess_tail
#>  occ_(Intercept) -0.0646 0.2792 -0.6259  0.4967    1     4268     2885
#>            occ_x -0.9257 0.3353 -1.6494 -0.3178    1     4557     2403
#>  det_(Intercept)  0.4965 0.1424  0.2276  0.7694    1     3412     2791
#>          p11_ref  0.6448 0.0315  0.5814  0.7035    1     3893     2873
#>   fp_(Intercept) -2.7320 0.2362 -3.2126 -2.2963    1     5341     3096
#>           fp_gap  3.2285 0.2818  2.7063  3.8093    1     3473     2585
#>          p10_ref  0.0625 0.0137  0.0387  0.0914    1     5341     3096
#>          p10_bar  0.0625 0.0137  0.0387  0.0914    1     5341     3096
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.002, min bulk/tail ESS 3412/2403
```

The simulation used an occupancy slope of -0.8 and a false-positive rate
of 0.05; both lie inside the 95% intervals. The rest of the analysis,
from convergence checks to predictions, is the same as in
[`vignette("crywolf")`](https://provetelab.org/crywolf/articles/crywolf.md).

## With real files

Real data differ from this example in a few ways that matter.

- The occasion table must list every recording, including those where
  BirdNET found nothing, with the site and any covariates. A list of the
  audio files on disk is a good start.
- An occasion can be shorter than a file. With one-minute occasions, add
  a column for the minute to both tables and pass
  `by = c("file", "minute")`.
- BirdNET tables hold every species; give `species` to
  [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
  to keep one.
- The threshold changes which occasions are detections, and with them
  the estimates.
  [`vignette("threshold")`](https://provetelab.org/crywolf/articles/threshold.md)
  shows how to check that.

## References

Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339. <https://doi.org/10.1890/14-1507.1>

Kahl, S., Wood, C. M., Eibl, M., & Klinck, H. (2021). BirdNET: A deep
learning solution for avian diversity monitoring. *Ecological
Informatics*, 61, 101236. <https://doi.org/10.1016/j.ecoinf.2021.101236>

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
<https://doi.org/10.1890/10-1396.1>
