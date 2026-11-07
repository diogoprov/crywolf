# Surveys of the American pika in five park units of the western United States

Sign of the American pika (*Ochotona princeps*) recorded on repeated
surveys of 771 plots in five park units of Oregon, Idaho and California,
from 2010 to 2021, with no surveys in 2015 and 2016. Plots are circles
of 12 m radius, searched for up to 30 minutes. Fresh scat, fresh
haypiles and calls can be false, whereas sightings and some fresh scat
are unambiguous, so the data fit a false-positive occupancy model with a
verified method.

## Usage

``` r
pika
```

## Format

A data frame with 5,126 rows, one per survey, and 29 columns:

- `plot`, `Plot_ID`:

  Plot, with and without the park code.

- `Park`:

  Park unit: Crater Lake (`CRLA`), Craters of the Moon (`CRMO`), Lava
  Beds (`LABE`), Lassen Volcanic (`LAVO`) or Newberry (`NEWB`, from
  2017).

- `Survey_ID`, `Location_ID`:

  Survey, and plot-year (4,012 values).

- `Field_Date`, `Year`, `Month`, `Yday`:

  Date of the survey (day/month/year, as text), and its year, month and
  day of the year.

- `Season`, `ssn`:

  Season of the survey, as text (`summer`, or `autumm` as spelled in the
  source) and as 0 or 1. The source describes `ssn` as the original
  survey (summer) or the resurvey (autumn), but 2,489 surveys are in
  autumn and 1,107 are resurveys by `resurv`.

- `FS`, `OS`, `FH`, `OH`:

  Number of fresh and old scat and fresh and old haypile observations.

- `PC`, `PS`:

  Number of calls and sightings.

- `FS_U`:

  Number of fresh scat observations classed as unambiguous. See Details.

- `hassign`:

  1 if the survey recorded any of the six signs above.

- `resurv`:

  1 for a resurvey, 0 for the original survey. These 1,107 are the
  second surveys of a plot-year; the seven third surveys have 0.

- `PerClass_Rock`, `PerClass_Bare`, `PerClass_Grass`, `PerClass_Shrub`,
  `PerClass_Tree`, `PerClass_Forb`:

  Cover class of rock, bare ground, grasses, shrubs, trees and forbs, as
  text. Goldman et al. (2023) give the method, but not the limits of
  each class or the meaning of `T`.

- `lava_complx`:

  Rock complexity, 1 to 3. Goldman et al. (2023) give the levels as low,
  intermediate and high. It is missing on all surveys at Lassen Volcanic
  and on 294 others, although a note to their Table 1 names Lava Beds.

- `tmaxmean`:

  Mean daily maximum summer temperature (°C), from gridMET.

- `dailysumppt`:

  Sum of daily precipitation over the previous water year (mm), from
  gridMET.

## Source

Goldman, M. R., Shinderman, M., Jeffress, M. R., Rodhouse, T. J., &
Shoemaker, K. T. (2023). Data from: Integrating multiple sign types to
improve occupancy estimation for inconspicuous species. Dryad.
[doi:10.5061/dryad.ttdz08m2n](https://doi.org/10.5061/dryad.ttdz08m2n) .
Licence CC0. Data collected by the OSU Cascades Human and Ecosystem
Resilience and Sustainability Lab with support from the National Park
Service; uses of them should attribute that source.

## Details

Rows and columns are those of `post_process.csv` in the public deposit,
with empty cells, which the source defines as no data available, as
`NA`. The script is `data-raw/pika.R` in the package repository.

Fresh scat is unambiguous when the surveyor's notes describe fresh
mucous holding it together, or when a fresh haypile was seen on the same
survey (Goldman et al., 2023). Those notes were kept consistently only
in 2019 to 2021. `FS_U` is `NA` on 682 surveys, all from 2010 to 2018,
all with fresh scat and none with a fresh haypile.
[`vignette("dynamic")`](https://provetelab.org/crywolf/articles/dynamic.md)
treats these `NA` as 0, so the scat of those surveys is ambiguous.

A verified method built from `FS_U` is not independent of the automated
occasions of the same survey. Every survey with both fresh scat and a
fresh haypile has `FS_U` above 0, so a haypile can turn scat into a
verified detection. The rule also differs between years, because of the
notes.

## References

Goldman, M. R., Shinderman, M., Jeffress, M. R., Rodhouse, T. J., &
Shoemaker, K. T. (2023). Integrating multiple sign types to improve
occupancy estimation for inconspicuous species. *Ecology and Evolution*,
13, e10019. [doi:10.1002/ece3.10019](https://doi.org/10.1002/ece3.10019)

## See also

[`vignette("dynamic")`](https://provetelab.org/crywolf/articles/dynamic.md)
fits
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md)
to these data.

## Examples

``` r
table(pika$Park, pika$Year)
#>       
#>        2010 2011 2012 2013 2014 2017 2018 2019 2020 2021
#>   CRLA  135  138  171  141  130    0  123  128    0   95
#>   CRMO  106  170  176  122  132  119  129  121  130  130
#>   LABE  151  101  130  114  150  126  130  130   89   68
#>   LAVO  117  143   75  124  123  128  129  127   49    0
#>   NEWB    0    0    0    0    0  130  133  133    0  130

# one automated occasion per ambiguous sign, and sightings and unambiguous
# fresh scat as the verified occasion, with years as seasons
u <- ifelse(is.na(pika$FS_U), 0, pika$FS_U)
year <- factor(pika$Year, levels = 2010:2021)
sign <- function(name, y) {
  data.frame(site = pika$plot, sign = name, y = as.integer(y), year = year)
}
auto <- rbind(sign("scat", pika$FS > 0 & u == 0), sign("haypile", pika$FH > 0),
              sign("call", pika$PC > 0))
verified <- data.frame(site = pika$plot, y = as.integer(pika$PS > 0 | u > 0),
                       year = year)
plots <- unique(pika[c("plot", "Park")])
d <- fp_data(auto, data.frame(site = plots$plot, park = plots$Park), verified,
             season = "year")
#> Warning: Automated effort is very unequal among sites (3 to 57 occasions); check that occasions were not lost or duplicated.
#> Warning: 2 season(s) have no data at any site (e.g. 2015, 2016); they enter only through the transitions of occupancy.
```
