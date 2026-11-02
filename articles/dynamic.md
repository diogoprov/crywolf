# Use from one night to the next

*[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md)
is experimental: its interface may change in a later version.*

The other articles on the `curiango` data (Piatti et al., 2026) treat
each site-night as an independent site. A recorder that picked up the
Common Pauraque last night is then no more likely to pick it up tonight
than one that has been silent for a month.
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md)
drops that assumption. It fits the dynamic occupancy model (MacKenzie et
al., 2003) with the false positives of the multiple-methods design
(Miller et al., 2013), with the nights of each site as seasons. This
article fits it to the `curiango` data and compares it with the model of
independent site-nights by blocked cross-validation.

## The model

Each site i is in use (z\_{it} = 1) or not (z\_{it} = 0) on night t. On
the first night, \Pr(z\_{i1} = 1) = \psi_1. On every later night, the
state follows a Markov chain:

\Pr(z\_{it} = 1 \mid z\_{i,t-1} = 0) = \gamma\_{it}, \qquad \Pr(z\_{it}
= 0 \mid z\_{i,t-1} = 1) = \epsilon\_{it},

where the colonization \gamma\_{it} and the extinction \epsilon\_{it}
have logit-linear predictors (`col` and `ext`). Within a night, the
observation model is that of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md): true
detections with probability p\_{11} on nights in use, false positives
with probability p\_{10} on every night, and verified recordings that
settle what was recorded. The states are summed out by the forward
algorithm, so the model has no latent discrete parameters for Stan to
sample.

With nights as seasons, “occupancy” means use of the site on a given
night. Colonization is the start of use after a night without it, and
extinction is the end of use after a night with it. Neither is a change
in the population of the site.

`fp_data(season = )` takes the column of the night. When the seasons are
dates, it fills in the dates without data, so that every transition
spans one night. A site without data on a night contributes to that
night only through the transitions. In `curiango`, the nights run from
25 November 2023 to 26 February 2025, with between 1 and 7 of the seven
sites recorded on each night.

## Fitting

The colonization and the extinction depend on the temperature and
relative humidity of the night. Transitions into nights without data
need these covariates too, so we use the mean over the sites recorded on
each date (`temp_d`, `rh_d`), passed as season covariates. Every date
has at least one site recorded. Detection depends on the soundscape
indices of each minute, and the false-positive rate on the temperature
and humidity at the site that night.

``` r

zs <- function(v) as.numeric(scale(v))
m <- curiango$minutes
nt <- curiango$nights
daily <- aggregate(cbind(temp, rh) ~ night, nt, mean)
daily <- data.frame(night = daily$night, temp_d = zs(daily$temp), rh_d = zs(daily$rh))
i <- match(m$site_night, nt$site_night)
auto <- data.frame(site = m$site, night = nt$night[i], time = m$time,
                   y = as.integer(m$conf_custom >= 0.40),
                   bio = zs(m$bio), aci = zs(m$aci),
                   temp = zs(nt$temp[i]), rh = zs(nt$rh[i]), recorder = m$site)
ver <- data.frame(site = curiango$verified$site,
                  night = nt$night[match(curiango$verified$site_night, nt$site_night)],
                  y = curiango$verified$y)
d <- fp_data(auto, data.frame(site = unique(auto$site)), ver, site = "site",
             recorder = "recorder", season = "night", season_covs = daily)

fit <- fp_occu_dyn(d, psi1 = ~ 1, col = ~ temp_d + rh_d, ext = ~ temp_d + rh_d,
                   det = ~ bio + aci, fp = ~ temp + rh, threads_per_chain = 2)
```

The fit took about 35 minutes with four chains of two threads. It
converged: the largest R-hat was 1.003, the smallest bulk effective
sample size 783, with 0 divergent transitions.

## Results

``` r

vars <- c("col_(Intercept)", "col_temp_d", "col_rh_d", "ext_(Intercept)",
          "ext_temp_d", "ext_rh_d", "p11_ref", "p10_ref")
est <- res$summary[res$summary$variable %in% vars, c("variable", "mean", "lower", "upper")]
est[c("mean", "lower", "upper")] <- signif(est[c("mean", "lower", "upper")], 3)
est
#> # A tibble: 8 × 4
#>   variable            mean    lower   upper
#>   <chr>              <dbl>    <dbl>   <dbl>
#> 1 p11_ref          0.0584   0.0249   0.13  
#> 2 p10_ref          0.00306  0.00264  0.0035
#> 3 col_(Intercept) -3.05    -3.34    -2.77  
#> 4 col_temp_d      -0.237   -0.486    0.0153
#> 5 col_rh_d        -0.761   -1.09    -0.414 
#> 6 ext_(Intercept) -1.23    -1.61    -0.879 
#> 7 ext_temp_d      -0.38    -0.804    0.0266
#> 8 ext_rh_d         0.402   -0.0147   0.854
```

On a night of average temperature and humidity, a site not in use starts
being used with probability 0.045. A site in use stops with probability
0.23, so a spell of use lasts about 4.4 nights on average. Humid nights
start fewer spells: the humidity coefficient of colonization is -0.76
(-1.09 to -0.41). Its coefficient of extinction, 0.40 (-0.01 to 0.85),
points to humid nights also ending spells, but its interval reaches
zero. The intervals of both temperature coefficients include zero.

A minute of a night in use yields a true detection with probability
0.058 (0.025 to 0.130) at the reference values. A minute yields a false
positive with probability 0.0031 (0.0026 to 0.0035).

The model also estimates, for each night, the fraction of the seven
sites in use (`psi_fs`, from the states sampled given all the data). On
nights when a site was not recorded, its state is inferred from the
nights around it.

``` r

fs <- res$summary[grepl("^psi_fs\\[", res$summary$variable), ]
fs$night <- as.Date(fs$season)
ggplot2::ggplot(fs, ggplot2::aes(night, mean)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#2DE2E6",
                       alpha = 0.4) +
  ggplot2::geom_line(linewidth = 0.3) +
  ggplot2::labs(x = NULL, y = "Fraction of sites in use") +
  ggplot2::theme_minimal()
```

![Fraction of the seven sites in use on each night (posterior mean and
95% interval of psi_fs).](dynamic_files/figure-html/psi-fs-1.png)

Fraction of the seven sites in use on each night (posterior mean and 95%
interval of psi_fs).

[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)
gives, for each site and night, the probability that the site was in use
given all its data:

``` r

s <- res$sites
rec <- paste(s$site, s$night) %in% paste(curiango$nights$site, curiango$nights$night)
s$recorded <- ifelse(rec, "recorded", "not recorded")
ggplot2::ggplot(s, ggplot2::aes(as.Date(night), site, fill = occupied, alpha = recorded)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
  ggplot2::scale_alpha_manual(values = c(recorded = 1, "not recorded" = 0.35)) +
  ggplot2::labs(x = NULL, y = NULL, fill = "Pr(in use)", alpha = NULL) +
  ggplot2::theme_minimal()
```

![Probability that each site was in use on each night given all the
data. Nights without recordings at a site are
paler.](dynamic_files/figure-html/use-1.png)

Probability that each site was in use on each night given all the data.
Nights without recordings at a site are paler.

### Posterior predictive check

The posterior predictive check compares the detections of each
site-night with those replicated from the fit. It gives p = 0: the
observed discrepancy (median 566) exceeds the replicated one (median
385) on every draw. The model of independent site-nights fails this
check by a similar margin
([`vignette("markov")`](https://provetelab.org/crywolf/articles/markov.md)).
The misfit therefore lies within nights, where the minutes are taken as
independent, and not in the dynamics across nights.
`fp_occu_dyn(markov = )` adds the Markov term of
[`vignette("markov")`](https://provetelab.org/crywolf/articles/markov.md)
for runs of detections within a night; it is not fitted here.

## Comparison with independent site-nights

LOO cannot compare the two models here. The log-likelihood of
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md)
has one value per site, because the nights of a site are not
independent. With seven sites, leaving one out would remove a seventh of
the data at a time.

Instead, we cross-validated on blocks of nights. The dates were cut into
weeks, and the weeks were dealt at random into ten folds, the same for
every site. The first and the last week stayed in every training set, so
that every fit spans the same nights. Each model was fitted ten times,
each time without one fold. For the dynamic model, the left-out nights
became nights without data, whose states are predicted from the
neighbouring nights. The model of independent site-nights used the same
detection, false-positive and recorder terms, with
`occ ~ temp_d + rh_d`. The score of each site and fold is the log of the
posterior predictive density of its left-out minutes and verified
recordings (elpd).

``` r

transform(cv, elpd = round(elpd, 1), se = round(se, 1))
#>       model    elpd    se
#> 1       dyn -6369.9 890.9
#> 2       ind -6561.4 910.1
#> 3 dyn - ind   191.5  29.0
```

The dynamic model predicts the left-out nights better, by 191.5 (SE
29.0) in elpd over 69 site and fold units. It does better in all ten
folds and at all seven sites:

``` r

by_site <- res$cv$by_site
by_site[-1] <- round(by_site[-1], 1)
names(by_site) <- c("site", "minutes left out", "dynamic", "independent", "difference")
by_site
#>   site minutes left out dynamic independent difference
#> 1  ARA            18479 -1463.1     -1494.1       31.0
#> 2  BAI            12302  -688.6      -712.9       24.3
#> 3  BEP             7615  -148.1      -156.2        8.1
#> 4  MIM            19370 -2059.1     -2122.2       63.0
#> 5  REF            12212  -972.2     -1009.5       37.3
#> 6  UEM            12460  -313.1      -328.5       15.4
#> 7  XAR            13438  -725.7      -738.0       12.3
```

Two caveats apply. The SE treats the site and fold units as independent,
although the units of a fold share a training fit, so it is likely too
small. And each left-out week has nights in the training data on both
sides, so the dynamic model interpolates as well as forecasts. All 20
training fits converged; one fit of the independent model was rerun with
another seed after a chain stuck in a second mode (see
`dev/dyn-kfold/README.md` in the source repository).

## When to use it

- With repeated nights or visits at the same sites, use on one night is
  informative about the next, and the dynamic model uses that
  information.
- Each season needs covariates for the transitions into it, including
  seasons without data. Choose covariates that exist on every night.
- Compare with the independent model on blocks of seasons, not by LOO,
  whenever the sites are few.

## References

MacKenzie, D. I., Nichols, J. D., Hines, J. E., Knutson, M. G., &
Franklin, A. B. (2003). Estimating site occupancy, colonization, and
local extinction when a species is detected imperfectly. *Ecology*,
84(8), 2200–2207. <https://doi.org/10.1890/02-3090>

Miller, D. A. W., Nichols, J. D., Gude, J. A., Rich, L. N., Podruzny, K.
M., Hines, J. E., & Mitchell, M. S. (2013). Determining occurrence
dynamics when false positives occur: Estimating the range dynamics of
wolves from public survey data. *PLoS ONE*, 8(6), e65808.
<https://doi.org/10.1371/journal.pone.0065808>

Piatti, L., de Sousa, D. L. H., Oliveira, B. S., de Melo, A. V., Sugai,
L. S. M., & Provete, D. B. (2026). Methodological choices influence
ecological inference in passive acoustic monitoring of a Neotropical
nightjar. *Journal of Avian Biology*, in press.
