# Reporting a crywolf analysis

A Bayesian occupancy analysis is reproducible and open to scrutiny only
if the paper says what the model was, how it was computed, how well the
computation went, and how much the conclusions depend on the priors. The
Bayesian Analysis Reporting Guidelines (BARG; Kruschke, 2021) list these
items, and Davis & Kay (2023) give advice on writing statistical methods
for ecologists. This article goes through the items for a crywolf
analysis and shows where each piece of information comes from in the
package. It ends with a template for the methods and results.

We use a simulated dataset so that the article runs anywhere; with real
data every step is the same.

``` r

library(crywolf)
sim <- simulate_fp(nsite = 200, nocc = 10, frac_verified = 0.3, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
             threads_per_chain = 1, seed = 1, refresh = 0, show_messages = FALSE)
```

## 1. The model and why

State the design (multiple methods, Miller et al., 2011, or observation
confirmation, Chambert et al., 2015), what an occasion is, which method
is automated and which is checked, and how the checked occasions were
chosen. Then write the model, with every submodel and its covariates.
[`print()`](https://rdrr.io/r/base/print.html) shows the formulas as
fitted:

``` r

m
#> Bayesian false-positive occupancy model (multiple-methods design)
#> 
#>   occupancy:       ~x 
#>   true detection:  ~1   
#>   false positive:  ~1  
#>   data:            200 sites, 2000 automated and 138 verified occasions
#> 
#>              variable    mean      sd   lower   upper rhat ess_bulk ess_tail
#>       occ_(Intercept)  0.0632 0.15216 -0.2369  0.3638    1     4136     3011
#>                 occ_x -0.5196 0.15798 -0.8399 -0.2173    1     4470     2496
#>       det_(Intercept)  0.3774 0.06691  0.2402  0.5057    1     4134     2467
#>  verified_(Intercept)  0.5939 0.24055  0.1305  1.0800    1     4163     2925
#>            p_verified  0.6424 0.05453  0.5326  0.7465    1     4163     2925
#>               p11_ref  0.5931 0.01614  0.5598  0.6238    1     4134     2467
#>        fp_(Intercept) -2.8435 0.15521 -3.1599 -2.5521    1     4416     2877
#>                fp_gap  3.2209 0.15815  2.9210  3.5493    1     3813     2912
#>               p10_ref  0.0556 0.00807  0.0407  0.0723    1     4416     2877
#>               p10_bar  0.0556 0.00807  0.0407  0.0723    1     4416     2877
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.001, min bulk/tail ESS 3813/2467
```

The equations of
[`vignette("crywolf")`](https://provetelab.org/crywolf/articles/crywolf.md)
can be copied into a paper. Say what the gap \Delta is for: it keeps the
true-detection probability above the false-positive rate at the
reference values of the covariates, which rules out the mirror-image
solution in which the two are swapped.

## 2. Priors, with a reason for each

List every prior.
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)
prints the ones used:

``` r

m$priors
#> <fp_priors> (logit scale)
#>   intercepts      ~ Normal(0, 1.5)
#>   slopes          ~ Normal(0, 2.5)
#>   fp_gap          ~ Normal(2, 1.5), truncated at 0
#>   sd_recorder     ~ half-Normal(0, 2)
#>   sd_fp_site      ~ half-Normal(0, 1)
#>   sd_det_site     ~ half-Normal(0, 1)
#>   sd among species ~ half-Normal(0, 1) (fp_occu_ms())
```

Give a reason for each, even a short one: the intercept and slope priors
are weakly informative on the logit scale, and the gap prior expresses
the belief that the classifier is right more often than not. If the gap
prior came from a reviewed sample
([`fp_gap_prior()`](https://provetelab.org/crywolf/reference/fp_gap_prior.md)),
give the counts and how the sample was drawn. A prior predictive check
([`fp_prior_predictive()`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md))
shows what the priors imply for data like yours, and is worth a
sentence.

## 3. Computation

Report the software and versions, the number of chains, warm-up and
sampling iterations, and the target acceptance rate. These are stored in
the fit:

``` r

c(crywolf = as.character(utils::packageVersion("crywolf")),
  cmdstanr = as.character(utils::packageVersion("cmdstanr")),
  CmdStan = cmdstanr::cmdstan_version())
#>      crywolf     cmdstanr      CmdStan 
#> "2.1.0.9000"      "0.9.0"     "2.40.0"
unlist(m$sampler)
#>            chains   parallel_chains threads_per_chain       iter_warmup 
#>              4.00              4.00              1.00           1000.00 
#>     iter_sampling       adapt_delta 
#>           1000.00              0.95
```

Then the diagnostics: the largest R-hat (Vehtari et al., 2021), the
smallest bulk and tail effective sample sizes, and the number of
divergent transitions.
[`summary()`](https://rdrr.io/r/base/summary.html) gives them per
parameter, and the fit keeps the extremes:

``` r

unlist(m$diagnostics)
#>   divergences max_treedepth      max_rhat  min_ess_bulk  min_ess_tail 
#>      0.000000      0.000000      1.001209   3813.305312   2466.593563
```

If the chains used more than one thread (`threads_per_chain`, `"auto"`
by default, recorded in `m$sampler`), say so: the posterior is the same,
but the draws cannot be reproduced exactly from the seed. Fit with
`threads_per_chain = 1` when they must be.

## 4. Model adequacy and comparison

A posterior predictive check tells whether the model reproduces the
features of the data that matter. Report the discrepancy measure and the
Bayesian p-value:

``` r

p <- fp_ppc(m, ndraws = 500, seed = 1)
p$p_value
#> [1] 0.028
```

Values near 0 or 1 flag a misfit. This model is the one that generated
the data, yet the p-value is low: a single small value calls for a look
at which sites contribute most to the discrepancy (`plot(p)`), not for
rejecting the model outright.

If models were compared, give the method (leave-one-site-out
cross-validation; Vehtari et al., 2017), the difference in expected log
predictive density with its standard error, and any sites with Pareto k
above 0.7, with how they were handled (`loo(reloo = TRUE)`).

## 5. The posterior

Report each quantity of interest with its posterior mean or median and a
credible interval, and say which interval (crywolf reports 95%
equal-tailed intervals). Report probabilities on the scale readers
understand, such as the false-positive rate (`p10_bar`) rather than its
logit.

``` r

s <- summary(m)$table
cols <- c("variable", "mean", "lower", "upper")
s[s$variable %in% c("occ_x", "p11_ref", "p10_bar"), cols]
#> # A tibble: 3 × 4
#>   variable    mean   lower   upper
#>   <chr>      <dbl>   <dbl>   <dbl>
#> 1 occ_x    -0.520  -0.840  -0.217 
#> 2 p11_ref   0.593   0.560   0.624 
#> 3 p10_bar   0.0556  0.0407  0.0723
```

## 6. Sensitivity to the priors

Refit with different priors and report whether the conclusions change.
For false-positive models the gap prior and the slope priors are the
natural candidates:

``` r

refit <- function(priors) {
  f <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
               priors = priors, threads_per_chain = 1, seed = 1, refresh = 0,
               show_messages = FALSE)
  t <- summary(f)$table
  t[t$variable %in% c("occ_x", "p10_bar"), cols]
}
sens <- list(
  default = s[s$variable %in% c("occ_x", "p10_bar"), cols],
  narrow_gap = refit(fp_priors(gap_mean = 2, gap_sd = 0.5)),
  tight_slopes = refit(fp_priors(slope_sd = 1))
)
out <- purrr::imap(sens, function(t, nm) data.frame(priors = nm, t)) |>
  purrr::list_rbind()
out[c("mean", "lower", "upper")] <- round(out[c("mean", "lower", "upper")], 3)
out
#>         priors variable   mean  lower  upper
#> 1      default    occ_x -0.520 -0.840 -0.217
#> 2      default  p10_bar  0.056  0.041  0.072
#> 3   narrow_gap    occ_x -0.520 -0.835 -0.213
#> 4   narrow_gap  p10_bar  0.060  0.045  0.078
#> 5 tight_slopes    occ_x -0.514 -0.822 -0.221
#> 6 tight_slopes  p10_bar  0.055  0.040  0.073
```

Here the posterior mean of the occupancy slope moves by at most 0.006
across the three sets of priors. When they do change, report it: it
means the data alone do not settle the question.

## 7. Data and code

Share the data and the script that fits the model, with the crywolf
version. The fit object is not needed: the script and the seed reproduce
it when one thread per chain is used.

## A template

The text below is filled in from the fit above; adapt it to the study.

> We fitted a multiple-methods false-positive occupancy model (Miller et
> al., 2011) in a Bayesian framework with the R package crywolf
> 2.1.0.9000, which uses Stan through CmdStan 2.40.0. Occupancy
> probability depended on the standardized covariate *x*; the
> probability of a true automated detection and the false-positive rate
> were constant. The false-positive intercept was the true-detection
> intercept minus a positive gap, with a Normal(2, 1.5) prior truncated
> at zero; intercepts had Normal(0, 1.5) priors and slopes Normal(0,
> 2.5) priors on the logit scale. We ran 4 chains of 1000 warm-up and
> 1000 sampling iterations. The largest R-hat was 1.001, the smallest
> bulk effective sample size was 3813, and there were 0 divergent
> transitions. A posterior predictive check based on the Freeman-Tukey
> discrepancy gave a Bayesian p-value of 0.03.
>
> Occupancy decreased with *x* (slope -0.52, 95% credible interval -0.84
> to -0.22). The mean false-positive rate per automated occasion was
> 0.056 (0.041 to 0.072). With a narrower gap prior or tighter slope
> priors, the posterior mean of the slope changed by at most 0.006.

## References

Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339. <https://doi.org/10.1890/14-1507.1>

Davis, A. J., & Kay, S. (2023). Writing statistical methods for
ecologists. *Ecosphere*, 14, e4539. <https://doi.org/10.1002/ecs2.4539>

Kruschke, J. K. (2021). Bayesian analysis reporting guidelines. *Nature
Human Behaviour*, 5, 1282–1291.
<https://doi.org/10.1038/s41562-021-01177-7>

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
<https://doi.org/10.1890/10-1396.1>

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27, 1413–1432.
<https://doi.org/10.1007/s11222-016-9696-4>

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC (with discussion). *Bayesian
Analysis*, 16, 667–718. <https://doi.org/10.1214/20-BA1221>
