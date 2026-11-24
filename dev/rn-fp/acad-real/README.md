# Realistic ACAD design and sensitivity to the prior on sigma_r

Request from the FrogNet project (24/11/2026). The questions:

- Can `fp_occu_rn()` recover the site effect on per-individual detection (σ_r) under the design of the ACAD data?
- How much does the prior on σ_r drive the answer?

## Design (`run_grid.R`)

- **Sites and occasions:** 21 sites with 11–19 automated occasions each. The occasion counts are one fixed draw.
- **Confirmation:** every occasion is confirmed (true and false observations known). There are no verified data.
- **Fixed parameters:**
  - r = 0.27.
  - log λ has a slope of 0.5 on a standardized covariate.
  - σ_λ = 0.7.
- **Cells:** ψ {0.97, 0.8, 0.3} × p10 {0.05, 0.15} × σ_r {0, 1}, with 100 datasets per cell.
- **Priors:** each dataset is fitted under four priors on σ_r:
  - half-Normal(0, 0.5);
  - half-Normal(0, 1), the package default;
  - half-Normal(0, 2);
  - half-Student-t(3, 0, 1), from `fp_occu_rn_t3.stan`, which changes only that line of the package program.
- **Fit:** 4 chains of 1000 + 1000 iterations, adapt_delta 0.95, `n_max` = 40 for ψ = 0.97 and 30 otherwise.
- **Run:** AWS spot (c6a.16xlarge, then c7a.16xlarge after one interruption), about 2.5 h, about US$2.2.
- **Summaries:** `summarise.R` writes `grid.rds`.

## Convergence

- **R-hat and ESS:** R-hat never exceeded 1.05. 1–9% of fits had a largest R-hat above 1.01, and bulk ESS of the key parameters was above 400 in more than 99% of fits.
- **Divergences:** these are the problem. Fits with at least one divergence:

  | ψ | Fits with ≥ 1 divergence | Fits with > 5 |
  |---|---|---|
  | 0.3 | 77% | 19% |
  | 0.8 | 34% | 3% |
  | 0.97 | 19% | 1% |

  The median is 0 divergences, and the 99th percentile is 17.
- **Likely cause:** a funnel in the site effects when few sites are occupied (about 7 of 21 at ψ = 0.3). This is my interpretation and was not checked directly.
- **Effect on the conclusions:** the metrics computed on fits without divergences differ little from those on all fits (table in the session log). Even so, the ψ = 0.3 cells rest on fits with divergences and should be read with care. A rerun at adapt_delta 0.99 would settle this.

## Results (all fits; signal = lower 5% bound of σ_r > 0.1)

- **p10:** recovered in every cell, with median bias < 0.005 and 90% coverage of 0.85–0.94.
- **True σ_r = 0:**
  - The posterior median of σ_r is 0.30–0.70 and follows the prior scale.
  - False signal is 0–1% under HN(0.5), 9–33% under HN(1), 20–43% under HN(2) and 11–28% under t3. It is highest at ψ = 0.97.
  - The Savage–Dickey BF01 is rarely above 3, except under HN(2) (18–34% of fits), where the prior puts little mass near 0.
  - The data almost never support σ_r = 0 against σ_r > 0.
- **True σ_r = 1:**
  - Coverage is 0.92–1.00 under HN(1), HN(2) and t3, but 0.51–0.78 under HN(0.5), which pulls the median down by 0.44–0.62.
  - Power is 38–40% at ψ = 0.3, 58–61% at ψ = 0.8 and 68–71% at ψ = 0.97 under HN(1) or t3. Under HN(2) it is 56–83%.
  - BF01 < 1/3 in 13–46% of fits.
- **σ_λ:** σ_λ (truth 0.7) is underestimated by 0.04–0.27, most when σ_r = 0. Across datasets, the posterior medians of σ_r and σ_λ are correlated (Spearman 0.6–0.7 when σ_r = 0). Both parameters seem to absorb the same site heterogeneity in detection frequency. This is an interpretation.

## Reading

- With 21 sites and confirmation only, the posterior of σ_r is mostly the prior.
- The design cannot tell σ_r = 0 from a moderate σ_r. The answer depends on the prior's scale:
  - HN(0.5) suppresses false signals but misses a real σ_r = 1.
  - HN(2) gives more power and more false signals.
- Compared with step 1 at the ACAD size (verified occasions, `../README.md`), confirmation-only data give more false signal at high ψ.
- An analysis of the ACAD data should report the posterior of σ_r under at least two priors, together with this simulation. A non-null σ_r should not be interpreted from a single fit.
