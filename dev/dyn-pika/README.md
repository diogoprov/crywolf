# Second data set for the dynamic model (criterion 8)

`run.R` fits `fp_occu_dyn()` to the American pika surveys of Goldman et al.
(2023, *Ecology and Evolution* 13: e10019, doi:10.1002/ece3.10019). Its header
gives the design. The data are the `pika` dataset, whose
`post_process.csv` comes unchanged from Dryad (doi:10.5061/dryad.ttdz08m2n, CC0)
and is in `data-raw/`. The data were collected by the OSU Cascades
Human and Ecosystem Resilience and Sustainability Lab with support from the
National Park Service, and uses of them should attribute that source.

The data have 771 plots in 5 park units and 5,126 surveys from 2010 to 2021. Each
fit takes about 4 minutes on 8 cores, so it runs locally:

```bash
Rscript dev/dyn-pika/run.R
```

## Run of 05/11/2026

All five fits have 0 divergences, R-hat of at most 1.004 and bulk ESS of at least
1,109.

| Requirement | Result |
|---|---|
| PPC, automated | p = 0.78 (`dyn`); p = 0.21 (`ind`) |
| PPC, verified | p = 0.16 (`dyn`) |
| Priors of `fp_gap` | means of 1, 2 and 3 and an SD of 3 change no reported estimate by more than 0.023 on its own scale |
| Independent seasons | `dyn` is ahead of `ind` by 178.7 elpd (SE 42.6), with no Pareto k above 0.7 |

Colonization is 0.048 and extinction 0.071 per year, and `fp_gap` is 3.14.
The article reports about 0.06 for both with a different model, with a false
detection rate for each sign. The comparison with it is only of order of
magnitude.

Two simplifications matter. The verified occasion shares its survey with the
automated ones, and unambiguous scat is a subset of fresh scat, so the two methods
are not fully independent given the state. And the sign enters `det` and `fp` as
a covariate, whereas the article gives each sign its own model.

## Sensitivity to the verified method

Goldman et al. class fresh scat as unambiguous when a fresh haypile was seen on the
same survey, so every survey with both has `FS_U` above 0. `sensitivity.R` refits
the default model with two other verified methods. All three fits have 0
divergences, R-hat of at most 1.003 and bulk ESS of at least 1,709. The base fit
reproduces `run.R`, with no mean differing by more than 0.007.

| Verified method | Verified detections | PPC automated | PPC verified | `psi_eq` | Colonization | Extinction |
|---|---|---|---|---|---|---|
| Sighting or `FS_U` above 0 (base) | 732 | 0.81 | 0.20 | 0.41 (0.34, 0.48) | 0.048 | 0.071 |
| Sighting, or `FS_U` above 0 with no fresh haypile | 380 | 0.00 | 0.00 | 0.42 (0.35, 0.49) | 0.038 | 0.053 |
| Sighting only | 186 | 0.00 | 0.11 | 0.34 (0.28, 0.41) | 0.035 | 0.068 |

The two alternatives fail the automated check, and the first of them also fails the
verified one, so their estimates only bound the effect and are not better fits. In
the base design the automated scat and haypile occasions never coincide, because
scat seen with a haypile moves to the verified occasion. In the alternatives they
coincide on 423 surveys, 2.5 to 3 times the number expected from their marginal
rates. Part of this excess may only reflect that occupied plots have more signs. Our
reading, not tested here, is that signs of one survey are more alike than the model
allows, and that the rule hides this from the check. The base PPC is therefore not
evidence that the signs are independent.

Goldman et al. treat calls as ambiguous, although the false calls they expect
come from pikas just outside the plot, not from other species.
