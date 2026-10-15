# crywolf: Bayesian false-positive occupancy models with calibration checks

Automated detectors, such as acoustic classifiers, return false
positives, and occupancy models that ignore them overestimate occupancy
and distort covariate effects. crywolf fits the multiple-methods
(site-confirmation) false-positive occupancy model of Miller et al.
(2011), in which a subset of detections is verified and therefore
certain, and asks whether the correction is calibrated for a given study
design.

The package has three parts:

- **Data and model.**
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) (or
  [`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md))
  assembles the detections, and
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
  fits the model in Stan through 'cmdstanr'. The false-positive rate can
  be constant (`fp = ~ 1`) or depend on site covariates (`fp = ~ w`),
  and the true-detection process accepts observation covariates and a
  recorder random effect. The observation-confirmation design of
  Chambert et al. (2015), in which an expert checks the detections of
  some automated occasions, is supported through the `confirmed`
  argument of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md).
  [`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md) runs
  a posterior predictive check, and
  [curiango](https://provetelab.org/crywolf/reference/curiango.md) holds
  a real dataset.

- **Simulation.**
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  generates data under homogeneous, unstructured, and structured
  ("confusion") false-positive regimes.

- **Calibration.**
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  and
  [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md)
  fit the model to many simulated datasets and report the coverage of
  the 95% credible intervals and the bias of the quantities of interest,
  so that users can see whether the correction holds for their design
  and how much verification it needs.

## Installation

crywolf is distributed through R-universe, with binaries for Windows,
macOS, and Linux. 'cmdstanr' comes from the Stan R-universe, and CmdStan
itself is installed once:

    install.packages("crywolf", repos = c("https://diogoprov.r-universe.dev",
                                          "https://stan-dev.r-universe.dev",
                                          getOption("repos")))
    cmdstanr::install_cmdstan()

The model is compiled the first time
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) needs
it, or ahead of time with
[`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md),
and the executable is kept in the user cache directory for later
sessions.

## Getting started

[`vignette("crywolf")`](https://provetelab.org/crywolf/articles/crywolf.md)
introduces the model and the workflow. Further articles cover simulation
([`vignette("simulation")`](https://provetelab.org/crywolf/articles/simulation.md)),
calibration and how much to verify
([`vignette("calibration")`](https://provetelab.org/crywolf/articles/calibration.md)),
a case study with real data
([`vignette("curiango")`](https://provetelab.org/crywolf/articles/curiango.md)),
the translation from 'unmarked'
([`vignette("unmarked")`](https://provetelab.org/crywolf/articles/unmarked.md)),
and the Stan program
([`vignette("stan_model")`](https://provetelab.org/crywolf/articles/stan_model.md)).
They are also on the package website, <https://provetelab.org/crywolf/>.

## Model

For site \\i\\, automated occasion \\k\\, and verified occasion \\m\\:
\$\$z_i \sim \mathrm{Bernoulli}(\psi_i), \quad \mathrm{logit}(\psi_i) =
X_i\beta\$\$ \$\$y^{aut}\_{ik} \sim \mathrm{Bernoulli}(z_i p\_{11,ik} +
(1 - z_i) p\_{10,i})\$\$ \$\$\mathrm{logit}(p\_{11,ik}) =
Z\_{ik}\alpha + u\_{j\[ik\]}, \quad \mathrm{logit}(p\_{10,i}) =
\alpha_1 - \Delta + W_i\gamma\$\$ \$\$y^{ver}\_{im} \sim
\mathrm{Bernoulli}(z_i r\_{11})\$\$ The latent state \\z\\ is
marginalized. The positive gap \\\Delta\\ keeps the true-detection
probability above the false-positive probability at the reference values
of the covariates, which prevents label switching. Priors are Normal on
the logit scale for intercepts and slopes, Normal truncated at zero for
\\\Delta\\, and half-Normal for the standard deviation of the recorder
effect;
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)
lists the defaults and changes them without recompiling.
[`vignette("stan_model")`](https://provetelab.org/crywolf/articles/stan_model.md)
walks through the Stan program.

## Computation

Each chain runs on one core by default. The likelihood is summed over
sites with Stan's `reduce_sum()`, so `threads_per_chain` in
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
splits each chain among several cores, which pays off with many sites or
occasions. Keep `parallel_chains * threads_per_chain` at or below the
number of physical cores.
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
fits many models, one chain per worker; compile the model with
[`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md)
before starting workers.

## Citation

`citation("crywolf")` gives the reference to the package. Every release
is archived on Zenodo;
[doi:10.5281/zenodo.23385045](https://doi.org/10.5281/zenodo.23385045)
always resolves to the latest version.

## References

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
[doi:10.1890/10-1396.1](https://doi.org/10.1890/10-1396.1)

Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339.
[doi:10.1890/14-1507.1](https://doi.org/10.1890/14-1507.1)

## See also

Useful links:

- <https://github.com/diogoprov/crywolf>

- <https://provetelab.org/crywolf/>

- Report bugs at <https://github.com/diogoprov/crywolf/issues>

## Author

**Maintainer**: Diogo B. Provete <dbprovete@gmail.com>
([ORCID](https://orcid.org/0000-0002-0097-0651))

Authors:

- Diogo B. Provete <dbprovete@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-0097-0651))

- Tessa Rhinehart ([ORCID](https://orcid.org/0000-0002-4352-3464))

- Liliana Piatti ([ORCID](https://orcid.org/0000-0001-6324-1723))

- Larissa Sayuri Moreira Sugai
  ([ORCID](https://orcid.org/0000-0002-6248-291X))
