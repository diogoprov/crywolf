#' crywolf: Bayesian false-positive occupancy models with calibration checks
#'
#' @description
#' Automated detectors, such as acoustic classifiers, return false positives,
#' and occupancy models that ignore them overestimate occupancy and distort
#' covariate effects. crywolf fits the multiple-methods (site-confirmation)
#' false-positive occupancy model of Miller et al. (2011), in which a subset of
#' detections is verified and therefore certain, and asks whether the
#' correction is calibrated for a given study design.
#'
#' The package has three parts:
#'
#' * **Data and model.** [fp_data()] (or [fp_data_wide()]) assembles the
#'   detections, and [fp_occu()] fits the model in Stan through 'cmdstanr'.
#'   The false-positive rate can be constant (`fp = ~ 1`) or depend on site
#'   covariates (`fp = ~ w`), and the true-detection process accepts
#'   observation covariates and a recorder random effect. The
#'   observation-confirmation design of Chambert et al. (2015), in which an
#'   expert checks the detections of some automated occasions, is supported
#'   through the `confirmed` argument of [fp_data()]. [fp_ppc()] runs a
#'   posterior predictive check, and [curiango] holds a real dataset.
#' * **Simulation.** [simulate_fp()] generates data under homogeneous,
#'   unstructured, and structured ("confusion") false-positive regimes.
#' * **Calibration.** [fp_calibrate()] and [fp_verification_curve()] fit the
#'   model to many simulated datasets and report the coverage of the 95%
#'   credible intervals and the bias of the quantities of interest, so that
#'   users can see whether the correction holds for their design and how much
#'   verification it needs.
#'
#' @section Getting started:
#' `vignette("crywolf")` introduces the model and the workflow. Further
#' articles cover simulation (`vignette("simulation")`), calibration and how
#' much to verify (`vignette("calibration")`), a case study with real data
#' (`vignette("curiango")`), the translation from 'unmarked'
#' (`vignette("unmarked")`), and the Stan program (`vignette("stan_model")`).
#' They are also on the package website, <https://provetelab.org/crywolf/>.
#'
#' @section Model:
#' For site \eqn{i}, automated occasion \eqn{k}, and verified occasion
#' \eqn{m}:
#' \deqn{z_i \sim \mathrm{Bernoulli}(\psi_i), \quad \mathrm{logit}(\psi_i) = X_i\beta}
#' \deqn{y^{aut}_{ik} \sim \mathrm{Bernoulli}(z_i p_{11,ik} + (1 - z_i) p_{10,i})}
#' \deqn{\mathrm{logit}(p_{11,ik}) = Z_{ik}\alpha + u_{j[ik]}, \quad
#'       \mathrm{logit}(p_{10,i}) = \alpha_1 - \Delta + W_i\gamma}
#' \deqn{y^{ver}_{im} \sim \mathrm{Bernoulli}(z_i r_{11})}
#' The latent state \eqn{z} is marginalized. The positive gap \eqn{\Delta}
#' keeps the true-detection probability above the false-positive probability
#' at the reference values of the covariates, which prevents label switching.
#'
#' @references
#' Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
#' Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
#' two types of observational error occur: non-detection and species
#' misidentification. *Ecology*, 92, 1422–1428. \doi{10.1890/10-1396.1}
#'
#' Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
#' positive detections in species occurrence data under different study
#' designs. *Ecology*, 96, 332–339. \doi{10.1890/14-1507.1}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats model.frame model.matrix terms delete.response quantile
#'   sd rnorm rbinom rbeta plogis qlogis setNames na.pass
#' @importFrom posterior as_draws_df
#' @importFrom loo loo
## usethis namespace: end
NULL
