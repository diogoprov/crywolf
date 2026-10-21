#' A prior for the gap from validated detections
#'
#' Turns the precision of the classifier, measured on reviewed detections,
#' into a prior for the gap \eqn{\Delta} between the true-detection and
#' false-positive intercepts, and returns it in [fp_priors()].
#'
#' @param n_true,n_false Numbers of reviewed detections confirmed as true and
#'   as false. They must be a random sample of the detections at the
#'   threshold of the analysis; detections picked for their high confidence
#'   overstate the precision.
#' @param p11 Probability that an occasion at an occupied site is a
#'   detection, at the reference values of the covariates. A single value,
#'   or a vector of draws that expresses its uncertainty.
#' @param psi Occupancy probability at the reference values of the
#'   covariates. A single value, or a vector of draws.
#' @param ndraws Number of draws used to propagate the uncertainty.
#' @param ... Further arguments for [fp_priors()].
#'
#' @details
#' In the model, an occasion at an occupied site is a detection with
#' probability \eqn{p_{11}}, which includes false alarms; an occasion at an
#' unoccupied site, with probability \eqn{p_{10}}. Assuming that false alarms
#' occur at the same rate \eqn{p_{10}} at occupied sites, independently of
#' true detections, the probability that an occasion holds a true detection
#' is \eqn{s = (p_{11} - p_{10}) / (1 - p_{10})}, and the share of automated
#' detections that are true is
#' \deqn{\pi = \frac{\psi s}{\psi p_{11} + (1 - \psi) p_{10}}.}
#' For given \eqn{\pi}, \eqn{p_{11}} and \eqn{\psi}, this has a single
#' solution \eqn{p_{10}} between 0 and \eqn{p_{11}}, the smaller root of a
#' quadratic. The function draws \eqn{\pi} from its
#' Beta(1 + `n_true`, 1 + `n_false`) posterior, solves for \eqn{p_{10}},
#' computes \eqn{\Delta = \mathrm{logit}(p_{11}) - \mathrm{logit}(p_{10})},
#' and uses the mean and standard deviation of the result as the parameters
#' of the Normal prior on \eqn{\Delta}. The model truncates that prior at
#' zero; when most of the derived values are far from zero, as is usual,
#' the truncation changes little. Guesses of `p11` and `psi` usually come
#' from a pilot fit or from earlier studies; give them as vectors when they
#' are uncertain, so that the prior widens accordingly.
#'
#' @return An `fp_priors` object; the derived values are in the attribute
#'   `"gap_draws"`.
#' @seealso [fp_priors()], [fp_read_validation()]
#' @export
#' @examples
#' # 46 of 50 randomly drawn detections were true; p11 and psi from a pilot
#' fp_gap_prior(n_true = 46, n_false = 4, p11 = 0.3, psi = 0.4)
#' # with uncertainty in p11
#' fp_gap_prior(46, 4, p11 = rbeta(1000, 30, 70), psi = 0.4)
fp_gap_prior <- function(n_true, n_false, p11, psi, ndraws = 4000, ...) {
  stopifnot(n_true >= 0, n_false >= 0, ndraws >= 2, all(p11 > 0 & p11 < 1),
            all(psi > 0 & psi < 1))
  draw <- function(v) if (length(v) == 1) rep(v, ndraws) else sample(v, ndraws, replace = TRUE)
  prec <- stats::rbeta(ndraws, 1 + n_true, 1 + n_false)
  p11 <- draw(p11)
  psi <- draw(psi)
  # pi * (psi * p11 + (1 - psi) * p10) * (1 - p10) = psi * (p11 - p10),
  # a quadratic in p10 with one root in (0, p11): the smaller one
  a <- prec * (1 - psi)
  b <- -(prec * (1 - psi) - prec * psi * p11 + psi)
  c0 <- psi * p11 * (1 - prec)
  p10 <- (-b - sqrt(b^2 - 4 * a * c0)) / (2 * a)
  gap <- stats::qlogis(p11) - stats::qlogis(p10)
  pr <- fp_priors(gap_mean = mean(gap), gap_sd = stats::sd(gap), ...)
  attr(pr, "gap_draws") <- gap
  pr
}
