#############################
# one proportion exact test #
#############################

#' Power Analysis for the Test of One Proportion (Exact Method)
#'
#' @description
#' Calculates power, sample size or effect size (only one can be NULL at a
#' time) for test of a proportion against a constant using the exact method.
#'
#' Formulas are validated using PASS documentation.
#'
#'
#' @param prob        probability of success under alternative.
#' @param req.sign    whether `prob` is smaller or larger than `null.prob`
#'                    (when minimum detectable prob is of interest).
#' @param null.prob   probability of success under null.
#' @param n           integer; sample size.
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided", "one.sided", or "two.one.sided". For
#'                    non-inferiority or superiority tests, add margin to the
#'                    null hypothesis value and use
#'                    \code{alternative = "one.sided"}.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @param utf         logical; whether the output should show Unicode
#'                    characters (if encoding allows for it). \code{FALSE} by
#'                    default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test ("exact").}
#'   \item{delta}{difference between `prob` and `null.prob`}
#'   \item{odds.ratio}{Odds-ratio \eqn{(prob / (1 - prob)) /
#'                     (null.prob / (1 - null.prob))}}
#'   \item{prob}{probability of success under alternative.}
#'   \item{null.prob}{probability of success under null.}
#'   \item{binom.alpha}{critical value(s).}
#'   \item{power}{statistical power \eqn{(1 - \beta)}.}
#'   \item{n}{sample size.}
#'
#' @references
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#' @examples
#' # power'
#' power.exact.oneprop(prob = 0.45, null.prob = 0.50,
#'                     alpha = 0.05, n = 500,
#'                     alternative = "one.sided")
#'
#' # sample size
#' power.exact.oneprop(prob = 0.45, null.prob = 0.50,
#'                     alpha = 0.05, power = 0.80,
#'                     alternative = "one.sided")
#'
#' @export power.exact.oneprop
power.exact.oneprop <- function(prob = NULL, req.sign = "+", null.prob = 0.50,
                                n = NULL, power = NULL, alpha = 0.05,
                                alternative = c("two.sided", "one.sided", "two.one.sided"),
                                verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  func.parms <- as.list(environment())

  if (!is.null(prob)) check.proportion(prob)
  null.prob <- check.margins(null.prob, check.proportion, alternative)
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = prob, n = n, power = power)

  if (requested == "n") {

    n <- power.z.oneprop(prob = prob, null.prob = null.prob, power = power, alpha = alpha,
                         alternative = alternative, ceil.n = TRUE, verbose = 0)$n

    # sort(...) generates a vector with 10, 30, 100, 300, ... < n
    for (m in sort(vapply(10 ^ seq(floor(log10(n))), function(p) p * c(1, 3), numeric(2)), decreasing = TRUE)) {
      if (n > m) {
        step.prob <- round(m / 10)
        achieved.power <- 0
        while (achieved.power < power) {
          achieved.power <- power.binom.test(size = n, prob = prob, null.prob = null.prob, alpha = alpha,
                                             alternative = alternative, plot = FALSE, verbose = 0)$power
          if (achieved.power < power) n <- n + step.prob
        } # while
        n <- n - ifelse(step.prob > 1, step.prob, 0)
      } # n > m
    }

  } else if (requested == "es") {

    prob <- power.binom.test(power = power, size = n, prob = NULL, req.sign = req.sign,
                             null.prob = null.prob, alpha = alpha, alternative = alternative,
                             plot = FALSE, verbose = 0, utf = FALSE)$prob

  }

  # calculate power (if requested == "power") or update it (if requested == "n")
  pwr.obj <- power.binom.test(size = n, prob = prob, null.prob = null.prob,
                              alpha = alpha, alternative = alternative,
                              plot = FALSE, verbose = 0)

  delta <- prob - null.prob
  odds.ratio <- (prob / (1 - prob)) /  (null.prob / (1 - null.prob))

  if (verbose > 0) {

    print.obj <- list(requested = requested,
                      tgt.effect = "prob",
                      test = "One Proportion",
                      alpha = pwr.obj$alpha,
                      alternative = alternative,
                      method = "exact",
                      prob = prob,
                      null.prob = null.prob,
                      delta = delta,
                      odds.ratio = odds.ratio,
                      size = n,
                      binom.alpha = pwr.obj$binom.alpha,
                      power = pwr.obj$power,
                      n = n)

    .print.pwrss.oneprop(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = "exact",
                           prob = prob,
                           null.prob = null.prob,
                           delta = delta,
                           odds.ratio = odds.ratio,
                           size = n,
                           binom.alpha = pwr.obj$binom.alpha,
                           alpha = pwr.obj$alpha,
                           power = pwr.obj$power,
                           n = n),
                      class = c("pwrss", "exact", "oneprop")))

} # power.exact.oneprop


#####################################
# one proportion approximate z test #
#####################################

#' Power Analysis for the Test of One Proportion (Normal Approximation Method)
#'
#' @description
#' Calculates power or sample size (only one can be NULL at a time) for test of
#' a proportion against a constant using the normal approximation method.
#'
#' Formulas are validated using the PASS documentation and G*Power.
#'
#' @details
#' * NB: The \code{pwrss.z.prop()} function is deprecated, but it will remain
#'   available as a wrapper for the \code{power.z.oneprop()} function during a
#'   transition period.
#'
#' @aliases power.z.oneprop pwrss.z.prop
#'
#'
#' @param prob        probability of success under alternative.
#' @param req.sign    whether `prob` is smaller or larger than `null.prob`
#'                    (when minimum detectable prob is of interest).
#' @param null.prob   probability of success under null.
#' @param n           integer; sample size.
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided", "one.sided", or "two.one.sided". For
#'                    non-inferiority or superiority tests, add margin to the
#'                    null hypothesis value and use
#'                    \code{alternative = "one.sided"}.
#' @param std.error   character; whether to calculate standard error using
#'                    "null" or "alternative" value. "null" by default.
#' @param arcsine     logical; whether arcsine transformation should be
#'                    applied. \code{FALSE} by default. Note that when
#'                    \code{arcsine = TRUE}, any specification to
#'                    \code{correct} and \code{std.error} will be ignored.
#' @param correct     logical; whether Yate's continuity correction should be
#'                    applied.
#' @param ceil.n      logical; whether sample size should be rounded up.
#'                    \code{TRUE} by default.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @param utf         logical; whether the output should show Unicode
#'                    characters (if encoding allows for it). \code{FALSE} by
#'                    default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test ("exact").}
#'   \item{delta}{difference between `prob` and `null.prob`}
#'   \item{odds.ratio}{Odds-ratio \eqn{(prob / (1 - prob)) /
#'                     (null.prob / (1 - null.prob))}}
#'   \item{mean}{mean of the alternative distribution.}
#'   \item{sd}{standard deviation of the alternative distribution.}
#'   \item{null.mean}{mean of the null distribution.}
#'   \item{null.sd}{standard deviation of the null distribution.}
#'   \item{z.alpha}{critical value(s).}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample size.}
#'
#' @references
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#' @examples
#' # power
#' power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05,
#'                 n = 500, alternative = "one.sided")
#'
#' # sample size
#' power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05,
#'                 power = 0.80, alternative = "one.sided")
#'
#' # effect size
#' power.z.oneprop(req.sign = "+", null.prob = 0.50, alpha = 0.05,
#'                 n = 500, power = 0.80, alternative = "one.sided")
#'
#' @export power.z.oneprop
power.z.oneprop <- function(prob = NULL, req.sign = "+", null.prob = 0.50,
                            n = NULL, power = NULL, alpha = 0.05,
                            alternative = c("two.sided", "one.sided", "two.one.sided"),
                            std.error = c("null", "alternative"),
                            arcsine = FALSE, correct = FALSE,
                            ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  std.error <- match.arg(std.error)
  func.parms <- as.list(environment())

  if (!is.null(prob)) check.proportion(prob)
  null.prob <- check.margins(null.prob, check.proportion, alternative)
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(arcsine, correct, ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = prob, n = n, power = power)

  if (alternative == "two.one.sided" && std.error == "null") {
    std.error <- "alternative"
    warning("`std.error` = \"null\" is ignored. Using \"alternative\" for equivalence or minimal effect testing.", call. = FALSE)
  }

  if (arcsine && correct) warning("Continuity correction does not apply to arcsine transformation approach.", call. = FALSE)

  pwr <- function(prob, null.prob, n, std.error, arcsine, correct, alpha, alternative) {

    if (arcsine) {

      var.num <- 1
      h <- switch(alternative,
                  `two.sided` = probs.to.h(prob, null.prob, FALSE)$h,
                  `one.sided` = probs.to.h(prob, null.prob, FALSE)$h,
                  `two.one.sided` = c(probs.to.h(prob, null.prob[1], FALSE)$h,
                                   probs.to.h(prob, null.prob[2], FALSE)$h))

    } else {

      var.num <- prob * (1 - prob)

      k <- c(0, 0)[c(TRUE, alternative == "two.one.sided")]

      if (correct) {
        k <- switch(alternative,
                    `one.sided` = ifelse(alternative == "one.sided" && prob > null.prob, -1, 1) / (2 * n),
                    `two.sided` = -1 / (2 * n),
                    `two.one.sided` =  c(-1 / (2 * n), 1 / (2 * n)))
      }

      h <- switch(alternative,
                  `two.sided` = abs(prob - null.prob) + k,
                  `one.sided` = prob - null.prob + k,
                  `two.one.sided` = c(prob - null.prob[1] + k[1], prob - null.prob[2] + k[2]))

    } # if arcsine

    null.sd <- ifelse(std.error == "null" && alternative %in% c("two.sided", "one.sided"),
                      sqrt((null.prob * (1 - null.prob)) / (prob * (1 - prob))), 1)

    lambda <- h / sqrt(var.num / n)

    if (alternative == "two.one.sided") {
      mean <- 0
      null.mean <- sort(-lambda)
    } else {
      mean <- ifelse(sign(prob - null.prob) == -1 && alternative == "two.sided", -lambda, lambda)
      null.mean <- 0
    }

    out.pwr <- power.z.test(mean = mean, sd = 1,
                            null.mean = null.mean, null.sd = null.sd,
                            alpha = alpha, alternative = alternative,
                            plot = FALSE, verbose = 0)

    out.pwr

  } # pwr()

  ss.no.correction <- function(prob, null.prob, power, std.error, arcsine, alpha, alternative) {

    beta <- 1 - power

    if (arcsine) {

      var.num <- 1
      h <- switch(alternative,
                  `two.sided` = probs.to.h(prob, null.prob, FALSE)$h,
                  `one.sided` = probs.to.h(prob, null.prob, FALSE)$h,
                  `two.one.sided` = c(probs.to.h(prob, null.prob[1], FALSE)$h,
                                      probs.to.h(prob, null.prob[2], FALSE)$h))

    } else {

      var.num <- prob * (1 - prob)

      h <- switch(alternative,
                  `two.sided` = abs(prob - null.prob),
                  `one.sided` = prob - null.prob,
                  `two.one.sided` = c(prob - null.prob[1], prob - null.prob[2]))

    } # if arcsine

    null.sd <- ifelse(std.error == "null", sqrt((null.prob * (1 - null.prob)) / (prob * (1 - prob))), 1)

    if (alternative == "two.sided") {

      M <- stats::qnorm(alpha / 2, sd = null.sd, lower.tail = FALSE) +
        stats::qnorm(beta, sd = null.sd, lower.tail = FALSE)
      n <- M ^ 2 * var.num / h ^ 2
    }

    if (alternative == "one.sided") {

      M <- stats::qnorm(alpha, sd = null.sd, lower.tail = FALSE) +
        stats::qnorm(beta, sd = null.sd, lower.tail = FALSE)
      n <- M ^ 2 * var.num / h ^ 2

    }

    if (alternative == "two.one.sided") {

      if (prob > min(null.prob) && prob < max(null.prob)) {
        # equivalence
        M <- stats::qnorm(alpha,     sd = null.sd, lower.tail = FALSE) +
             stats::qnorm(beta / 2,  sd = null.sd, lower.tail = FALSE)
      } else {
        # minimal effect
        M <- stats::qnorm(alpha / 2, sd = null.sd, lower.tail = FALSE) +
             stats::qnorm(beta,      sd = null.sd, lower.tail = FALSE)
      }
      n <- M ^ 2 * var.num / h ^ 2
      n <- max(n)

    }

    n

  } # ss.no.correction()

  ss <- function(prob, null.prob, power, std.error, arcsine, correct, alpha, alternative) {

    n.init <- ss.no.correction(prob = prob, null.prob = null.prob, power = power,
                               std.error = std.error, arcsine = arcsine,
                               alpha = alpha, alternative = alternative)

    n.init <- floor(n.init)

    stop <- FALSE
    sign.previous <- 0
    sign.switch <- 0
    while (isFALSE(stop)) {

      pwr.est <- pwr(prob = prob, null.prob = null.prob, n = n.init,
                     std.error = std.error, arcsine = arcsine, correct = correct,
                     alpha = alpha, alternative = alternative)$power

      sign.current <- sign(pwr.est - power)

      if (sign.current != 0 && sign.previous != 0 && sign.previous != sign.current)
        sign.switch <- sign.switch + 1

      n.init <- n.init + ifelse(sign.current > 0, -1, 1)
      sign.previous <- sign.current

      if (sign.switch >= 2) {
        if (sign.current > 0) n.init <- n.init + 1
        stop <- TRUE
      }

    } # while

    list(n = n.init)

  } # ss()

  min.pwr <- function(prob = NULL, n = NULL, power = NULL) {

    power - pwr(prob = prob, null.prob = null.prob, n = n, std.error = std.error, arcsine = arcsine,
                correct = correct, alpha = alpha, alternative = alternative)$power

  } # min.pwr() (for uniroot and optimize)

  if (requested == "n") {

    n <- ss(prob = prob, null.prob = null.prob, power = power, std.error = std.error, arcsine = arcsine,
            correct = correct, alpha = alpha, alternative = alternative)$n

    if (ceil.n) n <- ceiling(n)

  } else if (requested == "es") { # sample size

    if (check.null_sign(req.sign, alternative)) {

      int.lower <- c(min(null.prob) + 1e-6, mean(null.prob) - 1e-6)
      int.upper <- c(mean(null.prob) + 1e-6, max(null.prob) - 1e-6)
      prob.lower <- stats::optimize(f = function(prob) min.pwr(prob, n, power) ^ 2, interval = int.lower, tol = 1e-12)$minimum
      prob.upper <- stats::optimize(f = function(prob) min.pwr(prob, n, power) ^ 2, interval = int.upper, tol = 1e-12)$minimum

      prob <- mean(c(prob.lower, prob.upper))

      warn.txt <- ifelse(max(abs(c(min.pwr(prob.lower, n, power), min.pwr(prob.upper, n, power)))) < 1e-6,
                         sprintf("Target effect ranges from %.4f to %.4f within the null bounds.", prob.lower, prob.upper),
                         "The target power rate cannot be achieved within the null bounds.")
      warning(warn.txt, call. = FALSE)

    } else {

      val.rng <- get.interval(null.ncp = null.prob, distribution = "binom", alternative = alternative, req.sign = req.sign)
      prob <- stats::optimize(f = function(prob) min.pwr(prob, n, power) ^ 2, interval = val.rng, tol = 1e-12)$minimum

    }

  } # effect size

  # calculate power (if requested == "power") or update it (if requested == "n")
  pwr.obj <- pwr(prob = prob, null.prob = null.prob, n = n, std.error = std.error,
                 arcsine = arcsine, correct = correct, alpha = alpha, alternative = alternative)

  odds.ratio <- (prob / (1 - prob)) /  (null.prob / (1 - null.prob))

  if (verbose > 0) {

    print.obj <- list(requested = requested,
                      tgt.effect = "prob",
                      test = "One Proportion",
                      alpha = alpha,
                      alternative = alternative,
                      method = "z",
                      std.error = std.error,
                      arcsine = arcsine,
                      correct = correct,
                      prob = prob,
                      null.prob = null.prob,
                      delta = prob - null.prob,
                      odds.ratio = odds.ratio,
                      mean = pwr.obj$mean,
                      sd = pwr.obj$sd,
                      null.mean = pwr.obj$null.mean,
                      null.sd = pwr.obj$null.sd,
                      z.alpha = pwr.obj$z.alpha,
                      power = pwr.obj$power,
                      n = n)

    .print.pwrss.oneprop(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = "z",
                           prob = prob,
                           null.prob = null.prob,
                           delta = prob - null.prob,
                           odds.ratio = odds.ratio,
                           mean = pwr.obj$mean,
                           sd = pwr.obj$sd,
                           null.mean = pwr.obj$null.mean,
                           null.sd = pwr.obj$null.sd,
                           z.alpha = pwr.obj$z.alpha,
                           power = pwr.obj$power,
                           n = n),
                      class = c("pwrss", "z", "oneprop")))
} # power.z.oneprop

#' @export pwrss.z.prop
pwrss.z.prop <- function(p, p0 = 0.50, margin = 0, arcsin.trans = FALSE, alpha = 0.05,
                          alternative = c("not equal", "greater", "less",
                                          "equivalent", "non-inferior", "superior"),
                          n = NULL, power = NULL, verbose = TRUE) {

  alternative <- match.arg(alternative)
  verbose <- ensure.verbose(verbose)

  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"
  std.error <- ifelse(alternative == "two.one.sided", "alternative", "null")

  check.numeric(margin)
  check.logical(arcsin.trans)

  if (margin > 0.99 || margin < -0.99) stop("Provide a reasonable margin consistent with `p` - `p0`.", call. = FALSE)

  if (alternative == "two.one.sided") margin <- c(-margin, margin)
  null.prob <- p0 + margin

  prop.obj <-  power.z.oneprop(prob = p, null.prob = null.prob, arcsine = arcsin.trans,
                                      n = n, power = power, alpha = alpha,
                                      alternative = alternative, std.error = std.error,
                                      ceil.n = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.z.oneprop() function. \n")

  invisible(prop.obj)

} # pwrss.z.prop()


##############################
# two proportions exact test #
##############################

#' Power Analysis for Testing the Difference Between Two Proportions (Exact Method)
#'
#' @description
#' Calculates power or sample size (only one can be NULL at a time) for two
#' proportions using the exact method. The function is a wrapper for
#' `power.exact.mcnemar` (if `paired` == TRUE), or `power.exact.fisher` (if
#' `paired` == FALSE)
#'
#' Validated via G*Power and PASS documentation.
#'
#' @aliases power.exact.twoprops power.exact.twoprop
#'
#'
#' @param prob1       probability of success in the first group.
#' @param prob2       probability of success in the second group.
#' @param req.sign    whether estimated prob is smaller or larger than the
#'                    other (when minimum detectable prob is of interest).
#' @param n.ratio     sample size ratio (n1 / n2).
#' @param n2          integer; sample size for the second group.
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided" or "one.sided".
#' @param paired      logical; if \code{TRUE} samples are paired. \code{FALSE}
#'                    by default.
#' @param rho.paired  correlation between paired observations.
#' @param method      character; whether to use "approximate" or "exact"
#'                    method. Default is \code{"exact"} (only in the
#'                    \code{power.exact.twoprops()} function).
#' @param ceil.n      logical; \code{TRUE} rounds up sample size in each group.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @param utf         logical; whether the output should show Unicode
#'                    characters (if encoding allows for it). \code{FALSE} by
#'                    default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the test, which is "z" or "exact".}
#'   \item{delta}{difference between `prob1` and `prob2`}
#'   \item{odds.ratio}{Odds-ratio \eqn{(prob1 / (1 - prob1)) /
#'                     (prob2 / (1 - prob2))}}
#'   \item{size}{... (applies to paired proportions).}
#'   \item{prob}{... (applies to paired proportions).}
#'   \item{null.prob}{... (applies to paired proportions).}
#'   \item{binom.alpha}{critical value(s) (applies to paired proportions).}
#'   \item{mean}{mean of the alternative distribution.}
#'   \item{sd}{standard deviation of the alternative distribution.}
#'   \item{null.mean}{mean of the null distribution.}
#'   \item{null.sd}{standard deviation of the null distribution.}
#'   \item{alternative}{standard deviation of the null distribution.}
#'   \item{z.alpha}{critical value(s).}
#'   \item{power}{statistical power \eqn{(1 - \beta)}.}
#'   \item{n}{sample size in the form of c(n1, n2) (applies to independent
#'            proportions).}
#'   \item{n.total}{total sample size (applies to independent proportions).}
#'   \item{n.paired}{paired sample size (applies to paired proportions).}
#'
#' @references
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#' @examples
#'   # power
#'   power.exact.twoprops(prob1 = 0.70, prob2 = 0.60,
#'                        alpha = 0.05, n2 = 500,
#'                        alternative = "one.sided")
#'
#'   power.exact.twoprops(prob1 = 0.70, prob2 = 0.60,
#'                        alpha = 0.05, n2 = 500,
#'                        alternative = "one.sided",
#'                        paired = TRUE)
#'
#'   # sample size
#'   power.exact.twoprops(prob1 = 0.70, prob2 = 0.60,
#'                        alpha = 0.05, power = 0.80,
#'                        alternative = "one.sided")
#'
#'   power.exact.twoprops(prob1 = 0.70, prob2 = 0.60,
#'                        alpha = 0.05, power = 0.80,
#'                        alternative = "one.sided",
#'                        paired = TRUE)
#'
#' @export power.exact.twoprops
power.exact.twoprops <- function(prob1 = NULL, prob2 = NULL, req.sign = "+", n.ratio = 1, n2 = NULL,
                                 power = NULL, alpha = 0.05, alternative = c("two.sided", "one.sided"),
                                 paired = FALSE, rho.paired = 0.50, method = c("exact", "approximate"),
                                 ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  method <- match.arg(method)

  check.positive(n.ratio)
  if (!is.null(prob1)) check.proportion(prob1)
  if (!is.null(prob2)) check.proportion(prob2)
  if (!is.null(n2)) check.sample.size(n2)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(paired)
  check.correlation(rho.paired)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = list(prob1, prob2), n = n2, power = power)

  min.pwr.paired <- function(prob1, prob2, power) {

    jp <- joint.probs.2x2(prob1 = prob1, prob2 = prob2, rho = rho.paired, verbose = 0)

    power - power.exact.mcnemar(prob10 = jp$prob10, prob01 = jp$prob01, req.sign = req.sign, power = NULL, n.paired = n2,
                                alpha = alpha, alternative = alternative, method =  method, verbose = 0)$power

  } # min.pwr.paired (for stats::optimize)

  if (requested == "es" && paired) {

    if (is.null(prob1)) {

      val.rng <- sort(c(unname(prob.limits.paired(prob2 = prob2, rho = rho.paired)), prob2))[ifelse(check.pos_sign(req.sign), -1, -3)]
      prob1 <- stats::optimize(f = function(prob1) min.pwr.paired(prob1, prob2, power) ^ 2, interval = val.rng)$minimum

    } else {

      val.rng <- sort(c(unname(prob.limits.paired(prob1 = prob1, rho = rho.paired)), prob1))[ifelse(check.pos_sign(req.sign), -1, -3)]
      prob2 <- stats::optimize(f = function(prob2) min.pwr.paired(prob1, prob2, power) ^ 2, interval = val.rng)$minimum

    } # prob1 or prob2?

    power <- NULL # reset power (which is then updated / calculated below)

  }

  if (paired) {

    jp <- joint.probs.2x2(prob1 = prob1, prob2 = prob2, rho = rho.paired, verbose = 0)

    power.exact.mcnemar(prob10 = jp$prob10, prob01 = jp$prob01, power = power, n.paired = n2, alpha = alpha,
                        alternative = alternative, method =  method, ceil.n = ceil.n, verbose = verbose, utf = utf)

  } else {

    power.exact.fisher(prob1 = prob1, prob2 = prob2, req.sign = req.sign, power = power, n2 = n2, n.ratio = n.ratio,
                       alpha = alpha, alternative = alternative, method = method, ceil.n = ceil.n, verbose = verbose,
                       utf = utf)

  } # paired

} # power.exact.twoprops()

#' @rdname power.exact.twoprops
#' @export power.exact.twoprop
power.exact.twoprop <- power.exact.twoprops


########################
# two proportions test #
########################

#' Power Analysis for Testing the Difference Between Two Proportions (Normal
#' Approximation Method)
#'
#' @description
#' Calculates power or sample size (only one can be NULL at a time) for two
#' proportions using the normal approximation method.
#'
#' Validated via G*Power and PASS documentation.
#'
#' @details
#' * NB: The \code{pwrss.z.2props()} function is deprecated, but it will remain
#'   available as a wrapper for the \code{power.z.twoprops()} function during
#'   a transition period.
#'
#' @aliases power.z.twoprops power.z.twoprop pwrss.z.2props pwrss.z.2prop
#'
#'
#' @param prob1       probability of success in the first group.
#' @param prob2       probability of success in the second group.
#' @param req.sign    whether estimated prob is smaller or larger than the
#'                    other (when minimum detectable prob is of interest).
#' @param margin      ignorable \code{prob1} - \code{prob2} difference. For two
#'                    one-sided tests provide lower and upper margins in the
#'                    form of \code{c(lower, upper)}.
#' @param n.ratio     sample size ratio (n1 / n2).
#' @param n2          integer; sample size for the second group.
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided", "one.sided", or "two.one.sided".
#' @param arcsine     logical; whether arcsine transformation should be
#'                    applied. Note that this only applies to independent
#'                    proportions without continuity correction.
#' @param correct     logical; whether Yates' continuity correction should be
#'                    applied to the test statistic. Ignored for the paired
#'                    test.
#' @param paired      logical; if \code{TRUE} samples are paired. \code{FALSE}
#'                    by default.
#' @param rho.paired  correlation between paired observations.
#' @param std.error   character; whether to calculate standard error using
#'                    "pooled" or "unpooled" standard deviation. Ignored for
#'                    the paired test.
#' @param ceil.n      logical; \code{TRUE} rounds up sample size in each group.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @param utf         logical; whether the output should show Unicode
#'                    characters (if encoding allows for it). \code{FALSE} by
#'                    default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the test, which is "z" or "exact".}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{mean}{mean of the alternative distribution.}
#'   \item{sd}{standard deviation of the alternative distribution.}
#'   \item{null.mean}{mean of the null distribution.}
#'   \item{null.sd}{standard deviation of the null distribution.}
#'   \item{z.alpha}{critical value(s).}
#'   \item{n}{sample size in the form of c(n1, n2) (applies to independent
#'            proportions).}
#'   \item{n.total}{total sample size (applies to independent proportions).}
#'   \item{n.paired}{paired sample size (applies to paired proportions).}
#'
#' @references
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#' @examples
#'   # power
#'   power.z.twoprops(prob1 = 0.65, prob2 = 0.60,
#'                    alpha = 0.05, n2 = 500,
#'                    alternative = "one.sided")
#'
#'   # sample size
#'   power.z.twoprops(prob1 = 0.65, prob2 = 0.60,
#'                    alpha = 0.05, power = 0.80,
#'                    alternative = "one.sided")
#'
#' @export power.z.twoprops
power.z.twoprops <- function(prob1 = NULL, prob2 = NULL, req.sign = "+", margin = 0,
                             n.ratio = 1, n2 = NULL, power = NULL, alpha = 0.05,
                             alternative = c("two.sided", "one.sided", "two.one.sided"),
                             arcsine = FALSE, correct = FALSE, paired = FALSE, rho.paired = 0.50,
                             std.error = c("pooled", "unpooled"),
                             ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  std.error <- match.arg(std.error)
  func.parms <- as.list(environment())

  if (!is.null(prob1)) check.proportion(prob1)
  if (!is.null(prob2)) check.proportion(prob2)
  if (!is.null(n2)) check.sample.size(n2)
  if (!is.null(power)) check.power(power)
  check.positive(n.ratio)
  check.proportion(alpha)
  check.logical(arcsine, correct, paired)
  check.correlation(rho.paired)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = list(prob1, prob2), n = n2, power = power)

  if (!is.numeric(margin) || any(margin > 0.99) || any(margin < -0.99))
    stop("Provide a reasonable `margin` consistent with `prob1` - `prob2`.", call. = FALSE)
  if (alternative == "two.one.sided" && paired)
      stop("Two one-sided tests are currently not available for paired proportions.", call. = FALSE)
  if (alternative == "two.one.sided" && arcsine)
      stop("Arcsine transformation is currently not available for two one-sided tests.", call. = FALSE)
  if (alternative == "two.one.sided" && length(margin) != 2)
      stop("Provide margins in the form of margin = c(lower, upper).", call. = FALSE)
  if (arcsine && correct)
      stop("Continuity correction does not apply to arcsine transformation approach.", call. = FALSE)
  if (arcsine && paired)
      stop("Arcsine transformation is currently not available for paired proportions.", call. = FALSE)
  if (arcsine && margin != 0)
      stop("Arcsine transformation is currently not available for non-zero null.", call. = FALSE)
  if (correct && paired)
    stop("Continuity correction is currently not available for paired proportions.", call. = FALSE)
  if (requested %in% c("power", "n") && any(abs((prob1 - prob2) - margin) < 1e-6))
    stop("The value of margin should be different from the prob1 - prob2 difference.", call. = FALSE)
  if (requested == "es" && margin != 0)
    warning("`margin` argument is ignored.", call. = FALSE)

  pwr <- function(prob1, prob2, margin, n2, n.ratio, arcsine,
                  std.error, correct, alpha, alternative) {

    n1 <- n.ratio * n2

    if (arcsine) {

      stat <- probs.to.h(prob1, prob2, FALSE)$h
      null.stat <- 0

      stderr <- sqrt(1 / n1 + 1 / n2)
      null.sd <- 1

    } else {

      k <- 0

      if (correct) {
        # k <- -sign(prob1 - prob2)
        # k.margin <- -sign(margin)
        # stat <- prob1 - prob2 - 0 + (k / 2) * (1 / n1 + 1 / n2)
        # null.stat <- margin + (k.margin / 2) * (1 / n1 + 1 / n2)
        k <- -sign(prob1 - prob2 - margin)
        stat <- prob1 - prob2 - 0 + (k / 2) * (1 / n1 + 1 / n2)
      } else {
        stat <- prob1 - prob2
      }  # if correct

      null.stat <- margin
      stderr <- sqrt(prob1 * (1 - prob1) / n1 + prob2 * (1 - prob2) / n2)
      if (std.error == "pooled") {
        p.bar <- (prob1 * n1 + prob2 * n2) / (n1 + n2)
        stderr.pooled <- sqrt(p.bar * (1 - p.bar) * (1 / n1 + 1 / n2))
        null.sd <- stderr.pooled / stderr
      } else {
        null.sd <- 1
      }

    } # if arcsine

    lambda <- stat / stderr
    null.lambda <- null.stat / stderr

    power.z.test(mean = lambda, sd = 1, null.mean = null.lambda, null.sd = null.sd,
                 alpha = alpha, alternative =  alternative, plot = FALSE, verbose = 0)

  } # pwr()


  ss.nocorrection.unpooled <- function(prob1, prob2, margin, power, n.ratio,
                                      arcsine, alpha, alternative) {

    beta <- 1 - power
    null.sd <- 1

    if (arcsine) {

      stat <- probs.to.h(prob1, prob2, FALSE)$h

    } else {

      stat <- prob1 - prob2 - margin

    } # if arcsine

    if (alternative == "one.sided") {

      z.alpha <- stats::qnorm(alpha, sd = null.sd, lower.tail = FALSE)
      z.beta <- stats::qnorm(beta, sd = null.sd, lower.tail = FALSE)

      if (arcsine) {

        n2 <- ((z.alpha + z.beta) ^ 2 / stat ^ 2) * (n.ratio + 1) / n.ratio
        n1 <- n.ratio * n2

      } else {

        n1 <- (z.alpha + z.beta) ^ 2 * (prob1 * (1 - prob1) + prob2 * (1 - prob2) * n.ratio) / stat ^ 2
        n2 <- n1 / n.ratio

      } # if arcsine

    } else if (alternative == "two.sided") {

      z.alpha <- stats::qnorm(alpha / 2, sd = null.sd, lower.tail = FALSE)
      z.beta <- stats::qnorm(beta, sd = null.sd, lower.tail = FALSE)

      if (arcsine) {

        n2 <- ((z.alpha + z.beta) ^ 2 / stat ^ 2) * (n.ratio + 1) / n.ratio
        n1 <- n.ratio * n2

      } else {

        n1 <- (z.alpha + z.beta) ^ 2 * (prob1 * (1 - prob1) + prob2 * (1 - prob2) * n.ratio) / stat ^ 2
        n2 <- n1 / n.ratio

      } # if arcsine

    } else if (alternative == "two.one.sided") {

      if (prob1 - prob2 > min(margin) && prob1 - prob2 < max(margin)) {
        z.alpha <- stats::qnorm(alpha, sd = null.sd, lower.tail = FALSE)
        z.beta <- stats::qnorm(beta / 2, sd = null.sd, lower.tail = FALSE)
      } else {
        z.alpha <- stats::qnorm(alpha / 2, sd = null.sd, lower.tail = FALSE)
        z.beta <- stats::qnorm(beta, sd = null.sd, lower.tail = FALSE)
      }

# arcsine currently not permitted (l. 803)
#      if (arcsine) {

#        n2 <- ((z.alpha + z.beta) ^ 2 / stat ^ 2) * (n.ratio + 1) / n.ratio
#        n2 <- max(n2)
#        n1 <- n.ratio * n2

#      } else {

      n1 <- (z.alpha + z.beta) ^ 2 * (prob1 * (1 - prob1) + prob2 * (1 - prob2) * n.ratio) / stat ^ 2
      n1 <- max(n1)
      n2 <- n1 / n.ratio

#      } # if arcsine

    } # if alt

    n2 <- n1 / n.ratio

    list(n1 = n1, n2 = n2)

  } # ss.no.correction()

  ss <- function(prob1, prob2, margin, n.ratio, alpha, power, alternative, arcsine, std.error, correct) {

    n2 <- floor(ss.nocorrection.unpooled(prob1 = prob1, prob2 = prob2, margin = margin, power = power, n.ratio = n.ratio,
                                         alpha = alpha, arcsine = arcsine, alternative = alternative)$n2)

    sign.switch <- 0
    sign.previous <- 0
    while (sign.switch < 2) {

      pwr.est <- pwr(prob1 = prob1, prob2 = prob2, margin = margin, n2 = n2, n.ratio = n.ratio, arcsine = arcsine,
                     std.error = std.error, correct = correct, alpha = alpha, alternative = alternative)$power

      sign.current <- sign(pwr.est - power)
      sign.switch <- sign.switch + as.integer(sign.current != 0 && sign.previous != 0 && sign.previous != sign.current)
      sign.previous <- sign.current
      n2 <- n2 + ifelse(sign.current > 0, -1, 1)

    } # while

    n2 + as.integer(sign.current > 0)

  } # ss()


  if (paired) {

    if (requested == "es") {

      pwr.exact <- function(prob1, prob2, n.ratio, n2, power, alpha, alternative, paired, rho.paired, method) {

        jp <- joint.probs.2x2(prob1 = prob1, prob2 = prob2, rho = rho.paired, verbose = 0)

        power.exact.mcnemar(prob10 = jp$prob10, prob01 = jp$prob01,
                            power = power, n.paired = n2, alpha = alpha,
                            alternative = alternative, method =  "approx",
                            ceil.n = FALSE, verbose = 0, utf = FALSE)$power

      } # power.exact()

      val.rng <- c(0.0001, ifelse(is.null(prob1), prob2, prob1), 0.9999)[ifelse(check.pos_sign(req.sign), -1, -3)]

      prob.lim.obj <- prob.limits.paired(prob1 = prob1, prob2 = prob2, rho = rho.paired)

      if (is.null(prob1)) {

        prob1 <- stats::optimize(
          f = function(prob1) {
            (power - pwr.exact(prob1 = prob1, prob2 = prob2, n.ratio = n.ratio, n2 = n2, power = NULL, alpha = alpha,
                               alternative = alternative, paired = paired, rho.paired = rho.paired)) ^ 2
          },
          maximum = FALSE, lower = max(val.rng[1], prob.lim.obj[[1]]), upper = min(val.rng[2], prob.lim.obj[[2]]))$minimum

      } else {

        prob2 <- stats::optimize(
          f = function(prob2) {
            (power - pwr.exact(prob1 = prob1, prob2 = prob2, n.ratio = n.ratio, n2 = n2, power = NULL, alpha = alpha,
                               alternative = alternative, paired = paired, rho.paired = rho.paired)) ^ 2
          },
          maximum = FALSE, lower = max(val.rng[1], prob.lim.obj[[1]]), upper = min(val.rng[2], prob.lim.obj[[2]]))$minimum

      } # prob1 or prob2?

      power <- NULL

    } # effect size

    # update or estimate power or sample size
    jp <- joint.probs.2x2(prob1 = prob1, prob2 = prob2, rho = rho.paired, verbose = 0)

    pwr.obj <- power.exact.mcnemar(prob10 = jp$prob10, prob01 = jp$prob01, power = power, n.paired = n2, alpha = alpha,
                                   method =  "approx", alternative = alternative, ceil.n = ceil.n, verbose = 0)
    pwr.obj[c("prob1", "prob2")] <- c(prob1, prob2)

    if (verbose > 0) {

      print.obj <- list(requested = requested,
                        tgt.effect = "prob10",
                        test = "Paired Proportions",
                        alpha = pwr.obj$alpha,
                        alternative = pwr.obj$alternative,
                        method = pwr.obj$test,
                        prob10 = jp$prob10,
                        prob01 = jp$prob01,
                        delta = pwr.obj$delta,
                        odds.ratio = pwr.obj$odds.ratio,
                        size = pwr.obj$size,
                        prob = pwr.obj$prob,
                        null.prob = pwr.obj$null.prob,
                        binom.alpha = pwr.obj$binom.alpha,
                        mean = pwr.obj$mean,
                        sd = pwr.obj$sd,
                        null.mean = pwr.obj$null.mean,
                        null.sd = pwr.obj$null.sd,
                        z.alpha = pwr.obj$z.alpha,
                        power = pwr.obj$power,
                        n.paired = pwr.obj$n.paired)

      .print.pwrss.mcnemar(print.obj, verbose = verbose, utf = utf)

    } # verbose

    invisible(pwr.obj)

  } else { # paired?

    if (requested == "n") {

      n2 <- ss(prob1 = prob1, prob2 = prob2, margin = margin, power = power, n.ratio = n.ratio, arcsine = arcsine,
               std.error = std.error, correct = correct, alpha = alpha, alternative = alternative)

    } else if (requested == "es") {

      if  (check.null_sign(req.sign, alternative) && is.null(prob1)) {
        val.rng <- c(min(prob2 + margin), max(prob2 + margin))
      } else if (check.null_sign(req.sign, alternative) && is.null(prob2)) {
        val.rng <- c(min(prob1 - margin), max(prob1 - margin))
      } else if (check.pos_sign(req.sign)) {
        val.rng <- c(ifelse(is.null(prob1), max(prob2 + margin), max(prob1 - margin)), 0.9999)
      } else { # req.sign neither null nor positive
        val.rng <- c(0.0001, ifelse(is.null(prob1), min(prob2 + margin), min(prob1 - margin)))
      }

      if (is.null(prob1)) {

        prob1 <- stats::optimize(
          f = function(prob1) {
            (power - pwr(prob1 = prob1, prob2 = prob2, margin = margin, n2 = n2, n.ratio = n.ratio, arcsine = arcsine,
                         std.error = std.error, correct = correct, alpha = alpha, alternative = alternative)$power) ^ 2
          },
          maximum = FALSE, interval = val.rng)$minimum

      } else {

        prob2 <- stats::optimize(
          f = function(prob2) {
            (power - pwr(prob1 = prob1, prob2 = prob2, margin = margin, n2 = n2, n.ratio = n.ratio, arcsine = arcsine,
                         std.error = std.error, correct = correct, alpha = alpha, alternative = alternative)$power) ^ 2
          },
          maximum = FALSE, interval = val.rng)$minimum

      } # prob1 or prob2?

    } # n, es

    n1 <- ifelse(ceil.n, ceiling(n2 * n.ratio), n2 * n.ratio)
    n.ratio <- n1 / n2
    n.total <- n1 + n2

    # calculate power (if requested == "power") or update it (if requested == "n")
    pwr.obj <- pwr(prob1 = prob1, prob2 = prob2, margin = margin, n2 = n2, n.ratio = n.ratio, arcsine = arcsine,
                   std.error = std.error, correct = correct, alpha = alpha, alternative = alternative)

    delta <- prob1 - prob2
    odds.ratio <- (prob1 / (1 - prob1)) /  (prob2 / (1 - prob2))

    if (verbose > 0) {

      print.obj <- list(requested = requested,
                        tgt.effect = ifelse(is.null(func.parms[["prob2"]]), "prob2", "prob1"),
                        test = "Independent Proportions",
                        alpha = alpha,
                        alternative = alternative,
                        method = "z",
                        prob1 = prob1,
                        prob2 = prob2,
                        delta = delta,
                        margin = margin,
                        odds.ratio = odds.ratio,
                        mean = pwr.obj$mean,
                        sd = pwr.obj$sd,
                        null.mean = pwr.obj$null.mean,
                        null.sd = pwr.obj$null.sd,
                        z.alpha = pwr.obj$z.alpha,
                        power = pwr.obj$power,
                        n = c(n1 = n1, n2 = n2),
                        n.total = n.total)

      .print.pwrss.fisher(print.obj, verbose = verbose, utf = utf)

    } # verbose

    invisible(structure(list(parms = func.parms,
                             test = "z",
                             prob1 = prob1,
                             prob2 = prob2,
                             delta = delta,
                             odds.ratio = odds.ratio,
                             mean = pwr.obj$mean,
                             sd = pwr.obj$sd,
                             null.mean = pwr.obj$null.mean,
                             null.sd = pwr.obj$null.sd,
                             z.alpha = pwr.obj$z.alpha,
                             power = pwr.obj$power,
                             n = c(n1 = n1, n2 = n2),
                             n.total = n.total),
                        class = c("pwrss", "z", "twoprops")))

  } # if paired

} # power.z.twoprops()

#' @rdname power.z.twoprops
#' @export power.z.twoprop
power.z.twoprop <- power.z.twoprops

#' @export pwrss.z.2props
pwrss.z.2props <- function(p1, p2, margin = 0, arcsin.trans = FALSE,
                           kappa = 1, alpha = 0.05,
                           alternative = c("not equal", "greater", "less",
                                           "equivalent", "non-inferior", "superior"),
                           n2 = NULL, power = NULL, verbose = TRUE) {

  alternative <- match.arg(alternative)
  verbose <- ensure.verbose(verbose)

  check.proportion(p1, p2, alpha)
  check.positive(kappa)
  check.numeric(margin)
  check.logical(arcsin.trans)

  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"

  if (alternative == "two.one.sided") margin <- c(-margin, margin)

  twoprops.obj <- power.z.twoprops(prob1 = p1, prob2 = p2,
                                   margin = margin,
                                   n2 = n2, n.ratio = kappa,
                                   power = power, alpha = alpha,
                                   alternative = alternative,
                                   arcsine = arcsin.trans,
                                   correct = FALSE, paired = FALSE,
                                   std.error = "pooled",
                                   ceil.n = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.z.twoprops() function. \n")

  invisible(twoprops.obj)

} # pwrss.z.2props()

#' @export pwrss.z.2props
pwrss.z.2prop <- pwrss.z.2props
