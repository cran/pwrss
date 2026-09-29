#' Statistical Power for the Generic z-Test
#'
#' @description
#' Determines the power or the non-centrality parameter (mean) for the generic
#' z-Test with (optional) Type 1 and Type 2 error plots.
#'
#' @aliases power.z
#'
#' @param power       statistical power \eqn{(1-\beta)}; either `power`, or
#'                    `mean` needs to be NULL (and is then estimated).
#' @param mean        mean of the alternative; either `power`, or `mean` needs
#'                    to be NULL (and is then estimated).
#' @param sd          standard deviation of the alternative. Do not change this
#'                    value except when some sort of variance correction is
#'                    applied (e.g. as in logistic and Poisson regressions).
#' @param null.mean   mean of the null. When alternative = "two.one.sided", the
#'                    function expects two values in the form `c(lower,
#'                    upper)`. If a single value is provided, it is interpreted
#'                    as the absolute bound and automatically expanded to
#'                    `c(-value, +value)`.
#' @param null.sd     standard deviation of the null. Do not change this value
#'                    except when some sort of correction is applied.
#' @param req.sign    whether `mean` is expected to be greater '+1', less than
#'                    '-1', or within '0' the `null.mean` bounds; only relevant
#'                    if `mean` is to be estimated.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided", "one.sided", or "two.one.sided".
#'                    "two.one.sided" is used for equivalence and minimal
#'                    effect testing.
#' @param plot        logical; \code{FALSE} switches off Type 1 and Type 2
#'                    error plot. \code{TRUE} by default.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @param utf         logical; whether the output should show Unicode
#'                    characters (if encoding allows for it). \code{FALSE} by
#'                    default.
#'
#' @return
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{mean}{mean of the alternative distribution.}
#'   \item{sd}{standard deviation of the alternative distribution.}
#'   \item{null.mean}{mean of the null distribution.}
#'   \item{null.sd}{standard deviation of the null distribution.}
#'   \item{alpha}{type 1 error rate (user-specified).}
#'   \item{z.alpha}{critical value(s).}
#'   \item{beta}{type 2 error rate.}
#'   \item{type.s}{type S error rate (only for two-tailed test).}
#'   \item{type.m}{type M error rate (only for two-tailed test).}
#'
#' @examples
#' # two-sided
#' # power defined as the probability of observing test statistics greater than
#' # the positive critical value OR less than the negative critical value
#' power.z.test(mean = 1.96, alpha = 0.05, alternative = "two.sided")
#' power.z.test(power = 0.80, alpha = 0.05, alternative = "two.sided")
#'
#' # one-sided
#' # power is defined as the probability of observing a test statistic greater
#' # than the critical value
#' power.z.test(mean = 1.96, alpha = 0.05, alternative = "one.sided")
#' power.z.test(power = 0.80, alpha = 0.05, alternative = "one.sided")
#'
#' # equivalence
#' # power is defined as the probability of observing a test statistic greater
#' # than the upper critical value (for the lower bound) AND less than the
#' # lower critical value (for the upper bound)
#' power.z.test(mean = 0, null.mean = c(-2, 2), alpha = 0.05,
#'              alternative = "two.one.sided")
#' power.z.test(power = 0.80, req.sign = "0", null.mean = c(-2, 2),
#'              alpha = 0.05, alternative = "two.one.sided")
#'
#' # minimal effect testing
#' # power is defined as the probability of observing a test statistic greater
#' # than the upper critical value (for the upper bound) OR less than the lower
#' # critical value (for the lower bound).
#' power.z.test(mean = 2, null.mean = c(-1, 1), alpha = 0.05,
#'              alternative = "two.one.sided")
#' power.z.test(power = 0.80, req.sign = "+", null.mean = c(-1, 1),
#'              alpha = 0.05, alternative = "two.one.sided")
#'
#' @export power.z.test
power.z.test <- function(power = NULL, mean = NULL, sd = 1, null.mean = 0, null.sd = 1, req.sign = "+",
                         alpha = 0.05, alternative = c("two.sided", "one.sided", "two.one.sided"),
                         plot = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  if (!is.null(power)) check.power(power)
  if (!is.null(mean)) check.numeric(mean)
  check.positive(sd)
  null.mean <- check.margins(null.mean, check.numeric, alternative)
  check.nonnegative(null.sd)
  check.proportion(alpha)
  check.logical(plot, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = mean, n = NA, power = power)

  # calculate statistical power
  pwr <- function(mean, sd = 1, null.mean = 0, null.sd = 1, alpha = 0.05, alternative) {

    if (alternative == "two.sided") {

      z.alpha <- c(stats::qnorm(alpha / 2,  mean = null.mean,      sd = null.sd, lower.tail = TRUE),
                   stats::qnorm(alpha / 2,  mean = null.mean,      sd = null.sd, lower.tail = FALSE))
      power   <-   stats::pnorm(z.alpha[1], mean = mean,           sd = sd,      lower.tail = TRUE) +
                   stats::pnorm(z.alpha[2], mean = mean,           sd = sd,      lower.tail = FALSE)

      Phi.p <- stats::pnorm(q = max(z.alpha), mean = mean, sd = sd)
      Phi.m <- stats::pnorm(q = min(z.alpha), mean = mean, sd = sd)
      type.s <- min(Phi.m, 1 - Phi.p) / (Phi.m + 1 - Phi.p)

      phi.p <- stats::dnorm(x = max(z.alpha), mean = mean, sd = sd)
      phi.m <- stats::dnorm(x = min(z.alpha), mean = mean, sd = sd)
      Phi.p <- stats::pnorm(q = max(z.alpha), mean = mean, sd = sd)
      Phi.m <- stats::pnorm(q = min(z.alpha), mean = mean, sd = sd)
      type.m <- (sd ^ 2 * (phi.p + phi.m) + mean * (1 - Phi.p - Phi.m)) / (abs(mean) * (1 - Phi.p + Phi.m))

    } else if (alternative == "one.sided") {

      lower.tail <- mean < null.mean
      z.alpha <- stats::qnorm(alpha,        mean = null.mean,      sd = null.sd, lower.tail = lower.tail)
      power   <- stats::pnorm(z.alpha,      mean = mean,           sd = sd,      lower.tail = lower.tail)

      type.s <- 0
      type.m <- NA

    } else if (alternative == "two.one.sided" && (mean > min(null.mean) && mean < max(null.mean))) {  # equivalence test

      z.alpha <- c(stats::qnorm(alpha,      mean = min(null.mean), sd = null.sd, lower.tail = FALSE),
                   stats::qnorm(alpha,      mean = max(null.mean), sd = null.sd, lower.tail = TRUE))
      power   <-   stats::pnorm(z.alpha[2], mean = mean,           sd = sd,      lower.tail = TRUE) +
                   stats::pnorm(z.alpha[1], mean = mean,           sd = sd,      lower.tail = FALSE) - 1

      type.s <- NA
      type.m <- NA

    } else if (alternative == "two.one.sided" && (mean < min(null.mean) || mean > max(null.mean))) {  # minimum effect test

      z.alpha <- c(stats::qnorm(alpha / 2,  mean = min(null.mean), sd = null.sd, lower.tail = TRUE),
                   stats::qnorm(alpha / 2,  mean = max(null.mean), sd = null.sd, lower.tail = FALSE))
      power   <-   stats::pnorm(z.alpha[1], mean = mean,           sd = sd,      lower.tail = TRUE) +
                   stats::pnorm(z.alpha[2], mean = mean,           sd = sd,      lower.tail = FALSE)

      type.s <- NA
      type.m <- NA

    }

    power[power < 0] <- 0

    list(power = power, z.alpha = z.alpha, type.s = type.s, type.m = type.m)

  } # pwr()

  min.pwr <- function(mean, power) {

    power - pwr(mean = mean, sd = sd, null.mean = null.mean, null.sd = null.sd, alpha = alpha, alternative = alternative)$power

  } # min.pwr() (for uniroot and optimize)

  if (requested == "es") {

    if (check.null_sign(req.sign, alternative)) {

      lower.int <- c(min(null.mean), mean(null.mean))
      upper.int <- c(mean(null.mean), max(null.mean))
      mean.lower <- stats::optimize(f = function(mean) min.pwr(mean, power) ^ 2, interval = lower.int, tol = 1e-12)$minimum
      mean.upper <- stats::optimize(f = function(mean) min.pwr(mean, power) ^ 2, interval = upper.int, tol = 1e-12)$minimum
      mean <- mean(c(mean.lower, mean.upper))

      warn.txt <- ifelse(max(abs(c(min.pwr(mean.lower, power), min.pwr(mean.upper, power)))) < 1e-6,
                         sprintf("Target NCP ranges from %.4f to %.4f within the null bounds.", mean.lower, mean.upper),
                         "The target power rate cannot be achieved within the null bounds.")
      warning(warn.txt, call. = FALSE)

    } else {

      val.rng <- get.interval(null.ncp = null.mean, distribution = "z", alpha = alpha, alternative = alternative,
                              req.sign = req.sign, sd = sd)
      mean <- stats::optimize(f = function(mean) min.pwr(mean, power) ^ 2, interval = val.rng, tol = 1e-12)$minimum

    }

  }

  pwr.obj <- pwr(mean = mean, sd = sd, null.mean = null.mean, null.sd = null.sd, alpha = alpha, alternative = alternative)

  if (plot) {
    # if (sd != 1 || null.sd != 1)
    #   stop("Plotting is not available when the standard deviation of the standard normal distribution deviates from one.", call. = FALSE)

    suppressWarnings(.plot.t.t1t2(ncp = mean, null.ncp = null.mean, df = Inf, alpha = alpha, alternative = alternative))

  }

  if (verbose > 0) {

    print.obj <- list(test = "Generic z-Test",
                      requested = requested,
                      tgt.ncp = "mean",
                      mean = mean,
                      sd = sd,
                      null.mean = null.mean,
                      null.sd = null.sd,
                      alternative = alternative,
                      alpha = alpha,
                      z.alpha = pwr.obj$z.alpha,
                      power = pwr.obj$power)

    .print.pwrss.z(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(power = pwr.obj$power, mean = mean, sd = sd, null.mean = null.mean, null.sd = null.sd, alpha = alpha,
                           alternative = alternative, z.alpha = pwr.obj$z.alpha, beta = 1 - pwr.obj$power,
                           type.s = pwr.obj$type.s, type.m = pwr.obj$type.m),
                      class = c("pwrss", "generic", "z")))

} # end of power.z.test()

power.z <- power.z.test
