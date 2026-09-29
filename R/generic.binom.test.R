#' Power Analysis for the Generic Binomial Test
#'
#' @description
#' Calculates power or find Probability (Non-Centrality) for the generic
#' binomial test with (optional) Type 1 and Type 2 error plots.
#'
#'
#' @aliases power.binom
#'
#'
#' @param power       statistical power \eqn{(1 - \beta)}; either `power`,
#'                    `size`, or `prob` needs to be NULL (and is then
#'                    estimated).
#' @param size        number of trials (zero or more); either `power`, `size`,
#'                    or `prob` needs to be NULL (and is then estimated).
#' @param prob        probability of success on each trial under alternative;
#'                    either `power`, `size` or `prob` needs to be NULL (and is
#'                    then estimated).
#' @param null.prob   probability of success on each trial under null.
#' @param req.sign    whether 'prob' is expected to be greater '+1', less than
#'                    '-1', or within '0' the null.prob' bounds.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided", "one.sided", or "two.one.sided". For
#'                    non-inferiority or superiority tests, add or subtract the
#'                    margin from the null hypothesis value and use alternative
#'                    = "one.sided".
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
#'   \item{power}{statistical power \eqn{(1 - \beta)}.}
#'   \item{size}{number of trials (zero or more).}
#'   \item{prob}{probability of success on each trial under alternative.}
#'   \item{null.prob}{probability of success on each trial under null.}
#'   \item{alpha}{type 1 error rate.}
#'   \item{alternative}{direction or type of the hypothesis test.}
#'   \item{binom.alpha}{critical value(s).}
#'
#' @examples
#' # one-sided
#' power.binom.test(size = 200, prob = 0.6, null.prob = 0.5, alpha = 0.05,
#'                  alternative = "one.sided")
#' power.binom.test(power = 0.80, size = 200, req.sign = "+", null.prob = 0.5,
#'                  alpha = 0.05, alternative = "one.sided")
#'
#' # two-sided
#' power.binom.test(size = 200, prob = 0.4, null.prob = 0.5, alpha = 0.05,
#'                  alternative = "two.sided")
#' power.binom.test(power = 0.80, size = 200, req.sign = "+", null.prob = 0.5,
#'                  alpha = 0.05, alternative = "two.sided")
#'
#' # equivalence
#' power.binom.test(size = 200, prob = 0.5, null.prob = c(0.4, 0.6),
#'                  alpha = 0.05, alternative = "two.one.sided")
#' power.binom.test(power = 0.80, size = 200, req.sign = "0",
#'                  null.prob = c(0.4, 0.6), alpha = 0.05,
#'                  alternative = "two.one.sided")
#'
#' @export power.binom.test
power.binom.test <- function(power = NULL,
                             size = NULL,
                             prob = NULL,
                             null.prob = 0.5,
                             req.sign = "+",
                             alpha = 0.05,
                             alternative = c("two.sided", "one.sided", "two.one.sided"),
                             plot = TRUE,
                             verbose = 1,
                             utf = FALSE) {

  alternative <- match.arg(alternative)
  if (!is.null(power)) check.power(power)
  if (!is.null(size)) check.vector(size, check.size, 1)
  if (!is.null(prob)) check.proportion(prob)
  null.prob <- check.margins(null.prob, check.proportion, alternative)
  check.proportion(alpha)
  check.logical(plot, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = prob, n = size, power = power)

  pwr <- function(size = NULL, prob = NULL, null.prob = 0.5, alpha = 0.05, alternative) {

    # initialize variables to prevent "not found" errors
    approx.alpha <- NA
    binom.alpha  <- NA
    power        <- NA

    if (alternative == "two.sided") {

      q.low     <- stats::qbinom(alpha / 2, size, null.prob, lower.tail = TRUE)
      q.high    <- stats::qbinom(alpha / 2, size, null.prob, lower.tail = FALSE)
      prob.low  <- stats::pbinom(q.low,     size, null.prob, lower.tail = TRUE)
      prob.high <- stats::pbinom(q.high,    size, null.prob, lower.tail = FALSE)
      q.low[prob.low   > alpha / 2] <- q.low[prob.low   > alpha / 2] - 1
      q.high[prob.high > alpha / 2] <- q.high[prob.high > alpha / 2] + 1
      prob.low <- stats::pbinom(q.low, size, null.prob, lower.tail = TRUE)
      prob.high <- stats::pbinom(q.high, size, null.prob, lower.tail = FALSE)
      approx.alpha <-  prob.low + prob.high
      binom.alpha <- c(q.low, q.high)

      power <- stats::pbinom(q.low,  size, prob, lower.tail = TRUE) +
               stats::pbinom(q.high, size, prob, lower.tail = FALSE)

    } else if (alternative == "one.sided") {

      lower.tail <- prob < null.prob

      q          <- stats::qbinom(alpha, size, null.prob, lower.tail = lower.tail)
      prob.alpha <- stats::pbinom(q,     size, null.prob, lower.tail = lower.tail)
      q[prob.alpha > alpha] <- q[prob.alpha > alpha] + ifelse(lower.tail, -1, 1)
      approx.alpha <- stats::pbinom(q, size, null.prob, lower.tail = lower.tail)
      binom.alpha <- q

      power <- stats::pbinom(q, size, prob, lower.tail = lower.tail)

    } else if (alternative == "two.one.sided" && (prob > min(null.prob) && prob < max(null.prob))) {  # equivalence

      q.low     <- stats::qbinom(alpha,  size, prob = null.prob[1], lower.tail = FALSE)
      q.high    <- stats::qbinom(alpha,  size, prob = null.prob[2], lower.tail = TRUE)
      prob.low  <- stats::pbinom(q.low,  size, prob = null.prob[1], lower.tail = FALSE)
      prob.high <- stats::pbinom(q.high, size, prob = null.prob[2], lower.tail = TRUE)
      q.low[prob.low   > alpha] <- q.low[prob.low   > alpha] + 1
      q.high[prob.high > alpha] <- q.high[prob.high > alpha] - 1
      prob.low  <- stats::pbinom(q.low,  size, prob = null.prob[1], lower.tail = FALSE)
      prob.high <- stats::pbinom(q.high, size, prob = null.prob[2], lower.tail = TRUE)
      approx.alpha <-  (prob.low + prob.high) / 2 # average
      binom.alpha  <- c(q.low, q.high)

      power <- stats::pbinom(q.low, size, prob, lower.tail = FALSE) +
               stats::pbinom(q.high, size, prob, lower.tail = TRUE) - 1

    } else if (alternative == "two.one.sided" && (prob < min(null.prob) || prob > max(null.prob))) {  # minimal effect

      q.low     <- stats::qbinom(alpha / 2, size, prob = null.prob[1], lower.tail = TRUE)
      q.high    <- stats::qbinom(alpha / 2, size, prob = null.prob[2], lower.tail = FALSE)
      prob.low  <- stats::pbinom(q.low,     size, prob = null.prob[1], lower.tail = TRUE)
      prob.high <- stats::pbinom(q.high,    size, prob = null.prob[2], lower.tail = FALSE)
      q.low[prob.low   > alpha / 2] <- q.low[prob.low   > alpha / 2] - 1
      q.high[prob.high > alpha / 2] <- q.high[prob.high > alpha / 2] + 1
      prob.low  <- stats::pbinom(q.low,     size, prob = null.prob[1], lower.tail = TRUE)
      prob.high <- stats::pbinom(q.high,    size, prob = null.prob[2], lower.tail = FALSE)
      approx.alpha <- prob.low + prob.high
      binom.alpha  <- c(q.low, q.high)

      power <- stats::pbinom(q.low,  size, prob = prob, lower.tail = TRUE) +
               stats::pbinom(q.high, size, prob = prob, lower.tail = FALSE)

    }

    power[power < 0] <- 0

    list(power = power, approx.alpha = approx.alpha, binom.alpha = binom.alpha)

  } # pwr()

  min.pwr <- function(prob = NULL, size = NULL, power = NULL) {

    power - pwr(size = size, prob = prob, null.prob = null.prob, alpha = alpha, alternative = alternative)$power

  } # min.pwr() (for uniroot and optimize)

  if (requested == "n") {

    size <- round(stats::uniroot(f = function(size) min.pwr(prob, round(size), power), interval = c(1, 1e10))$root)

  } else if (requested == "es") {

    if (check.null_sign(req.sign, alternative)) {

      int.lower <- c(min(null.prob) + 1e-6, mean(null.prob))
      int.upper <- c(mean(null.prob), max(null.prob) - 1e-6)
      prob.lower <- stats::optimize(f = function(prob) min.pwr(prob, size, power) ^ 2, interval = int.lower, tol = 1e-12)$minimum
      prob.upper <- stats::optimize(f = function(prob) min.pwr(prob, size, power) ^ 2, interval = int.upper, tol = 1e-12)$minimum
      prob <- mean(c(prob.lower, prob.upper))

      warn.txt <- ifelse(max(abs(c(min.pwr(prob.lower, size, power), min.pwr(prob.upper, size, power)))) < 1e-6,
                         sprintf("Target NCP ranges from %.4f to %.4f within the null bounds.", prob.lower, prob.upper),
                         "The target power rate cannot be achieved within the null bounds.")
      warning(warn.txt, call. = FALSE)

    } else {

      val.rng <- get.interval(null.ncp = null.prob, distribution = "binom", alternative = alternative, req.sign = req.sign)
      prob <- stats::optimize(f = function(prob) min.pwr(prob, size, power) ^ 2, interval = val.rng, tol = 1e-12)$minimum

    }

  }

  # calculate power (if requested == "power") or update it (if requested == "n" or "es")
  pwr.obj <- pwr(size = size, prob = prob, null.prob = null.prob, alpha = alpha, alternative = alternative)

  if (plot)
    suppressWarnings(.plot.binom.t1t2(size = size, prob = prob, null.prob = null.prob, alpha = pwr.obj$approx.alpha,
                                      alternative = alternative))

  if (verbose > 0) {

    print.obj <- list(test = "Generic Binomial Test",
                      requested = requested,
                      tgt.ncp = "prob",
                      size = size,
                      prob = prob,
                      null.prob = null.prob,
                      alpha = pwr.obj$approx.alpha,
                      alternative = alternative,
                      binom.alpha = pwr.obj$binom.alpha,
                      power = pwr.obj$power)

    .print.pwrss.binom(print.obj, verbose = verbose, utf = utf)

  } # end of verbose

  invisible(structure(list(power = pwr.obj$power,
                           size = size,
                           prob = prob,
                           null.prob = null.prob,
                           alpha = pwr.obj$approx.alpha,
                           alternative = alternative,
                           binom.alpha = pwr.obj$binom.alpha),
                      class = c("pwrss", "generic", "binom")))

} # power.binom.test()

power.binom <- power.binom.test
