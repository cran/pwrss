# type 1 and type 2 error plots are not available for this function

#' Power Analysis for Fisher's Exact Test (Independent Proportions)
#'
#' @description
#' Calculates power or sample size for Fisher's exact test on independent
#' binary outcomes. Approximate and exact methods are available.
#'
#' Validated using the PASS documentation and G*Power.
#'
#' @aliases power.exact.fisher power.exact.twoprops.fisher
#'
#'
#' @param prob1       probability of success in the first group.
#' @param prob2       probability of success in the second group.
#' @param req.sign    whether estimated prob is smaller or larger than the other
#'                    (when minimum detectable prob is of interest).
#' @param n2          integer; sample size for the second group.
#' @param n.ratio     n1 / n2 ratio.
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided" or "one.sided".
#' @param method      character; method used for power calculation. "exact"
#'                    specifies Fisher's exact test, while "approximate" refers
#'                    to the Z-Test based on the normal approximation.
#' @param ceil.n      logical; if \code{TRUE} rounds up sample size in each
#'                    group.
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
#'   \item{test}{type of the test, which is "exact" or "z".}
#'   \item{odds.ratio}{odds ratio.}
#'   \item{mean}{mean of the alternative distribution.}
#'   \item{sd}{standard deviation of the alternative distribution.}
#'   \item{null.mean}{mean of the null distribution.}
#'   \item{null.sd}{standard deviation of the null distribution.}
#'   \item{z.alpha}{critical value(s).}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample sizes for the first and second groups, specified as
#'            c(n1, n2).}
#'   \item{n.total}{total sample size, which is sum of cell frequencies in the
#'                  2 x 2 table (f11 + f10 + f01 + f00), or number of rows in a
#'                  data frame with group variable stacked.}
#'
#' @references
#'   Bennett, B. M., & Hsu, P. (1960). On the power function of the exact test
#'   for the 2 x 2 contingency table. *Biometrika, 47*(3/4), 393-398.
#'   https://doi.org/10.2307/2333309
#'
#'   Fisher, R. A. (1935). The logic of inductive inference. *Journal of the
#'   Royal Statistical Society, 98*(1), 39-82.
#'   https://doi.org/10.2307/2342435
#'
#' @examples
#' # example data for a randomized controlled trial
#' # subject  group    success
#' # <int>    <dbl>      <dbl>
#' #   1        1          1
#' #   2        0          1
#' #   3        1          0
#' #   4        0          1
#' #   5        1          1
#' #   ...     ...        ...
#' #   100      0          0
#'
#' # prob1 = mean(success | group = 1)
#' # prob2 = mean(success | group = 0)
#'
#' # post-hoc exact power
#' power.exact.fisher(prob1 = 0.60, prob2 = 0.40, n2 = 50)
#'
#' # we may have 2 x 2 joint probs such as
#' # -------------------------------------
#' #             | group (1) | group (0) |
#' # -------------------------------------
#' # success (1) |  0.24    |   0.36     |
#' # -------------------------------------
#' # success (0) |  0.16    |   0.24     |
#' # -------------------------------------
#'
#' # convert joint probs to marginal probs
#' marginal.probs.2x2(prob11 = 0.24, prob10 = 0.36,
#'                    prob01 = 0.16, prob00 = 0.24)
#'
#' @export power.exact.fisher
power.exact.fisher <- function(prob1 = NULL, prob2 = NULL, req.sign = "+",
                               n.ratio = 1, n2 = NULL,
                               power = NULL, alpha = 0.05,
                               alternative = c("two.sided", "one.sided"),
                               method = c("exact", "approximate"),
                               ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  method <- match.arg(method)
  func.parms <- as.list(environment())

  if (!is.null(prob1)) check.proportion(prob1)
  if (!is.null(prob2)) check.proportion(prob2)
  if (!is.null(n2))    check.sample.size(n2)
  if (!is.null(power)) check.power(power)
  check.positive(n.ratio)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = list(prob1, prob2), n = n2, power = power)

  pwr.approx <- function(prob1, prob2, n2, n.ratio, alpha, alternative,
                         pooled.stderr = FALSE, correct.continuity = FALSE) {

    n1 <- n.ratio * n2
    # Gpower
    # sigma0 <- sqrt(((n1 * (1 - prob1) + n2 * (1 - prob2)) / (n1 * n2)) * ((n1 * prob1 + n2 * prob2) / (n1 + n2)))
    # stderr <- (1 / sigma0) * sqrt((prob1 * (1 - prob1)) / n1 + (prob2 * (1 - prob2)) / n2)
    # delta <- (1 / sigma0) * (prob1 - prob2 - (k / 2) * (1 / n1 + 1 / n2)) # w/ continuity correction

#    if (pooled.stderr) {

#      p.bar <- (n1 * prob1 + n2 * prob2) / (n1 + n2)
#      stderr <- sqrt(p.bar * (1 - p.bar) * (1 / n1 + 1 / n2))

#    } else {

      stderr <- sqrt((prob1 * (1 - prob1)) / n1 + (prob2 * (1 - prob2)) / n2)

#    }

#    if (correct.continuity) {

#      if (alternative %in% c("not equal", "two.sided")) {
#        k <- c(-1, 1)
#      } else {
#        k <- ifelse(prob1 < prob2, -1, 1)
#      }
#      delta <- (prob1 - prob2 - (k / 2) * (1 / n1 + 1 / n2))

#    } else {

      delta <- prob1 - prob2

#    }

    power.z.test(mean = delta / stderr, sd = 1, null.mean = 0, null.sd = 1,
                 alpha = alpha, alternative = alternative, plot = FALSE, verbose = 0)

  } # pwr.approx


  ss.approx <- function(prob1, prob2, power, n.ratio,
                        alpha, alternative,
                        pooled.stderr = FALSE) {

    n2 <- stats::uniroot(function(n2) {
      power - pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                         alpha = alpha, alternative = alternative,
                         pooled.stderr = pooled.stderr)$power
    }, interval = c(2, 1e10))$root

    n2

  } # ss.approx()


  pwr.exact <- function(prob1, prob2, n2, n.ratio, alpha, alternative) {

    if (n2 > 2000) {

      power <- pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                      alpha = alpha, alternative = alternative)$power

      n.total <- 2000 + ceiling(n.ratio * 2000)

      stop(sprintf("Consider `method` = 'approximate' for total sample size > %d.", n.total), call. = FALSE)

    } else {

      eps <- 1e-10

      n1 <- round(n.ratio * n2)
      x1.seq <- rep(0:n1, times = n2 + 1)
      x2.seq <- rep(0:n2, each = n1 + 1)
      m <- x1.seq + x2.seq
      n.total <- n1 + n2
      k <- n1

      if (alternative == "one.sided") {

        joint.probs <- stats::dbinom(x1.seq, n1, prob1) * stats::dbinom(x2.seq, n2, prob2)
        p.values <- stats::phyper(x1.seq - as.integer(prob1 >= prob2), m, n.total - m, k, lower.tail = prob1 < prob2)
        reject <- !is.na(p.values) & p.values <= alpha
        power <- sum(joint.probs[reject])

      } else if (alternative == "two.sided") {

        joint.probs <- stats::dbinom(x1.seq, n1, prob1) * stats::dbinom(x2.seq, n2, prob2)
        valid <- joint.probs > 0
        x1.seq <- x1.seq[valid]
        x2.seq <- x2.seq[valid]
        joint.probs <- joint.probs[valid]
        m.seq <- x1.seq + x2.seq

        m.unique <- unique(m.seq)
        pmf.lookup <- lapply(m.unique, function(m) {
          support <- max(0, m - n2):min(n1, m)
          probs <- stats::dhyper(support, m, n.total - m, n1)
          list(support = support, probs = probs)
        })
        names(pmf.lookup) <- as.character(m.unique)

        passed <- vapply(seq_along(x1.seq), function(i) {
          x1 <- x1.seq[i]
          m  <- m.seq[i]
          key <- as.character(m)
          lookup <- pmf.lookup[[key]]
          p.obs <- stats::dhyper(x1, m, n.total - m, n1)
          p.value <- sum(lookup$probs[lookup$probs <= (p.obs + eps)])
          !is.na(p.value) && p.value <= alpha
        }, logical(1))

        power <- sum(joint.probs[passed])

      }

    } # if n2 > 2000

    list(power = power, mean = NA, sd = NA, null.mean = NA, null.sd = NA, z.alpha = NA) # to match return from pwr.approx

  } # pwr.exact()


  ss.exact <- function(prob1, prob2, power, n.ratio, alpha, alternative) {

    n2 <- ss.approx(prob1 = prob1, prob2 = prob2, power = power, n.ratio = n.ratio,
                    alpha = alpha, alternative = alternative)
    n2 <- ceiling(n2)

    if (n2 > 500) {
      n.total <- 500 + ceiling(n.ratio * 500)
      stop("Consider `method` = 'approximate' for total sample size > ", n.total, call. = FALSE)
    }

    # power is approached in successively smaller steps (pwr.exact is time-consuming to calculate):
    # 30 (if n2 > 300), 10 (if n2 > 100), 3 (if n2 > 30), and 1 (in any case; for the exact n2)
    for (step in c(30, 10, 3, 1)[c(n2 >= 300, n2 >= 100, n2 >= 30, TRUE)]) {
      achieved.power <- 0
      while (achieved.power < power) {
        achieved.power <- pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                                    alpha = alpha, alternative = alternative)$power
        if (achieved.power < power) n2 <- n2 + step
      } # while
      # step is subtracted from n2 (for the larger steps), so that n2 can be approached with smaller steps
      n2 <- n2 - ifelse(step > 1, step, 0)
    }

    n2
  } #  ss.exact()


  # method
  if (method == "exact") {

    if (requested == "n") {

      n2 <- ss.exact(prob1 = prob1, prob2 = prob2, power = power, n.ratio = n.ratio,
                     alpha = alpha, alternative = alternative)
      n1 <- n.ratio * n2

      if (ceil.n) {
        n1 <- ceiling(n1)
        n2 <- ceiling(n2)
      }

    } else if (requested == "power") {

      n1 <- ifelse(ceil.n, ceiling(n.ratio * n2), n.ratio * n2)

    } else if (requested == "es") {

      if (is.null(prob1)) {

        prob1 <- stats::optimize(
          f = function(prob1) {
            (power - pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                               alpha = alpha, alternative = alternative)$power) ^ 2
          },
          maximum = FALSE,
          lower = ifelse(check.pos_sign(req.sign),  prob2, 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 0.9999,  prob2))$minimum

      } else {

        prob2 <- stats::optimize(
          f = function(prob2) {
            (power - pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                               alpha = alpha, alternative = alternative)$power) ^ 2
          },
          maximum = FALSE,
          lower = ifelse(check.pos_sign(req.sign),  prob1, 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 0.9999,  prob1))$minimum

      } # prob1 or prob2?

      n1 <- ifelse(ceil.n, ceiling(n.ratio * n2), n.ratio * n2)

    } # n, power, es?

    n.total <- n1 + n2

    # calculate power (if requested == "power") or update it (if requested == "n")
    pwr.obj <- pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                         alpha = alpha, alternative = alternative)

  } else if (method == "approximate") {

    if (requested == "n") {

      n2 <- ss.approx(prob1 = prob1, prob2 = prob2, power = power, n.ratio = n.ratio,
                      alpha = alpha, alternative = alternative)
      n2 <- ifelse(ceil.n, ceiling(n2), n2)

    } else if (requested == "es") {

      if (is.null(prob1)) {

        prob1 <- stats::optimize(
          f = function(prob1) {
            (power - pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                               alpha = alpha, alternative = alternative)$power) ^ 2
          },
          maximum = FALSE,
          lower = ifelse(check.pos_sign(req.sign),  prob2, 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 0.9999,  prob2))$minimum

      } else {

        prob2 <- stats::optimize(
          f = function(prob2) {
            (power - pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                               alpha = alpha, alternative = alternative)$power) ^ 2
          },
          maximum = FALSE,
          lower = ifelse(check.pos_sign(req.sign),  prob1, 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 0.9999,  prob1))$minimum

      } # prob1 or prob2?

    } # n, power, es?

    n1 <- ifelse(ceil.n, ceiling(n.ratio * n2), n.ratio * n2)
    n.total <- n1 + n2

    # calculate power (if requested == "power") or update it (if requested == "n")
    pwr.obj <- pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                          alpha = alpha, alternative = alternative)

  }  # method

  if (verbose > 0) {

    print.obj <- list(requested = requested,
                      tgt.effect = ifelse(is.null(func.parms[["prob2"]]), "prob2", "prob1"),
                      test = "Independent Proportions",
                      alpha = alpha,
                      alternative = alternative,
                      method = ifelse(method == "exact", "exact", "z"),
                      prob1 = prob1,
                      prob2 = prob2,
                      delta = prob1 - prob2,
                      margin = 0,
                      odds.ratio = (prob1 / (1 - prob1)) /  (prob2 / (1 - prob2)),
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
                           test = ifelse(method == "exact", "exact", "z"),
                           prob1 = prob1,
                           prob2 = prob2,
                           delta = prob1 - prob2,
                           odds.ratio = (prob1 / (1 - prob1)) /  (prob2 / (1 - prob2)),
                           mean = pwr.obj$mean,
                           sd = pwr.obj$sd,
                           null.mean = pwr.obj$null.mean,
                           null.sd = pwr.obj$null.sd,
                           alternative = alternative,
                           z.alpha = pwr.obj$z.alpha,
                           power = pwr.obj$power,
                           n = c(n1 = n1, n2 = n2),
                           n.total = n.total),
                      class = c("pwrss", ifelse(method == "exact", "exact", "z"), ifelse(method == "exact", "fisher", "twoprops"))))

} # power.exact.fisher()

#' @rdname power.exact.fisher
#' @export power.exact.twoprops.fisher
power.exact.twoprops.fisher <- power.exact.fisher
