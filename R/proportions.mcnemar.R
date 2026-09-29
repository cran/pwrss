#' Power Analysis for McNemar's Exact Test (Paired Proportions)
#'
#' @description
#' Calculates power or sample size for McNemar's test on paired binary
#' outcomes. Approximate and exact methods are available (check references for
#' details).
#'
#' Validated using the PASS documentation and G*Power.
#'
#' @aliases power.exact.mcnemar power.exact.twoprops.mcnemar
#'
#'
#' @param prob10      (joint) probability of success in case (or after) but
#'                    failure in matched control (or before). 'prob10' and
#'                    'prob01' are known as discordant probs.
#' @param prob01      (joint) probability of failure in case (or after) but
#'                    success in matched control (or before). prob10' and
#'                    'prob01' are known as discordant probs.
#' @param req.sign    whether estimated prob is smaller or larger than the other
#'                    (when minimum detectable prob is of interest).
#' @param n.paired    number of pairs, which is sum of cell frequencies in the
#'                    2 x 2 table (f11 + f10 + f01 + f00), or number of rows in
#'                    a data frame with matched variables 'case' and 'control'
#'                    or 'after' and 'before'.
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param method      character; the method used for power calculation. "exact"
#'                    specifies Fisher's exact test, while "approximate" refers
#'                    to the z-test based on the normal approximation.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided" or "one.sided".
#' @param ceil.n      logical; if \code{TRUE} rounds up sample size in each
#'                    cell. This procedure assumes symmetry for concordant
#'                    probs, which are `p11` and `p00`). Thus results may
#'                    differ from other software by a few units. To match
#'                    results set `ceil.n` to FALSE.
#' @param verbose     \code{1} by default (returns test, hypotheses, and
#'                    results), if \code{2} a more detailed output is given
#'                    (plus key parameters and definitions), if \code{0} no
#'                    output is printed on the console.
#' @param utf         logical; whether the output should show Unicode
#'                    characters (if encoding allows for it). \code{FALSE}
#'                    by default.
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
#'   \item{n.paired}{paired sample size, which is sum of cell frequencies in
#'                   the 2 x 2 table (f11 + f10 + f01 + f00), or number of rows
#'                   in a data frame with variables 'case' and 'control' or
#'                   'after' and 'before'.}
#'
#' @references
#'   Bennett, B. M., & Underwood, R. E. (1970). 283. Note: On McNemar's Test
#'   for the 2 * 2 Table and Its Power Function. *Biometrics, 26*(2), 339-343.
#'   https://doi.org/10.2307/2529083
#'
#'   Connor, R. J. (1987). Sample size for testing differences in proportions
#'   for the paired-sample design. *Biometrics, 43*(1), 207-211.
#'   https://doi.org/10.2307/2531961
#'
#'   Duffy, S. W. (1984). Asymptotic and exact power for the McNemar test and
#'   its analogue with R controls per case. *Biometrics, 40*(4) 1005-1015.
#'   https://doi.org/10.2307/2531151
#'
#'   McNemar, Q. (1947). Note on the sampling error of the difference between
#'   correlated proportions or percentages. *Psychometrika, 12*(2), 153-157.
#'   https://doi.org/10.1007/BF02295996
#'
#'   Miettinen, O. S. (1968). The matched pairs design in the case of
#'   all-or-none responses. *Biometrics, 24*(2), 339-352.
#'   https://doi.org/10.2307/2528039
#'
#' @examples
#' # example data for a matched case-control design
#' # subject  case     control
#' # <int>    <dbl>    <dbl>
#' #   1        1        1
#' #   2        0        1
#' #   3        1        0
#' #   4        0        1
#' #   5        1        1
#' #   ...     ...      ...
#' #   100      0        0
#'
#' # example data for a before-after design
#' # subject  before   after
#' # <int>    <dbl>    <dbl>
#' #   1        1        1
#' #   2        0        1
#' #   3        1        0
#' #   4        0        1
#' #   5        1        1
#' #   ...     ...      ...
#' #   100      0        0
#'
#' # convert to a 2 x 2 frequency table
#' freqs <- matrix(c(30, 10, 20, 40), nrow = 2, ncol = 2)
#' colnames(freqs) <- c("control_1", "control_0")
#' rownames(freqs) <- c("case_1", "case_0")
#' freqs
#'
#' # convert to a 2 x 2 proportion table
#' props <- freqs / sum(freqs)
#' props
#'
#' # discordant pairs (0 and 1, or 1 and 0) in 'props' matrix
#' # are the sample estimates of prob01 and prob10
#'
#' # post-hoc exact power
#' power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10,
#'                     n.paired = 100, alpha = 0.05,
#'                     alternative = "two.sided", method = "exact")
#'
#' # required sample size for exact test
#' # assuming prob01 and prob10 are population parameters
#' power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10,
#'                     power = 0.80, alpha = 0.05,
#'                     alternative = "two.sided", method = "exact")
#'
#' # we may not have 2 x 2 joint probs
#' # convert marginal probs to joint probs
#' joint.probs.2x2(prob1 = 0.55, # mean of case group (or after)
#'                 prob2 = 0.45, # mean of matched control group (or before)
#'                 # correlation between matched case-control or before-after
#'                 rho = 0.4141414
#' )
#'
#' @export power.exact.mcnemar
power.exact.mcnemar <- function(prob10 = NULL, prob01 = NULL,
                                req.sign = "+", n.paired = NULL,
                                power = NULL,  alpha = 0.05,
                                alternative = c("two.sided", "one.sided"),
                                method = c("exact", "approximate"),
                                ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  method <- match.arg(method)
  func.parms <- as.list(environment())

  if (!is.null(prob10)) check.proportion(prob10)
  if (!is.null(prob01)) check.proportion(prob01)
  if (!is.null(n.paired)) check.sample.size(n.paired)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = list(prob01, prob10), n = n.paired, power = power)

  pwr.exact <- function(prob10, prob01, n.paired, alpha, alternative) {

    OR <- prob10 / prob01
    OR2prob <- c(1 / (1 + OR), OR / (1 + OR))

    prob <- if (alternative == "two.sided") OR2prob[1] else if (prob10 < prob01) min(OR2prob) else max(OR2prob)

    prod1 <- stats::dbinom(x = seq(0, ceiling(n.paired)), size = ceiling(n.paired), prob = prob01 + prob10)
    prod2 <- power.binom.test(prob = prob, null.prob = 0.50, size = seq(0, ceiling(n.paired)), alpha = alpha,
                              alternative = alternative, plot = FALSE, verbose = 0)$power

    list(power = sum(prod1 * prod2), mean = NA, sd = NA, null.mean = NA, null.sd = NA, z.alpha = NA)

  } # pwr.exact()

  ss.exact <- function(prob10, prob01, power, alpha, alternative) {

    n.paired <- ss.approx(prob10, prob01, power, alpha, alternative)
    achieved.power <- 0

    # sort(...) generates a vector with 10, 30, 100, 300, ... < n.paired
    if (n.paired < 1e5) {
      for (m in sort(vapply(10 ^ seq(floor(log10(n.paired))), function(p) p * c(1, 3), numeric(2)), decreasing = TRUE)) {
        if (n.paired > m) {
          step.prob <- round(m / 10)
          achieved.power <- 0
          while (achieved.power < power) {
            achieved.power <- pwr.exact(prob10 = prob10, prob01 = prob01, n.paired = n.paired,
                                        alpha = alpha, alternative = alternative)$power
            if (achieved.power < power) n.paired <- n.paired + step.prob
          } # while
          n.paired <- n.paired - ifelse(step.prob > 1, step.prob, 0)
        } # n.paired > m
      }
    }

    if (n.paired > 1e5) stop("Sample size exceeds 100,000. Please check the assumptions.", call. = FALSE)

    n.paired

  } #  ss.exact()

  pwr.approx <- function(prob10, prob01, n.paired, alpha, alternative) {

    OR <- prob10 / prob01
    PD <- prob10 + prob01

    # Machin, Campbell, Fayers, and Pinol (1997)
    z.alpha <- stats::qnorm(1 - (alpha / ifelse(alternative == "two.sided", 2, 1)), mean = 0, sd = 1, lower.tail = TRUE)
    z.beta <- (sqrt((OR - 1) ^ 2 * PD * n.paired) - z.alpha * (1 + OR)) / sqrt((OR + 1) ^ 2 - (OR - 1) ^ 2 * PD)
    power <- stats::pnorm(z.beta, mean = 0, sd = 1, lower.tail = TRUE)

    mean <- z.alpha + z.beta

    if (OR < 1) {
      mean <- -mean
      if (alternative == "one.sided") z.alpha <- -z.alpha
    }
    if (alternative == "two.sided") z.alpha <- c(-z.alpha, z.alpha)

    list(power = power, mean = mean, sd = 1, null.mean = 0, null.sd = 1, z.alpha = z.alpha)

  } # pwr.approx()

  ss.approx <- function(prob10, prob01, power, alpha, alternative) {

    OR <- prob10 / prob01
    PD <- prob10 + prob01
    beta <- 1 - power

    # Machin, Campbell, Fayers, and Pinol (1997)
    if (alternative == "two.sided") alpha <- alpha / 2
    z.alpha <- stats::qnorm(1 - alpha, mean = 0, sd = 1, lower.tail = TRUE)
    z.beta <- stats::qnorm(1 - beta, mean = 0, sd = 1, lower.tail = TRUE)
    n.paired <- (z.alpha * (1 + OR) + z.beta * sqrt((OR + 1) ^ 2 - (OR - 1) ^ 2 * PD)) ^ 2 / ((OR - 1) ^ 2 * PD)

    n.paired

  } # ss.approx()

  # for reasonable round-up
  prob11 <- (1 - (prob10 + prob01)) / 2
  prob00 <- (1 - (prob10 + prob01)) / 2

  # method
  if (method == "exact") {

    if (requested == "n") {

      n.paired <- ss.exact(prob10 = prob10, prob01 = prob01, power = power, alpha = alpha, alternative = alternative)
      n.paired <- ifelse(ceil.n, sum(ceiling(n.paired * c(prob11, prob10, prob01, prob00))), n.paired)

    } else if (requested == "es") {

      if (is.null(prob10)) {

        prob10 <- stats::optimize(
          f = function(prob10) {
            pwr.est <- pwr.exact(prob10 = prob10, prob01 = prob01, n.paired = n.paired,
                                 alpha = alpha, alternative = alternative)$power
            if (pwr.est == 1) 1 + max(prob10, prob01) else (power - pwr.est) ^ 2
          },
          lower = ifelse(check.pos_sign(req.sign),     prob01, 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 1 - prob01, min(prob01, 1 - prob01)))$minimum

      } else {

        prob01 <- stats::optimize(
          f = function(prob01) {
            pwr.est <- pwr.exact(prob10 = prob10, prob01 = prob01, n.paired = n.paired,
                                 alpha = alpha, alternative = alternative)$power
            if (pwr.est == 1) 1 + max(prob10, prob01) else (power - pwr.est) ^ 2
          },
          lower = ifelse(check.pos_sign(req.sign),     prob10, 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 1 - prob10, min(prob10, 1 - prob10)))$minimum

      } # prob10 or prob01?

    } # effect size

    # calculate power (if requested == "power") or update it (if requested == "n")
    pwr.obj <- pwr.exact(prob10 = prob10, prob01 = prob01, n.paired = n.paired, alpha = alpha, alternative = alternative)

  } else if (method == "approximate") {

    if (requested == "n") {

      n.paired <- ss.approx(prob10 = prob10, prob01 = prob01, power = power, alpha = alpha, alternative = alternative)
      n.paired <- ifelse(ceil.n, sum(ceiling(n.paired * c(prob11, prob10, prob01, prob00))), n.paired)

    } else if (requested == "es") {

      if (is.null(prob10)) {

        prob10 <- stats::optimize(
          f = function(prob10) {
            pwr.est <- pwr.approx(prob10 = prob10, prob01 = prob01, n.paired = n.paired,
                                  alpha = alpha, alternative = alternative)$power
            if (pwr.est == 1) 1 + max(prob10, prob01) else (power - pwr.est) ^ 2
          },
          lower = ifelse(check.pos_sign(req.sign),     prob01,                 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 1 - prob01, min(prob01, 1 - prob01)),
          tol = 1e-12)$minimum

      } else {

        prob01 <- stats::optimize(
          f = function(prob01) {
            pwr.est <- pwr.approx(prob10 = prob10, prob01 = prob01, n.paired = n.paired,
                                  alpha = alpha, alternative = alternative)$power
            if (pwr.est == 1) 1 + max(prob10, prob01) else (power - pwr.est) ^ 2
          },
          lower = ifelse(check.pos_sign(req.sign),     prob10,                 0.0001),
          upper = ifelse(check.pos_sign(req.sign), 1 - prob10, min(prob10, 1 - prob10)),
          tol = 1e-12)$minimum

      } # prob10 or prob01?

    } # effect size

    # calculate power (if requested == "power") or update it (if requested == "n")
    pwr.obj <- pwr.approx(prob10 = prob10, prob01 = prob01, n.paired = n.paired, alpha = alpha, alternative = alternative)

  }  # method

  # efficient non-centrality parameter for the chi-squared dist
  # f10 <- prob10 * n.paired
  # f01 <- prob01 * n.paired
  # ncp <- (f10 - f01) ^ 2 / (f10 + f01)
  # df <- 1

  # critical values for the binomial approach
  prob <- prob10 / (prob10 + prob01)
  null.prob <- 0.50
  size <- ceiling(sum(n.paired * c(prob10, prob01)))
  if (alternative == "one.sided") {
    if (prob < null.prob) {
      # less
      q.binom.alpha <- stats::qbinom(alpha, size = size, prob = null.prob, lower.tail = TRUE)
      p.binom.alpha <- stats::pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = TRUE)
      if (p.binom.alpha > alpha) q.binom.alpha <- q.binom.alpha - 1
      if (method == "exact") {
        approx.alpha <- stats::pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = TRUE)
      } else {
        approx.alpha <- alpha
      }
    } else {
      # greater
      q.binom.alpha <- stats::qbinom(alpha, size = size, prob = null.prob, lower.tail = FALSE)
      p.binom.alpha <- stats::pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = FALSE)
      if (p.binom.alpha > alpha) q.binom.alpha <- q.binom.alpha + 1
      if (method == "exact") {
        approx.alpha <- stats::pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = FALSE)
      } else {
        approx.alpha <- alpha
      }
    }
  } else if (alternative == "two.sided") {
    q.binom.lower <- stats::qbinom(alpha / 2, size = size, prob = null.prob, lower.tail = TRUE)
    p.binom.lower <- stats::pbinom(q.binom.lower, size = size, prob = null.prob, lower.tail = TRUE)
    if (p.binom.lower > alpha / 2) q.binom.lower <- q.binom.lower - 1
    p.binom.lower <- stats::pbinom(q.binom.lower, size = size, prob = null.prob, lower.tail = TRUE)

    q.binom.upper <- stats::qbinom(alpha / 2, size = size, prob = null.prob, lower.tail = FALSE)
    p.binom.upper <- stats::pbinom(q.binom.upper, size = size, prob = null.prob, lower.tail = FALSE)
    if (p.binom.upper > alpha / 2) q.binom.upper <- q.binom.upper + 1
    p.binom.upper <- stats::pbinom(q.binom.upper, size = size, prob = null.prob, lower.tail = FALSE)

    q.binom.alpha <- c(q.binom.lower, q.binom.upper)
    if (method == "exact") {
      approx.alpha <- p.binom.upper + p.binom.lower
    } else {
      approx.alpha <- alpha
    }
  }

  if (verbose > 0) {

    print.obj <- list(requested = requested,
                      tgt.effect = ifelse(is.null(func.parms[["prob01"]]), "prob01", "prob10"),
                      test = "Paired Proportions",
                      alpha = approx.alpha,
                      alternative = alternative,
                      method = ifelse(method == "exact", "exact", "z"),
                      prob10 = prob10,
                      prob01 = prob01,
                      delta = prob10 - prob01,
                      odds.ratio = prob10 / prob01,
                      size = size,
                      prob = prob10 / (prob10 + prob01),
                      null.prob = 0.50,
                      binom.alpha = q.binom.alpha,
                      mean = pwr.obj$mean,
                      sd = pwr.obj$sd,
                      null.mean = pwr.obj$null.mean,
                      null.sd = pwr.obj$null.sd,
                      z.alpha = pwr.obj$z.alpha,
                      power = pwr.obj$power,
                      n.paired = n.paired)

    .print.pwrss.mcnemar(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = ifelse(method == "exact", "exact", "z"),
                           prob10 = prob10,
                           prob01 = prob01,
                           delta = prob10 - prob01,
                           odds.ratio = prob10 / prob01,
                           size = size,
                           prob = prob10 / (prob10 + prob01),
                           null.prob = 0.50,
                           binom.alpha = q.binom.alpha,
                           mean = pwr.obj$mean,
                           sd = pwr.obj$sd,
                           null.mean = pwr.obj$null.mean,
                           null.sd = pwr.obj$null.sd,
                           z.alpha = pwr.obj$z.alpha,
                           alpha = approx.alpha,
                           power = pwr.obj$power,
                           n.paired = n.paired),
                      class = c("pwrss", ifelse(method == "exact", "exact", "z"), ifelse(method == "exact", "mcnemar", "twoprops"))))

} # power.exact.mcnemar()

#' @rdname power.exact.mcnemar
#' @export power.exact.twoprops.mcnemar
power.exact.twoprops.mcnemar <- power.exact.mcnemar
