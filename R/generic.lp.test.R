#' Statistical Power for the Lambda-Prime Distribution
#'
#' @description
#' Determines the power, the non-centrality parameter, or the degrees of
#' freedom for the lambda-prime distribution with (optional) Type 1 and Type 2
#' error plots.
#'
#' @aliases power.lp
#'
#'
#' @param power       statistical power \eqn{(1 - \beta)}; either `power`,
#'                    `ncp` or `df` needs to be NULL (and is then estimated).
#' @param ncp         non-centrality parameter for the alternative; either
#'                    `power`, `ncp` or `df` needs to be NULL (and is then
#'                    estimated).
#' @param null.ncp    non-centrality parameter for the null. When alternative =
#'                    "two.one.sided", the function expects two values in the
#'                    form `c(lower, upper)`. If a single value is provided, it
#'                    is interpreted as the absolute bound and automatically
#'                    expanded to `c(-value, +value)`.
#' @param req.sign    whether `ncp` is expected to be greater '+1', less than
#'                    '-1', or within '0' the `null.ncp` bounds; only relevant
#'                    if `ncp` is to be estimated.
#' @param df          degrees of freedom; either `power`, `ncp` or `df` needs
#'                    to be NULL (and is then estimated).
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
#'   \item{power}{statistical power \eqn{(1 - \beta)}.}
#'   \item{ncp}{non-centrality parameter under alternative.}
#'   \item{null.ncp}{non-centrality parameter under null.}
#'   \item{df}{degrees of freedom.}
#'   \item{alpha}{type 1 error rate (user-specified).}
#'   \item{alternative}{the direction or type of the hypothesis test.}
#'   \item{t.alpha}{critical value(s).}
#'   \item{beta}{type 2 error rate.}
#'   \item{type.s}{type S error rate (only for two-tailed test).}
#'   \item{type.m}{type M error rate (only for two-tailed test).}
#'
#' @examples
#' # two-sided
#' # power defined as the probability of observing test statistics greater
#' # than the positive critical value OR less than the negative critical value
#' power.lp.test(ncp = 1.960, df = 100, alpha = 0.05,
#'               alternative = "two.sided", plot = FALSE)
#' power.lp.test(power = 0.800, df = 100, alpha = 0.05,
#'               alternative = "two.sided", plot = FALSE)
#'
#' # the two examples below estimate the df's based upon the first example
#' # (revealing a power of 0.498; df = 94.11) and the second example (revealing
#' # a ncp of 2.825; df = 101.06)
#' power.lp.test(ncp = 1.960, power = 0.498, alpha = 0.05,
#'               alternative = "two.sided", plot = FALSE)
#' power.lp.test(ncp = 2.825, power = 0.800, alpha = 0.05,
#'               alternative = "two.sided", plot = FALSE)
#'
#' # one-sided
#' # power is defined as the probability of observing a test statistic greater
#' # than the critical value
#' power.lp.test(ncp = 1.960, df = 100, alpha = 0.05, alternative = "one.sided")
#' power.lp.test(power = 0.800, df = 100, alpha = 0.05, alternative = "one.sided")
#' # the two examples below estimate the df's based upon the first example
#' # (revealing a power of 0.6207; df = 100.323) and the second example (revealing
#' # a ncp of 2.506; df = 99.12)
#' power.lp.test(ncp = 1.960, power = 0.6207, alpha = 0.05,
#'               alternative = "one.sided", plot = FALSE)
#' power.lp.test(ncp = 2.506, power = 0.8000, alpha = 0.05,
#'               alternative = "one.sided", plot = FALSE)
#'
#' # equivalence
#' # power is defined as the probability of observing a test statistic greater
#' # than the upper critical value (for the lower bound) AND less than the
#' # lower critical value (for the upper bound)
#' power.lp.test(ncp = 0, null.ncp = c(-3, 3), df = 100, alpha = 0.05,
#'               alternative = "two.one.sided", plot = FALSE)
#' power.lp.test(power = 0.80, req.sign = "0", null.ncp = c(-3, 3),
#'               df = 100, alpha = 0.05, alternative = "two.one.sided", plot = FALSE)
#' # adjust the power based upon what is returned from the example above in
#' # order to get a valid estimate of the df's (100.321; power = 0.8 -> 58.911)
#' power.lp.test(ncp = 0, power = 0.8103, req.sign = "0", null.ncp = c(-3, 3),
#'               alpha = 0.05, alternative = "two.one.sided", plot = FALSE)
#'
#' # minimal effect testing
#' # power is defined as the probability of observing a test statistic greater
#' # than the upper critical value (for the upper bound) OR less than the lower
#' # critical value (for the lower bound).
#' power.lp.test(ncp = 2, null.ncp = c(-1, 1), df = 100, alpha = 0.05,
#'               alternative = "two.one.sided", plot = FALSE)
#' power.lp.test(power = 0.80, req.sign = "+", null.ncp = c(-1, 1),
#'               df = 100, alpha = 0.05, alternative = "two.one.sided")
#' # the first example (ncp = 2) reveals insufficient power (0.169), hence
#' # use the ncp returned from the example above for estimating the df's
#' power.lp.test(ncp = 3.844, power = 0.8, req.sign = "+", null.ncp = c(-3, 3),
#'               alpha = 0.05, alternative = "two.one.sided", plot = FALSE)
#'
#' @export power.lp.test
power.lp.test <- function(power = NULL, ncp = NULL, req.sign = "+", null.ncp = 0,
                          df = NULL, alpha = 0.05,
                          alternative = c("two.sided", "one.sided", "two.one.sided"),
                          plot = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  if (!is.null(power)) check.power(power)
  if (!is.null(ncp)) check.numeric(ncp)
  null.ncp <- check.margins(null.ncp, check.numeric, alternative)
  if (!is.null(df)) check.positive(df)
  check.proportion(alpha)
  check.logical(plot, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = ncp, n = df, power = power)

  if (!is.null(df) && df < 3)
    stop("`df` can not be smaller than 3.", call. = FALSE)
  if (!is.null(ncp) && !is.null(null.ncp) && (ncp > 35 || any(null.ncp > 35)))
    warning("Consider using a z-test. Lambda-prime distribution with a large non-centrality parameter can be unreliable.", call. = FALSE)

  # calculate statistical power
  pwr <- function(ncp = NULL, null.ncp = 0, df, alpha = 0.05, alternative) {

    if (alternative == "two.sided") {

      t.alpha <- c(qlambdap(p = alpha / 2, df = df, ncp = 0, lower.tail = TRUE),
                   qlambdap(p = alpha / 2, df = df, ncp = 0, lower.tail = FALSE))
      power <- 1 - plambdap(q = t.alpha[2], df = df, ncp = abs(ncp)) +
                   plambdap(q = t.alpha[1], df = df, ncp = abs(ncp))

      Phi.p <- stats::pt(q = max(t.alpha), df = df, ncp = ncp)
      Phi.m <- stats::pt(q = min(t.alpha), df = df, ncp = ncp)
      type.s <- min(Phi.m, 1 - Phi.p) / (Phi.m + 1 - Phi.p)

      type.m <- suppressMessages({
        bounds <- qlambdap(c(1e-10, 1 - 1e-10), df = df, ncp = ncp)
        integrand <- function(t) abs(t) * dlambdap(t, df = df, ncp = ncp)
        numerator <- stats::integrate(integrand, min(bounds), min(t.alpha))$value +
                     stats::integrate(integrand, max(t.alpha), max(bounds))$value
        denominator  <- abs(ncp) * (plambdap(min(t.alpha), df = df, ncp = ncp) +
                                    plambdap(max(t.alpha), df = df, ncp = ncp, lower.tail = FALSE))
        numerator / denominator
      })

    } else if (alternative == "one.sided") {

      lower.tail <- ncp < null.ncp
      t.alpha <- qlambdap(p = alpha,   df = df, ncp = null.ncp, lower.tail = lower.tail)
      power   <- plambdap(q = t.alpha, df = df, ncp = ncp,      lower.tail = lower.tail)

      type.s <- 0
      type.m <- NA

    } else if (alternative == "two.one.sided" && (ncp > min(null.ncp) && ncp < max(null.ncp))) {  # equivalence test

      t.alpha.left  <- qlambdap(p = alpha,     df = df, ncp = min(null.ncp), lower.tail = FALSE)
      t.alpha.right <- qlambdap(p = 1 - alpha, df = df, ncp = max(null.ncp), lower.tail = FALSE)
      power <- plambdap(q = t.alpha.right, df = df, ncp = ncp) -
               plambdap(q = t.alpha.left,  df = df, ncp = ncp)

      t.alpha <- c(t.alpha.left,  t.alpha.right)

      type.s <- NA
      type.m <- NA

    } else if (alternative == "two.one.sided" && (ncp < min(null.ncp) || ncp > max(null.ncp))) {  # minimum effect test

      t.alpha.right <- qlambdap(p = alpha / 2, df = df, ncp = max(null.ncp), lower.tail = FALSE)
      t.alpha.left  <- qlambdap(p = alpha / 2, df = df, ncp = min(null.ncp), lower.tail = TRUE)
      power <- plambdap(q = t.alpha.right, df = df, ncp = ncp, lower.tail = FALSE) +
               plambdap(q = t.alpha.left,  df = df, ncp = ncp, lower.tail = TRUE)

      t.alpha <- c(t.alpha.left,  t.alpha.right)

      type.s <- NA
      type.m <- NA

    }

    power[power < 0] <- 0

    list(power = power, t.alpha = t.alpha, type.s = type.s, type.m = type.m)

  } # pwr()

  min.pwr <- function(ncp, df, power) {

    power - pwr(ncp = ncp, null.ncp = null.ncp, df = df, alpha = alpha, alternative = alternative)$power

  } # min.pwr() (for uniroot and optimize)

  if (requested == "es") {

    if (check.null_sign(req.sign, alternative)) {

      lower.int <- c(min(null.ncp), mean(null.ncp))
      upper.int <- c(mean(null.ncp), max(null.ncp))
      ncp.lower <- stats::optimize(f = function(ncp) min.pwr(ncp, df, power) ^ 2, interval = lower.int, tol = 1e-12)$minimum
      ncp.upper <- stats::optimize(f = function(ncp) min.pwr(ncp, df, power) ^ 2, interval = upper.int, tol = 1e-12)$minimum
      ncp <- mean(c(ncp.lower, ncp.upper))

      warn.txt <- ifelse(max(abs(c(min.pwr(ncp.lower, df, power), min.pwr(ncp.upper, df, power)))) < 1e-6,
                         sprintf("Target NCP ranges from %.4f to %.4f within the null bounds.", ncp.lower, ncp.upper),
                         "The target power rate cannot be achieved within the null bounds.")
      warning(warn.txt, call. = FALSE)

    } else {

      val.rng <- get.interval(null.ncp = null.ncp, distribution = "lp", alpha = alpha, alternative = alternative,
                              req.sign = req.sign, df = df)
      ncp <- stats::optimize(f = function(ncp) min.pwr(ncp, df, power) ^ 2, interval = val.rng, tol = 1e-12)$minimum

    }

  } else if (requested == "n") {

    df <- stats::optimize(f = function(df) min.pwr(ncp, df, power) ^ 2, interval = c(1, 1e10))$minimum

  }

  pwr.obj <- pwr(ncp = ncp, null.ncp = null.ncp, df = df, alpha = alpha, alternative = alternative)

  if (plot)
    .plot.lp.t1t2(ncp = ncp, null.ncp = null.ncp, df = df, alpha = alpha, alternative = alternative)

  if (verbose > 0) {

    print.obj <- list(test = "Generic Lambda-Prime Distribution",
                      requested = requested,
                      tgt.ncp = "lambda",
                      lambda = ncp,
                      null.lambda = null.ncp,
                      df = df,
                      alpha = alpha,
                      alternative = alternative,
                      t.alpha = pwr.obj$t.alpha,
                      power = pwr.obj$power)

    .print.pwrss.t(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(power = pwr.obj$power,
                           ncp = ncp,
                           null.ncp = null.ncp,
                           df = df,
                           alpha = alpha,
                           alternative = alternative,
                           t.alpha = pwr.obj$t.alpha,
                           beta = 1 - pwr.obj$power,
                           type.s = pwr.obj$type.s,
                           type.m = pwr.obj$type.m),
                      class = c("pwrss", "generic", "lp")))

} # end of power.lp.test()

power.lp <- power.lp.test
