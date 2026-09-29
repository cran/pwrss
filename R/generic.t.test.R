#' Statistical Power for the Generic t-Test
#'
#' @description
#' Determines the power, the non-centrality parameter, or the degrees of
#' freedom for the generic t-Test with (optional) Type 1 and Type 2 error
#' plots.
#'
#' @aliases power.t
#'
#'
#' @param power       statistical power \eqn{(1 - \beta)}; either `power`,
#'                    `ncp` or `df` needs to be NULL (and is then estimated).
#' @param ncp         non-centrality parameter for the alternative; ; either
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
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{ncp}{non-centrality parameter under alternative.}
#'   \item{null.ncp}{non-centrality parameter under null.}
#'   \item{df}{degrees of freedom.}
#'   \item{alpha}{type 1 error rate (user-specified).}
#'   \item{t.alpha}{critical value(s).}
#'   \item{beta}{type 2 error rate.}
#'   \item{type.s}{type S error rate (only for two-tailed test).}
#'   \item{type.m}{type M error rate (only for two-tailed test).}
#'
#' @examples
#' # two-sided
#' # power defined as the probability of observing test statistics greater
#' # than the positive critical value OR less than the negative critical value
#' power.t.test(ncp = 1.96, df = 100, alpha = 0.05, alternative = "two.sided")
#' power.t.test(power = 0.80, df = 100, alpha = 0.05, alternative = "two.sided")
#'
#' # one-sided
#' # power is defined as the probability of observing a test statistic greater
#' # than the critical value
#' power.t.test(ncp = 1.96, df = 100, alpha = 0.05, alternative = "one.sided")
#' power.t.test(power = 0.80, df = 100, alpha = 0.05, alternative = "one.sided")
#'
#' # equivalence
#' # power is defined as the probability of observing a test statistic greater
#' # than the upper critical value (for the lower bound) AND less than the
#' # lower critical value (for the upper bound)
#' power.t.test(ncp = 0, null.ncp = c(-2, 2), df = 100, alpha = 0.05,
#'              alternative = "two.one.sided")
#' power.t.test(power = 0.80, req.sign = "0", null.ncp = c(-2, 2),
#'              df = 100, alpha = 0.05, alternative = "two.one.sided")
#'
#' # minimal effect testing
#' # power is defined as the probability of observing a test statistic greater
#' # than the upper critical value (for the upper bound) OR less than the lower
#' # critical value (for the lower bound).
#' power.t.test(ncp = 2, null.ncp = c(-1, 1), df = 100, alpha = 0.05,
#'              alternative = "two.one.sided")
#' power.t.test(power = 0.80, req.sign = "+", null.ncp = c(-1, 1),
#'              df = 100, alpha = 0.05, alternative = "two.one.sided")
#'
#' @export power.t.test
power.t.test <- function(power = NULL, ncp = NULL, req.sign = "+", null.ncp = 0,
                         df = NULL, alpha = 0.05,
                         alternative = c("two.sided", "one.sided", "two.one.sided"),
                         plot = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  if (!is.null(power)) check.power(power)
  if (!is.null(ncp)) check.numeric(ncp)
  null.ncp <- check.margins(null.ncp, check.numeric, alternative)
  if (!is.null(df) && df != Inf) check.positive(df)
  check.proportion(alpha)
  check.logical(plot, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = ncp, n = df, power = power)

  if (!is.null(df) && df < 3)
    stop("Degrees of freedom can not be smaller than 3.", call. = FALSE)

  # calculate statistical power
  pwr <- function(ncp, null.ncp = 0, df, alpha = 0.05, alternative) {

    if (alternative == "two.sided") {

      t.alpha <- c(stats::qt(alpha / 2,  df = df, ncp = null.ncp,      lower.tail = TRUE),
                   stats::qt(alpha / 2,  df = df, ncp = null.ncp,      lower.tail = FALSE))
      power   <-   stats::pt(t.alpha[1], df = df, ncp = ncp,           lower.tail = TRUE) +
                   stats::pt(t.alpha[2], df = df, ncp = ncp,           lower.tail = FALSE)

      Phi.p <- suppressWarnings(stats::pt(q = max(t.alpha), df = df, ncp = ncp))
      Phi.m <-                  stats::pt(q = min(t.alpha), df = df, ncp = ncp)
      type.s <- min(Phi.m, 1 - Phi.p) / (Phi.m + 1 - Phi.p)

      type.m <- try(suppressWarnings({
        bounds <- stats::qt(c(1e-10, 1 - 1e-10), df = df, ncp = ncp)
        integrand <- function(t) abs(t) * stats::dt(t, df = df, ncp = ncp)
        numerator <- stats::integrate(integrand, min(bounds), min(t.alpha))$value +
                     stats::integrate(integrand, max(t.alpha), max(bounds))$value
        denominator  <- abs(ncp) * (stats::pt(min(t.alpha), df = df, ncp = ncp, lower.tail = TRUE) +
                                    stats::pt(max(t.alpha), df = df, ncp = ncp, lower.tail = FALSE))
        numerator / denominator
      }), silent = TRUE)
      if (inherits(type.m, "try-error")) type.m <- NA

    } else if (alternative == "one.sided") {

      lower.tail <- ncp < null.ncp
      t.alpha <- suppressWarnings(stats::qt(alpha,        df = df, ncp = null.ncp,      lower.tail = lower.tail))
      power   <- suppressWarnings(stats::pt(t.alpha,      df = df, ncp = ncp,           lower.tail = lower.tail))

      type.s <- 0
      type.m <- NA

    } else if (alternative == "two.one.sided" && (ncp > min(null.ncp) && ncp < max(null.ncp))) {  # equivalence test

      t.alpha <- suppressWarnings(c(stats::qt(alpha,      df = df, ncp = min(null.ncp), lower.tail = FALSE),
                                    stats::qt(alpha,      df = df, ncp = max(null.ncp), lower.tail = TRUE)))
      power   <-   stats::pt(t.alpha[2], df = df, ncp = ncp,           lower.tail = TRUE) +
                   stats::pt(t.alpha[1], df = df, ncp = ncp,           lower.tail = FALSE) - 1

      type.s <- NA
      type.m <- NA

    } else if (alternative == "two.one.sided" && (ncp < min(null.ncp) || ncp > max(null.ncp))) {  # minimum effect test

      t.alpha <- suppressWarnings(c(stats::qt(alpha / 2,  df = df, ncp = min(null.ncp), lower.tail = TRUE),
                                    stats::qt(alpha / 2,  df = df, ncp = max(null.ncp), lower.tail = FALSE)))
      power   <-   suppressWarnings(stats::pt(t.alpha[1], df = df, ncp = ncp,           lower.tail = TRUE) +
                                    stats::pt(t.alpha[2], df = df, ncp = ncp,           lower.tail = FALSE))

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

      val.rng <- get.interval(null.ncp = null.ncp, distribution = "t", alpha = alpha, alternative = alternative,
                              req.sign = req.sign, df = df)
      ncp <- stats::optimize(f = function(ncp) min.pwr(ncp, df, power) ^ 2, interval = val.rng, tol = 1e-12)$minimum

    }

  } else if (requested == "n") {

    df <- stats::optimize(f = function(df) min.pwr(ncp, df, power) ^ 2, interval = c(1, 1e10))$minimum

  }

  pwr.obj <- pwr(ncp = ncp, null.ncp = null.ncp, df = df, alpha = alpha, alternative = alternative)

  if (plot)
    suppressWarnings(.plot.t.t1t2(ncp = ncp, null.ncp = null.ncp, df = df, alpha = alpha, alternative = alternative))

  if (verbose > 0) {

    print.obj <- list(test = "Generic t-Test",
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
                      class = c("pwrss", "generic", "t")))

} # end of power.t.test()

power.t <- power.t.test
