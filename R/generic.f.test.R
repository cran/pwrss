#' Statistical Power for the Generic F-Test
#'
#' @description
#' Determines the power or the non-centrality parameter for the generic F-Test
#' with (optional) Type 1 and Type 2 error plots.
#'
#' @aliases power.f
#'
#'
#' @param power    statistical power \eqn{(1 - \beta)}; either `power` or `ncp`
#'                 needs to be NULL (and is then estimated).
#' @param ncp      non-centrality parameter for the alternative; either `power`
#'                  or `ncp` needs to be NULL (and is then estimated).
#' @param null.ncp non-centrality parameter for the null.
#' @param alpha    type 1 error rate, defined as the probability of incorrectly
#'                 rejecting a true null hypothesis, denoted as \eqn{\alpha}.
#' @param df1      integer; numerator degrees of freedom.
#' @param df2      integer; denominator degrees of freedom.
#' @param plot     logical; \code{FALSE} switches off Type 1 and Type 2 error
#'                 plot. \code{TRUE} by default.
#' @param verbose  \code{1} by default (returns test, hypotheses, and results),
#'                 if \code{2} a more detailed output is given (plus key
#'                 parameters and definitions), if \code{0} no output is printed
#'                 on the console.
#' @param utf      logical; whether the output should show Unicode characters
#'                 (if encoding allows for it). \code{FALSE} by default.
#'
#' @return
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{ncp}{non-centrality parameter under alternative.}
#'   \item{null.ncp}{non-centrality parameter under null.}
#'   \item{df1}{numerator degrees of freedom.}
#'   \item{df2}{denominator degrees of freedom.}
#'   \item{f.alpha}{critical value(s).}
#'
#' @examples
#' # power is defined as the probability of observing a test statistics greater
#' # than the critical value
#' power.f.test(ncp = 1, df1 = 4, df2 = 100, alpha = 0.05)
#' power.f.test(power = 0.80, df1 = 4, df2 = 100, alpha = 0.05)
#'
#' @export power.f.test
power.f.test <- function(power = NULL, ncp = NULL, null.ncp = 0, df1, df2,
                         alpha = 0.05, plot = TRUE, verbose = 1, utf = FALSE) {

  if (!is.null(power)) check.power(power)
  if (!is.null(ncp)) check.positive(ncp)
  check.nonnegative(null.ncp)
  check.positive(df1, df2)
  check.proportion(alpha)
  check.logical(plot, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = ncp, n = NA, power = power)

  if (!is.null(ncp) && ncp < null.ncp)
    stop("`ncp` should be greater than or equal to `null.ncp`.", call. = FALSE)
  if (is.null(df1) || df1 < 1)
    stop("`df1` can not be NULL, and needs to be at least 1.", call. = FALSE)
  if (is.null(df2) || df2 < 3)
    stop("`df2` can not be NULL, and needs to be at least 3.", call. = FALSE)

  pwr <- function(ncp, null.ncp = 0, df1, df2, alpha = 0.05) {

    f.alpha <- stats::qf(alpha, df1 = df1, df2 = df2, ncp = null.ncp, lower.tail = FALSE)
    power   <- stats::pf(f.alpha, df1 = df1, df2 = df2, ncp = ncp, lower.tail = FALSE)

    list(power = power, f.alpha = f.alpha)

  } # pwr()

  min.pwr <- function(ncp, power) {
    power - pwr(ncp = ncp, null.ncp = null.ncp, df1 = df1, df2 = df2, alpha = alpha)$power
  } # min.pwr()

  if (requested == "es") {

    max.thresh <- stats::qf(1 - 1e-12, ncp = null.ncp, df1 = df1, df2 = df2)
    # adjust max.thresh: if min.pwr < 0, max.thresh is multiplied with 10 / 9, if min.pwr > 0 it is multiplied with 9 / 10
    for (chk.sign in c(-1, 1)) {
      if (sign(min.pwr(ncp = max.thresh, power = power)) == chk.sign) {
        while (sign(min.pwr(ncp = max.thresh, power = power)) == chk.sign) {
          max.thresh <- max.thresh * ifelse(chk.sign == -1, 9 / 10, 10 / 9)
        }
        max.thresh <- max.thresh * ifelse(chk.sign == -1, 10 / 9, 9 / 10)
      }
    }

    ncp <- stats::optimize(f = function(ncp) min.pwr(ncp, power) ^ 2, interval = c(0, max.thresh), tol = 1e-12)$minimum

  }

  pwr.obj <- pwr(ncp = ncp, null.ncp = null.ncp, df1 = df1, df2 = df2, alpha = alpha)

  if (plot)
    suppressWarnings(.plot.f.t1t2(ncp = ncp, null.ncp = null.ncp, df1 = df1, df2 = df2, alpha = alpha))

  if (verbose > 0) {

    print.obj <- list(test = "Generic F-Test",
                      requested = requested,
                      tgt.ncp = "lambda",
                      lambda = ncp,
                      null.lambda = null.ncp,
                      df1 = df1,
                      df2 = df2,
                      n.pres = df1 + df2 + 1,
                      power = pwr.obj$power,
                      alpha = alpha,
                      f.alpha = pwr.obj$f.alpha)

    .print.pwrss.f(print.obj, verbose = verbose, utf = utf)

  } # end of verbose

  invisible(structure(list(power = pwr.obj$power,
                           ncp = ncp,
                           null.ncp = null.ncp,
                           df1 = df1,
                           df2 = df2,
                           alpha = alpha,
                           f.alpha = pwr.obj$f.alpha),
                      class = c("pwrss", "generic", "f")))

} # end of power.f.test()

power.f <- power.f.test
