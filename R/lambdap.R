#' Distribution functions for the Lambda prime / non-central Lambda distribution
#'
#' @name lambdap
#'
#' @description
#' Density, distribution function, quantile function and random generation for
#' the lambda prime distribution.
#'
#' @details
#' These functions compute LeCoutre's Lambda prime \eqn{\Lambda'} distribution
#' with df degrees of freedom (denoted df or \eqn{\nu}) and a non-centrality
#' parameter (denoted ncp or t, and being the observed t-statistic). It is a
#' continuous probability distribution that frequently arises in the sampling
#' distribution of confidence limits for a normal mean and for inferences
#' regarding signal-to-noise or standardized effect sizes. The distribution is
#' generally asymmetric, and its shape adapts based on its parameters. When the
#' non-centrality parameter (t) is zero, or if the degrees of freedom grow
#' large (\eqn{\chi^2_{df} / df \to 0}), it reduces / converges to the standard
#' normal distribution.
#' Formally: \deqn{\Lambda'_{df}(t) = z + t \sqrt{\chi^2_{df} / df}}
#' The non-central t distribution is the non-centrality parameter \eqn{\Lambda}
#' plus the standard normal z distribution, all divided by the square root of
#' the usual chi-square distribution divided by the degrees of freedom:
#' \deqn{t'_{df}(\Lambda) = (\Lambda + z) / \sqrt{\chi^2_{df} / df}}
#' A \eqn{\Lambda'} distributed random variable can be viewed as a confidence
#' level on a non-central t (with the confidence intervals being computed as
#' percent points of the \eqn{\Lambda'} distribution).
#'
#' @usage
#' dlambdap(x, df, ncp, log = FALSE)
#' plambdap(q, df, ncp, lower.tail = TRUE, log.p = FALSE)
#' qlambdap(p, df, ncp, lower.tail = TRUE, log.p = FALSE)
#' rlambdap(n, df, ncp)
#'
#' @param x,q         vector of quantiles
#' @param p           vector of probabilities
#' @param n           number of observations
#' @param df          the degrees of freedom of the distribution
#' @param ncp         the non-centrality parameter (t) of the distribution
#' @param lower.tail  logical; if `TRUE` (default), probabilities are
#'                    \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}
#' @param log,log.p   logical; if `TRUE`, probabilities / densities are given
#'                    as logarithms
#'
#' @return \code{dlambdap} gives the density, \code{plambdap} gives the
#' distribution function (probabilities), \code{qlambdap} gives the quantile
#' function, and \code{rlambdap} generates a random vector with lambda prime
#' distributed values.
#'
#' @aliases dlambdap plambdap qlambdap rlambdap
#'
#' @seealso
#' t distribution functions: [stats::dt()], [stats::pt()], [stats::qt()], and
#' [stats::rt()].
#'
#' @examples
#' set.seed(1)
#' dlambdap(11.1, df = 9, ncp = 10) # 0.1294471
#' plambdap(11.1, df = 9, ncp = 10) # 0.7134134
#' qlambdap(0.01, df = 9, ncp = 10) # 4.245347
#' rv <- rlambdap(100, df = 50, ncp = 2)
#' mean(rv) # 2.077029
#' pv <- plambdap(rv,  df = 50, ncp = 2)
#' summary(pv)
#' #     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' #  0.01715 0.31579 0.51642 0.52259 0.74839 0.99703
#' qv <- qlambdap(pv,  df = 50, ncp = 2)
#' summary(qv)
#' #     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#' #  -0.1677  1.5007  2.0319  2.0770  2.6726  4.7966
#' # absolute difference between the original random vector and
#' # the quantile vector calculated from the probabilities of the
#' # original random vector (< 1e-12)
#' max(abs(qv - rv)) # 0.0000000000002498002
#'
#' @references
#' LeCoutre, B. (2007). Another look at confidence intervals for the noncentral
#' t distribution. Journal of Modern Applied Statistical Methods, 6(1),
#' 107–116. https://doi.org/10.22237/jmasm/1177992600
#'
#' @export plambdap
plambdap <- function(q, df, ncp, lower.tail = TRUE, log.p = FALSE) {
    suppressWarnings(stats::pt(q = ncp, df = df, ncp = q, lower.tail = !lower.tail, log.p = log.p))
}

.qlambdap <- function(p, df, ncp, lower.tail = TRUE, log.p = FALSE) {
    # handle border-cases
    if (log.p) {
        if (p == -Inf) return(ifelse(lower.tail, -Inf, Inf))
        if (p == 0)    return(ifelse(lower.tail,  Inf, -Inf))
        if (p  > 0)    return(NA)
    } else {
        if (p == 1) return(ifelse(lower.tail,  Inf, -Inf))
        if (p == 0) return(ifelse(lower.tail, -Inf,  Inf))
        if (p < 0 || p > 1) return(NA)
    }

    if (lower.tail) {
        zerof <- function(q)     plambdap(q, df = df, ncp = ncp, lower.tail = lower.tail, log.p = log.p) - p
    } else {
        zerof <- function(q) p - plambdap(q, df = df, ncp = ncp, lower.tail = lower.tail, log.p = log.p)
    }
    zmax <- 2 * max(stats::qnorm(p, lower.tail = TRUE,  log.p = log.p),
                    stats::qnorm(p, lower.tail = FALSE, log.p = log.p))
    flim <- c(min(-1, ncp - zmax), max(1, ncp + zmax))
    for (i in seq(2)) while (sign(zerof(flim[i])) == ifelse(i == 1, 1, -1) && abs(flim[i]) < 1e8) flim[i] <- 2 * flim[i]

    stats::uniroot(zerof, flim, tol = 1e-12)$root
}

#' @rdname lambdap
#' @export qlambdap
qlambdap <- Vectorize(.qlambdap, vectorize.args = c("p", "df", "ncp"), SIMPLIFY = TRUE)

.dlambdap <- function(x, df, ncp, log = FALSE) {
    gradf <- function(q) plambdap(q, df, ncp, lower.tail = TRUE, log.p = log) # used 3 and 18 lines below

    x <- as.numeric(x)
    d <- stats::numericDeriv(quote(gradf(x)), "x")
    if (attr(d, "gradient") >  0) {
        as.numeric(attr(d, "gradient"))
    } else {
        # if no gradient is found, this can either be due to the value being positioned at the outer ends
        # of the tails of the distribution (either being smaller than the tol or larger than 1 - tol);
        # here 0 is the appropriate return value
        tol <- sqrt(.Machine$double.eps)
        if (d < tol || d > (1 - tol)) {
            0
        # the alternative where no gradient is found is around the peak, where tol is to low to
        # derive a value for density; this happens very rarely; to return a value, several parameter
        # settings for eps are tried out, and the one of those which reveals the most consistent
        # results is selected, and the median of these results is selected
        } else {
            epsf <- function(e) as.numeric(attr(stats::numericDeriv(quote(gradf(x)), "x", eps = e), "gradient"))
            epsv <- vapply(seq(0, 3), function(e) vapply(seq(tol, 10 ^ -e, length.out = 1e3), epsf, numeric(1)), numeric(1e3))
            epsv <- epsv[, colSums(epsv == 0) < 1e2, drop = FALSE]  # ensure sufficient valid values
            epsv <- epsv[, which.min(apply(epsv, 2, stats::sd))]    # select column with the highest consistency
            stats::median(epsv[epsv != 0 & !is.na(epsv)])
        }
    }
}

#' @rdname lambdap
#' @export dlambdap
dlambdap <- Vectorize(.dlambdap, vectorize.args = c("x", "df", "ncp"), SIMPLIFY = TRUE)

#' @rdname lambdap
#' @export rlambdap
rlambdap <- function(n, df, ncp) {
    stats::rnorm(n) + ncp * sqrt(stats::rchisq(n, df = df) / df)
}
