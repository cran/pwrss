#############################################################
#################### non-parametric tests ###################
#### one sample, paired samples, and independent samples ####
#############################################################

#' Power Analysis for Non-parametric Rank-Based Tests (One-Sample, Independent,
#' and Paired Designs)
#'
#' @description
#' Calculates power, sample size or effect size (only one can be NULL at a
#' time) for non-parametric rank-based tests. The following tests and designs
#' are available:
#' * Wilcoxon Signed-Rank Test (One Sample)
#' * Wilcoxon Rank-Sum or Mann-Whitney U Test (Independent Samples)
#' * Wilcoxon Matched-Pairs Signed-Rank Test (Paired Samples)
#'
#' Formulas are validated using G*Power and tables in the PASS documentation.
#' However, we adopt the rounding convention used by G*Power.
#'
#' @details
#' * Use \code{means.to.d()} to convert raw means and standard deviations to
#'   Cohen's d, and \code{d.to.cles()} to convert Cohen's d to the probability
#'   of superiority. Note that this interpretation is appropriate only when the
#'   underlying distribution is approximately normal and the two groups have
#'   similar population variances.
#' * NB: \code{pwrss.np.2means()} function is depreciated and no longer
#'   supported. \code{pwrss.np.2groups()} will remain available for some time.
#' * Note that R has a partial matching feature which allows you to specify
#'   shortened versions of arguments, such as \code{alt} instead of
#'   \code{alternative}, or \code{dist} instead of \code{distribution}.
#'
#' @aliases power.np.wilcox pwrss.np.2groups pwrss.np.2means
#'
#'
#' @param d            Cohen's d or Hedges' g.
#' @param null.d       Cohen's d or Hedges' g under null, typically 0 (zero).
#' @param margin       margin - ignorable \code{d} - \code{null.d} difference.
#' @param req.sign     whether `d` is smaller or larger than `null.d` (when
#'                     minimum detectable prob is of interest).
#' @param n.ratio      \code{n1 / n2} ratio (applies to independent samples
#'                     only)
#' @param n2           integer; sample size in the second group (or for the
#'                     single group in paired samples or one-sample)
#' @param power        statistical power, defined as the probability of
#'                     correctly rejecting a false null hypothesis, denoted as
#'                     \eqn{1 - \beta}.
#' @param alpha        type 1 error rate, defined as the probability of
#'                     incorrectly rejecting a true null hypothesis, denoted as
#'                     \eqn{\alpha}.
#' @param design       character; "independent" (default), "one.sample", or
#'                     "paired".
#' @param alternative  character; the direction or type of the hypothesis test:
#'                     "two.sided", "one.sided", or "two.one.sided".
#' @param distribution character; parent distribution: "normal", "uniform",
#'                     "double.exponential", "laplace", or "logistic".
#' @param method       character; non-parametric approach: "guenther" (default)
#'                     or "noether"
#' @param ceil.n       logical; whether sample size should be rounded up.
#'                     \code{TRUE} by default.
#' @param verbose      \code{1} by default (returns test, hypotheses, and
#'                     results), if \code{2} a more detailed output is given
#'                     (plus key parameters and definitions), if \code{0} no
#'                     output is printed on the console.
#' @param utf          logical; whether the output should show Unicode
#'                     characters (if encoding allows for it). \code{FALSE} by
#'                     default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test (Z- or T-Test).}
#'   \item{df}{degrees of freedom (applies when method = 'guenther').}
#'   \item{ncp}{non-centrality parameter for the alternative (applies when
#'              method = 'guenther').}
#'   \item{null.ncp}{non-centrality parameter for the null (applies when
#'                   method = 'guenther').}
#'   \item{t.alpha}{critical value(s) (applies when method = 'guenther').}
#'   \item{mean}{mean of the alternative (applies when method = 'noether').}
#'   \item{null.mean}{mean of the null (applies when method = 'noether').}
#'   \item{sd}{standard deviation of the alternative (applies when method =
#'             'noether').}
#'   \item{null.sd}{standard deviation of the null (applies when method =
#'         'noether').}
#'   \item{z.alpha}{critical value(s) (applies when method = 'noether').}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample size (`n` or `c(n1, n2)` depending on the design.}
#'
#' @references
#'   Al-Sunduqchi, M. S. (1990). *Determining the appropriate sample size for
#'   inferences based on the Wilcoxon statistics* \[Unpublished doctoral
#'   dissertation\]. University of Wyoming - Laramie
#'
#'   Chow, S. C., Shao, J., Wang, H., and Lokhnygina, Y. (2018). *Sample size
#'   calculations in clinical research* (3rd ed.). Taylor & Francis/CRC.
#'
#'   Lehmann, E. (1975). *Nonparameterics: Statistical methods based on ranks*.
#'   McGraw-Hill.
#'
#'   Noether, G. E. (1987). Sample size determination for some common
#'   nonparametric tests. *Journal of the American Statistical Association,
#'   82*(1), 645-647.
#'
#'   Ruscio, J. (2008). A probability-based measure of effect size: Robustness
#'   to base rates and other factors. *Psychological Methods, 13*(1), 19-30.
#'
#'   Ruscio, J., & Mullen, T. (2012). Confidence intervals for the probability
#'   of superiority effect size measure and the area under a receiver operating
#'   characteristic curve. *Multivariate Behavioral Research, 47*(2), 201-223.
#'
#'   Zhao, Y.D., Rahardja, D., & Qu, Y. (2008). Sample size calculation for the
#'   Wilcoxon-Mann-Whitney test adjusting for ties. *Statistics in Medicine,
#'   27*(3), 462-468.
#'
#' @examples
#' # see `?pwrss::power.t.student` for further examples
#'
#' # Mann-Whitney U or Wilcoxon rank-sum test
#' # (a.k.a Wilcoxon-Mann-Whitney test) for independent samples
#'
#' ## difference between group 1 and group 2 is not equal to zero
#' ## estimated difference is Cohen'd = 0.25
#' power.np.wilcoxon(d = 0.25,
#'                 power = 0.80)
#'
#' ## difference between group 1 and group 2 is greater than zero
#' ## estimated difference is Cohen'd = 0.25
#' power.np.wilcoxon(d = 0.25,
#'                 power = 0.80,
#'                 alternative = "one.sided")
#'
#' ## mean of group 1 is practically not smaller than mean of group 2
#' ## estimated difference is Cohen'd = 0.10 and can be as small as -0.05
#' power.np.wilcoxon(d = 0.10,
#'                 margin = -0.05,
#'                 power = 0.80,
#'                 alternative = "one.sided")
#'
#' ## mean of group 1 is practically greater than mean of group 2
#' ## estimated difference is Cohen'd = 0.10 and can be as small as 0.05
#' power.np.wilcoxon(d = 0.10,
#'                 margin = 0.05,
#'                 power = 0.80,
#'                 alternative = "one.sided")
#'
#' ## mean of group 1 is practically same as mean of group 2
#' ## estimated difference is Cohen'd = 0
#' ## and can be as small as -0.05 and as high as 0.05
#' power.np.wilcoxon(d = 0,
#'                 margin = c(-0.05, 0.05),
#'                 power = 0.80,
#'                 alternative = "two.one.sided")
#'
#'
#' # Wilcoxon signed-rank test for matched pairs (dependent samples)
#'
#' ## difference between time 1 and time 2 is not equal to zero
#' ## estimated difference between time 1 and time 2 is Cohen'd = -0.25
#' power.np.wilcoxon(d = -0.25,
#'                 power = 0.80,
#'                 design = "paired")
#'
#' ## difference between time 1 and time 2 is greater than zero
#' ## estimated difference between time 1 and time 2 is Cohen'd = -0.25
#' power.np.wilcoxon(d = -0.25,
#'                 power = 0.80,
#'                 design = "paired",
#'                 alternative = "one.sided")
#'
#' ## mean of time 1 is practically not smaller than mean of time 2
#' ## estimated difference is Cohen'd = -0.10 and can be as small as 0.05
#' power.np.wilcoxon(d = -0.10,
#'                 margin = 0.05,
#'                 power = 0.80,
#'                 design = "paired",
#'                 alternative = "one.sided")
#'
#' ## mean of time 1 is practically greater than mean of time 2
#' ## estimated difference is Cohen'd = -0.10 and can be as small as -0.05
#' power.np.wilcoxon(d = -0.10,
#'                 margin = -0.05,
#'                 power = 0.80,
#'                 design = "paired",
#'                 alternative = "one.sided")
#'
#' ## mean of time 1 is practically same as mean of time 2
#' ## estimated difference is Cohen'd = 0
#' ## and can be as small as -0.05 and as high as 0.05
#' power.np.wilcoxon(d = 0,
#'                 margin = c(-0.05, 0.05),
#'                 power = 0.80,
#'                 design = "paired",
#'                 alternative = "two.one.sided")
#'
#' @export power.np.wilcoxon
power.np.wilcoxon <- function(d = NULL, null.d = 0, margin = 0, req.sign = "+",
                              n.ratio = 1, n2 = NULL, power = NULL, alpha = 0.05,
                              alternative = c("two.sided", "one.sided", "two.one.sided"),
                              design = c("independent", "paired", "one.sample"),
                              distribution = c("normal", "uniform", "double.exponential", "laplace", "logistic"),
                              method = c("guenther", "noether"),
                              ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  distribution <- match.arg(distribution)
  method <- match.arg(method)
  design <- match.arg(design)
  func.parms <- as.list(environment())

  if (!is.null(d)) check.numeric(d)
  check.numeric(null.d)
  margin <- check.margins(margin, check.numeric, alternative)
  check.positive(n.ratio)
  if (!is.null(n2)) check.sample.size(n2)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = d, n = n2, power = power)

  if (any(abs(margin) > 10))
      stop("Possibly incorrect value for `margin` (should be within -10 ... 10).", call. = FALSE)

  if (method == "noether" && design != "independent")
    stop("Specify `method` = \"guenther\" to request Wilcoxon signed-rank test for matched pairs.", call. = FALSE)

  # if (requested == "es" && alternative == "two.one.sided")
  #  stop("Determining the effect size is not possible if `alternative` is \"two.one.sided\".", call. = FALSE)

  propss <- n.ratio / (n.ratio + 1)

  pwr.wilcox <- function(d, null.d, margin, n2, n.ratio, alpha, design, alternative, method) {

    n1 <- n2 * n.ratio

    if (method == "noether") {

      prob <- d.to.cles(d = d, design = design, verbose = 0)$cles
      null.prob <- d.to.cles(d = null.d, design = design, verbose = 0)$cles

      if (alternative == "two.one.sided") {
        ignorable.prob.lower <- d.to.cles(d = min(margin) + null.d, design = design, verbose = 0)$cles
        ignorable.prob.upper <- d.to.cles(d = max(margin) + null.d, design = design, verbose = 0)$cles
        margin.prob <-  c(ignorable.prob.lower, ignorable.prob.upper) - null.prob
      } else {
        margin.prob <- d.to.cles(d = margin + null.d, design = design, verbose = 0)$cles - null.prob
      }

      lambda      <- sqrt(n2 + n2 * n.ratio) * sqrt(12 * propss * (1 - propss)) * (prob - null.prob)
      null.lambda <- sqrt(n2 + n2 * n.ratio) * sqrt(12 * propss * (1 - propss)) * (margin.prob)

      power.z.test(mean = lambda, null.mean = null.lambda, alpha = alpha, alternative = alternative,
                   plot = FALSE, verbose = 0)[c("power", "z.alpha", "mean", "null.mean")]

    } else if (method == "guenther") {

      if (design == "independent") {
        df <- n1 + n2 - 2
        lambda      <- (d - null.d) / sqrt(1 / n1 + 1 / n2)
        null.lambda <- (margin)     / sqrt(1 / n1 + 1 / n2)
      } else {
        df <- n2 - 1
        lambda      <- (d - null.d) / sqrt(1 / n2)
        null.lambda <- (margin)     / sqrt(1 / n2)
      }

      power.t.test(ncp = lambda, null.ncp = null.lambda, df = df, alpha = alpha, alternative = alternative,
                   plot = FALSE, verbose = 0)[c("power", "t.alpha", "ncp", "null.ncp", "df")]

    }

  } # pwr.wilcox

  # wilcoxon adjustment for guenther method
  if (method == "guenther") {
    w <- switch(distribution, `uniform` = 1, `double.exponential` = 2 / 3, `laplace` = 2 / 3,
                              `logistic` = 9 / pi ^ 2, `normal` = pi / 3)
  } else {
    w <- 1
  }

  min.pwr.wilcox <- function(d, n2, power) {
    power - pwr.wilcox(d = d, null.d = null.d, margin = margin, n2 = n2 / w, n.ratio = n.ratio, alpha = alpha,
                       design = design, alternative = alternative, method = method)$power
  } # min.pwr (for stats::uniroot)

  # calculate sample size
  if (requested == "n") {

    H1_H0.min <- 0.001
    lambda.max <- 4

    if (method == "noether") {
      n2.max <- (lambda.max / (sqrt(12 * propss * (1 - propss)) * (H1_H0.min))) ^ 2 / (1 + n.ratio)
    } else {
      n2.max <- ifelse(design == "independent", (1 + 1 / n.ratio) / (H1_H0.min / lambda.max) ^ 2, (lambda.max / H1_H0.min) ^ 2)
    }
    val.rng <- c(4, n2.max) # power.t.test requires at least df = 3, hence the minumum N is 4

    # estimation, using wilcoxon adjustment (* w)
    n2 <- try(stats::uniroot(function(n2) min.pwr.wilcox(d, n2 * w, power), interval = val.rng)$root * w, silent = TRUE)
    if (inherits(n2, "try-error") || n2 == 1e10) stop("Design is not feasible.", call. = FALSE)

    n2 <- ifelse(ceil.n, ceiling(n2), n2)

  } else if (requested == "es") {

    if (alternative == "two.one.sided" && check.null_sign(req.sign, alternative)) {

      lower.int <- c(min(margin) + null.d, mean(margin) + null.d) + c(+1e-7, 0)
      upper.int <- c(mean(margin) + null.d, max(margin) + null.d) + c(0, -1e-7)
      d.lower <- stats::optimize(f = function(d) min.pwr.wilcox(d, n2, power) ^ 2, interval = lower.int, tol = 1e-12)$minimum
      d.upper <- stats::optimize(f = function(d) min.pwr.wilcox(d, n2, power) ^ 2, interval = upper.int, tol = 1e-12)$minimum

      d <- mean(c(d.lower, d.upper))

      pwr.lower <- pwr.wilcox(d = d.lower, null.d = null.d, margin = margin, n2 = n2 / w, n.ratio = n.ratio, alpha = alpha,
                              design = design, alternative = alternative, method = method)$power
      pwr.upper <- pwr.wilcox(d = d.upper, null.d = null.d, margin = margin, n2 = n2 / w, n.ratio = n.ratio, alpha = alpha,
                              design = design, alternative = alternative, method = method)$power

      if (round(pwr.lower, 3) >= power && round(pwr.upper, 3) >= power) {

        warning("Target effect ranges from ", round(d.lower, 4), " to ", round(d.upper, 4), " within the null bounds.",
                call. = FALSE)

      } else {

        warning("The target power rate cannot be achieved within the null bounds.", call. = FALSE)

      }

    } else {

      if (check.pos_sign(req.sign)) {
        d.int <- c(max(margin) + null.d, 10) + c(+1e-7, -1e-7)
      } else {
        d.int <- c(-10, min(margin) + null.d) + c(+1e-7, -1e-7)
      }

      d <- try(stats::uniroot(f = function(d) min.pwr.wilcox(d, n2, power), interval = d.int, tol = 1e-12)$root, silent = TRUE)
      if (inherits(d, "try-error")) stop("Design is not feasible.", call. = FALSE)

    } # two.one.sided?

  } # ss or es

  pwr.obj <- pwr.wilcox(d = d, null.d = null.d, margin = margin,
                        n2 = n2 / w, n.ratio = n.ratio, alpha = alpha,
                        design = design, alternative = alternative,
                        method = method)

  n1 <- ifelse(ceil.n, ceiling(n2 * n.ratio), n2 * n.ratio)
  if (design == "independent") n <- c(n1 = n1, n2 = n2) else n <- n2

  if (method == "guenther") {
    list.out <- list(d = d,
                     n = n,
                     power = pwr.obj$power,
                     t.alpha = pwr.obj$t.alpha,
                     ncp = pwr.obj$ncp,
                     null.ncp = pwr.obj$null.ncp,
                     df = pwr.obj$df)
  } else if (method == "noether") {
    list.out <- list(d = d,
                     n = n,
                     power = pwr.obj$power,
                     z.alpha = pwr.obj$z.alpha,
                     mean = pwr.obj$mean,
                     sd = 1,
                     null.mean = pwr.obj$null.mean,
                     null.sd = 1,
                     df = Inf)
  }

  if (verbose > 0) {

    print.obj <- c(list(requested = requested, test = fmt_test_wilcoxon(design),
                        design = design, method = method, dist = distribution,
                        null.d = null.d, margin = margin, # d is in list.out
                        alpha = alpha, alternative = alternative),
                   list.out)

    .print.pwrss.wilcoxon(print.obj, verbose = verbose, utf = utf)

  }

  invisible(structure(c(list(parms = func.parms,
                             test = ifelse(method == "noether", "z", "t")),
                        list.out),
                      class = c("pwrss", "np", "wilcoxon", ifelse(method == "noether", "z", "t"))))

} # end of power.np.wilcoxon()

#' @rdname power.np.wilcoxon
#' @export power.np.wilcox
power.np.wilcox <- power.np.wilcoxon

#' @export pwrss.np.2groups
pwrss.np.2groups <- function(mu1 = 0.20, mu2 = 0,
                             sd1 = ifelse(paired, sqrt(1 / (2 * (1 - paired.r))), 1),
                             sd2 = sd1, margin = 0, alpha = 0.05, paired = FALSE,
                             paired.r = 0.50, kappa = 1, n2 = NULL, power = NULL,
                             alternative = c("not equal", "greater", "less",
                                             "non-inferior", "superior", "equivalent"),
                             distribution = c("normal", "uniform", "double.exponential",
                                              "laplace", "logistic"),
                             method = c("guenther", "noether"),
                             verbose = TRUE) {

  method <- match.arg(method)
  alternative <- match.arg(alternative)
  distribution <- match.arg(distribution)
  verbose <- ensure.verbose(verbose)

  null.d <- 0
  d <- means.to.d(mu1 = mu1, mu2 = mu2,
                  sd1 = sd1, sd2 = sd2,
                  n2 = 1e4, n.ratio = kappa,
                  paired = paired,
                  rho.paired = paired.r,
                  verbose = 0)$d

  if (alternative == "equivalent") {
    margin <- c(min(-margin, margin), max(-margin, margin))
    margin.lower <- means.to.d(mu1 = margin[1], mu2 = 0,
                               sd1 = sd1, sd2 = sd2,
                               n2 = 1e10, n.ratio = kappa,
                               paired = paired,
                               rho.paired = paired.r,
                               verbose = 0)$d
    margin.upper <- means.to.d(mu1 = margin[2], mu2 = 0,
                               sd1 = sd1, sd2 = sd2,
                               n2 = 1e10, n.ratio = kappa,
                               paired = paired,
                               rho.paired = paired.r,
                               verbose = 0)$d
    margin <- c(margin.lower, margin.upper)
  } else {
    margin <- means.to.d(mu1 = margin, mu2 = 0,
                         sd1 = sd1, sd2 = sd2,
                         n2 = 1e10, n.ratio = kappa,
                         paired = paired,
                         rho.paired = paired.r,
                         verbose = 0)$d
  }


  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"
  if (alternative == "not equal") alternative <- "two.sided"

  design <- ifelse(paired, "paired", "independent")

  wilcox.obj <- power.np.wilcoxon(d = d,
                                  null.d = null.d,
                                  margin = margin,
                                  alpha = alpha,
                                  n2 = n2, n.ratio = kappa,
                                  power = power,
                                  design = design,
                                  alternative = alternative,
                                  distribution = distribution,
                                  method = method,
                                  ceil.n = TRUE,
                                  verbose = verbose)

  # cat("This function will be removed in the future. \n Please use `power.np.wilcoxon()`. \n")

  invisible(wilcox.obj)

} # end of pwrss.np.2groups()

fmt_test_wilcoxon <- function(design) {
  switch(design,
         `independent` = "Wilcoxon Rank-Sum Test (Independent Samples) \n(Wilcoxon-Mann-Whitney or Mann-Whitney U Test)",
         `paired`      = "Wilcoxon Signed-Rank Test (Paired Samples)",
         `one.sample`  = "Wilcoxon Signed-Rank Test (One Sample)")
}

#' @export pwrss.np.2means
pwrss.np.2means <- function(...) {
  stop("This function is no longer available. Please use `power.np.wilcoxon()`.", call. = FALSE)
}
