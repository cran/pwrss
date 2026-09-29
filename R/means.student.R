# one can specify Cohen's d or Hedges' g
# because small sample correction applies to the effect size not its standard error

#' Power Analysis for Student's t-Test
#'
#' @description
#' Calculates power, sample size or effect size (only one can be NULL at a
#' time) for Student's t-Test.
#'
#' In contrast to previous versions, users can now specify whether their claims
#' will be based on raw score mean difference with P-values or standardized
#' mean difference with confidence intervals. While results typically differ by
#' only a few units, these distinctions can be particularly consequential in
#' studies with small sample sizes or high-risk interventions.
#'
#' Formulas are validated using Monte Carlo simulations (see Bulus, 2024),
#' G*Power, and tables in the PASS documentation. One key difference between
#' PASS and \code{pwrss} lies in how they handle non-inferiority and
#' superiority tests-that is, one-sided tests defined by a negligible effect
#' margin (implemented as of this version). PASS shifts the test statistic so
#' that the null hypothesis assumes a zero effect, treating the negligible
#' margin as part of the alternative hypothesis. As a result, the test
#' statistic is evaluated against a central distribution. In contrast,
#' \code{pwrss} treats the negligible effect as the true null value, and the
#' test statistic is evaluated under a non-central distribution. This leads to
#' slight differences up to third decimal place. To get the same results,
#' reflect the margin in \code{null.d} and specify \code{margin = 0}.
#'
#' Equivalence tests are implemented in line with Bulus and Polat (2023), Chow
#' et al. (2018) and Lakens (2017).
#'
#' @details
#' * Use \code{means.to.d()} to convert raw means and standard deviations to
#'   Cohen's d, and \code{d.to.cles()} to convert Cohen's d to the probability
#'   of superiority. Note that this interpretation is appropriate only when the
#'   underlying distribution is approximately normal and the two groups have
#'   similar population variances.
#' * NB: The functions \code{pwrss.z.mean()} and \code{pwrss.z.2means()} are no
#'   longer supported. The \code{pwrss.t.mean()} and \code{pwrss.t.2means()}
#'   functions are deprecated, but they will remain available as wrappers for
#'   \code{power.t.student()} or \code{power.t.welch()} functions during a
#'   transition period.
#'
#' @param d           Cohen's d or Hedges' g.
#' @param null.d      Cohen's d or Hedges' g under null, typically 0 (zero).
#' @param margin      margin - ignorable \code{d} - \code{null.d} difference.
#' @param req.sign    whether `d` is smaller or larger than `null.d` (when
#'                    minimum detectable prob is of interest).
#' @param n.ratio     \code{n1 / n2} ratio (applies to independent samples
#'                    only)
#' @param n2          integer; sample size in the second group (or for the
#'                    single group in paired samples or one-sample).
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided", "one.sided", or "two.one.sided". For
#'                    non-inferiority or superiority tests, add or subtract the
#'                    margin from the null hypothesis value and use
#'                    \code{alternative = "one.sided"}.
#' @param claim.basis character; "md.pval" when claims are based on raw mean
#'                    differences and p-values, "smd.ci" when claims are based
#'                    on standardized mean differences and confidence
#'                    intervals.
#' @param design      character; "independent", "paired" or "one.sample".
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
#'   \item{test}{type of the statistical test (T-Test).}
#'   \item{df}{degrees of freedom.}
#'   \item{ncp}{non-centrality parameter for the alternative.}
#'   \item{null.ncp}{non-centrality parameter for the null.}
#'   \item{t.alpha}{critical value(s).}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample size (`n` or `c(n1, n2)` depending on the design.}
#'
#' @references
#'   Bulus, M. (2024). *Robust standard errors and confidence intervals for
#'   standardized mean differences*. https://doi.org/10.31219/osf.io/k6mbs
#'
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#'   Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). *Sample size
#'   calculations in clinical research* (3rd ed.). Taylor & Francis/CRC.
#'
#'   Cohen, J. (1988). *Statistical power analysis for the behavioral sciences*
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#'   Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
#'   correlations, and meta-analyses. *Social psychological and personality
#'   science, 8*(4), 355-362. https://doi.org/10.1177/1948550617697177
#'
#'
#' @aliases power.t.student pwrss.t.2means pwrss.z.2means pwrss.t.mean pwrss.z.mean
#'
#' @examples
#'
#' #######################
#' # Independent Samples #
#' #######################
#'
#' ## difference between group 1 and group 2 is not equal to zero
#' ## targeting minimal difference of Cohen'd = 0.20
#' ## non-parametric
#' power.np.wilcoxon(d = 0.20,
#'                   power = 0.80,
#'                   alternative = "two.sided",
#'                   design = "independent")
#'
#' ## parametric
#' power.t.student(d = 0.20,
#'                 power = 0.80,
#'                 alternative = "two.sided",
#'                 design = "independent")
#'
#' ## when sample size ratio and group variances differ
#' power.t.welch(d = 0.20,
#'               n.ratio = 2,
#'               var.ratio = 2,
#'               power = 0.80,
#'               alternative = "two.sided")
#'
#'
#' ## difference between group 1 and group 2 is greater than zero
#' ## targeting minimal difference of Cohen'd = 0.20
#' ## non-parametric
#' power.np.wilcoxon(d = 0.20,
#'                   power = 0.80,
#'                   alternative = "one.sided",
#'                   design = "independent")
#'
#' ## parametric
#' power.t.student(d = 0.20,
#'                 power = 0.80,
#'                 alternative = "one.sided",
#'                 design = "independent")
#'
#' ## when sample size ratio and group variances differ
#' power.t.welch(d = 0.20,
#'               n.ratio = 2,
#'               var.ratio = 2,
#'               power = 0.80,
#'               alternative = "one.sided")
#'
#'
#' ## mean of group 1 is practically not smaller than mean of group 2
#' ## targeting minimal difference of Cohen'd = 0.20 and can be as small as -0.05
#' ## non-parametric
#' power.np.wilcoxon(d = 0.20,
#'                   margin = -0.05,
#'                   power = 0.80,
#'                   alternative = "one.sided",
#'                   design = "independent")
#'
#' ## parametric
#' power.t.student(d = 0.20,
#'                 margin = -0.05,
#'                 power = 0.80,
#'                 alternative = "one.sided",
#'                 design = "independent")
#'
#' ## when sample size ratio and group variances differ
#' power.t.welch(d = 0.20,
#'               margin = -0.05,
#'               n.ratio = 2,
#'               var.ratio = 2,
#'               power = 0.80,
#'               alternative = "one.sided")
#'
#'
#' ## mean of group 1 is practically greater than mean of group 2
#' ## targeting minimal difference of Cohen'd = 0.20 and can be as small as 0.05
#' ## non-parametric
#' power.np.wilcoxon(d = 0.20,
#'                   margin = 0.05,
#'                   power = 0.80,
#'                   alternative = "one.sided",
#'                   design = "independent")
#'
#' ## parametric
#' power.t.student(d = 0.20,
#'                 margin = 0.05,
#'                 power = 0.80,
#'                 alternative = "one.sided",
#'                 design = "independent")
#'
#' ## when sample size ratio and group variances differ
#' power.t.welch(d = 0.20,
#'               margin = 0.05,
#'               n.ratio = 2,
#'               var.ratio = 2,
#'               power = 0.80,
#'               alternative = "one.sided")
#'
#'
#' ## mean of group 1 is practically same as mean of group 2
#' ## targeting minimal difference of Cohen'd = 0
#' ## and can be as small as -0.05 or as high as 0.05
#' ## non-parametric
#' power.np.wilcoxon(d = 0,
#'                   margin = c(-0.05, 0.05),
#'                   power = 0.80,
#'                   alternative = "two.one.sided",
#'                   design = "independent")
#'
#' ## parametric
#' power.t.student(d = 0,
#'                 margin = c(-0.05, 0.05),
#'                 power = 0.80,
#'                 alternative = "two.one.sided",
#'                 design = "independent")
#'
#' ## when sample size ratio and group variances differ
#' power.t.welch(d = 0,
#'               margin = c(-0.05, 0.05),
#'               n.ratio = 2,
#'               var.ratio = 2,
#'               power = 0.80,
#'               alternative = "two.one.sided")
#'
#'
#' ##################
#' # Paired Samples #
#' ##################
#'
#' ## difference between time 1 and time 2 is not equal to zero
#' ## targeting minimal difference of Cohen'd = -0.20
#' ## non-parametric
#' power.np.wilcoxon(d = -0.20,
#'                   power = 0.80,
#'                   alternative = "two.sided",
#'                   design = "paired")
#'
#' ## parametric
#' power.t.student(d = -0.20,
#'                 power = 0.80,
#'                 alternative = "two.sided",
#'                 design = "paired")
#'
#' ## difference between time 1 and time 2 is less than zero
#' ## targeting minimal difference of Cohen'd = -0.20
#' ## non-parametric
#' power.np.wilcoxon(d = -0.20,
#'                   power = 0.80,
#'                   alternative = "one.sided",
#'                   design = "paired")
#'
#' ## parametric
#' power.t.student(d = -0.20,
#'                 power = 0.80,
#'                 alternative = "one.sided",
#'                 design = "paired")
#'
#' ## mean of time 1 is practically not greater than mean of time 2
#' ## targeting minimal difference of Cohen'd = -0.20 and can be as small as 0.05
#' ## non-parametric
#' ## non-parametric
#' power.np.wilcoxon(d = 0.20,
#'                   margin = 0.05,
#'                   power = 0.80,
#'                   alternative = "one.sided",
#'                   design = "paired")
#'
#' ## parametric
#' power.t.student(d = 0.20,
#'                 margin = 0.05,
#'                 power = 0.80,
#'                 alternative = "one.sided",
#'                 design = "paired")
#'
#' ## mean of time 1 is practically greater than mean of time 2
#' ## targeting minimal difference of Cohen'd = -0.20 and can be as small as -0.05
#' ## non-parametric
#' power.np.wilcoxon(d = 0.20,
#'                   margin = -0.05,
#'                   power = 0.80,
#'                   alternative = "one.sided",
#'                   design = "paired")
#'
#' ## parametric
#' power.t.student(d = 0.20,
#'                 margin = -0.05,
#'                 power = 0.80,
#'                 alternative = "one.sided",
#'                 design = "paired")
#'
#'
#' ## mean of time 1 is practically same as mean of time 2
#' ## targeting minimal difference of Cohen'd = 0
#' ## and can be as small as -0.05 or as high as 0.05
#' ## non-parametric
#' ## non-parametric
#' power.np.wilcoxon(d = 0,
#'                   margin = c(-0.05, 0.05),
#'                   power = 0.80,
#'                   alternative = "two.one.sided",
#'                   design = "paired")
#'
#' ## parametric
#' power.t.student(d = 0,
#'                 margin = c(-0.05, 0.05),
#'                 power = 0.80,
#'                 alternative = "two.one.sided",
#'                 design = "paired")
#'
#' @export power.t.student
power.t.student <- function(d = NULL, null.d = 0, margin = 0, req.sign = "+",
                            n2 = NULL, n.ratio = 1, power = NULL, alpha = 0.05,
                            alternative = c("two.sided", "one.sided", "two.one.sided"),
                            design = c("independent", "paired", "one.sample"),
                            claim.basis = c("md.pval", "smd.ci"),
                            ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  design <- match.arg(design)
  claim.basis <- match.arg(claim.basis)
  func.parms <- as.list(environment())

  if (!is.null(d)) check.numeric(d)
  check.numeric(null.d)
  margin <- check.margins(margin, check.numeric, alternative)
  if (!is.null(n2)) check.sample.size(n2)
  check.positive(n.ratio)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = d, n = n2, power = power)

  pwr.student <- function(d, null.d, margin, n2, n.ratio,
                          alpha, alternative,
                          design, claim.basis) {
    n1 <- n2 * n.ratio

    if (design == "independent") {
      df <- n1 + n2 - 2
      var.d <- (n1 + n2) / (n1 * n2) + ifelse(claim.basis == "md.pval", 0, d ^ 2 / (2 * df))
    } else {
      df <- n2 - 1
      var.d <- 1 / n2 + ifelse(claim.basis == "md.pval", 0, d ^ 2 / (2 * n2)) # not df, like above?
    }

    se.d <- sqrt(var.d)
    lambda <- (d - null.d) / se.d
    null.lambda <- margin / se.d

    suppressWarnings(power.t.test(ncp = lambda, null.ncp = null.lambda, df = df, alpha = alpha, alternative = alternative,
                                  plot = FALSE, verbose = 0))[c("power", "t.alpha", "ncp", "null.ncp", "df")]

  } # pwr.student()

  min.pwr.student <- function(d, n2, power) {
    power - pwr.student(d = d, null.d = null.d, margin = margin, n2 = n2, n.ratio = n.ratio, alpha = alpha,
                        alternative = alternative, design = design, claim.basis = claim.basis)$power
  } # min.pwr.student (for uniroot)

  if (requested == "n") {

    n2 <- try(stats::uniroot(function(n2) min.pwr.student(d, n2, power), interval = c(4, 1e10))$root, silent = TRUE)
    if (inherits(n2, "try-error") || n2 == 1e10) stop("Design is not feasible. Possibly very large effect size.", call. = FALSE)

    n2 <- ifelse(ceil.n, ceiling(n2), n2)

  } else if (requested == "es") {

    if (alternative == "two.one.sided" && check.null_sign(req.sign, alternative)) {

      lower.int <- c(min(margin) + null.d, mean(margin) + null.d) + c(+1e-7, 0)
      upper.int <- c(mean(margin) + null.d, max(margin) + null.d) + c(0, -1e-7)
      d.lower <- stats::optimize(f = function(d) min.pwr.student(d, n2, power) ^ 2, interval = lower.int, tol = 1e-12)$minimum
      d.upper <- stats::optimize(f = function(d) min.pwr.student(d, n2, power) ^ 2, interval = upper.int, tol = 1e-12)$minimum

      d <- mean(c(d.lower, d.upper))

      pwr.lower <- pwr.student(d = d.lower, null.d = null.d, margin = margin, n2 = n2, n.ratio = n.ratio, alpha = alpha,
                               alternative = alternative, design = design, claim.basis = claim.basis)$power
      pwr.upper <- pwr.student(d = d.upper, null.d = null.d, margin = margin, n2 = n2, n.ratio = n.ratio, alpha = alpha,
                               alternative = alternative, design = design, claim.basis = claim.basis)$power

      if (round(pwr.lower, 3) >= power && round(pwr.upper, 3) >= power) {

        warning("Target effect ranges from ", round(d.lower, 4), " to ", round(d.upper, 4), " within the null bounds.", call. = FALSE)

      } else {

        warning("The target power rate cannot be achieved within the null bounds.", call. = FALSE)

      }

    } else {

      if (check.pos_sign(req.sign)) {
        d.int <- c(max(margin) + null.d, 10) + c(+1e-7, -1e-7)
      } else {
        d.int <- c(-10, min(margin) + null.d) + c(+1e-7, -1e-7)
      }

      d <- try(stats::uniroot(f = function(d) min.pwr.student(d, n2, power), interval = d.int, tol = 1e-12)$root, silent = TRUE)
      if (inherits(d, "try-error")) stop("Design is not feasible.", call. = FALSE)

    } # two.one.sided?

  } # ss or es?

  n1 <- ifelse(ceil.n, ceiling(n.ratio * n2), n.ratio * n2)
  if (design == "independent") n <- c(n1 = n1, n2 = n2) else n <- n2

  # calculate power (if requested == "power") or update it (if requested == "n" / "es")
  pwr.obj <- pwr.student(d = d, null.d = null.d, margin = margin, n2 = n2, n.ratio = n.ratio, alpha = alpha,
                         alternative = alternative, design = design, claim.basis = claim.basis)

  if (verbose > 0) {

    print.obj <- list(test = fmt_test_student(design),
                      requested = requested,
                      d = d,
                      null.d = null.d,
                      margin = margin,
                      alpha = alpha,
                      t.alpha = pwr.obj$t.alpha,
                      alternative = alternative,
                      n = n,
                      df = pwr.obj$df,
                      ncp = pwr.obj$ncp,
                      null.ncp = pwr.obj$null.ncp,
                      power = pwr.obj$power)

    .print.pwrss.ttest(print.obj, verbose = verbose, utf = utf)

  }

  invisible(structure(list(parms = func.parms,
                           test = "t",
                           df = pwr.obj$df,
                           ncp = pwr.obj$ncp,
                           null.ncp = pwr.obj$null.ncp,
                           t.alpha = pwr.obj$t.alpha,
                           d = d,
                           power = pwr.obj$power,
                           n = n,
                           n.total = sum(n)),
                      class = c("pwrss", "t", "student")))

} # power.t.student()

# provide a table of conversions from d to d.hc2
# mention that it can not be interpreted same as d
# when both var.ratio and n.ratio deviate from 1
# minor deviations are OK

#' Power Analysis for Welch's t-Test
#'
#' @description
#' Calculates power, sample size or effect size (only one can be NULL at a
#' time) for Welch's t-Tests. Welch's T-Test implementation relies on formulas
#' proposed by Bulus (2024).
#'
#' In contrast to previous versions, users can now specify whether their claims
#' will be based on raw score mean difference with p-values or standardized
#' mean difference with confidence intervals. While results typically differ by
#' only a few units, these distinctions can be particularly consequential in
#' studies with small sample sizes or high-risk interventions.
#'
#' Formulas are validated using Monte Carlo simulations (see Bulus, 2024),
#' G*Power, and tables in the PASS documentation. One key difference between
#' PASS and \code{pwrss} lies in how they handle non-inferiority and
#' superiority tests-that is, one-sided tests defined by a negligible effect
#' margin (implemented as of this version). PASS shifts the test statistic so
#' that the null hypothesis assumes a zero effect, treating the negligible
#' margin as part of the alternative hypothesis. As a result, the test
#' statistic is evaluated against a central distribution. In contrast,
#' \code{pwrss} treats the negligible effect as the true null value, and the
#' test statistic is evaluated under a non-central distribution. This leads to
#' slight differences up to third decimal place. To get the same results,
#' reflect the margin in \code{null.d} and specify \code{margin = 0}.
#'
#' Equivalence tests are implemented in line with Bulus and Polat (2023), Chow
#' et al. (2018) and Lakens (2017).
#'
#' @details
#' * Use \code{means.to.d()} to convert raw means and standard deviations to
#'   Cohen's d, and \code{d.to.cles()} to convert Cohen's d to the probability
#'   of superiority. Note that this interpretation is appropriate only when the
#'   underlying distribution is approximately normal and the two groups have
#'   similar population variances.
#' * NB: The functions \code{pwrss.z.mean()} and \code{pwrss.z.2means()} are no
#'   longer supported. The \code{pwrss.t.mean()} and \code{pwrss.t.2means()}
#'   functions are deprecated, but they will remain available as wrappers for
#'   \code{power.t.student()} or \code{power.t.welch()} during a transition
#'   period.
#'
#'
#' @param d           Cohen's d or Hedges' g.
#' @param null.d      Cohen's d or Hedges' g under null, typically 0 (zero).
#' @param margin      margin - ignorable \code{d} - \code{null.d} difference.
#' @param req.sign    whether `d` is smaller or larger than `null.d` (when
#'                    minimum detectable prob is of interest).
#' @param var.ratio   variance ratio in the form of sd1 ^ 2 / sd2 ^ 2.
#' @param n.ratio     \code{n1 / n2} ratio (applies to independent samples
#'                    only)
#' @param n2          integer; sample size in the second group (or for the
#'                    single group in paired samples or one-sample).
#' @param power       statistical power, defined as the probability of
#'                    correctly rejecting a false null hypothesis, denoted as
#'                    \eqn{1 - \beta}.
#' @param alpha       type 1 error rate, defined as the probability of
#'                    incorrectly rejecting a true null hypothesis, denoted as
#'                    \eqn{\alpha}.
#' @param alternative character; the direction or type of the hypothesis test:
#'                    "two.sided", "one.sided", or "two.one.sided". For
#'                    non-inferiority or superiority tests, add or subtract the
#'                    margin from the null hypothesis value and use
#'                    \code{alternative = "one.sided"}.
#' @param claim.basis character; "md.pval" when claims are based on raw mean
#'                    differences and p-values, "smd.ci" when claims are based
#'                    on standardized mean differences and confidence
#'                    intervals.
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
#'   \item{test}{type of the statistical test (T-Test).}
#'   \item{df}{degrees of freedom.}
#'   \item{ncp}{non-centrality parameter for the alternative.}
#'   \item{null.ncp}{non-centrality parameter for the null.}
#'   \item{t.alpha}{critical value(s).}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample size (`n` or `c(n1, n2)`.}
#'
#' @references
#'   Bulus, M. (2024). *Robust standard errors and confidence intervals for
#'   standardized mean differences*. https://doi.org/10.31219/osf.io/k6mbs
#'
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#'   Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). *Sample size
#'   calculations in clinical research* (3rd ed.). Taylor & Francis/CRC.
#'
#'   Cohen, J. (1988). *Statistical power analysis for the behavioral sciences*
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#'   Lakens, D. (2017). Equivalence tests: A practical primer for t tests,
#'   correlations, and meta-analyses. *Social psychological and personality
#'   science, 8*(4), 355-362. https://doi.org/10.1177/1948550617697177
#'
#' @examples
#' # see `?pwrss::power.t.student` for examples
#'
#' @export power.t.welch
power.t.welch <- function(d = NULL, null.d = 0, margin = 0, req.sign = "+",
                          var.ratio = 1, n.ratio = 1, n2 = NULL,
                          power = NULL, alpha = 0.05,
                          alternative = c("two.sided", "one.sided", "two.one.sided"),
                          claim.basis = c("md.pval", "smd.ci"),
                          ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  claim.basis <- match.arg(claim.basis)
  func.parms <- as.list(environment())

  if (!is.null(d)) check.numeric(d)
  check.numeric(null.d, var.ratio)
  margin <- check.margins(margin, check.numeric, alternative)
  check.positive(n.ratio)
  if (!is.null(n2)) check.sample.size(n2)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = d, n = n2, power = power)

  pwr.welch <- function(d, null.d, margin, var.ratio, n2, n.ratio, alpha, alternative, claim.basis) {

    n1 <- n.ratio * n2

    # variance ratio constraint
    sd2  <- sqrt((n2 * (n.ratio + 1) - 2) / (n2 - 1 + var.ratio * (n2 * n.ratio - 1)))
    sd1 <- sqrt(var.ratio * sd2 ^ 2)

    # Welch df
    df <- (sd1 ^ 2 / n1 + sd2 ^ 2 / n2) ^ 2 / (sd1 ^ 4 / (n1 ^ 2 * (n1 - 1)) + sd2 ^ 4 / (n2 ^ 2 * (n2 - 1)))

    var.delta <- 1 / n1 * sd1 ^ 2 + 1 / n2 * sd2 ^ 2
    std <- sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
    var.std <- ((n1 - 1) * sd1 ^ 4 + (n2 - 1) * sd2 ^ 4) / ((2 * (n1 - 1) * sd1 ^ 2 + 2 * (n2 - 1) * sd2 ^ 2) * (n1 + n2 - 2))

    if (claim.basis == "md.pval") {
      var.d.hc2 <- var.delta
    } else {
      var.d.hc2 <- var.delta / std ^ 2 + d ^ 2 * var.std / std ^ 2
    }

    se.d <- sqrt(var.d.hc2)
    lambda <- (d - null.d) / se.d
    null.lambda <- margin / se.d

    suppressWarnings(power.t.test(ncp = lambda, null.ncp = null.lambda, df = df, alpha = alpha, alternative = alternative,
                                  plot = FALSE, verbose = 0))[c("power", "t.alpha", "ncp", "null.ncp", "df")]

  } # pwr.welch()

  min.pwr.welch <- function(d, n2, power) {
    power - pwr.welch(d = d, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                      alpha = alpha, alternative = alternative, claim.basis = claim.basis)$power
  } # min.pwr.welch (for uniroot)

  if (requested == "n") {

    n2 <- try(stats::uniroot(function(n2) min.pwr.welch(d, n2, power), interval = c(3, 1e10))$root, silent = TRUE)
    if (inherits(n2, "try-error") || n2 == 1e10) stop("Design is not feasible. Possibly very large effect size.", call. = FALSE)

    n2 <- ifelse(ceil.n, ceiling(n2), n2)

  } else if (requested == "es") {

    if (alternative == "two.one.sided" && check.null_sign(req.sign, alternative)) {

      lower.int <- c(min(margin) + null.d, mean(margin) + null.d) + c(+1e-7, 0)
      upper.int <- c(mean(margin) + null.d, max(margin) + null.d) + c(0, -1e-7)
      d.lower <- stats::optimize(f = function(d) min.pwr.welch(d, n2, power) ^ 2, interval = lower.int, tol = 1e-12)$minimum
      d.upper <- stats::optimize(f = function(d) min.pwr.welch(d, n2, power) ^ 2, interval = upper.int, tol = 1e-12)$minimum

      d <- mean(c(d.lower, d.upper))

      pwr.lower <- pwr.welch(d = d.lower, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                             alpha = alpha, alternative = alternative, claim.basis = claim.basis)$power
      pwr.upper <- pwr.welch(d = d.upper, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                             alpha = alpha, alternative = alternative, claim.basis = claim.basis)$power

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

      d <- try(stats::uniroot(f = function(d) min.pwr.welch(d, n2, power), interval = d.int, tol = 1e-12)$root, silent = TRUE)
      if (inherits(d, "try-error")) stop("Design is not feasible.", call. = FALSE)

    } # two.one.sided?

  } # ss or es?

  n1 <- ifelse(ceil.n, ceiling(n.ratio * n2), n.ratio * n2)
  n <- c(n1 = n1, n2 = n2)

  # calculate power (if requested == "power") or update it (if requested == "n")
  pwr.obj <- pwr.welch(d = d, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                       alpha = alpha, alternative = alternative, claim.basis = claim.basis)

  if (verbose > 0) {

    print.obj <- list(requested = requested,
                      test = "Welch's T-Test (Independent Samples)",
                      d = d, null.d = null.d,
                      margin = margin,
                      alpha = alpha,
                      t.alpha = pwr.obj$t.alpha,
                      alternative = alternative,
                      n = n,
                      df = pwr.obj$df,
                      ncp = pwr.obj$ncp,
                      null.ncp = pwr.obj$null.ncp,
                      power = pwr.obj$power)

    .print.pwrss.ttest(print.obj, verbose = verbose, utf = utf)

  }

  invisible(structure(list(parms = func.parms,
                           test = "t",
                           df = pwr.obj$df,
                           ncp = pwr.obj$ncp,
                           null.ncp = pwr.obj$null.ncp,
                           t.alpha = pwr.obj$t.alpha,
                           d = d,
                           power = pwr.obj$power,
                           n = n,
                           n.total = sum(n)),
                      class = c("pwrss", "t", "welch")))

} # power.t.welch()




###################
# one mean t test #
###################

#' @export pwrss.t.mean
pwrss.t.mean <- function(mu, sd = 1, mu0 = 0, margin = 0, alpha = 0.05,
                         alternative = c("not equal", "greater", "less",
                                         "equivalent", "non-inferior", "superior"),
                         n = NULL, power = NULL, verbose = TRUE) {

  alternative <- match.arg(alternative)
  verbose <- ensure.verbose(verbose)

  check.positive(sd)
  check.numeric(mu, mu0, margin)
  if (!is.null(power)) check.power(power)
  if (!is.null(n)) check.sample.size(n)

  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"

  d <- (mu - mu0) / sd

  if (alternative == "two.one.sided") {

    margin.d.lower <- - margin / sd
    margin.d.upper <- + margin / sd
    margin.d <- c(margin.d.lower, margin.d.upper)

  } else if (alternative == "two.sided") {

    if (margin != 0) warning("Margin is forced to be 0 for the 'two.sided' test.", call. = FALSE)
    margin.d <- 0

  } else {

    margin.d <- margin / sd

  }

  student.obj <- power.t.student(d = d, margin = margin.d,
                                 n2 = n, power = power, alpha = alpha,
                                 alternative = alternative,
                                 design = "one.sample",
                                 claim.basis = "md.pval",
                                 ceil.n = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.t.student() function. \n")

  invisible(student.obj)

} # pwrss.t.mean()

#' @export pwrss.t.2means
pwrss.t.2means <- function(mu1, mu2 = 0, margin = 0,
                            sd1 = ifelse(paired, sqrt(1 / (2 * (1 - paired.r))), 1),
                            sd2 = sd1, kappa = 1, paired = FALSE, paired.r = 0.50,
                            alpha = 0.05, welch.df = TRUE,
                            alternative = c("not equal", "greater", "less",
                                            "equivalent", "non-inferior", "superior"),
                            n2 = NULL, power = NULL, verbose = TRUE) {

  alternative <- match.arg(alternative)
  verbose <- ensure.verbose(verbose)

  if (isFALSE(welch.df)) warning("Forcing welch.df = TRUE.", call. = FALSE)

  check.positive(sd1, sd2)
  check.correlation(paired.r)
  check.numeric(mu1, mu2, margin)
  if (!is.null(power)) check.power(power)
  if (!is.null(n2)) check.sample.size(n2)

  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"

  calculate.d <- function(mu1, mu2, sd1, sd2, n2,
                     n.ratio, paired, rho.paired) {
    if (paired) {
      psd <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * sd1 * sd2 * rho.paired)
    } else {
      n1 <- n.ratio * n2
      psd <- sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
    }

    (mu1 - mu2) / psd
  }

  # population d
  d <- calculate.d(mu1 = mu1, mu2 = mu2,
                   sd1 = sd1, sd2 = sd2,
                   n2 = 1e10, n.ratio = kappa,
                   paired = paired, rho.paired = paired.r)

  if (alternative == "two.one.sided") {

    margin.d.lower <- calculate.d(mu1 = -margin, mu2 = 0,
                             sd1 = sd1, sd2 = sd2,
                             n2 = 1e5, n.ratio = kappa,
                             paired = paired, rho.paired = paired.r)
    margin.d.upper <- calculate.d(mu1 = margin, mu2 = 0,
                             sd1 = sd1, sd2 = sd2,
                             n2 = 1e5, n.ratio = kappa,
                             paired = paired, rho.paired = paired.r)
    margin.d <- c(margin.d.lower, margin.d.upper)

  } else if (alternative == "two.sided") {

    if (margin != 0) warning("Margin is forced to be 0 for the 'two.sided' test.", call. = FALSE)
    margin.d <- 0

  } else {

    margin.d <- calculate.d(mu1 = margin, mu2 = 0,
                            sd1 = sd1, sd2 = sd2,
                            n2 = 1e5, n.ratio = kappa,
                            paired = paired, rho.paired = paired.r)

  }

  if (paired) {
    t.obj <- power.t.student(d = d, margin = margin.d,
                             n2 = n2, n.ratio = kappa,
                             power = power, alpha = alpha,
                             alternative = alternative,
                             design = "paired",
                             claim.basis = "md.pval",
                             ceil.n = TRUE, verbose = verbose)
  } else {
    t.obj <- power.t.welch(d = d, margin = margin.d,
                           var.ratio = sd1 ^ 2 / sd2 ^ 2,
                           n2 = n2, n.ratio = kappa,
                           power = power, alpha = alpha,
                           alternative = alternative,
                           claim.basis = "md.pval",
                           ceil.n = TRUE, verbose = verbose)
  }

  # cat("This function will be removed in the future. \n Please use power.t.student() or power.t.welch() function. \n")

  invisible(t.obj)

} # pwrss.t.2means()

fmt_test_student <- function(design) {
  sprintf("Student's T-Test (%s)",
          switch(design, `independent` = "Independent Samples", `paired` = "Paired Samples", `one.sample` = "One Sample"))
}

# defunct
#' @export pwrss.z.mean
pwrss.z.mean   <- function(...) {
  stop("This function is no longer available. Please use `power.t.student()`.", call. = FALSE)
}
#' @export pwrss.z.2means
pwrss.z.2means <- function(...) {
  stop("This function is no longer available. Please use `power.t.student()` or `power.t.welch()`.", call. = FALSE)
}
