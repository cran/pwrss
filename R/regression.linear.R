############################
# linear regression f test #
############################

#' Power Analysis for Linear Regression: R-squared or R-squared Change (F-Test)
#'
#' @description
#' Calculates power or sample size (only one can be NULL at a time) to test
#' R-squared deviation from 0 (zero) in linear regression or to test R-squared
#' change between two linear regression models. The test of R-squared change is
#' often used to evaluate incremental contribution of a set of predictors in
#' hierarchical linear regression.
#'
#' Formulas are validated using Monte Carlo simulation, G*Power, and tables in
#' the PASS documentation.
#'
#' @details
#' * NB: The \code{pwrss.f.regression} function and its alias
#'   \code{pwrss.f.reg()} are deprecated, but they will remain available as a
#'   wrapper for the \code{power.f.regression()} function during a transition
#'   period.
#'
#' @aliases pwrss.f.reg power.f.reg pwrss.f.regression power.f.regression
#'
#'
#' @param r.squared.change R-squared (or R-squared change).
#' @param margin           margin - ignorable R-squared (or R-squared change).
#' @param k.total          integer; total number of predictors.
#' @param k.tested         integer; number of predictors in the subset of
#'                         interest. By default \code{k.tested = k.total},
#'                         which implies that one is interested in the
#'                         contribution of all predictors, and tests whether
#'                         R-squared value is different from 0 (zero).
#' @param n                integer; sample size.
#' @param power            statistical power, defined as the probability of
#'                         correctly rejecting a false null hypothesis, denoted
#'                         as \eqn{1 - \beta}.
#' @param alpha            type 1 error rate, defined as the probability of
#'                         incorrectly rejecting a true null hypothesis,
#'                         denoted as \eqn{\alpha}.
#' @param ceil.n           logical; whether sample size should be rounded up.
#'                         \code{TRUE} by default.
#' @param verbose          \code{1} by default (returns test, hypotheses, and
#'                         results), if \code{2} a more detailed output is
#'                         given (plus key parameters and definitions), if
#'                         \code{0} no output is printed on the console.
#' @param utf              logical; whether the output should show Unicode
#'                         characters (if encoding allows for it).
#'                         \code{FALSE} by default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test (F-Test).}
#'   \item{df1}{numerator degrees of freedom.}
#'   \item{df2}{denominator degrees of freedom.}
#'   \item{ncp}{non-centrality parameter for the alternative.}
#'   \item{null.ncp}{non-centrality parameter for the null.}
#'   \item{f.alpha}{critical value.}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample size.}
#'
#' @references
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#'   Cohen, J. (1988). Statistical power analysis for the behavioral sciences
#'   (2nd ed.). Lawrence Erlbaum Associates.
#'
#' @examples
#'
#' # in the outcome (R-squared = 0.15).
#' power.f.regression(r.squared = 0.15,
#'                    k.total = 3, # total number of predictors
#'                    power = 0.80)
#'
#' # adding two more variables will increase R-squared
#' # from 0.15 (with 3 predictors) to 0.25 (with 3 + 2 predictors)
#' power.f.regression(r.squared.change = 0.10, # R-squared change
#'                    k.total = 5, # total number of predictors
#'                    k.tested = 2, # predictors to be tested
#'                    power = 0.80)
#'
#' @export power.f.regression
power.f.regression <- function(r.squared.change = NULL,
                               margin = 0,
                               k.total,
                               k.tested = k.total,
                               n = NULL, power = NULL, alpha = 0.05,
                               ceil.n = TRUE, verbose = 1, utf = FALSE) {

  func.parms <- as.list(environment())

  if (!is.null(r.squared.change)) check.proportion(r.squared.change)
  check.numeric(margin)
  check.positive(k.total, k.tested)
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = r.squared.change, n = n, power = power)

  if (k.tested >  k.total)
    stop("`k.tested` can not be greater than `k.total`.", call. = FALSE)
  if (!is.null(r.squared.change) && (r.squared.change == 0 || r.squared.change == 1))
    stop("Value for `r.squared.change` must be larger than 0 and smaller than 1.", call. = FALSE)

  pwr.f.reg <- function(r.squared.change, margin, k.total, k.tested, n, alpha) {

    df1 <- k.tested
    df2 <- n - k.total - 1

    f.squared <- rsq.to.f(r.squared.full = r.squared.change, verbose = 0)$f.squared
    null.f.squared <- rsq.to.f(r.squared.full = margin, verbose = 0)$f.squared

    lambda <- f.squared * n
    null.lambda <- null.f.squared * n


    f.alpha <- stats::qf(alpha, df1 = df1, df2 = df2, ncp = null.lambda, lower.tail = FALSE)
    power <- stats::pf(f.alpha, df1 = df1, df2 = df2, ncp = lambda, lower.tail = FALSE)

    list(power = power, lambda = lambda, null.lambda = null.lambda,
         df1 = df1, df2 = df2, f.alpha = f.alpha)

  } # pwr.f.reg()

  min.pwr <- function(r.squared.change, n, power) {

    power - suppressWarnings(pwr.f.reg(r.squared.change = r.squared.change, margin = margin, k.total = k.total,
                                       k.tested = k.tested, n = n, alpha = alpha))$power

  } # min.pwr (for uniroot)

  if (requested == "n") {

    n <- try(stats::uniroot(function(n) min.pwr(r.squared.change, n, power), interval = c(k.total + 2, 1e10))$root, silent = TRUE)
    if (inherits(n, "try-error") || n == 1e10) stop("Design is not feasible.", call. = FALSE)

    if (ceil.n) n <- ceiling(n)

  } else if (requested == "es") {

    r.squared.change <- try(stats::uniroot(function(r.squared.change) min.pwr(r.squared.change, n, power),
                                           interval = c(0, 1 - 1e-8), tol = 1e-12)$root,
                            silent = TRUE)
    if (inherits(r.squared.change, "try-error")) stop("Design is not feasible.", call. = FALSE)

  }

  # calculate power (if requested == "power") or update it (if requested == "n" / "es")
  pwr.obj <- pwr.f.reg(r.squared.change = r.squared.change, margin = margin,
                       k.total = k.total, k.tested = k.tested, n = n, alpha = alpha)

  if (verbose > 0) {

    print.obj <- list(test = paste0(ifelse(k.tested == k.total, "", "Hierarchical "), "Linear Regression (F-Test)"),
                      requested = requested,
                      r.squared.change = r.squared.change,
                      k.total = k.total,
                      k.tested = k.tested,
                      margin = margin,
                      n = n,
                      df1 = pwr.obj$df1,
                      df2 = pwr.obj$df2,
                      ncp = pwr.obj$lambda,
                      null.ncp = pwr.obj$null.lambda,
                      f.alpha = pwr.obj$f.alpha,
                      alpha = alpha,
                      power = pwr.obj$power)

    .print.pwrss.f.regression(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = "F",
                           df1 = pwr.obj$df1,
                           df2 = pwr.obj$df2,
                           ncp = pwr.obj$lambda,
                           null.ncp = pwr.obj$null.lambda,
                           f.alpha = pwr.obj$f.alpha,
                           r.squared.change = r.squared.change,
                           power = pwr.obj$power,
                           n = n),
                      class = c("pwrss", "f", "regression")))

} # pwrss.f.regression()

#' @rdname power.f.regression
#' @export power.f.reg
power.f.reg <- power.f.regression

#' @export pwrss.f.regression
pwrss.f.regression <- function(r2 = NULL, f2 = NULL,
                        k = 1, m = k, alpha = 0.05,
                        n = NULL, power = NULL, verbose = TRUE) {

  verbose <- ensure.verbose(verbose)

  if (all(check.not_null(r2, f2))) {
    stop("Effect size conflict for the alternative. Specify only either `r2` or `f2`.", call. = FALSE)
  } else if (check.not_null(f2)) {
    r2 <- f.to.rsq(f = sqrt(f2), verbose = 0)$r.squared.full
  }

  pwrss.f.reg.obj <- power.f.regression(r.squared.change = r2, margin = 0,
                                        k.total = k, k.tested = m,
                                        n = n, power = power, alpha = alpha,
                                        ceil.n = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.f.regression() function. \n")

  invisible(pwrss.f.reg.obj)

} # pwrss.f.regression

#' @export pwrss.f.reg
pwrss.f.reg <- pwrss.f.regression


#####################
# linear regression #
#####################

# if the predictor is binary provide sd.predictor = sqrt(p * (1 - p))
# p = proportion of subjects in treatment group
# use defaults if beta is standardized

#' Power Analysis for Linear Regression: Single Coefficient (T-Test)
#'
#' @description
#' Calculates power, sample size or effect size (only one can be NULL at a
#' time) to test a single coefficient in multiple linear regression. By
#' default, the predictor is assumed to be continuous. However, one can
#' calculate power, sample size or effect size for a binary predictor (such
#' as treatment and control groups in an experimental design) by specifying
#' \code{sd.predictor = sqrt(p * (1 - p))} where \code{p} is the proportion
#' of subjects in one of the groups. The sample size in each group would be
#' \code{n * p} and \code{n * (1 - p)}.
#'
#' Minimal effect and equivalence tests are implemented in line with Hodges and
#' Lehmann (1954), Kim and Robinson (2019), Phillips (1990), and Dupont and
#' Plummer (1998).
#'
#' Formulas are validated using Monte Carlo simulation, G*Power, tables in the
#' PASS documentation, and tables in Bulus (2021).
#'
#' @details
#' * When it is requested to calculate the effect size (by giving both their
#'   parameters `n` and `power`), `r.squared` must be empty and `beta` can be
#'   empty (`NULL`). By default, a linear regression with one predictor is
#'   assumed. If `beta`, `sd.predictor` and `sd.outcome` are given (i.e., in
#'   cases where it is not requested to calculate the effect size), `r.squared`
#'   is calculated as follows:
#'   \code{r.squared = (beta * sd.predictor / sd.outcome) ^ 2}. If the given
#'   `beta` results in an `r.squared` below this value, a warning will be
#'   issued. To calculate `beta` from `r.squared`, the following formula can be
#'   used: \code{beta = (sqrt(r.squared) * sd.outcome / sd.predictor)}.
#' * To consider other covariates in the model provide a value greater than the
#'   default value for \code{r.squared} (see above) along with the argument
#'   \code{k.total > 1}. However, in such case, the above formula for
#'   calculating `beta` can not be used.
#' * \code{power.t.regression()}, \code{pwrss.t.regression()},
#'   \code{power.t.reg()} and \code{pwrss.t.reg()} are the same functions.
#' * NB: The \code{pwrss.z.regression()} function and its alias
#'   \code{pwrss.z.reg()} are deprecated, but they will remain available as a
#'   wrapper for the \code{power.t.regression()} function during a transition
#'   period.
#'
#'
#' @aliases power.t.regression pwrss.t.regression pwrss.z.regression
#' @aliases pwrss.t.reg pwrss.z.reg power.t.reg
#'
#' @param beta         regression coefficient. One can use standardized
#'                     regression coefficient, but should keep
#'                     \code{sd.predictor = 1} and \code{sd.outcome = 1} or
#'                     leave them out as they are default specifications.
#' @param null.beta    regression coefficient under null hypothesis (typically
#'                     zero). One can use standardized regression coefficient,
#'                     but should keep \code{sd.predictor = 1} and
#'                     \code{sd.outcome = 1} or leave them out as they are
#'                     default specifications.
#' @param margin       margin - ignorable \code{beta} - \code{null.beta}
#'                     difference.
#' @param req.sign     sign of the regression coefficient `beta` (when minimum
#'                     detectable effect or `beta` is of interest).
#' @param sd.predictor standard deviation of the predictor. For a binary
#'                     predictor, \code{sd.predictor = sqrt(p * (1 - p))} where
#'                     \code{p} is the proportion of subjects in one of the
#'                     groups.
#' @param sd.outcome   standard deviation of the outcome.
#' @param k.total      integer; total number of predictors, including the
#'                     predictor of interest.
#' @param r.squared    model R-squared. Please see also Details below.
#' @param n            integer; sample size.
#' @param power        statistical power, defined as the probability of
#'                     correctly rejecting a false null hypothesis, denoted as
#'                     \eqn{1 - \beta}.
#' @param alpha        type 1 error rate, defined as the probability of
#'                     incorrectly rejecting a true null hypothesis, denoted as
#'                     \eqn{\alpha}.
#' @param alternative  character; the direction or type of the hypothesis test:
#'                     "two.sided", "one.sided", or "two.one.sided".
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
#'   \item{test}{type of the statistical test (T-Test).}
#'   \item{df}{degrees of freedom.}
#'   \item{ncp}{non-centrality parameter for the alternative.}
#'   \item{null.ncp}{non-centrality parameter for the null.}
#'   \item{t.alpha}{critical value(s).}
#'   \item{r.squared}{model R-squared.}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample size.}
#'
#' @references
#'   Bulus, M. (2021). Sample size determination and optimal design of
#'   randomized / non-equivalent pretest-post-test control-group designs.
#'   *Adiyaman University Journal of Educational Sciences, 11*(1), 48-69.
#'   https://doi.org/10.17984/adyuebd.941434
#'
#'   Hodges Jr, J. L., & Lehmann, E. L. (1954). Testing the approximate
#'   validity of statistical hypotheses. *Journal of the Royal Statistical
#'   Society Series B: Statistical Methodology, 16*(2), 261-268.
#'   https://doi.org/10.1111/j.2517-6161.1954.tb00169.x
#'
#'   Kim, J. H., & Robinson, A. P. (2019). Interval-based hypothesis testing
#'   and its applications to economics and finance. *Econometrics, 7*(2), 21.
#'   https://doi.org/10.1111/10.3390/econometrics7020021
#'
#'   Phillips, K. F. (1990). Power of the two one-sided tests procedure in
#'   bioequivalence. *Journal of Pharmacokinetics and Biopharmaceutics, 18*(2),
#'   137-144. https://doi.org/10.1007/bf01063556
#'
#'   Dupont, W. D., and Plummer, W. D. (1998). Power and sample size
#'   calculations for studies involving linear regression. *Controlled Clinical
#'   Trials, 19*(6), 589-601.
#'   https://doi.org/10.1007/10.1016/s0197-2456(98)00037-3
#'
#' @examples
#'
#' # continuous predictor x (and 4 covariates)
#' power.t.regression(beta = 0.20,
#'             k.total = 5,
#'             r.squared = 0.30,
#'             power = 0.80)
#'
#' # binary predictor x (and 4 covariates)
#' p <- 0.50 # proportion of subjects in one group
#' power.t.regression(beta = 0.20,
#'             sd.predictor = sqrt(p * (1 - p)),
#'             k.total = 5,
#'             r.squared = 0.30,
#'             power = 0.80)
#'
#' # non-inferiority test with binary predictor x (and 4 covariates)
#' p <- 0.50 # proportion of subjects in one group
#' power.t.regression(beta = 0.20, # Cohen's d
#'             margin = -0.05, # non-inferiority margin in Cohen's d unit
#'             alternative = "one.sided",
#'             sd.predictor = sqrt(p*(1-p)),
#'             k.total = 5,
#'             r.squared = 0.30,
#'             power = 0.80)
#'
#' # superiority test with binary predictor x (and 4 covariates)
#' p <- 0.50 # proportion of subjects in one group
#' power.t.regression(beta = 0.20, # Cohen's d
#'             margin = 0.05, # superiority margin in Cohen's d unit
#'             alternative = "one.sided",
#'             sd.predictor = sqrt(p * (1 - p)),
#'             k.total = 5,
#'             r.squared = 0.30,
#'             power = 0.80)
#'
#' # equivalence test with binary predictor x (and 4 covariates)
#' p <- 0.50 # proportion of subjects in one group
#' power.t.regression(beta = 0, # Cohen's d
#'             margin = c(-0.05, 0.05), # equivalence bounds in Cohen's d unit
#'             alternative = "two.one.sided",
#'             sd.predictor = sqrt(p * (1 - p)),
#'             k.total = 5,
#'             r.squared = 0.30,
#'             power = 0.80)
#'
#' @export power.t.regression
power.t.regression <- function(beta = NULL, null.beta = 0, margin = 0, req.sign = "+",
                               sd.predictor = 1, sd.outcome = 1,
                               r.squared = NULL,
                               k.total = 1,
                               n = NULL, power = NULL, alpha = 0.05,
                               alternative = c("two.sided", "one.sided", "two.one.sided"),
                               ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  func.parms <- as.list(environment())

  if (!is.null(beta)) check.numeric(beta)
  check.numeric(null.beta)
  margin <- check.margins(margin, check.numeric, alternative)
  check.positive(sd.predictor, sd.outcome)
  if (!is.null(r.squared)) check.proportion(r.squared)
  check.positive(k.total)
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)

  if (is.null(r.squared) && !is.null(beta))
    r.squared <- (beta * sd.predictor / sd.outcome) ^ 2
  if (!is.null(r.squared) && !is.null(beta) && r.squared > 0 && r.squared < (beta * sd.predictor / sd.outcome) ^ 2)
    warning("`r.squared` is possibly larger.", call. = FALSE)
# if (is.null(r.squared) && k.total > 1)
#  warning("When requesting to calculate the effect size, `r.squared` is calculated assuming only one predictor. ",
#          "With several predictors, `beta` should not be calculated using the formula under Details in the help ",
#          "for this function.", call. = FALSE)

  # NB: Needs more careful consideration, how the different options (beta, k.total) will affect r.squared
  requested <- get.requested(es = beta, n = n, power = power)

  pwr.t.reg <- function(beta, null.beta, margin, sd.outcome, sd.predictor, n, k.total, r.squared, alpha, alternative) {

    df <- n - k.total - 1
    lambda <- (beta - null.beta)  / ((sd.outcome / sd.predictor) * sqrt((1 - r.squared) / n))
    null.lambda <- margin / ((sd.outcome / sd.predictor) * sqrt((1 - r.squared) / n))

    pwr.obj <- power.t.test(ncp = lambda, null.ncp = null.lambda, df = df, alpha = alpha,
                            alternative = alternative, plot = FALSE, verbose = 0)

    list(power = pwr.obj$power, lambda = lambda, null.lambda = null.lambda, df = df, t.alpha = pwr.obj$t.alpha)

  } # pwr.t.reg()

  min.pwr.t.reg <- function(beta, n, power) {

    # beta2min <- ifelse(!is.null(beta), beta, sqrt(r.squared) * sd.outcome / sd.predictor)

    power - suppressWarnings(pwr.t.reg(beta = beta, null.beta = null.beta, margin = margin, sd.outcome = sd.outcome,
                                       sd.predictor = sd.predictor, n = n, k.total = k.total, r.squared = r.squared,
                                       alpha = alpha, alternative =  alternative))$power

  } # min.pwr (for uniroot)

  if (requested == "n") {

    n <- try(stats::uniroot(function(n) min.pwr.t.reg(beta, n, power),
                            interval = c(k.total + 4, 1e10))$root,
             silent = TRUE)
    if (inherits(n, "try-error") || n == 1e10) stop("Design is not feasible.", call. = FALSE)

    if (ceil.n) n <- ceiling(n)

  } else if (requested == "es") {

    if (is.null(r.squared)) {
      r.squared <- 0
      warning("Setting r.squared = 0 to calculate the minimum detectable effect.", call. = FALSE)
    }

    if (alternative == "two.one.sided" && check.null_sign(req.sign, alternative)) {

      lower.int <- c(min(margin), mean(margin)) + c(+1e-7, 0)
      upper.int <- c(mean(margin), max(margin)) + c(0, -1e-7)
      beta.lower <- suppressWarnings(stats::optimize(f = function(beta) min.pwr.t.reg(beta, n, power) ^ 2, interval = lower.int, tol = 1e-12))$minimum
      beta.upper <- suppressWarnings(stats::optimize(f = function(beta) min.pwr.t.reg(beta, n, power) ^ 2, interval = upper.int, tol = 1e-12))$minimum

      beta <- mean(c(beta.lower, beta.upper))

      pwr.lower <- suppressWarnings(pwr.t.reg(beta = beta.lower, null.beta = null.beta, margin = margin, sd.outcome = sd.outcome,
                                              sd.predictor = sd.predictor, n = n, k.total = k.total, r.squared = r.squared,
                                              alpha = alpha, alternative =  alternative))$power
      pwr.upper <- suppressWarnings(pwr.t.reg(beta = beta.upper, null.beta = null.beta, margin = margin, sd.outcome = sd.outcome,
                                              sd.predictor = sd.predictor, n = n, k.total = k.total, r.squared = r.squared,
                                              alpha = alpha, alternative =  alternative))$power

      if (round(pwr.lower, 3) >= power && round(pwr.upper, 3) >= power) {

        warning("Target effect ranges from ", round(beta.lower, 4), " to ", round(beta.upper, 4), " within the null bounds.",
                call. = FALSE)

      } else {

        warning("The target power rate cannot be achieved within the null bounds.", call. = FALSE)

      }

    } else {

      if (check.pos_sign(req.sign)) {
        beta.int <- c(max(margin), 1e10) + c(+1e-7, -1e-7)
      } else {
        beta.int <- c(-1e10, min(margin)) + c(+1e-7, -1e-7)
      }

      beta <- try(stats::uniroot(f = function(beta) min.pwr.t.reg(beta, n, power), interval = beta.int, tol = 1e-12)$root, silent = TRUE)
      if (inherits(beta, "try-error")) stop("Design is not feasible.", call. = FALSE)

    } # two.one.sided?

    # r.squared <- try(stats::uniroot(function(r.squared) min.pwr(r.squared, beta, n, power),
    #                                interval = c(1e-8, 1 - 1e-8), tol = 1e-12)$root,
    #                 silent = TRUE)
    # if (inherits(r.squared, "try-error")) stop("Design is not feasible.", call. = FALSE)

    # beta <- ifelse(!is.null(beta), beta, sqrt(r.squared) * sd.outcome / sd.predictor)

  }

  # calculate power (if requested == "power") or update it (if requested == "n" / "es")
  pwr.obj <- pwr.t.reg(beta = beta, null.beta = null.beta, margin = margin,
                       sd.outcome = sd.outcome, sd.predictor = sd.predictor,
                       n = n, k.total = k.total, r.squared = r.squared,
                       alpha = alpha, alternative =  alternative)

  std.beta <- beta * (sd.predictor / sd.outcome)
  std.null.beta <- null.beta * (sd.predictor / sd.outcome)
  std.margin <- margin * (sd.predictor / sd.outcome)

  if (verbose > 0) {

    print.obj <- list(test = "Linear Regression Coefficient (T-Test)",
                      requested = requested,
                      tgt.effect = "std.beta",
                      alternative = alternative,
                      std.beta = std.beta,
                      std.null.beta = std.null.beta,
                      std.margin = std.margin,
                      margin = margin,
                      n = n,
                      df = pwr.obj$df,
                      ncp = pwr.obj$lambda,
                      null.ncp = pwr.obj$null.lambda,
                      alpha = alpha,
                      t.alpha = pwr.obj$t.alpha,
                      r.squared = r.squared,
                      power = pwr.obj$power)

    .print.pwrss.t.regression(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = "t",
                           std.beta = std.beta,
                           std.null.beta = std.null.beta,
                           std.margin = std.margin,
                           df = pwr.obj$df,
                           t.alpha = pwr.obj$t.alpha,
                           ncp = pwr.obj$lambda,
                           null.ncp = pwr.obj$null.lambda,
                           r.squared = r.squared,
                           power = pwr.obj$power,
                           n = n),
                      class = c("pwrss", "t", "regression")))
} # power.t.regression()

#' @rdname power.t.regression
#' @export power.t.reg
power.t.reg <- power.t.regression

#' @export pwrss.t.regression
pwrss.t.regression <- function(beta1 = 0.25, beta0 = 0, margin = 0,
                               sdx = 1, sdy = 1,
                               k = 1, r2 = (beta1 * sdx / sdy) ^ 2,
                               alpha = 0.05, n = NULL, power = NULL,
                               alternative = c("not equal", "less", "greater",
                                               "non-inferior", "superior", "equivalent"),
                               verbose = TRUE) {

  verbose <- ensure.verbose(verbose)
  alternative <- match.arg(alternative)
  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  if (alternative == "equivalent") {
    alternative <- "two.one.sided"
    if (length(margin) == 1) margin <- c(min(-margin, margin), max(-margin, margin))
  }

  pwrss.t.reg.obj <- power.t.regression(beta = beta1, null.beta = beta0, margin = margin,
                                        sd.predictor = sdx, sd.outcome = sdy,
                                        r.squared = r2, k.total = k,
                                        n = n, power = power, alpha = alpha,
                                        alternative = alternative,
                                        ceil.n = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.t.regression() function. \n")

  invisible(pwrss.t.reg.obj)

} # pwrss.t.regression

#' @export pwrss.t.reg
pwrss.t.reg <- pwrss.t.regression


# defunct
#' @export pwrss.z.regression
#' @export pwrss.z.reg
pwrss.z.regression <- pwrss.z.reg <- function(...) {
  stop("This function is no longer available. Please use `power.t.regression()`.", call. = FALSE)
}
