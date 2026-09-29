# f.squared <- eta.squared / (1 - eta.squared)
# eta.squared <- f.squared / (1 + f.squared)

# helper functions -----------------------------------------------------------------------------------------------------
# common function to calculate df2 for between-subject AN(C)OVAs
calc.df2 <- function(n.total, n.groups, k.covariates) n.total - (n.groups + k.covariates)

fmt_test_ancova <- function(n.way, cov) {
  paste0(c("One", "Two", "Three")[n.way], "-Way Analysis of ", ifelse(cov > 0, "Cov", "V"), "ariance (F-Test)")
}

check_var.ratio <- function(sd.vector, n.vector, alpha = 0.05) {
  sd.vector.squared <- sd.vector ^ 2
  var.ratio <- max(sd.vector.squared) / min(sd.vector.squared)
  n.max <- n.vector[sd.vector.squared == max(sd.vector.squared)][1]
  n.min <- n.vector[sd.vector.squared == min(sd.vector.squared)][1]

  if (var.ratio <= stats::qf(alpha,     df1 = n.max - 1, df2 = n.min - 1, lower.tail = TRUE) ||
      var.ratio >= stats::qf(1 - alpha, df1 = n.max - 1, df2 = n.min - 1, lower.tail = TRUE))
    warning("Interpretation of results may no longer be valid when variances differ beyond sampling error.", call. = FALSE)

  invisible(NULL)
}

# ----------------------------------------------------------------------------------------------------------------------
#' Power Analysis for One-, Two-, Three-Way ANOVA/ANCOVA Using Effect Size
#' (F-Test)
#'
#' @description
#' Calculates power, sample size or effect size for one-way, two-way, or
#' three-way ANOVA/ANCOVA. Set \code{k.covariates = 0} for ANOVA, and
#' \code{k.covariates > 0} for ANCOVA. Note that in the latter, the effect
#' size (\code{eta.squared} should be obtained from the relevant ANCOVA model,
#'  which is already adjusted for the explanatory power of covariates (thus,
#' an additional R-squared argument is not required as an input).
#'
#' Formulas are validated using G*Power and tables in the PASS documentation.
#'
#' @details
#' Note that R has a partial matching feature which allows you to specify
#' shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead
#' of \code{mu.vector}, or such as \code{k} or \code{k.covariates} instead of
#' \code{k.covariates}.
#'
#'
#' @aliases pwrss.f.ancova
#'
#' @param eta.squared      (partial) eta-squared for the alternative.
#' @param null.eta.squared (partial) eta-squared for the null.
#' @param factor.levels    integer; number of levels or groups in each factor.
#'                         For example, for two factors each having two levels
#'                         or groups use e.g. c(2, 2), for three factors each
#'                         having two levels (groups) use e.g. c(2, 2, 2).
#' @param target.effect    character; determine the main effect or interaction
#'                         that is of interest, e.g., in a three-way-design, it
#'                         is possible to define "A" (main effect of the first
#'                         factor), "B:C" (interaction of factor two and three)
#'                         or "A:B:C" (the three-way interaction of all
#'                         factors); if target is not used, the three-way
#'                         interaction is the default.
#' @param k.covariates     integer; number of covariates in an ANCOVA model
#' @param n.total          integer; total sample size
#' @param power            statistical power, defined as the probability of
#'                         correctly rejecting a false null hypothesis,
#'                         denoted as \eqn{1 - \beta}.
#' @param alpha            type 1 error rate, defined as the probability of
#'                         incorrectly rejecting a true null hypothesis,
#'                         denoted as \eqn{\alpha}.
#' @param ceil.n           logical; if \code{FALSE} sample size in each cell
#'                         is not rounded up.
#' @param verbose          \code{1} by default (returns test, hypotheses, and
#'                         results), if \code{2} a more detailed output is
#'                         given (plus key parameters and definitions), if
#'                         \code{0} no output is printed on the console.
#' @param utf              logical; whether the output should show Unicode
#'                         characters (if encoding allows for it). \code{FALSE}
#'                         by default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test (F-Test).}
#'   \item{eta.squared}{(partial) eta-squared for the alternative.}
#'   \item{df1}{numerator degrees of freedom.}
#'   \item{df2}{denominator degrees of freedom.}
#'   \item{ncp}{non-centrality parameter for the alternative.}
#'   \item{null.ncp}{non-centrality parameter for the null.}
#'   \item{f.alpha}{critical value.}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n.total}{total sample size.}
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
#' #############################################
#' #              one-way ANOVA                #
#' #############################################
#'
#' # Cohen's d = 0.50 between treatment and control
#' # translating into eta-squared = 0.059
#'
#' # estimate sample size using ANOVA approach
#' power.f.ancova(eta.squared = 0.059,
#'                factor.levels = 2,
#'                power = .80, alpha = 0.05)
#'
#' # estimate sample size using regression approach(F-Test)
#' power.f.regression(r.squared = 0.059,
#'                    k.total = 1,
#'                    power = 0.80, alpha = 0.05)
#'
#' # estimate sample size using regression approach (t-Test)
#' p <- 0.50 # proportion of sample in treatment (allocation rate)
#' power.t.regression(beta = 0.50, r.squared = 0,
#'                    k.total = 1,
#'                    sd.predictor = sqrt(p * (1 - p)),
#'                    power = 0.80, alpha = 0.05)
#'
#' # estimate sample size using t test approach
#' power.t.student(d = 0.50, power = 0.80, alpha = 0.05)
#'
#' #############################################
#' #              two-way ANOVA                #
#' #############################################
#'
#' # a researcher is expecting a partial eta-squared = 0.03
#' # for interaction of treatment (Factor A) with
#' # gender consisting of two levels (Factor B)
#'
#' power.f.ancova(eta.squared = 0.03,
#'                factor.levels = c(2,2),
#'                power = 0.80, alpha = 0.05)
#'
#' # estimate sample size using regression approach (F test)
#' # one dummy for treatment, one dummy for gender, and their interaction (k = 3)
#' # partial eta-squared is equivalent to the increase in R-squared by adding
#' # only the interaction term (m = 1)
#' power.f.regression(r.squared = 0.03,
#'                    k.total = 3, k.test = 1,
#'                    power = 0.80, alpha = 0.05)
#'
#' #############################################
#' #              one-way ANCOVA               #
#' #############################################
#'
#' # a researcher is expecting an adjusted difference of
#' # Cohen's d = 0.45 between treatment and control after
#' # controllling for the pretest (k.covariates = 1)
#' # translating into eta-squared = 0.048
#'
#' power.f.ancova(eta.squared = 0.048,
#'                factor.levels = 2,
#'                k.covariates = 1,
#'                power = .80, alpha = 0.05)
#'
#' #############################################
#' #              two-way ANCOVA               #
#' #############################################
#'
#' # a researcher is expecting an adjusted partial eta-squared = 0.02
#' # for interaction of treatment (Factor A) with
#' # gender consisting of two levels (Factor B)
#'
#' power.f.ancova(eta.squared = 0.02,
#'                factor.levels = c(2,2),
#'                k.covariates = 1,
#'                power = .80, alpha = 0.05)
#'
#' @export power.f.ancova
power.f.ancova <- function(eta.squared = NULL,
                           null.eta.squared = 0,
                           factor.levels = 2,
                           target.effect = NULL,
                           k.covariates = 0,
                           n.total = NULL,
                           power = NULL,
                           alpha = 0.05,
                           ceil.n = TRUE,
                           verbose = 1,
                           utf = FALSE) {

  func.parms <- as.list(environment())

  if (!is.null(eta.squared)) check.nonnegative(eta.squared)
  check.nonnegative(null.eta.squared, k.covariates)
  check.vector(factor.levels, check.factor.level, min.length = 1)
  if (!is.null(n.total)) check.sample.size(n.total)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = eta.squared, n = n.total, power = power)

  if (!is.null(eta.squared)) f.squared <- eta.squared / (1 - eta.squared)
  null.f.squared <- null.eta.squared / (1 - null.eta.squared)

  n.way <- length(factor.levels)
  if (n.way > 3)
    stop("More than three-way ANOVA or ANCOVA is not allowed at the moment.", call. = FALSE)

  n.groups <- prod(factor.levels)
  fac.letters <- LETTERS[seq_along(factor.levels)]
  effect <- paste0(fac.letters, "(", factor.levels, ")", collapse = ":")

  if (is.null(target.effect)) {
    df1 <- prod(factor.levels - 1)
  } else {
    if (all(strsplit(target.effect, ":", fixed = TRUE)[[1]] %in% fac.letters)) {
      fac.select <- fac.letters %in% strsplit(target.effect, ":", fixed = TRUE)[[1]]
      df1 <- prod(factor.levels[fac.select] - 1)
      effect <- paste(target.effect, "from", effect)
    } else {
      stop("Invalid specification of `target.effect`. It must be either a single letter \"A\", \"B\" or \"C\" ",
           "(depending on the length of `factor.levels`), assessing a main effect, or a combination of these ",
           "letters separated by \":\", e.g., \"A:B\", assessing an interaction.", call. = FALSE)
    }
  }

  pwr <- function(df1, n.total, n.groups, k.covariates, f.squared, null.f.squared, alpha) {
    df2 <- calc.df2(n.total, n.groups, k.covariates)
    lambda <- f.squared * n.total
    null.lambda <- null.f.squared * n.total
    f.alpha <- stats::qf(alpha, df1 = df1, df2 = df2, ncp = null.lambda, lower.tail = FALSE)
    power <- stats::pf(f.alpha, df1 = df1, df2 = df2, ncp = lambda, lower.tail = FALSE)

    list(power = power, lambda = lambda, null.lambda = null.lambda, f.alpha = f.alpha)
  } # pwr

  min.pwr <- function(f.squared, n.total, power) {
    power - pwr(df1, n.total, n.groups, k.covariates, f.squared, null.f.squared, alpha)$power
  } # min.pwr (for uniroot)

  if (requested == "n") {

    n.total <- try(stats::uniroot(function(n.total) min.pwr(f.squared, n.total, power),
                                  interval = c(n.groups + k.covariates + 2, 1e10))$root,
                   silent = TRUE)
    if (inherits(n.total, "try-error") || n.total == 1e10) stop("Design is not feasible.", call. = FALSE)

    if (ceil.n) n.total <- ceiling(n.total / n.groups) * n.groups

  } else if (requested == "es") {

    f.squared <- try(stats::uniroot(function(f.squared) min.pwr(f.squared, n.total, power),
                                    interval = c(0, 1), tol = 1e-12)$root,
                     silent = TRUE)
    if (inherits(f.squared, "try-error")) stop("Design is not feasible.", call. = FALSE)

    eta.squared <- f.squared / (1 + f.squared)

  }

  # calculate power (if requested == "power") or update it (if requested == "n" / "es")
  pwr.obj <- pwr(df1 = df1, n.total = n.total, n.groups = n.groups, k.covariates = k.covariates,
                 f.squared = f.squared, null.f.squared = null.f.squared, alpha = alpha)

  df2 <- calc.df2(n.total, n.groups, k.covariates)

  if (verbose > 0) {

    print.obj <- list(test = fmt_test_ancova(n.way, k.covariates),
                      effect = effect,
                      eta.squared = eta.squared,
                      n.total = n.total,
                      n.way = n.way,
                      requested = requested,
                      factor.levels = factor.levels,
                      power = pwr.obj$power,
                      ncp = pwr.obj$lambda,
                      null.ncp = pwr.obj$null.lambda,
                      alpha = alpha,
                      f.alpha = pwr.obj$f.alpha,
                      df1 = df1,
                      df2 = df2)

    .print.pwrss.ancova(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = "F",
                           effect = effect,
                           eta.squared = eta.squared,
                           df1 = df1,
                           df2 = df2,
                           ncp = pwr.obj$lambda,
                           null.ncp = pwr.obj$null.lambda,
                           f.alpha = pwr.obj$f.alpha,
                           power = pwr.obj$power,
                           n.total = n.total),
                      class = c("pwrss", "f", "ancova")))

} # end of power.f.ancova()

#' @export pwrss.f.ancova
pwrss.f.ancova <- function(eta2 = NULL, f2 = NULL,
                           n.way = length(n.levels),
                           n.levels = 2, n.covariates = 0, alpha = 0.05,
                           n = NULL, power = NULL, verbose = TRUE) {

  verbose <- ensure.verbose(verbose)

  if (all(check.not_null(f2, eta2))) {
    stop("Effect size conflict for the alternative. Specify only either `eta2` or `f2`.", call. = FALSE)
  } else if (check.not_null(f2)) {
    eta2 <- f2 / (1 + f2)
  }
  # eta2 doesn't need conversion, and falls back to the default if neither f2 nor eta2 is given explicitly

  ancova.obj <- power.f.ancova(eta.squared = eta2,
                               null.eta.squared = 0,
                               factor.levels = n.levels,
                               k.covariates = n.covariates,
                               n.total = n,
                               power = power,
                               alpha = alpha,
                               ceil.n = TRUE,
                               verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.f.ancova() function. \n")

  invisible(ancova.obj)

} # pwrss.f.ancova


#' Power Analysis for One-Way ANOVA/ANCOVA Using Means and Standard Deviations
#' (F test)
#'
#' @description
#' Calculates power, sample size or effect size for one-way ANOVA/ANCOVA. Set
#' \code{k.covariates = 0} for one-way ANOVA (without any pretest or covariate
#' adjustment). Set \code{k.covariates > 0} in combination with
#' \code{r.squared > 0} for one-way ANCOVA (with pretest or covariate
#' adjustment).
#'
#' Formulas are validated using the PASS documentation.
#'
#' @details
#' Note that R has a partial matching feature which allows you to specify
#' shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead
#' of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of
#' \code{k.covariates}.
#'
#' @param mu.vector     vector of adjusted means (or estimated marginal means)
#'                      for each level of a factor.
#' @param sd.vector     vector of unadjusted standard deviations for each level
#'                      of a factor.
#' @param n.vector      vector of sample sizes for each level of a factor.
#' @param p.vector      vector of proportion of total sample size in each level
#'                      of a factor. These proportions should sum to one.
#' @param factor.levels integer; number of levels or groups in each factor. For
#'                      example, for two factors each having two levels or
#'                      groups use e.g. c(2, 2), for three factors each having
#'                      two levels or groups use e.g. c(2, 2, 2)
#' @param r.squared     explanatory power of covariates (R-squared) in the
#'                      ANCOVA model. The default is \code{r.squared = 0},
#'                      which means an ANOVA model would be of interest.
#' @param k.covariates  integer; number of covariates in the ANCOVA model. The
#'                      default is \code{k.covariates = 0}, which means an
#'                      ANOVA model would be of interest.
#' @param power         statistical power, defined as the probability of
#'                      correctly rejecting a false null hypothesis, denoted as
#'                      \eqn{1 - \beta}.
#' @param alpha         type 1 error rate, defined as the probability of
#'                      incorrectly rejecting a true null hypothesis, denoted
#'                      as \eqn{\alpha}.
#' @param ceil.n        logical; whether sample size should be rounded up.
#'                      \code{TRUE} by default.
#' @param verbose       \code{1} by default (returns test, hypotheses, and
#'                      results), if \code{2} a more detailed output is given
#'                      (plus key parameters and definitions), if \code{0} no
#'                      output is printed on the console.
#' @param utf           logical; whether the output should show Unicode
#'                      characters (if encoding allows for it). \code{FALSE}
#'                      by default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test (F-Test).}
#'   \item{df1}{numerator degrees of freedom.}
#'   \item{df2}{denominator degrees of freedom.}
#'   \item{ncp}{non-centrality parameter under alternative.}
#'   \item{null.ncp}{non-centrality parameter under null.}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n.total}{total sample size.}
#'
#' @references
#'   Keppel, G., & Wickens, T. D. (2004). Design and analysis: A researcher's
#'   handbook (4th ed.). Pearson.
#'
#' @examples
#'
#' # required sample size to detect a mean difference of
#' # Cohen's d = 0.50 for a one-way two-group design
#' power.f.ancova.keppel(mu.vector = c(0.50, 0), # marginal means
#'                       sd.vector = c(1, 1), # unadjusted standard deviations
#'                       n.vector = NULL, # sample size (will be calculated)
#'                       p.vector = c(0.50, 0.50), # balanced allocation
#'                       k.covariates = 1, # number of covariates
#'                       r.squared = 0.50, # explanatory power of covariates
#'                       alpha = 0.05, # Type 1 error rate
#'                       power = .80)
#'
#' # effect size approach
#' power.f.ancova(eta.squared = 0.111, # effect size that is already adjusted for covariates
#'                factor.levels = 2, # one-way ANCOVA with two levels (groups)
#'                k.covariates = 1, # number of covariates
#'                alpha = 0.05, # Type 1 error rate
#'                power = .80)
#'
#' # regression approach
#' p <- 0.50
#' power.t.regression(beta = 0.50,
#'                    sd.predictor = sqrt(p * (1 - p)),
#'                    sd.outcome = 1,
#'                    k.total = 1,
#'                    r.squared = 0.50,
#'                    n = NULL, power = 0.80)
#'
#' @export power.f.ancova.keppel
power.f.ancova.keppel <- function(mu.vector,
                                  sd.vector,
                                  n.vector = NULL,
                                  p.vector = NULL,
                                  factor.levels = NULL,
                                  r.squared = 0,
                                  k.covariates = 0,
                                  power = NULL,
                                  alpha = 0.05,
                                  ceil.n = TRUE,
                                  verbose = 1,
                                  utf = FALSE) {

  func.parms <- as.list(environment())

  # value and consistency checks
  if (!is.vector(mu.vector) || !is.numeric(mu.vector))
    stop("Provide a vector of means (`mu.vector`) with its length equal to number of groups.", call. = FALSE)
  check.vector(mu.vector, check.numeric)
  check.vector(sd.vector, check.positive)
  if (!is.null(n.vector)) check.vector(n.vector, check.sample.size)
  if (!is.null(p.vector)) check.vector(p.vector, check.proportion)
  check.same.lengths(mu.vector, sd.vector, n.vector, p.vector)
  if (is.null(factor.levels)) factor.levels <- length(mu.vector)
  if (length(factor.levels) > 1)
    stop("Factorial designs are not allowed in Keppel's approach.", call. = FALSE)
  check.factor.level(factor.levels)
  if (length(mu.vector) != factor.levels)
    stop("Length of the vector of means (`mu.vector`) does not match number of levels.", call. = FALSE)
  if (r.squared > 1 || r.squared < 0 || !is.numeric(r.squared) || length(r.squared) != 1)
    stop("R-squared (explanatory power of covariates) takes a value between 0 and 1.", call. = FALSE)
  check.nonnegative(k.covariates)
  if (r.squared > 0 && k.covariates < 1)
    stop("Explanatory power of covariates is expected to be non-zero when number of covariates is non-zero.", call. = FALSE)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = NA, n = n.vector, power = power) # effect size can not be calculated

  if (requested == "n" && is.null(p.vector))
    stop("`p.vector` can not be NULL when sample size is requested", call. = FALSE)
  if (requested == "n" && round(sum(p.vector), 5) != 1)
    stop("The elements of the `p.vector` should sum to 1", call. = FALSE)

  pwr.keppel <- function(mu.vector, sd.vector, n.vector, k.covariates, r.squared, alpha, factor.levels) {

    n.total <- sum(n.vector)
    mu_bar <- sum(n.vector * mu.vector) / n.total

    sigma2_pooled <- sum((n.vector - 1) * sd.vector ^ 2) / (n.total - length(mu.vector))
    sigma2_between <- sum(n.vector * (mu.vector - mu_bar) ^ 2) / n.total
    sigma2_error <- sigma2_pooled * (1 - r.squared)

    f.squared <- sigma2_between / sigma2_error

    df1 <- prod(factor.levels - 1)
    df2 <- calc.df2(n.total, length(mu.vector), k.covariates)
    lambda <- f.squared * n.total

    f.alpha <- stats::qf(alpha, df1 = df1, df2 = df2, lower.tail = FALSE)
    power <- stats::pf(f.alpha, df1 = df1, df2 = df2, ncp = lambda, lower.tail = FALSE)

    list(power = power, f.squared = f.squared, df1 = df1, df2 = df2, lambda = lambda, f.alpha = f.alpha)

  }

  if (requested == "n") {

    n.total <- stats::uniroot(function(n.total) {
      power - pwr.keppel(mu.vector = mu.vector, sd.vector = sd.vector, n.vector = n.total * p.vector,
                         k.covariates = k.covariates, r.squared = r.squared, alpha = alpha,
                         factor.levels = factor.levels)$power
    }, interval = c(length(mu.vector) + k.covariates + 2, 1e10))$root

    n.vector <- n.total * p.vector
    if (ceil.n) n.vector <- ceiling(n.vector)

  }

  # calculate power (if requested == "power") or update it (if requested == "n")
  pwr.obj <- pwr.keppel(mu.vector = mu.vector, sd.vector = sd.vector, n.vector = n.vector, k.covariates = k.covariates,
                        r.squared = r.squared, alpha = alpha, factor.levels = factor.levels)

  n.total <- sum(n.vector)
  eta.squared <- pwr.obj$f.squared / (1 + pwr.obj$f.squared)
  effect <- paste0("A", "(", factor.levels, ")")
  n.way <- length(factor.levels)

  check_var.ratio(sd.vector, n.vector)

  if (verbose > 0) {

    print.obj <- list(test = fmt_test_ancova(1, k.covariates),
                      requested = requested,
                      effect = effect,
                      n.total = n.total,
                      n.way = n.way,
                      factor.levels = factor.levels,
                      eta.squared = eta.squared,
                      power = pwr.obj$power,
                      ncp = pwr.obj$lambda,
                      null.ncp = 0,
                      alpha = alpha,
                      f.alpha = pwr.obj$f.alpha,
                      df1 = pwr.obj$df1,
                      df2 = pwr.obj$df2)

    .print.pwrss.ancova(print.obj, verbose = verbose, utf = utf)

  } #verbose

  invisible(structure(list(parms = func.parms,
                           test = "F",
                           effect = effect,
                           eta.squared = eta.squared,
                           f = sqrt(pwr.obj$f.squared),
                           df1 = pwr.obj$df1,
                           df2 = pwr.obj$df2,
                           ncp = pwr.obj$lambda,
                           null.ncp = 0,
                           f.alpha = pwr.obj$f.alpha,
                           power = pwr.obj$power,
                           n.vector = n.vector,
                           n.total = n.total),
                      class = c("pwrss", "f", "ancova", "keppel")))

} # power.f.ancova.keppel


# default base functions stats::contr.treatment() and stats::contr.sum()
# provides coding for the design matrix (dummy, effect, etc.)
# some further manipulations are needed to get contrasts
# that is, add intercept and take inverse
# https://rpubs.com/timflutre/tuto_contrasts
# https://stats.oarc.ucla.edu/spss/faq/coding-systems-for-categorical-variables-in-regression-analysis/

#' Factorial Contrasts
#'
#' @description
#' Helper function to construct the default contrast coefficients for various
#' coding schemes.
#'
#' Validated using \code{lm()} and \code{aov()} functions.
#'
#' @details
#' Note that R has a partial matching feature which allows you to specify
#' shortened versions of arguments, such as \code{factor} instead of
#' \code{factor.levels}, or such as \code{cod} or \code{coding} instead of
#' \code{coding.scheme}.
#'
#'
#' @param factor.levels integer; Number of levels or groups in each factor. For
#'                      example, for two factors each having two levels or
#'                      groups use e.g. c(2, 2), for three factors each having
#'                      two levels or groups use e.g. c(2, 2, 2).
#' @param coding.scheme character vector; Coding scheme for each factor. "sum"
#'                      for deviation or effect coding, "treatment" for dummy
#'                      coding, "helmert" for Helmert type of coding, and
#'                      "poly" for polynomial coding. Each factor can have
#'                      their own coding scheme. If a single character value is
#'                      provided, it will be copied to other factors.
#' @param base          integer vector; Specifies which group is considered the
#'                      baseline group. Ignored for coding schemes other than
#'                      "treatment".
#' @param intercept     logical; \code{FALSE} by default. If \code{TRUE}
#'                      contrast matrix includes the intercept.
#' @param verbose       \code{1} by default. If \code{0} no output is printed
#'                      on the console.
#'
#' @return
#'   \item{factor.levels}{Number of levels (or groups) in each factor}
#'   \item{factor.data}{Unique combination of factor levels}
#'   \item{model.matrix}{Model (design) matrix based on unique combination of
#'                       factor levels}
#'   \item{contrast.matrix}{Contrast matrix}
#'
#' @examples
#' ###################################################
#' ############### dummy coding scheme ###############
#' ####################################################
#'
#' # one factor w/ 3 levels
#' contrast.object <- factorial.contrasts(factor.levels = 3,
#'                                        coding = "treatment")
#' # model matrix
#' contrast.object$model.matrix
#'
#' # contrast matrix
#' contrast.object$contrast.matrix
#'
#' #######################################################
#' ###############  deviation coding scheme ##############
#' #######################################################
#'
#' # especially useful for factorial designs
#' # two factors w/ 2 and 3 levels, respectively
#' contrast.object <- factorial.contrasts(factor.levels = c(2, 3),
#'                                        coding = "sum")
#'
#' # model matrix
#' contrast.object$model.matrix
#'
#' # contrast matrix
#' contrast.object$contrast.matrix
#'
#'
#' ######################################################
#' ###############  Helmert coding scheme ###############
#' ######################################################
#'
#' # one factor w/ 3 levels
#' contrast.object <- factorial.contrasts(factor.levels = 3,
#'                                        coding = "helmert")
#'
#' # model matrix
#' contrast.object$model.matrix
#'
#' # contrast matrix
#' contrast.object$contrast.matrix
#'
#' #########################################################
#' ###############  polynomial coding scheme ###############
#' #########################################################
#'
#' # one factor w/ 3 levels
#' contrast.object <- factorial.contrasts(factor.levels = 3,
#'                                        coding = "poly")
#'
#' # model matrix
#' contrast.object$model.matrix
#'
#' # contrast matrix
#' contrast.object$contrast.matrix
#'
#' @export factorial.contrasts
factorial.contrasts <- function(factor.levels = c(3, 2),
                                coding.scheme = rep("deviation", length(factor.levels)),
                                base = factor.levels, # only used with dummy or treatment coding
                                intercept = FALSE,
                                verbose = 1) {

  verbose <- ensure.verbose(verbose)

  if (length(coding.scheme) > length(factor.levels)) {

    coding.scheme <- coding.scheme[seq_along(factor.levels)]
    if (verbose >= 0)
      message(sprintf("Provide as many coding schemes as number of factors. Using the first %d.", length(factor.levels)))

  } # coding.scheme

  if (length(factor.levels) > 1) {

    if (length(coding.scheme) == 1 && verbose >= 0) {

      coding.scheme <- rep(coding.scheme, length(factor.levels))
      if (verbose >= 0)
        message("Assuming the same coding scheme applies to the other factor(s)")

    }

    if (length(base) == 1) {

      base <- rep(base, length(factor.levels))

    } else if (length(base) > length(factor.levels)) {

      base <- base[seq_along(factor.levels)]

    } # base

  } # factor.levels


  temp.contrast.list <- rep(list(NULL), length(factor.levels))

  for (i in seq_along(factor.levels)) {

    if (coding.scheme[i] %in% c("dummy", "treatment")) {

      temp.contrast.list[[i]] <- stats::contr.treatment(n = factor.levels[i], base = base[i])

    } else if (coding.scheme[i] %in% c("deviation", "effect", "sum")) {

      temp.contrast.list[[i]] <- stats::contr.sum(n = factor.levels[i])

    } else if (coding.scheme[i] == "helmert") {

      temp.contrast.list[[i]] <- stats::contr.helmert(n = factor.levels[i])

    } else if (coding.scheme[i] %in% c("poly", "polynomial")) {

      temp.contrast.list[[i]] <- stats::contr.poly(n = factor.levels[i])

    } else {

      stop("Contrast type \"", coding.scheme[i], "\" not supported at the moment.", call. = FALSE)

    }

  } # for loop



  if (length(factor.levels) == 1) {

    factor.data <- data.frame(A = gl(factor.levels[1], 1))

    contrasts.list <- list(A = temp.contrast.list[[1]])

    model.mat <- stats::model.matrix(~ A, factor.data, contrasts.arg = contrasts.list)

    contrast.mat <- solve(model.mat)

    col.names <- expand.grid(
      A = paste0("A", seq_len(factor.levels[1]))
    )

    means.cell.names <- apply(col.names, 1, paste, collapse = ":")
    colnames(contrast.mat) <- means.cell.names


  } else if (length(factor.levels) == 2) {

    factor.data <- data.frame(A = gl(factor.levels[1], factor.levels[2]),
                              B = gl(factor.levels[2], 1, prod(factor.levels)))

    contrasts.list <- list(A = temp.contrast.list[[1]],
                           B = temp.contrast.list[[2]])

    model.mat <- stats::model.matrix(~ A + B + A:B, factor.data,
                                 contrasts.arg = contrasts.list)

    contrast.mat <- solve(model.mat)

    col.names <- expand.grid(
      B = paste0("B", seq_len(factor.levels[2])),
      A = paste0("A", seq_len(factor.levels[1]))
    )
    col.names <- col.names[, c("A", "B")]
    means.cell.names <- apply(col.names, 1, paste, collapse = ":")

    colnames(contrast.mat) <- means.cell.names


  } else if (length(factor.levels) == 3) {

    factor.data <- data.frame(A = gl(factor.levels[1], prod(factor.levels[2:3])),
                              B = gl(factor.levels[2], factor.levels[3], prod(factor.levels)),
                              C = gl(factor.levels[3], 1, prod(factor.levels)))


    contrasts.list <- list(A = temp.contrast.list[[1]],
                           B = temp.contrast.list[[2]],
                           C = temp.contrast.list[[3]])

    model.mat <- stats::model.matrix(~ A + B + C + A:B + A:C + B:C + A:B:C, factor.data,
                                 contrasts.arg = contrasts.list)

    contrast.mat <- solve(model.mat)

    col.names <- expand.grid(
      C = paste0("C", seq_len(factor.levels[3])),
      B = paste0("B", seq_len(factor.levels[2])),
      A = paste0("A", seq_len(factor.levels[1]))
    )
    col.names <- col.names[, c("A", "B", "C")]
    means.cell.names <- apply(col.names, 1, paste, collapse = ":")
    colnames(contrast.mat) <- means.cell.names

  } else {

    stop("This version supports only up to three-way interactions.", call. = FALSE)

  }

  if (length(factor.levels) > 1 && verbose >= 0)
    message("Elements of `mu.vector`, `sd.vector`, `n.vector` or `p.vector` should follow this specific order:\n",
            paste(means.cell.names, " "), "\n")

  if (isFALSE(intercept)) contrast.mat <- contrast.mat[-1, ]

  if (verbose > 0) {
    print(round(contrast.mat, 3))
  }

  invisible(list(factor.levels = factor.levels,
                 factor.data = factor.data,
                 model.matrix = model.mat,
                 contrast.matrix = contrast.mat))

} # factorial.contrasts()


#' Power Analysis for One-, Two-, Three-Way ANCOVA Using Means, Standard
#' Deviations, and (Optionally) Contrasts (F test)
#'
#' @description
#' Calculates power, sample size or effect size for one-, two-, three-way
#' ANCOVA. For factorial designs, use the argument \code{factor.levels} but
#' note that unique combination of levels (cells in this case) should follow a
#' specific order for the test of interaction. The order of marginal means and
#' standard deviations is printed as a warning message.
#'
#' Formulas are validated using examples and tables in Shieh (2020).
#'
#' @details
#' Note that R has a partial matching feature which allows you to specify
#' shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead
#' of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of
#' \code{k.covariates}.
#'
#'
#' @param mu.vector       vector; adjusted means (or estimated marginal means)
#'                        for each level of a factor.
#' @param sd.vector       vector; unadjusted standard deviations for each level
#'                        of a factor. If a pooled standard deviation is
#'                        provided, repeat its value to match the number of
#'                        group means. A warning will be issued if group
#'                        standard deviations differ substantially beyond what
#'                        is expected due to sampling error.
#' @param n.vector        vector; sample sizes for each level of a factor.
#' @param p.vector        vector; proportion of total sample size in each level
#'                        of a factor. These proportions should sum to one.
#' @param factor.levels   integer; number of levels or groups in each factor.
#'                        For example, for two factors each having two levels
#'                        or groups use e.g. c(2, 2), for three factors each
#'                        having two levels or groups use e.g. c(2, 2, 2).
#' @param r.squared       explanatory power of covariates (R-squared) in the
#'                        ANCOVA model.
#' @param k.covariates    integer; number of covariates in the ANCOVA model.
#' @param contrast.matrix vector or matrix; contrasts should not be confused
#'                        with the model (design) matrix. Rows of contrast
#'                        matrix indicate independent vector of contrasts
#'                        summing to zero. The default contrast matrix is
#'                        constructed using deviation coding scheme (a.k.a.
#'                        effect coding). Columns in the contrast matrix
#'                        indicate number of levels or groups (or cells in
#'                        factorial designs).
#' @param power           statistical power, defined as the probability of
#'                        correctly rejecting a false null hypothesis, denoted
#'                        as \eqn{1 - \beta}.
#' @param alpha           type 1 error rate, defined as the probability of
#'                        incorrectly rejecting a true null hypothesis, denoted
#'                        as \eqn{\alpha}.
#' @param ceil.n          logical; \code{TRUE} by default. If \code{FALSE}
#'                        sample sizes in each cell are NOT rounded up.
#' @param verbose         \code{1} by default (returns test, hypotheses, and
#'                        results), if \code{2} a more detailed output is given
#'                        (plus key parameters and definitions), if \code{0} no
#'                        output is printed on the console.
#' @param utf             logical; whether the output should show Unicode
#'                        characters (if encoding allows for it). \code{FALSE}
#'                        by default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test (F-Test)}
#'   \item{eta.squared}{(partial) eta-squared.}
#'   \item{f}{Cohen's f statistic.}
#'   \item{df1}{numerator degrees of freedom.}
#'   \item{df2}{denominator degrees of freedom.}
#'   \item{ncp}{non-centrality parameter for the alternative.}
#'   \item{null.ncp}{non-centrality parameter for the null.}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n.total}{total sample size.}
#'
#' @references
#'   Shieh, G. (2020). Power analysis and sample size planning in ANCOVA
#'   designs. *Psychometrika, 85*(1), 101-120.
#'   https://doi.org/10.1007/s11336-019-09692-3
#'
#' @examples
#'
#' ###################################################################
#' ##########################  main effect  ##########################
#' ###################################################################
#'
#' # power for one-way ANCOVA (two levels or groups)
#' power.f.ancova.shieh(mu.vector = c(0.20, 0), # marginal means
#'                      sd.vector = c(1, 1), # unadjusted standard deviations
#'                      n.vector = c(150, 150), # sample sizes
#'                      r.squared = 0.50, # proportion of variance explained by covariates
#'                      k.covariates = 1, # number of covariates
#'                      alpha = 0.05)
#'
#'
#' # sample size for one-way ANCOVA (two levels or groups)
#' power.f.ancova.shieh(mu.vector = c(0.20, 0), # marginal means
#'                      sd.vector = c(1, 1), # unadjusted standard deviations
#'                      p.vector = c(0.50, 0.50), # allocation, should sum to 1
#'                      r.squared = 0.50,
#'                      k.covariates = 1,
#'                      alpha = 0.05,
#'                      power = 0.80)
#'
#' ###################################################################
#' #######################  interaction effect  ######################
#' ###################################################################
#'
#' # sample size for two-way ANCOVA (2 x 2)
#' power.f.ancova.shieh(mu.vector = c(0.20, 0.25, 0.15, 0.05), # marginal means
#'                      sd.vector = c(1, 1, 1, 1), # unadjusted standard deviations
#'                      p.vector = c(0.25, 0.25, 0.25, 0.25), # allocation, should sum to 1
#'                      factor.levels = c(2, 2), # 2 by 2 factorial design
#'                      r.squared = 0.50,
#'                      k.covariates = 1,
#'                      alpha = 0.05,
#'                      power = 0.80)
#' # Elements of `mu.vector`, `sd.vector`, `n.vector` or `p.vector` should follow this specific order:
#' #  A1:B1  A1:B2  A2:B1  A2:B2
#'
#' ###################################################################
#' #######################  planned contrasts  #######################
#' ###################################################################
#'
#' #########################
#' ## dummy coding scheme ##
#' #########################
#'
#' contrast.object <- factorial.contrasts(factor.levels = 3, # one factor w/ 3 levels
#'                                        coding = "treatment") # use dummy coding scheme
#'
#' # get contrast matrix from the contrast object
#' contrast.matrix <- contrast.object$contrast.matrix
#'
#' # calculate sample size given design characteristics
#' ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
#'                                       sd.vector = c(1, 1, 1), # unadjusted standard deviations
#'                                       p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
#'                                       contrast.matrix = contrast.matrix,
#'                                       r.squared = 0.50,
#'                                       k.covariates = 1,
#'                                       alpha = 0.05,
#'                                       power = 0.80)
#'
#' # power of planned contrasts, adjusted for alpha level
#' power.t.contrasts(ancova.design, adjust.alpha = "fdr")
#'
#' ###########################
#' ## Helmert coding scheme ##
#' ###########################
#'
#' contrast.object <- factorial.contrasts(factor.levels = 3, # one factor w/ 4 levels
#'                                        coding = "helmert") # use helmert coding scheme
#'
#' # get contrast matrix from the contrast object
#' contrast.matrix <- contrast.object$contrast.matrix
#'
#' # calculate sample size given design characteristics
#' ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
#'                                       sd.vector = c(1, 1, 1), # unadjusted standard deviations
#'                                       p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
#'                                       contrast.matrix = contrast.matrix,
#'                                       r.squared = 0.50,
#'                                       k.covariates = 1,
#'                                       alpha = 0.05,
#'                                       power = 0.80)
#'
#' # power of planned contrasts
#' power.t.contrasts(ancova.design)
#'
#' ##############################
#' ## polynomial coding scheme ##
#' ##############################
#'
#' contrast.object <- factorial.contrasts(factor.levels = 3, # one factor w/ 4 levels
#'                                        coding = "poly") # use polynomial coding scheme
#'
#' # get contrast matrix from the contrast object
#' contrast.matrix <- contrast.object$contrast.matrix
#'
#' # calculate sample size given design characteristics
#' ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
#'                                       sd.vector = c(1, 1, 1), # unadjusted standard deviations
#'                                       p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
#'                                       contrast.matrix = contrast.matrix,
#'                                       r.squared = 0.50,
#'                                       k.covariates = 1,
#'                                       alpha = 0.05,
#'                                       power = 0.80)
#'
#' # power of the planned contrasts
#' power.t.contrasts(ancova.design)
#'
#' ######################
#' ## custom contrasts ##
#' ######################
#'
#' # custom contrasts
#' contrast.matrix <- rbind(
#'   cbind(A1 = 1, A2 = -0.50, A3 = -0.50),
#'   cbind(A1 = 0.50, A2 = 0.50, A3 = -1)
#' )
#' # labels are not required for custom contrasts,
#' # but they make it easier to understand power.t.contrasts() output
#'
#' # calculate sample size given design characteristics
#' ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
#'                                       sd.vector = c(1, 1, 1), # unadjusted standard deviations
#'                                       p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
#'                                       contrast.matrix = contrast.matrix,
#'                                       r.squared = 0.50,
#'                                       k.covariates = 1,
#'                                       alpha = 0.05,
#'                                       power = 0.80)
#'
#' # power of the planned contrasts
#' power.t.contrasts(ancova.design)
#'
#' @export power.f.ancova.shieh
power.f.ancova.shieh <- function(mu.vector,
                                 sd.vector,
                                 n.vector = NULL,
                                 p.vector = NULL,
                                 factor.levels = NULL,
                                 r.squared = 0,
                                 k.covariates = 1,
                                 contrast.matrix = NULL,
                                 power = NULL,
                                 alpha = 0.05,
                                 ceil.n = TRUE,
                                 verbose = 1,
                                 utf = FALSE) {

  func.parms <- as.list(environment())

  # value and consistency checks
  if (!is.vector(mu.vector) || !is.numeric(mu.vector))
    stop("Provide a vector of means (`mu.vector`) with its length equal to the number of groups.", call. = FALSE)
  check.vector(mu.vector, check.numeric)
  check.vector(sd.vector, check.positive)
  if (!is.null(n.vector)) check.vector(n.vector, check.sample.size)
  if (!is.null(p.vector)) check.vector(p.vector, check.proportion)
  check.same.lengths(mu.vector, sd.vector, n.vector, p.vector)
  if (is.null(factor.levels)) factor.levels <- length(mu.vector)
  check.vector(factor.levels, check.factor.level, min.length = 1)
  if (length(mu.vector) != prod(factor.levels))
    stop("Provide a vector of means (`mu.vector`) with its length equal to the the product of `factor.levels`.", call. = FALSE)
  if (r.squared > 1 || r.squared < 0 || !is.numeric(r.squared) || length(r.squared) != 1)
    stop("R-squared (explanatory power of covariates) takes a value between 0 and 1.", call. = FALSE)
  check.positive(k.covariates)
  if (!is.null(power)) check.power(power)
  if (alpha > 1 || alpha < 0 || !is.numeric(alpha) || length(alpha) != 1)
    stop("Type 1 error rate (alpha) takes a value between 0 and 1.", call. = FALSE)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = NA, n = n.vector, power = power) # calculation of effect size not possible

  n.way <- length(factor.levels)
  if (n.way > 3)
    stop("More than three-way ANCOVA is not allowed at the moment.", call. = FALSE)

  if (is.null(contrast.matrix)) {

    contrast.matrix <- factorial.contrasts(factor.levels = factor.levels,
                                           intercept = FALSE,
                                           verbose = min(verbose, 0))$contrast.matrix
    if (is.vector(contrast.matrix)) contrast.matrix <- t(as.matrix(contrast.matrix))
    contrast.matrix <- utils::tail(contrast.matrix, n = prod(factor.levels - 1))

  } else {

    if (is.vector(contrast.matrix))
      contrast.matrix <- t(as.matrix(contrast.matrix))
    if (!is.matrix(contrast.matrix))
      stop("Contrast coefficients are not provided in the form of a matrix.", call. = FALSE)
    if (dim(contrast.matrix)[2] != length(mu.vector))
      stop("The number of columns in the contrast matrix should match number of groups.", call. = FALSE)
    if (dim(contrast.matrix)[1] > length(mu.vector) - 1)
      stop("The number of rows in the contrast matrix should be less than or equal to number of groups minus one.", call. = FALSE)
  }

  pwr.shieh <- function(mu.vector, sd.vector, n.vector, k.covariates, r.squared, alpha, contrast.matrix, calculate.lambda = FALSE) {

    n.total <- sum(n.vector)
    sigma2_pooled <- sum((n.vector - 1) * sd.vector ^ 2) / (n.total - length(mu.vector))
    sigma2_error <- sigma2_pooled * (1 - r.squared)

    Q.mat <- diag(n.total / n.vector)
    mean.mat <- matrix(mu.vector, length(mu.vector), 1)
    gamma2 <- (t(contrast.matrix %*% mean.mat) %*% solve(contrast.matrix %*% Q.mat %*%
               t(contrast.matrix)) %*% (contrast.matrix %*% mean.mat)) / sigma2_error
    gamma2 <- as.numeric(gamma2)

    df1 <- nrow(contrast.matrix)
    df2 <- calc.df2(n.total, length(mu.vector), k.covariates)

    f.alpha <- stats::qf(1 - alpha, df1 = df1, df2 = df2)

    if (k.covariates > 1) {

      shape1 <- (df2 + 1) / 2
      shape2 <- k.covariates / 2
      mean.beta <- shape1 / (shape1 + shape2)
      sd.beta <- sqrt((shape1 * shape2) / ((shape1 + shape2) ^ 2 * (shape1 + shape2 + 1)))
      lower.beta <- max(0, mean.beta - 10 * sd.beta)

      integrand <- function(x) {
        stats::dbeta(x, shape1 = (df2 + 1) / 2, shape2 = k.covariates / 2) *
        stats::pf(f.alpha, df1, df2, n.total * gamma2 * x, lower.tail = FALSE)
      }
      power <- stats::integrate(integrand, lower = lower.beta, upper = 1)$value

      lambda <- ifelse(calculate.lambda, n.total * gamma2 * (df2 + 1) / (df2 + 1 + k.covariates), NA)

    } else if (k.covariates == 1) {

      integrand <- function(x) {
        stats::dt(x, (df2 + 1)) * stats::pf(f.alpha, df1, df2, n.total * gamma2 / (1 + (k.covariates / (df2 + 1)) * x ^ 2), lower.tail = FALSE)
      }
      power <- stats::integrate(integrand, lower = -10, upper = 10)$value

      if (calculate.lambda) {
        integrand <- function(x) {
          stats::dt(x, df2) * n.total * gamma2 / (1 + (k.covariates / df2) * x ^ 2)
        }
        lambda <- stats::integrate(integrand, lower = -Inf, upper = Inf)$value
      } else {
        lambda <- NA
      }

    }

    list(power = power, df1 = df1, df2 = df2, lambda = lambda, f.alpha = f.alpha)

  }


  if (requested == "n") {

    if (is.null(p.vector)) stop("`p.vector` can not be NULL when sample size is requested.", call. = FALSE)
    if (round(sum(p.vector), 5) != 1) stop("The elements of the `p.vector` should sum to 1.", call. = FALSE)

    n.total <- stats::uniroot(function(n.total) {
      power -  pwr.shieh(mu.vector = mu.vector, sd.vector = sd.vector, n.vector = n.total * p.vector,
                         k.covariates = k.covariates, r.squared =  r.squared, alpha = alpha,
                         contrast.matrix = contrast.matrix)$power
    }, interval = c(length(mu.vector) + k.covariates + 1, 1e10))$root
    n.vector <- n.total * p.vector

    if (ceil.n) n.vector <- ceiling(n.vector)

  }

  # calculate power (if requested == "power") or update it (if requested == "n")
  pwr.obj <- pwr.shieh(mu.vector = mu.vector, sd.vector = sd.vector,
                       n.vector = n.vector, k.covariates = k.covariates,
                       r.squared =  r.squared, alpha = alpha, contrast.matrix = contrast.matrix,
                       calculate.lambda = TRUE)

  n.total <- sum(n.vector)
  effect <- paste0(LETTERS[seq_along(factor.levels)], "(", factor.levels, ")", collapse = ":")
  f.squared <- pwr.obj$lambda / n.total
  eta.squared <- f.squared / (1 + f.squared)

  check_var.ratio(sd.vector, n.vector)

  if (verbose > 0) {

    print.obj <- list(test = fmt_test_ancova(n.way, k.covariates),
                      requested = requested,
                      effect = effect,
                      n.total = n.total,
                      factor.levels = factor.levels,
                      eta.squared = eta.squared,
                      power = pwr.obj$power,
                      ncp = pwr.obj$lambda,
                      null.ncp = 0,
                      alpha = alpha,
                      f.alpha = pwr.obj$f.alpha,
                      df1 = pwr.obj$df1,
                      df2 = pwr.obj$df2)

    .print.pwrss.ancova(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = "F",
                           effect = effect,
                           eta.squared = eta.squared,
                           f = sqrt(f.squared),
                           df1 = pwr.obj$df1,
                           df2 = pwr.obj$df2,
                           ncp = pwr.obj$lambda,
                           null.ncp = 0,
                           f.alpha = pwr.obj$f.alpha,
                           power = pwr.obj$power,
                           n.vector = n.vector,
                           n.total = n.total),
                      class = c("pwrss", "f", "ancova", "shieh")))


} # end of power.f.ancova.shieh()

#' Power Analysis for One-, Two-, Three-Way ANCOVA Contrasts and Multiple
#' Comparisons (T-Tests)
#'
#' @description
#' Calculates power or sample size for a single one-, two-, three-Way ANCOVA
#' contrast.
#'
#' Formulas are validated using examples and tables in Shieh (2017).
#'
#' @details
#' Note that R has a partial matching feature which allows you to specify
#' shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead
#' of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of
#' \code{k.covariates}.
#'
#'
#' @aliases pwrss.t.contrast
#'
#' @param mu.vector       vector; adjusted means (or estimated marginal means)
#'                        for each level of a factor.
#' @param sd.vector       vector; unadjusted standard deviations for each level
#'                        of a factor.
#' @param n.vector        vector; sample sizes for each level of a factor.
#' @param p.vector        vector; proportion of total sample size in each level
#'                        of a factor. These proportions should sum to one.
#' @param r.squared       explanatory power of covariates (R-squared) in the
#'                        ANCOVA model.
#' @param k.covariates    Number of covariates in the ANCOVA model.
#' @param contrast.vector vector; a single contrast in the form of a vector
#'                        with as many elements as number of levels or groups
#'                        (or cells in factorial designs). Ignored when 'x' is
#'                        specified.
#' @param power           statistical power, defined as the probability of
#'                        correctly rejecting a false null hypothesis, denoted
#'                        as \eqn{1 - \beta}.
#' @param alpha           type 1 error rate, defined as the probability of
#'                        incorrectly rejecting a true null hypothesis, denoted
#'                        as \eqn{\alpha}.
#' @param tukey.kramer    logical; \code{FALSE} by default. If \code{TRUE}
#'                        adjustments will be made to control Type 1 error.
#' @param ceil.n          logical; \code{TRUE} by default. If \code{FALSE}
#'                        sample sizes in each cell are NOT rounded up.
#' @param verbose         \code{1} by default (returns test, hypotheses, and
#'                        results), if \code{2} a more detailed output is given
#'                        (plus key parameters and definitions), if \code{0} no
#'                        output is printed on the console.
#' @param utf             logical; whether the output should show Unicode
#'                        characters (if encoding allows for it). \code{FALSE}
#'                        by default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test (T-Test).}
#'   \item{psi}{contrast-weighted mean difference.}
#'   \item{d}{contrast-weighted standardized mean difference.}
#'   \item{df}{degrees of freedom.}
#'   \item{t.alpha}{critical values.}
#'   \item{ncp}{non-centrality parameter for the alternative.}
#'   \item{null.ncp}{non-centrality parameter for the null.}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n.vector}{sample sizes for each level of a factor.}
#'   \item{n.total}{total sample size.}
#'
#' @references
#'   Shieh, G. (2017). Power and sample size calculations for contrast analysis
#'   in ANCOVA. *Multivariate Behavioral Research, 52*(1), 1-11.
#'   https://doi.org/10.1080/00273171.2016.1219841
#'
#' @examples
#' # dummy coding example (uses the first contrast from a three-level- / two-contrasts-design)
#' contrast.object <- factorial.contrasts(factor.levels = 3, coding = "treatment", verbose = 0)
#' contrast.vector <- contrast.object[["contrast.matrix"]][1, ]
#' power.t.contrast(mu.vector = c(0.15, 0.30, 0.20),
#'                  sd.vector = c(1,    1,    1),
#'                  p.vector  = c(1/3,  1/3,  1/3),
#'                  r.squared = 0.50, k.covariates = 1,
#'                  contrast.vector = contrast.vector,
#'                  power = 0.80, alpha = 0.05)
#'
#' @export power.t.contrast
power.t.contrast <- function(mu.vector, sd.vector,
                             n.vector = NULL, p.vector = NULL,
                             contrast.vector,
                             r.squared = 0, k.covariates = 1,
                             power = NULL,
                             alpha = 0.05, tukey.kramer = FALSE,
                             ceil.n = TRUE, verbose = 1, utf = FALSE) {

  func.parms <- as.list(environment())

  # value checks
  check.vector(mu.vector, check.numeric)
  check.vector(sd.vector, check.positive)
  if (!is.null(n.vector)) check.vector(n.vector, check.sample.size)
  if (!is.null(p.vector)) check.vector(p.vector, check.proportion)
  check.same.lengths(mu.vector, sd.vector, n.vector, p.vector)
  if (r.squared > 1 || r.squared < 0 || !is.numeric(r.squared) || length(r.squared) != 1)
    stop("R-squared (explanatory power of covariates) takes a value between 0 and 1.", call. = FALSE)
  check.nonnegative(k.covariates)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(tukey.kramer, ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = NA, n = n.vector, power = power) # calculation of effect size not possible

  if (is.matrix(contrast.vector) && dim(contrast.vector)[1] > 1)
    stop("The number of rows in the contrast matrix should be one.", call. = FALSE)
  if (!is.matrix(contrast.vector) && !is.vector(contrast.vector))
    stop("`contrast.vector` must be either a vector or a matrix.", call. = FALSE)

  if (is.matrix(contrast.vector))
    contrast.vector <- contrast.vector[1, ]

  if (length(contrast.vector) != length(mu.vector))
    stop("The number of columns / elements in the contrast vector should match number of groups.", call. = FALSE)

  if (requested == "n" && is.null(p.vector))
    stop("The `p.vector` can not be NULL when the sample size is requested.", call. = FALSE)
  if (requested == "n" && abs(sum(p.vector) - 1) > 1e-6)
    stop("The elements of the `p.vector` should sum to 1.", call. = FALSE)

  pwr.contrast <- function(mu.vector, sd.vector, n.vector, k.covariates,
                           r.squared, alpha, contrast.vector, tukey.kramer,
                           calculate.lambda = FALSE) {

    n.total <- sum(n.vector)
    sigma2_pooled <- sum((n.vector - 1) * sd.vector ^ 2) / (n.total - length(mu.vector))
    sigma2_error <- sigma2_pooled * (1 - r.squared)

    psi <- sum(contrast.vector * mu.vector)
    df  <- n.total - length(mu.vector) - k.covariates

    if (tukey.kramer == 1) {

      t.alpha <- stats::qtukey(1 - alpha, length(mu.vector), df) / sqrt(2)

      } else {

      t.alpha <- stats::qt(1 - alpha / 2, df)

    }

    a <- sum(contrast.vector ^ 2 / n.vector)
    # delta <- psi / sqrt(a * sigma2_error)
    d <- psi / sqrt(ifelse(k.covariates == 0, sigma2_pooled, sigma2_error))

    if (k.covariates == 0) {

      # Shieh (2023, p. 3/18)
      # Shieh G (2023) Assessing standardized
      # contrast effects in ANCOVA: Confidence intervals,
      # precision evaluations, and sample size
      # requirements. PLoS ONE 18(2): e0282161.
      # https://doi.org/10.1371/journal.pone.0282161

      ss_wts <- sum((contrast.vector ^ 2) / n.vector)
      se_psi <- sqrt(sigma2_error * ss_wts)

      power <- power.t.test(ncp = psi / se_psi, df = df, plot = FALSE, verbose = 0)$power

      if (calculate.lambda)
        lambda <- psi / sigma2_pooled
      else
        lambda <- NA

    } else if (k.covariates == 1) {

      ncp.x <- function(x) psi / sqrt((1 + x ^ 2 / (df + 1)) * sigma2_error * a)
      integrand <- function(x) {
        stats::dt(x, df + 1) * (stats::pt(-t.alpha, df, ncp.x(x)) + stats::pt(+t.alpha, df, ncp.x(x), lower.tail = FALSE))
      }

      power <- stats::integrate(integrand, -10, 10)$value
      if (calculate.lambda)
        lambda <- stats::integrate(function(x) stats::dt(x, df + 1) * ncp.x(x), -10, 10)$value
      else
        lambda <- NA

    } else if (k.covariates > 1) {

      shape1 <- (df + 1) / 2
      shape2 <- k.covariates / 2
      mean.beta <- shape1 / (shape1 + shape2)
      sd.beta <- sqrt((shape1 * shape2) / ((shape1 + shape2) ^ 2 * (shape1 + shape2 + 1)))
      lower.beta <- max(0, mean.beta - 10 * sd.beta)
      ncp.x <- function(x) sqrt(x) * psi / sqrt(sigma2_error * a)

      integrand <- function(x) {
        stats::dbeta(x, (df + 1) / 2, k.covariates / 2) *
          (stats::pt(-t.alpha, df, ncp.x(x)) + stats::pt(+t.alpha, df, ncp.x(x), lower.tail = FALSE))
      }
      power <- stats::integrate(integrand, lower = lower.beta, upper = 1)$value

      if (calculate.lambda) {
        lambda <- stats::integrate(function(x) stats::dbeta(x, (df + 1) / 2, k.covariates / 2) * ncp.x(x), lower.beta, 1)$value
      } else {
        lambda <- NA
      }

    }

    list(power = power, df = df, lambda = lambda, t.alpha = t.alpha, psi = psi, d = d)

  } # pwr.contrast()

  if (requested == "n") {

    if (sum(contrast.vector * mu.vector) == 0) {

      n.total <- as.integer(.Machine$integer.max)
      warning("Using infinity (maximum integer number as defined in R) for `n.total` because `psi` = 0.", call. = FALSE)

    } else {

      n.total <- try(stats::uniroot(function(n.total) {
        power - pwr.contrast(mu.vector = mu.vector, sd.vector = sd.vector,
                             n.vector = n.total * p.vector, k.covariates = k.covariates,
                             r.squared = r.squared, alpha = alpha, contrast.vector = contrast.vector,
                             tukey.kramer = tukey.kramer, calculate.lambda = FALSE)$power
      }, interval = c(length(mu.vector) + k.covariates + 1, 1e10))$root, silent = TRUE)
      if (inherits(n.total, "try-error") || n.total == 1e10) stop("Design not feasible.", call. = FALSE)

    }

    n.vector <- n.total * p.vector
    if (ceil.n) n.vector <- ceiling(n.vector)

  }

  # calculate power (if requested == "power") or update it (if requested == "n")
  pwr.obj <- pwr.contrast(mu.vector = mu.vector, sd.vector = sd.vector, n.vector = n.vector,
                          k.covariates = k.covariates, r.squared = r.squared, alpha = alpha,
                          contrast.vector = contrast.vector, tukey.kramer = tukey.kramer, calculate.lambda = TRUE)

  n.total <- sum(n.vector)

  if (verbose > 0) {

    print.obj <- list(test = "Single Contrast Analysis (T-Test)",
                      requested = requested,
                      psi = pwr.obj$psi,
                      d = pwr.obj$d,
                      n.total = n.total,
                      power = pwr.obj$power,
                      ncp = pwr.obj$lambda,
                      null.ncp = 0,
                      alpha = alpha,
                      t.alpha = pwr.obj$t.alpha * c(-1, 1),
                      df = pwr.obj$df)

    .print.pwrss.contrast(print.obj, verbose = verbose, utf = utf)

  }

  invisible(structure(list(parms = func.parms,
                           test = "t",
                           psi = pwr.obj$psi,
                           d = pwr.obj$d,
                           df = pwr.obj$df,
                           t.alpha = pwr.obj$t.alpha * c(-1, 1),
                           ncp = pwr.obj$lambda,
                           null.ncp = 0,
                           power = pwr.obj$power,
                           n.vector = n.vector,
                           n.total = n.total),
                      class = c("pwrss", "t", "contrast")))


} # end of power.t.contrast()


adjust.alpha <- function(n, alpha = 0.05,
                         method = c("bonferroni", "holm", "hochberg",
                                    "hommel", "BH", "BY", "fdr", "none")) {

  check.proportion(alpha)
  check.sample.size(n)
  method <- match.arg(method)

  p.adj <- stats::uniroot(function(p) {
    alpha - stats::p.adjust(p = p, method = method, n = n)
  }, interval = c(0, 1), tol = 1e-12)$root

  p.adj

}


#' Power Analysis for One-, Two-, Three-Way ANCOVA Contrasts and Multiple
#' Comparisons (T-Tests)
#'
#' @description
#' Calculates power or sample size for one-, two-, three-Way ANCOVA contrasts
#' and multiple comparisons. The \code{power.t.contrasts()} function permits
#' to test multiple contrasts (multiple comparisons) and also allows adjustment
#' to alpha due to multiple testing. Furthermore, \code{power.t.contrasts()}
#' accepts an object returned from the \code{power.f.ancova.shieh()} function
#' for convenience. Beware that, in this case, all other arguments are ignored
#' except \code{alpha} and \code{adjust.alpha}.
#'
#' Formulas are validated using examples and tables in Shieh (2017).
#'
#' @details
#' Note that R has a partial matching feature which allows you to specify
#' shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead
#' of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of
#' \code{k.covariates}.
#'
#'
#' @aliases pwrss.t.contrasts
#'
#' @param x               object; an object returned from the
#'                        \code{power.f.ancova.shieh()} function.
#' @param mu.vector       vector; adjusted means (or estimated marginal means)
#'                        for each level of a factor. Ignored when 'x' is
#'                        specified.
#' @param sd.vector       vector; unadjusted standard deviations for each level
#'                        of a factor. Ignored when 'x' is specified.
#' @param n.vector        vector; sample sizes for each level of a factor.
#'                        Ignored when 'x' is specified.
#' @param p.vector        vector; proportion of total sample size in each level
#'                        of a factor. These proportions should sum to one.
#'                        Ignored when 'x' is specified.
#' @param r.squared       explanatory power of covariates (R-squared) in the
#'                        ANCOVA model. Ignored when 'x' is specified.
#' @param k.covariates    Number of covariates in the ANCOVA model. Ignored
#'                        when 'x' is specified.
#' @param contrast.matrix vector or matrix; contrasts should not be confused
#'                        with the model (design) matrix. Rows of contrast
#'                        matrix indicate independent vector of contrasts
#'                        summing to zero. The default contrast matrix is
#'                        constructed using deviation coding scheme (a.k.a.
#'                        effect coding). Columns in the contrast matrix
#'                        indicate number of levels or groups (or cells in
#'                        factorial designs). Ignored when 'x' is specified.
#' @param power           statistical power, defined as the probability of
#'                        correctly rejecting a false null hypothesis, denoted
#'                        as \eqn{1 - \beta}. Ignored when 'x' is specified.
#' @param alpha           type 1 error rate, defined as the probability of
#'                        incorrectly rejecting a true null hypothesis, denoted
#'                        as \eqn{\alpha}. Note that this should be specified
#'                        even if 'x' is specified. The 'alpha' value within
#'                        the 'x' object pertains to the omnibus test, NOT the
#'                        test of contrasts.
#' @param adjust.alpha    character; one of the methods in c("none", "tukey",
#'                        "bonferroni", "holm", "hochberg", "hommel", "BH",
#'                        "BY", "fdr") to control Type 1 error. See
#'                        \code{?stats::p.adjust}.
#' @param ceil.n          logical; \code{TRUE} by default. If \code{FALSE}
#'                        sample sizes in each cell are NOT rounded up.
#' @param verbose         \code{1} by default (returns test, hypotheses, and
#'                        results), if \code{2} a more detailed output is given
#'                        (plus key parameters and definitions), if \code{0} no
#'                        output is printed on the console.
#' @param utf             logical; whether the output should show Unicode
#'                        characters (if encoding allows for it). \code{FALSE}
#'                        by default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test (T-Test).}
#'   \item{contrast}{contrast number (one contrast per line).}
#'   \item{comparison}{which factor levels are compared (one contrast per line).}
#'   \item{psi}{contrast-weighted mean difference (one contrast per line).}
#'   \item{d}{contrast-weighted standardized mean difference (one contrast per line).}
#'   \item{ncp}{non-centrality parameter for the alternative (one contrast per line).}
#'   \item{df}{degrees of freedom (one contrast per line).}
#'   \item{t.alpha}{critical values (one contrast per line).}
#'   \item{n.total}{total sample size (one contrast per line).}
#'   \item{power}{statistical power \eqn{(1-\beta)} (one contrast per line).}
#'
#' @references
#'   Shieh, G. (2017). Power and sample size calculations for contrast analysis
#'   in ANCOVA. *Multivariate Behavioral Research, 52*(1), 1-11.
#'   https://doi.org/10.1080/00273171.2016.1219841
#'
#' @examples
#' # see `?pwrss::power.f.ancova.shieh` for further examples
#'
#' # dummy coding example
#' contrast.object <- factorial.contrasts(factor.levels = 3, coding = "treatment", verbose = 0)
#' contrast.matrix <- contrast.object[["contrast.matrix"]]
#' power.t.contrasts(mu.vector = c(0.15, 0.30, 0.20),
#'                   sd.vector = c(1,    1,    1),
#'                   p.vector  = c(1/3,  1/3,  1/3),
#'                   r.squared = 0.50, k.covariates = 1,
#'                   contrast.matrix = contrast.matrix,
#'                   power = 0.80,
#'                   alpha = 0.05, adjust.alpha = "holm")
#'
#' @export power.t.contrasts
power.t.contrasts <- function(x = NULL,
                              mu.vector = NULL, sd.vector = NULL,
                              n.vector = NULL, p.vector = NULL,
                              r.squared = 0, k.covariates = 1,
                              contrast.matrix = NULL,
                              power = NULL,
                              alpha = 0.05,
                              adjust.alpha = c("none", "tukey", "bonferroni",
                                               "holm", "hochberg", "hommel",
                                               "BH", "BY", "fdr"),
                              ceil.n = TRUE, verbose = 1, utf = FALSE) {

  if (!is.null(x)) {

    if (inherits(x, c("pwrss", "f", "ancova", "shieh"))) {

      # transfer pwrss.f.ancova.shieh object into the input parameters and remove the object
      mu.vector <- x$parms$mu.vector
      sd.vector <- x$parms$sd.vector
      n.vector <- x$n.vector
      p.vector <- x$parms$p.vector
      r.squared <- x$parms$r.squared
      k.covariates <- x$parms$k.covariates
      contrast.matrix <- x$parms$contrast.matrix
      ceil.n <- x$parms$ceil.n
      power <- NULL

    } else {

      stop("This function only works with an object of type `pwrss`, `f`, `ancova`, and `shieh`.", call. = FALSE)

    }

  } else {

    # value checks
    check.vector(mu.vector, check.numeric)
    check.vector(sd.vector, check.positive)
    if (!is.null(n.vector)) check.vector(n.vector, check.sample.size)
    if (!is.null(p.vector)) check.vector(p.vector, check.proportion)
    check.proportion(r.squared)
    check.nonnegative(k.covariates)
    if (!is.null(power)) check.power(power)
    check.proportion(alpha)
    check.logical(ceil.n)

  } # if data is null

  rm(x)
  adjust.alpha <- match.arg(adjust.alpha)
  func.parms <- as.list(environment())
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = NA, n = n.vector, power = power) # calculation of effect size not possible

  if (is.vector(contrast.matrix))
    contrast.matrix <- t(as.matrix(contrast.matrix))
  levels <- colnames(contrast.matrix)

  if (adjust.alpha == "tukey") {

    tukey.kramer <- TRUE

  } else {

    tukey.kramer <- FALSE

    if (nrow(contrast.matrix) > 1)
      alpha <- adjust.alpha(alpha = alpha, n = nrow(contrast.matrix),  method = adjust.alpha)

  } # adjust.alpha

  comparison <- NULL
  power.out <- NULL
  for (i in seq_len(nrow(contrast.matrix))) {

    contrast.vector <- contrast.matrix[i, ]

    # on some machines, contr.poly() is not working as expected (it may contain very small values in places that should
    # be 0) - thus, the code below checks at which position the contrast vector contains positive or negative value
    # above / below neg.eps / -.neg.eps (OBS: neg.obs is slightly smaller than .eps and thus used in this comparison)
    comparison <- c(comparison, sprintf("%s <=> %s",
                                        paste(levels[which(contrast.vector >=  .Machine$double.neg.eps)], collapse = " "),
                                        paste(levels[which(contrast.vector <= -.Machine$double.neg.eps)], collapse = " ")))

    pwr.t.contr.obj <- power.t.contrast(mu.vector = mu.vector,
                                        sd.vector = sd.vector,
                                        n.vector = n.vector,
                                        p.vector = p.vector,
                                        contrast.vector = contrast.vector,
                                        r.squared = r.squared,
                                        k.covariates = k.covariates,
                                        alpha = alpha, power = power,
                                        tukey.kramer = tukey.kramer,
                                        ceil.n = ceil.n,
                                        verbose = 0)

    power.out.i <- cbind(psi = pwr.t.contr.obj$psi,
                         d = pwr.t.contr.obj$d,
                         ncp = pwr.t.contr.obj$ncp,
                         df = pwr.t.contr.obj$df,
                         t.alpha = max(pwr.t.contr.obj$t.alpha),
                         n.total = pwr.t.contr.obj$n.total,
                         power = pwr.t.contr.obj$power)

    power.out <- rbind(power.out, power.out.i)

  }

  rownames(comparison) <- NULL
  contrast.number <- seq_len(nrow(contrast.matrix))

  if (verbose > 0) {

    power.out <- as.data.frame(power.out)
    print.data <- data.frame(contr = contrast.number, comparison = comparison,
                             psi = round(power.out$psi, 3),
                             d = round(power.out$d, 3),
                             ncp = round(power.out$ncp, 3),
                             n.total = round(power.out$n.total, 3),
                             power = round(power.out$power, 3))

    print.obj <- list(test = "Multiple Contrast Analyses (T-Tests)",
                      requested = requested,
                      alpha = alpha,
                      adjust.alpha = adjust.alpha,
                      null.ncp = 0,
                      data = print.data)

    .print.pwrss.contrasts(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(c(list(parms = func.parms,
                             test = "t"),
                        data.frame(contrast = contrast.number,
                                   comparison = comparison,
                                   power.out)),
                      class = c("pwrss", "t", "contrasts")))

} # power.t.contrasts
