#' Power Analysis for Mixed-Effects Analysis of Variance (F-Test)
#'
#' @description
#' Calculates power, sample size or effect size for mixed-effects ANOVA design
#' with two factors (between and within). When there is only one group observed
#' over time, this design is often referred to as repeated-measures ANOVA.
#'
#' Formulas are validated using G*Power and tables in the PASS documentation.
#'
#' @details
#' * NB: The \code{pwrss.f.rmanova()} function is deprecated and will no longer
#'   be supported, but it will remain available as a wrapper for the
#'   \code{power.f.mixed.anova()} function during a transition period.
#'
#' @aliases power.f.mixed.anova pwrss.f.rmanova
#'
#'
#' @param eta.squared      (partial) eta-squared for the alternative.
#' @param null.eta.squared (partial) eta-squared for the null.
#' @param rho.within       Correlation between repeated measures. For example,
#'                         for pretest/post-test designs, this is the
#'                         correlation between pretest and post-test scores
#'                         regardless of group membership. The default is 0.50.
#'                         If \code{eta.squared} is already adjusted for this
#'                         correlation specify \code{rho.within = NA}.
#' @param factor.levels    vector; integer; length of two representing the
#'                         number of levels for groups and measures. For
#'                         example, in randomized controlled trials with two
#'                         arms (treatment and control) where pre-test,
#'                         post-test, and follow-up test are administered,
#'                         this would be represented as c(2, 3).
#' @param factor.type      vector; character; length of two indicating the
#'                         order of between-subject and within-subject
#'                         factors. By default, the first value represents the
#'                         between-subject factor and the second value
#'                         represents the within-subject factor. This argument
#'                         is rarely needed, except when unsure which element
#'                         in 'factor.levels' represent between-subject or
#'                         within-subject factors. Therefore, specify the
#'                         'factor.levels' accordingly.
#' @param epsilon          non-sphericity correction factor, default is 1
#'                         (means no violation of sphericity). Lower bound for
#'                         this argument is
#'                         \code{epsilon = 1 / (factor.levels[2] - 1)}.
#' @param n.total          integer; total sample size.
#' @param power            statistical power, defined as the probability of
#'                         correctly rejecting a false null hypothesis, denoted
#'                         as \eqn{1 - \beta}.
#' @param alpha            type 1 error rate, defined as the probability of
#'                         incorrectly rejecting a true null hypothesis,
#'                         denoted as \eqn{\alpha}.
#' @param effect           character; the effect of interest: "between",
#'                         "within", or "interaction".
#' @param ceil.n           logical; \code{TRUE} by default. If \code{FALSE}
#'                         sample size in each group is NOT rounded up.
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
#'   \item{ncp}{non-centrality parameter under alternative.}
#'   \item{null.ncp}{non-centrality parameter under null.}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n.total}{total sample size.}
#'
#' @references
#'   Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc
#'   analizi \[Statistical power analysis with pwrss R package\]. *Ahi Evran
#'   Universitesi Kirsehir Egitim Fakultesi Dergisi, 24*(3), 2207-2328.
#'   https://doi.org/10.29299/kefad.1209913
#'
#' @examples
#' ######################################################
#' # pretest-post-test design with treatment group only  #
#' ######################################################
#'
#' # a researcher is expecting a difference of Cohen's d = 0.30
#' # between post-test and pretest score translating into
#' # Eta-squared = 0.022
#'
#' # adjust effect size for correlation with 'rho.within'
#' power.f.mixed.anova(eta.squared = 0.022,
#'                     factor.levels = c(1, 2), # 1 between 2 within
#'                     rho.within = 0.50,
#'                     effect = "within",
#'                     power = 0.80, alpha = 0.05)
#'
#' # if effect size is already adjusted for correlation
#' # use 'rho.within = NA'
#' power.f.mixed.anova(eta.squared = 0.08255,
#'                     factor.levels = c(1, 2), # 1 between 2 within
#'                     rho.within = NA,
#'                     effect = "within",
#'                     power = 0.80, alpha = 0.05)
#'
#' ##########################################################
#' # post-test only design with treatment and control groups #
#' ##########################################################
#'
#' # a researcher is expecting a difference of Cohen's d = 0.50
#' # on the post-test score between treatment and control groups
#' # translating into Eta-squared = 0.059
#' power.f.mixed.anova(eta.squared = 0.059,
#'                     factor.levels = c(2, 1),  # 2 between 1 within
#'                     effect = "between",
#'                     power = 0.80, alpha = 0.05)
#'
#'
#' #############################################################
#' # pretest-post-test design with treatment and control groups #
#' #############################################################
#'
#' # a researcher is expecting a difference of Cohen's d = 0.40
#' # on the post-test score between treatment and control groups
#' # after controlling for the pretest translating into
#' # partial Eta-squared = 0.038
#' power.f.mixed.anova(eta.squared = 0.038,
#'                     factor.levels = c(2, 2),  # 2 between 2 within
#'                     rho.within = 0.50,
#'                     effect = "between",
#'                     power = 0.80, alpha = 0.05)
#'
#' # a researcher is expecting an interaction effect
#' # (between groups and time) of Eta-squared = 0.01
#' power.f.mixed.anova(eta.squared = 0.01,
#'                     factor.levels = c(2, 2),  # 2 between 2 within
#'                     rho.within = 0.50,
#'                     effect = "interaction",
#'                     power = 0.80, alpha = 0.05)
#'
#' # a researcher is expecting an interaction effect
#' # (between groups and time) of Eta-squared = 0.01
#' power.f.mixed.anova(eta.squared = 0.01,
#'                     factor.levels = c(2, 2),  # 2 between 2 within
#'                     rho.within = 0.50,
#'                     effect = "within",
#'                     power = 0.80, alpha = 0.05)
#'
#' @export power.f.mixed.anova
power.f.mixed.anova <- function(eta.squared = NULL,
                                null.eta.squared = 0,
                                factor.levels = c(2, 2),
                                factor.type = c("between", "within"),
                                rho.within = 0.50,
                                epsilon = 1,
                                n.total = NULL,
                                power = NULL, alpha = 0.05,
                                effect = c("between", "within", "interaction"),
                                ceil.n = TRUE, verbose = 1, utf = FALSE) {

  effect <- match.arg(effect)
  func.parms <- as.list(environment())

  if (!is.null(eta.squared)) check.nonnegative(eta.squared)
  check.nonnegative(null.eta.squared)
  if (!is.null(n.total)) check.sample.size(n.total)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = eta.squared, n = n.total, power = power)

  if (!all(c("between", "within") %in% factor.type))
    stop("The `factor.type` argument must be specified as either c('between', 'within') or c('within', 'between'), ",
         "indicating the order in which the corresponding values in `factor.levels` are interpreted - specifically, ",
         "which factor is treated as between-subjects and which as within-subjects.", call. = FALSE)
  if (length(factor.levels) != 2 || length(factor.type) != 2)
    stop("Exactly two factors are allowed in this procedure.", call. = FALSE)

  n.levels.between <- factor.levels[which(tolower(factor.type) == "between")]
  n.levels.within <- factor.levels[which(tolower(factor.type) == "within")]
  if (n.levels.within > 1 && epsilon <  1 / (n.levels.within - 1))
    stop("Incorrect value for the non-sphericity correction factor (`epsilon`).", call. = FALSE)

  if (is.na(rho.within)) {
    if (requested %in% c("power", "n")) {
      warning("Assuming that `eta.squared` and `null.eta.squared` are already adjusted for within-subject correlation.", call. = FALSE)
      correct4rho <- 1
    } else if (requested == "es") {
      stop("When the effect size shall be estimated, `rho.within` needs to be defined.", call. = FALSE)
    }
  } else {
    correct4rho <- (n.levels.within / (1 + ifelse(effect == "between", n.levels.within - 1, -1) * rho.within))
  }

  f.squared      <- eta.squared      / (1 - eta.squared) *      correct4rho
  null.f.squared <- null.eta.squared / (1 - null.eta.squared) * correct4rho

  pwr.mixed <- function(f.squared, null.f.squared, n.total,
                        n.levels.between, n.levels.within, epsilon, alpha, effect) {

    df1 <- prod(c(n.levels.between - 1, (n.levels.within - 1) * epsilon)[c(effect != "within", effect != "between")])
    df2 <- (n.total - n.levels.between) * ifelse(effect == "between", 1, (n.levels.within - 1) * epsilon)

    if (df1 < 1 || df2 < 1) stop("Design is not feasible", call. = FALSE)
    lambda <- f.squared * n.total * epsilon
    null.lambda <- null.f.squared * n.total * epsilon
    f.alpha <- stats::qf(alpha, df1 = df1, df2 = df2, ncp = null.lambda, lower.tail = FALSE)
    power <- stats::pf(f.alpha, df1 = df1, df2 = df2, ncp = lambda, lower.tail = FALSE)

    list(power = power, lambda = lambda, null.lambda = null.lambda,
         df1 = df1, df2 = df2, f.alpha = f.alpha)

  } # pwr.mixed()

  min.pwr <- function(f.squared, n.total, power) {
    power - pwr.mixed(f.squared, null.f.squared, n.total, n.levels.between, n.levels.within, epsilon, alpha, effect)$power
  } # min.pwr (for uniroot)

  if (requested == "n") {

    n.min <- n.levels.between + 1
    n.total <- try(stats::uniroot(function(n.total) min.pwr(f.squared, n.total, power), interval = c(n.min, 1e10))$root,
                   silent = TRUE)
    if (inherits(n.total, "try-error") || n.total == 1e+10) stop("Design is not feasible.", call. = FALSE)

  } else if (requested == "es") {

    f.squared <- try(stats::uniroot(function(f.squared) min.pwr(f.squared, n.total, power), interval = c(0, 1), tol = 1e-12)$root,
                     silent = TRUE)
    if (inherits(f.squared, "try-error")) stop("Design is not feasible.", call. = FALSE)

    eta.squared <- (f.squared / correct4rho) / (1 + (f.squared / correct4rho))

  }

  if (ceil.n) n.total <- ceiling(n.total / n.levels.between) * n.levels.between

  # calculate power (if requested == "power") or update it (if requested == "n" / "es")
  pwr.obj <- pwr.mixed(f.squared = f.squared, null.f.squared = null.f.squared,
                       n.total = n.total, n.levels.between = n.levels.between, n.levels.within = n.levels.within,
                       epsilon = epsilon, alpha = alpha, effect = effect)

  effect_bw <- fmt_effect_bw(effect, n.levels.between, n.levels.within)

  if (verbose > 0) {

    print.obj <- list(test = fmt_test_anovamxd(n.levels.between, n.levels.within),
                      requested = requested,
                      effect = effect_bw,
                      n.total = n.total,
                      factor.levels = factor.levels,
                      ncp = pwr.obj$lambda,
                      null.ncp = pwr.obj$null.lambda,
                      eta.squared = eta.squared,
                      power = pwr.obj$power,
                      alpha = alpha,
                      f.alpha = pwr.obj$f.alpha,
                      df1 = pwr.obj$df1,
                      df2 = pwr.obj$df2)

    .print.pwrss.ancova(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = "F",
                           effect = effect_bw,
                           df1 = pwr.obj$df1,
                           df2 = pwr.obj$df2,
                           ncp = pwr.obj$lambda,
                           null.ncp = pwr.obj$null.lambda,
                           f.alpha = pwr.obj$f.alpha,
                           eta.squared = eta.squared,
                           power = pwr.obj$power,
                           n.total = n.total),
                      class = c("pwrss", "f", "anova_mixed")))

} # power.f.anova.mixed

#' @export pwrss.f.rmanova
pwrss.f.rmanova <- function(eta2 = NULL, f2 = NULL,
                            corr.rm = 0.50, n.levels = 2, n.rm = 2,
                            epsilon = 1, alpha = 0.05,
                            type = c("between", "within", "interaction"),
                            n = NULL, power = NULL, verbose = TRUE) {

  type <- match.arg(type)
  verbose <- ensure.verbose(verbose)

  if (all(check.not_null(f2, eta2))) {
    stop("Effect size conflict for the alternative. Specify only either `eta2` or `f2`.", call. = FALSE)
  } else if (check.not_null(f2)) {
    eta2 <- f2 / (1 + f2)
  }

  mixed.anova.obj <- power.f.mixed.anova(eta.squared = eta2, n.total = n,
                                         factor.levels = c(n.levels, n.rm),
                                         factor.type = c("between", "within"),
                                         rho.within = corr.rm, epsilon = epsilon,
                                         power = power, alpha = alpha,
                                         effect = type, ceil.n = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.f.mixed.anova() function. \n")

  invisible(mixed.anova.obj)

} # pwrss.f.rmanova()

fmt_test_anovamxd <- function(n.levels.between, n.levels.within) {
  paste0(rep("Mixed-Effects ",     n.levels.between >  1 && n.levels.within > 1),
         rep("Repeated Measures ", n.levels.between == 1 && n.levels.within > 1),
         "Analysis of Variance (F-Test)")
}

fmt_effect_bw <- function(effect, n.levels.between, n.levels.within) {
  effects <- c(sprintf("B(%d)", n.levels.between), sprintf("W(%d)", n.levels.within))
  paste(effects[sort(c(1, 2), decreasing = (effect == "within"))], collapse = ifelse(effect == "interaction", ":", "|"))
}
