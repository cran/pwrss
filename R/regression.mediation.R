####################
# mediation z test #
####################
# 'cp = 0' by default, implying complete mediation (it increases explanatory power of the covariate)
# use `r.squared.mediator` and `r.squared.outcome` to adjust standard error for other predictors in
# mediation and outcome model

#' Power Analysis for Indirect Effects in a Mediation Model (Z, Joint, and
#' Monte Carlo Tests)
#'
#' @description
#' Calculates power or sample size (only one can be NULL at a time) to test
#' indirect effects in a mediation model (Z-Test, Joint Test, and Monte Carlo
#' Interval Test). One can consider explanatory power of the covariates in the
#' mediator and outcome model via specifying R-squared values accordingly.
#' \code{power.z.mediation()} and \code{power.z.med()} are the same functions.
#'
#' Formulas are validated using Monte Carlo simulations.
#'
#' @details
#' * NOTE: The function \code{pwrss.z.mediation()} (or its alias
#'   \code{pwrss.z.med()}) are no longer supported. However, they will remain
#'   available as wrappers for the \code{power.z.mediation} function.
#'
#' @aliases power.z.mediation pwrss.z.mediation power.z.med pwrss.z.med
#'
#'
#' @param beta.a              regression coefficient for X -> M path. One can
#'                            use standardized regression coefficient, but
#'                            should keep \code{sd.predictor = 1} and
#'                            \code{sd.mediator = 1} or leave them out as they
#'                            are default specifications.
#' @param beta.b              regression coefficient for M -> Y path. One can
#'                            use standardized regression coefficient, but
#'                            should keep \code{sd.mediator = 1} and
#'                            \code{sd.outcome = 1} or leave them out as they
#'                            are default specifications.
#' @param ab.ratio            'beta.a' / 'beta.b' ratio (can be negative) when minimum
#'                            detectable effect is of interest (either beta.a or beta.b, or both).
#' @param req.sign            Sign of the indirect effect (beta.a * beta.b product), when minimum
#'                            detectable effect is of interest (either beta.a or beta.b, or both)).
#' @param beta.cp             regression coefficient for X -> Y path (the
#'                            direct path). One can use standardized regression
#'                            coefficient, but should keep
#'                            \code{sd.predictor = 1} and \code{sd.outcome = 1}
#'                            or leave them out as they are default
#'                            specifications.
#' @param sd.predictor        standard deviation of the predictor (X). For a
#'                            binary predictor,
#'                            \code{sd.predictor = sqrt(p * (1 - p))} where
#'                            \code{p} is the proportion of subjects in one of
#'                            the groups.
#' @param sd.mediator         standard deviation of the mediator (M).
#' @param sd.outcome          standard deviation of the outcome (Y).
#' @param r.squared.mediator  R-squared value for the mediator model (M ~ X).
#'                            The default is \code{r.squared.mediator =
#'                            beta.a ^ 2 * sd.predictor ^ 2 / sd.mediator ^ 2}
#'                            assuming that X is the only predictor. Thus, an
#'                            \code{r.squared.mediator} below this value will
#'                            throw a warning. To consider other covariates in
#'                            the mediator model provide a value greater than
#'                            the default.
#' @param r.squared.outcome   R-squared value for the outcome model
#'                            (Y ~ M + X). The default is
#'                            \code{r.squared.outcome = (beta.b ^ 2 *
#'                            sd.mediator ^ 2 + beta.cp ^ 2 * sd.predictor ^ 2)
#'                            / sd.outcome ^ 2} assuming that M and X are the
#'                            only predictors. Thus, an
#'                            \code{r.squared.outcome} below this value will
#'                            throw a warning. To consider other covariates in
#'                            the outcome model provide a value greater than
#'                            the default.
#' @param n                   integer; sample size.
#' @param power               statistical power, defined as the probability of
#'                            correctly rejecting a false null hypothesis,
#'                            denoted as \eqn{1 - \beta}.
#' @param alpha               type 1 error rate, defined as the probability of
#'                            incorrectly rejecting a true null hypothesis,
#'                            denoted as \eqn{\alpha}.
#' @param alternative         character; the direction or type of the
#'                            hypothesis test: "two.sided" or "one.sided".
#' @param method              character; "sobel", "aroian", "goodman", "joint"
#'                            or "monte.carlo". "joint" and "monte.carlo"
#'                            methods can not be used for sample size
#'                            calculation.
#' @param n.simulation        integer; number of replications (applies when
#'                            method = "monte.carlo").
#' @param n.draws             integer; number of draws from the distribution of
#'                            the path coefficients for each replication
#'                            (applies when method = "monte.carlo").
#' @param ceil.n              logical; whether sample size should be rounded
#'                            up. \code{TRUE} by default.
#' @param verbose             \code{1} by default (returns test, hypotheses,
#'                            and results), if \code{2} a more detailed output
#'                            is given (plus key parameters and definitions), if
#'                            \code{0} no output is printed on the console.
#' @param utf                 logical; whether the output should show Unicode
#'                            characters (if encoding allows for it).
#'                            \code{FALSE} by default.
#'
#' @return
#'   \item{parms}{list of parameters used in calculation.}
#'   \item{test}{type of the statistical test ("Z-Test", "Joint Test", or
#'               "Monte Carlo Interval Test").}
#'   \item{mean}{mean of the alternative distribution.}
#'   \item{sd}{standard deviation of the alternative distribution.}
#'   \item{null.mean}{mean of the null distribution.}
#'   \item{null.sd}{standard deviation of the null distribution.}
#'   \item{z.alpha}{critical value(s).}
#'   \item{power}{statistical power \eqn{(1-\beta)}.}
#'   \item{n}{sample size.}
#'
#' @references
#'   Aroian, L. A. (1947). The probability function of the product of
#'   two normally distributed variables. *Annals of Mathematical Statistics,
#'   18*(2), 265-271.
#'
#'   Goodman, L. A. (1960). On the exact variance of products. *Journal of the
#'   American Statistical Association, 55*(292), 708-713.
#'
#'   MacKinnon, D. P., & Dwyer, J. H. (1993). Estimating mediated effects in
#'   prevention studies. *Evaluation Review, 17*(2), 144-158.
#'
#'   MacKinnon, D. P., Warsi, G., & Dwyer, J. H. (1995). A simulation study of
#'   mediated effect measures. *Multivariate Behavioral Research, 30*(1),
#'   41-62.
#'
#'   Preacher, K. J., & Hayes, A. F. (2004). SPSS and SAS procedures for
#'   estimating indirect effects in simple mediation models. *Behavior Research
#'   Methods, Instruments, & Computers, 36*, 717-731.
#'
#'   Preacher, K. J., & Hayes, A. F. (2008). Asymptotic and resampling strategies
#'   for assessing and comparing indirect effects in multiple mediator models.
#'   *Behavior Research Methods, 40*, 879-891.
#'
#'   Sobel, M. E. (1982). Asymptotic intervals for indirect effects in structural
#'   equations models. In S. Leinhart (Ed.), *Sociological methodology 1982* (pp.
#'   290-312). Jossey-Bass.
#'
#' @examples
#' # with standardized coefficients
#'
#' ## statistical power
#' power.z.mediation(beta.a = 0.25,
#'             beta.b = 0.25,
#'             beta.cp = 0.10,
#'             n = 200)
#'
#' ## minimum required sample size
#' power.z.mediation(beta.a = 0.25,
#'             beta.b = 0.25,
#'             beta.cp = 0.10,
#'             power = 0.80)
#'
#' ## adjust for covariates in the outcome model
#' power.z.mediation(beta.a = 0.25,
#'             beta.b = 0.25,
#'             beta.cp = 0.10,
#'             r.squared.outcome = 0.50,
#'             power = 0.80)
#'
#' # with binary predictor X such as treatment/control variable
#' # in this case standardized coefficients for path a and cp would be Cohen's d values
#'
#' ## statistical power
#' p <- 0.50 # proportion of subjects in one group
#' power.z.mediation(beta.a = 0.40,
#'             beta.b = 0.25,
#'             beta.cp = 0.10,
#'             sd.predictor = sqrt(p*(1-p)),
#'             n = 200)
#'
#' ## minimum required sample size
#' power.z.mediation(beta.a = 0.40,
#'             beta.b = 0.25,
#'             beta.cp = 0.10,
#'             sd.predictor = sqrt(p*(1-p)),
#'             power = 0.80)
#'
#' ## adjust for covariates in the outcome model
#' power.z.mediation(beta.a = 0.40,
#'             beta.b = 0.25, beta.cp = 0.10,
#'             r.squared.outcome = 0.50,
#'             sd.predictor = sqrt(p*(1-p)),
#'             power = 0.80)
#'
#' @export power.z.mediation
power.z.mediation  <- function(beta.a = NULL, beta.b = NULL, ab.ratio = 1, req.sign = "+", beta.cp = 0,
                               sd.predictor = 1, sd.mediator = 1, sd.outcome = 1,
                               r.squared.mediator = NULL, r.squared.outcome = NULL,
                               n = NULL, power = NULL, alpha = 0.05,
                               alternative = c("two.sided", "one.sided"),
                               method = c("sobel", "aroian", "goodman", "joint", "monte.carlo"),
                               n.simulation = 1000, n.draws = 1000,
                               ceil.n = TRUE, verbose = 1, utf = FALSE) {

  alternative <- match.arg(alternative)
  method <- match.arg(method)
  func.parms <- as.list(environment())

  if (!is.null(beta.a)) check.numeric(beta.a)
  if (!is.null(beta.b)) check.numeric(beta.b)
  check.numeric(ab.ratio)
  beta.cp <- ifelse(is.null(beta.cp), 0, beta.cp)
  check.numeric(beta.cp)
  check.nonnegative(sd.predictor, sd.mediator, sd.outcome)
  if (!is.null(r.squared.mediator)) {
    check.proportion(r.squared.mediator)
    if (!is.null(beta.a))
      if (r.squared.mediator < beta.a ^ 2 * sd.predictor ^ 2 / sd.mediator ^ 2)
        warning("Specified `r.squared.mediator` is smaller than the base `r.squared.mediator`.", call. = FALSE)
  }
  if (!is.null(r.squared.outcome)) {
    check.proportion(r.squared.outcome)
    if (r.squared.outcome == 0 && "beta.cp" %in% names(match.call()))
      warning("Ignoring any specification to `beta.cp`.", call. = FALSE)
    if (!is.null(beta.b))
      if (r.squared.outcome < (beta.b ^ 2 * sd.mediator ^ 2 + beta.cp ^ 2 * sd.predictor ^ 2) / sd.outcome ^ 2)
        warning("Specified `r.squared.outcome` is smaller than the base `r.squared.outcome`.", call. = FALSE)
  }
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.power(power)
  check.proportion(alpha)
  check.sample.size(n.simulation, n.draws)
  check.logical(ceil.n, utf)
  verbose <- ensure.verbose(verbose)
  requested <- get.requested(es = NA, n = n, power = power)
# requested <- get.requested(es = list(beta.a, beta.b), n = n, power = power)

  pwr <- function(beta.a, beta.b, sd.predictor, sd.mediator, sd.outcome, r.squared.mediator, r.squared.outcome,
                  n, method) {

    se.beta.a <- sqrt((1 / n) * (sd.mediator ^ 2) * (1 - r.squared.mediator) / (sd.predictor ^ 2))
    se.beta.b <- sqrt((1 / n) * (sd.outcome ^ 2) * (1 - r.squared.outcome) / ((sd.mediator ^ 2) * (1 - r.squared.mediator)))

    if (method == "sobel") {

      sobel.se <- sqrt(beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2)
      lambda <- (beta.a * beta.b) / sobel.se

    } else if (method == "aroian") {

      aroian.se <- sqrt(beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2 + se.beta.a ^ 2 * se.beta.b ^ 2)
      lambda <- (beta.a * beta.b) / aroian.se

    } else if (method == "goodman") {

      goodman.var <- beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2 - se.beta.a ^ 2 * se.beta.b ^ 2
      if (goodman.var <= 0) stop("Design is not feasible for Goodman's Z-Test", call. = FALSE)
      goodman.se <- sqrt(goodman.var)
      lambda <- (beta.a * beta.b) / goodman.se

    } else if (method == "joint") {

      power.beta.a <- power.z.test(mean = beta.a / se.beta.a, alpha = alpha, alternative = alternative,
                                   plot = FALSE, verbose = 0)$power
      power.beta.b <- power.z.test(mean = beta.b / se.beta.b, alpha = alpha, alternative = alternative,
                                   plot = FALSE, verbose = 0)$power
      power <- power.beta.a * power.beta.b

    } else if (method == "monte.carlo") {

      beta.a <- abs(beta.a)
      beta.b <- abs(beta.b)
      reject <- numeric(0)
      for (i in 1:n.simulation) {
        beta.a.star <- stats::rnorm(1, beta.a, se.beta.a)
        beta.b.star <- stats::rnorm(1, beta.b, se.beta.b)
        reject <- c(reject, stats::quantile(stats::rnorm(n.draws, beta.a.star, se.beta.a) * stats::rnorm(n.draws, beta.b.star, se.beta.b),
                                            probs = ifelse(alternative == "two.sided", alpha / 2, alpha), na.rm = TRUE) > 0)
      }
      power <- mean(reject)

    }

    if (method %in% c("sobel", "aroian", "goodman")) {
      power.z.test(mean = lambda, alpha = alpha, alternative = alternative, plot = FALSE, verbose = 0)
    } else {
      list(alternative = alternative, mean = NA, sd = NA, null.mean = NA, alpha = alpha, z.alpha = NA, power = power)
    }

  } # pwr()

  min.pwr <- function(beta.a, beta.b, n, power) {
    if (is.null(r.squared.mediator)) r.squared.mediator <-  beta.a ^ 2 * sd.predictor ^ 2 / sd.mediator ^ 2
    if (is.null(r.squared.outcome))  r.squared.outcome  <- (beta.b ^ 2 * sd.mediator ^ 2 + beta.cp ^ 2 * sd.predictor ^ 2) / sd.outcome ^ 2

    power - pwr(beta.a = beta.a, beta.b = beta.b,
                sd.predictor = sd.predictor, sd.mediator = sd.mediator, sd.outcome = sd.outcome,
                r.squared.mediator = r.squared.mediator, r.squared.outcome = r.squared.outcome,
                n = n, method = method)$power

  } # min.pwr (for uniroot and optimize)

  if (requested == "n") {

    if (method %in% c("joint", "monte.carlo"))
      stop("Sample size calculation not supported by this method.", call. = FALSE)

    n <- stats::uniroot(function(n) min.pwr(beta.a, beta.b, n, power), interval = c(10, 1e10))$root

#  } else if (requested == "es") {

#    if (method %in% c("joint", "monte.carlo"))
#      stop("Minimum detectable effect calculation not supported by this method.", call. = FALSE)

#    if (is.null(r.squared.mediator) || is.null(r.squared.outcome))
#      stop("R-squared values cannot be NULL when minimum detectable effect is of interest. Needed for identification.", call. = FALSE)

#    # lim.a <- sqrt(1 - r.squared.mediator) * (sd.mediator / sd.predictor)
#    # lim.b <- sqrt(1 - r.squared.outcome) * (sd.outcome / sd.mediator)
#    lim.a <- abs(ab.ratio) * (sd.mediator / sd.predictor)
#    lim.b <- abs(1 /  ab.ratio) * (sd.outcome / sd.mediator)
#    pos.sign <- check.pos_sign(req.sign)

#    if (is.null(beta.a) && !is.null(beta.b)) {

#      val.rng <- c(-lim.a, 0, lim.a)[ifelse((sign(beta.b) > 0) == pos.sign, -1, -3)]

#      beta.a <- try(stats::optimize(f = function(beta.a) { min.pwr(beta.a, beta.b, n, power) ^ 2; }, interval = val.rng)$minimum)
#      if (inherits(beta.a, "try-error"))
#        stop(sprintf("Design is not feasible. Try `req.sign` = \"%s\".", ifelse(pos.sign, "-", "+")), call. = FALSE)

#    } else if (!is.null(beta.a) && is.null(beta.b)) {

#      val.rng <- c(-lim.b, 0, lim.b)[ifelse((sign(beta.a) > 0) == pos.sign, -1, -3)]

#      beta.b <- try(stats::optimize(f = function(beta.b) min.pwr(beta.a, beta.b, n, power) ^ 2, interval = val.rng)$minimum)
#      if (inherits(beta.b, "try-error"))
#        stop(sprintf("Design is not feasible. Try `req.sign` = \"%s\".", ifelse(pos.sign, "-", "+")), call. = FALSE)

#    } else if (is.null(beta.a) && is.null(beta.b)) {

#      if (pos.sign != (sign(ab.ratio) > 0)) {
#        warning("`req.sign` is positive but 'ab.ratio' is negative. Ignoring the sign of 'ab.ratio' but not its value.", call. = FALSE)
#        ab.ratio <- -ab.ratio
#      }

#      val.rng <- c(-lim.b, 0, lim.b)[ifelse((sign(ab.ratio) > 0) == pos.sign, -1, -3)]

#      beta.b <-  try(stats::optimize(f = function(beta.b) min.pwr(beta.b * ab.ratio, beta.b, n, power) ^ 2, interval = val.rng)$minimum)
#      if (inherits(beta.b, "try-error"))
#        stop(sprintf("Design is not feasible. Try `req.sign` = \"%s\".", ifelse(pos.sign, "-", "+")), call. = FALSE)
#      beta.a <- beta.b * ab.ratio

#    } # beta.a and beta.b

  } # sample size or effect size

  if (ceil.n) n <- ceiling(n)

  if (is.null(r.squared.mediator)) r.squared.mediator <-  beta.a ^ 2 * sd.predictor ^ 2 / sd.mediator ^ 2
  if (is.null(r.squared.outcome))  r.squared.outcome  <- (beta.b ^ 2 * sd.mediator ^ 2 + beta.cp ^ 2 * sd.predictor ^ 2)  / sd.outcome ^ 2

  # calculate power (if requested == "power") or update it (if requested == "n" or "es")
  pwr.obj <- pwr(beta.a = beta.a, beta.b = beta.b,
                 sd.predictor = sd.predictor, sd.mediator = sd.mediator, sd.outcome = sd.outcome,
                 r.squared.mediator = r.squared.mediator, r.squared.outcome = r.squared.outcome,
                 n = n, method = method)

  beta.indirect <- beta.a * beta.b
  std.beta.a <- beta.a * sd.predictor / sd.mediator
  std.beta.b <- beta.b * sd.mediator / sd.outcome
  std.beta.cp <- beta.cp * sd.predictor / sd.outcome
  std.beta.indirect <-  std.beta.a * std.beta.b

  if (verbose > 0) {

    print.obj <- list(test = "Indirect Effect in a Mediation Model",
                      method = method,
                      requested = requested,
                      std.beta.a = std.beta.a,
                      std.beta.b = std.beta.b,
                      std.beta.cp = std.beta.cp,
                      std.beta.indirect = std.beta.indirect,
                      mean = pwr.obj$mean,
                      sd = 1,
                      null.mean = pwr.obj$null.mean,
                      null.sd = 1,
                      alpha = alpha,
                      alternative = alternative,
                      z.alpha = pwr.obj$z.alpha,
                      n = n,
                      power = pwr.obj$power)

    .print.pwrss.med(print.obj, verbose = verbose, utf = utf)

  } # verbose

  invisible(structure(list(parms = func.parms,
                           test = switch(method, `joint` = "joint", `monte.carlo` = "monte.carlo", "z"),
                           beta.a = beta.a,
                           beta.b = beta.b,
                           beta.indirect = beta.indirect,
                           beta.cp = beta.cp,
                           r.squared.mediator = r.squared.mediator,
                           r.squared.outcome = r.squared.outcome,
                           std.beta.a = std.beta.a,
                           std.beta.b = std.beta.b,
                           std.beta.cp = std.beta.cp,
                           std.beta.indirect = std.beta.indirect,
                           mean = pwr.obj$mean,
                           sd = 1,
                           null.mean = pwr.obj$null.mean,
                           null.sd = 1,
                           z.alpha = pwr.obj$z.alpha,
                           power = pwr.obj$power,
                           n = n),
                      class = c("pwrss", "z", "mediation")))

} # end of power.z.mediation()

#' @rdname power.z.mediation
#' @export power.z.med
power.z.med <- power.z.mediation

#' @export pwrss.z.mediation
pwrss.z.mediation  <- function(a, b, cp = 0,
                               sdx = 1, sdm = 1, sdy = 1,
                               r2m.x = a ^ 2 * sdx ^ 2 / sdm ^ 2,
                               r2y.mx = (b ^ 2 * sdm ^ 2 + cp ^ 2 * sdx ^ 2) / sdy ^ 2,
                               n = NULL, power = NULL,
                               alpha = 0.05, alternative = c("not equal", "less", "greater"),
                               mc = TRUE, nsims = 1000, ndraws = 1000, verbose = TRUE) {

  alternative <- match.arg(alternative)
  if (alternative %in% c("less", "greater")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  func.parms <- as.list(environment())
  verbose <- ensure.verbose(verbose)

  if (is.null(power)) {

    sobel.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                   sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                   r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                   n = n, power = power, alpha = alpha,
                                   alternative = alternative,
                                   method = "sobel",
                                   n.simulation = nsims,
                                   n.draws = ndraws,
                                   ceil.n = TRUE,
                                   verbose = 0)

    aroian.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                    sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                    r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                    n = n, power = power, alpha = alpha,
                                    alternative = alternative,
                                    method = "aroian",
                                    n.simulation = nsims,
                                    n.draws = ndraws,
                                    ceil.n = TRUE,
                                    verbose = 0)

    goodman.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                     sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                     r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                     n = n, power = power, alpha = alpha,
                                     alternative = alternative,
                                     method = "goodman",
                                     n.simulation = nsims,
                                     n.draws = ndraws,
                                     ceil.n = TRUE,
                                     verbose = 0)

    joint.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                   sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                   r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                   n = n, power = power, alpha = alpha,
                                   alternative = alternative,
                                   method = "joint",
                                   n.simulation = nsims,
                                   n.draws = ndraws,
                                   ceil.n = TRUE,
                                   verbose = 0)

    monte.carlo.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                         sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                         r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                         n = n, power = power, alpha = alpha,
                                         alternative = alternative,
                                         method = "monte.carlo",
                                         n.simulation = nsims,
                                         n.draws = ndraws,
                                         ceil.n = TRUE,
                                         verbose = 0)

    power.out <- data.frame(rbind(
      c(sobel.obj$mean, sobel.obj$n, sobel.obj$power),
      c(aroian.obj$mean, aroian.obj$n, aroian.obj$power),
      c(goodman.obj$mean, goodman.obj$n, goodman.obj$power),
      c(joint.obj$mean, joint.obj$n, joint.obj$power),
      c(monte.carlo.obj$mean, monte.carlo.obj$n, monte.carlo.obj$power)
    ))

    colnames(power.out) <- c("non-centrality", "n", "power")
    power.out$method <- c("Sobel", "Aroian", "Goodman", "Joint", "Monte Carlo")

  } else if (is.null(n)) {

    sobel.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                   sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                   r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                   n = n, power = power, alpha = alpha,
                                   alternative = alternative,
                                   method = "sobel",
                                   n.simulation = nsims,
                                   n.draws = ndraws,
                                   ceil.n = TRUE,
                                   verbose = 0)

    aroian.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                    sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                    r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                    n = n, power = power, alpha = alpha,
                                    alternative = alternative,
                                    method = "aroian",
                                    n.simulation = nsims,
                                    n.draws = ndraws,
                                    ceil.n = TRUE,
                                    verbose = 0)

    goodman.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                     sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                     r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                     n = n, power = power, alpha = alpha,
                                     alternative = alternative,
                                     method = "goodman",
                                     n.simulation = nsims,
                                     n.draws = ndraws,
                                     ceil.n = TRUE,
                                     verbose = 0)

    power.out <- data.frame(rbind(
      c(sobel.obj$mean, sobel.obj$n, sobel.obj$power),
      c(aroian.obj$mean, aroian.obj$n, aroian.obj$power),
      c(goodman.obj$mean, goodman.obj$n, goodman.obj$power),
      c(NA, NA, NA),
      c(NA, NA, NA)
    ))

    colnames(power.out) <- c("non-centrality", "n", "power")
    power.out$method <- c("Sobel", "Aroian", "Goodman", "Joint", "Monte Carlo")

  }


  if (is.null(power)) {
    ncp.vec <- c(sobel = power.out$mean[1], aroian = power.out$mean[2], goodman = power.out$mean[3], joint = power.out$mean[4],
                 mc = power.out$mean[5])
    n.vec <- c(sobel = power.out$n[1], aroian = power.out$n[2], goodman = power.out$n[3], joint = power.out$n[4], mc = power.out$n[5])
    power.vec <- c(sobel = power.out$power[1], aroian = power.out$power[2], goodman = power.out$power[3], joint = power.out$power[4],
                   mc = power.out$power[5])
  } else if (is.null(n)) {
    ncp.vec <- c(sobel = power.out$mean[1], aroian = power.out$mean[2], goodman = power.out$mean[3], joint = power.out$mean[4],
                 mc = power.out$mean[5])
    n.vec <- c(sobel = power.out$n[1], aroian = power.out$n[2], goodman = power.out$n[3], joint = power.out$n[4], mc = power.out$n[5])
    power.vec <- c(sobel = power.out$power[1], aroian = power.out$power[2], goodman = power.out$power[3], joint = power.out$power[4],
                   mc = power.out$power[5])
  }

  if (verbose > 0) {
    cat(" Indirect Effect in a Mediation Model",
        "\n ====================================\n",
        sep = "")

    ifelse(is.null(power),
           print(power.out[,    c("method", "non-centrality", "n", "power")], row.names = FALSE),
           print(power.out[1:3, c("method", "non-centrality", "n", "power")], row.names = FALSE))

    cat(" ------------------------------------\n",
        " Type 1 error rate: ", round(alpha, 3), "\n\n",
        sep = "")

    # cat2("beta.indirect := beta.a * beta.b \n", "green")

  } # verbose

# cat("This function will be removed in the future. \n Please use power.z.mediation() function. \n")

invisible(structure(list(parms = func.parms,
                         test = c("z", "joint", "monte.carlo"),
                         ncp = ncp.vec,
                         power = power.vec,
                         n = n.vec),
                    class = c("pwrss", "z", "med", "defunct")))

} # pwrss.z.mediation()

#' @export pwrss.z.med
pwrss.z.med  <- pwrss.z.mediation
