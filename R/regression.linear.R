rsq.to.f <- function(r.squared.full, r.squared.reduced = 0, verbose = TRUE) {

  check.nonnegative(r.squared.full, r.squared.reduced)
  if (r.squared.full < r.squared.reduced) stop("Expecting `r.squared.full` > `r.squared.reduced`.", call. = FALSE)

  f.squared <- (r.squared.full - r.squared.reduced) / (1 - r.squared.full)

  if (verbose) print(c(f.squared = f.squared, f = sqrt(f.squared),
                      r.squared.full = r.squared.full,
                      r.squared.reduced = r.squared.reduced))

  invisible(list(f.squared = f.squared, f = sqrt(f.squared),
                 r.squared.full = r.squared.full,
                 r.squared.reduced = r.squared.reduced))
} # rsq.to.f

f.to.rsq <- function(f, r.squared.full = NULL, verbose = TRUE) {

  check.nonnegative(f)

  f.squared <- f ^ 2
  if (is.null(r.squared.full)) {
    r.squared.full <- f.squared / (1 + f.squared)
    r.squared.reduced <- 0
  } else {
    change.r.squared <- f.squared * (1 - r.squared.full)
    r.squared.reduced <- r.squared.full - change.r.squared
  }

  if (verbose) print(c(f.squared = f.squared, f = sqrt(f.squared),
                      r.squared.full = r.squared.full,
                      r.squared.reduced = r.squared.reduced))

  invisible(list(f.squared = f.squared, f = sqrt(f.squared),
                 r.squared.full = r.squared.full,
                 r.squared.reduced = r.squared.reduced))
} # f.to.rsq



############################
# linear regression f test #
############################

power.f.regression <- function(r.squared.change = NULL,
                               margin = 0,
                               k.total,
                               k.tested = k.total,
                               n = NULL, power = NULL, alpha = 0.05,
                               ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  check.proportion(alpha)
  check.logical(ceiling)
  check.sample.size(k.total, k.tested)
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.proportion(power)
  if (is.null(n) && is.null(power)) stop("`n` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n) && !is.null(power)) stop("Exactly one of the `n` or `power` should be `NULL`.", call. = FALSE)
  if (k.tested >  k.total) stop("'m.tested' cannot be greater than 'k.total'", call. = FALSE)
  if (!is.numeric(r.squared.change) || r.squared.change > 1 || r.squared.change < 0) stop("Incorrect value for `r.squared.change`.", call. = FALSE)

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  pwr.f.reg <- function(r.squared.change, margin, k.total, k.tested, n, alpha) {

    df1 <- k.tested
    df2 <- n - k.total - 1

    f.squared <- rsq.to.f(r.squared.full = r.squared.change, verbose = FALSE)$f.squared
    null.f.squared <- rsq.to.f(r.squared.full = margin, verbose = FALSE)$f.squared

    lambda <- f.squared * n
    null.lambda <- null.f.squared * n


    f.alpha <- qf(alpha, df1 = df1, df2 = df2, ncp = null.lambda, lower.tail = FALSE)
    power <- pf(f.alpha, df1 = df1, df2 = df2, ncp = lambda, lower.tail = FALSE)

    list(power = power, lambda = lambda, null.lambda = null.lambda,
         df1 = df1, df2 = df2, f.alpha = f.alpha)

  } # pwr

  ss.f.reg <- function(r.squared.change, margin, k.total, k.tested, power, alpha) {

    n <- try(silent = TRUE,
             suppressWarnings({
               uniroot(function(n) {
                 power - pwr.f.reg(r.squared.change = r.squared.change, margin = margin,
                                   k.total = k.total, k.tested = k.tested,
                                   n = n, alpha = alpha)$power
               }, interval = c(k.total + 2, 1e10))$root
             }) # supressWarnings
    ) # try

    if (inherits(n, "try-error") || n == 1e10) stop("Design is not feasible.", call. = FALSE)

    n

  } # ss.f.reg()


  if (is.null(power)) {

    pwr.obj <- pwr.f.reg(r.squared.change = r.squared.change, margin = margin,
                         k.total = k.total, k.tested = k.tested,
                         n = n, alpha = alpha)
    power <- pwr.obj$power
    df1 <- pwr.obj$df1
    df2 <- pwr.obj$df2
    lambda <- pwr.obj$lambda
    null.lambda <- pwr.obj$null.lambda
    f.alpha <- pwr.obj$f.alpha

  } # pwr

  if (is.null(n)) {

    n <- ss.f.reg(r.squared.change = r.squared.change, margin = margin,
                  k.total = k.total, k.tested =  k.tested,
                  power = power, alpha = alpha)

    if (ceiling) {
      n <- ceiling(n)
    }

    pwr.obj <- pwr.f.reg(r.squared.change = r.squared.change, margin = margin,
                         k.total = k.total, k.tested = k.tested,
                         n = n, alpha = alpha)
    power <- pwr.obj$power
    df1 <- pwr.obj$df1
    df2 <- pwr.obj$df2
    lambda <- pwr.obj$lambda
    null.lambda <- pwr.obj$null.lambda
    f.alpha <- pwr.obj$f.alpha

  } # ss

  # verbose check
  if (is.logical(verbose)) {
    ifelse(isTRUE(verbose),
           verbose <- 1,
           verbose <- 0)
  } else if (is.numeric(verbose)) {
    if (length(verbose) == 1 && verbose %% 1 == 0) {
      ifelse(verbose %in% c(0, 1, 2),
             verbose <- verbose,
             verbose <- 1)
    }
  } else {
    verbose <- 1
  } # verbose

  if (verbose != 0) {

    ifelse(k.tested == k.total,
           test <- "Linear Regression (F-Test)",
           test <- "Hierarchical Linear Regression (F-Test)")

    print.obj <- list(requested = requested,
                      test = test,
                      r.squared.change = r.squared.change,
                      k.total = k.total,
                      k.tested = k.tested,
                      margin = margin,
                      n = n,
                      df1 = df1,
                      df2 = df2,
                      ncp.alternative = lambda,
                      ncp.null = null.lambda,
                      f.alpha = f.alpha,
                      alpha = alpha,
                      power = power)

    if (pretty) {
      .print.pwrss.f.regression(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.f.regression(print.obj, verbose = verbose)
    }

  } # verbose

  invisible(structure(list(parms = list(r.squared.change = r.squared.change,
                                        margin = margin,
                                        k.total = k.total, k.tested = k.tested,
                                        alpha = alpha,
                                        ceiling = ceiling, verbose = verbose,
                                        pretty = pretty),
                           test = "F",
                           df1 = df1,
                           df2 = df2,
                           ncp = lambda,
                           null.ncp = null.lambda,
                           f.alpha = f.alpha,
                           power = power,
                           n = n),
                      class = c("pwrss", "f", "regression")))


} # pwrss.f.regression()
power.f.reg <- power.f.regression


pwrss.f.regression <- function(r2 = 0.10, f2 = r2 / (1 - r2),
                        k = 1, m = k, alpha = 0.05,
                        n = NULL, power = NULL, verbose = TRUE) {

  user.parms <- as.list(match.call())
  names.user.parms <- names(user.parms)

  if (all(c("r2", "f2") %in% names.user.parms)) {
    stop("Specify either `r2` or `f2`.", call. = FALSE)
  } else {
    if ("f2" %in% names.user.parms) r2 <- f.to.rsq(f = sqrt(f2))$r.squared.full
  }

  pwrss.f.reg.obj <- power.f.regression(r.squared.change = r2, margin = 0,
                                        k.total = k, k.tested = m,
                                        n = n, power = power, alpha = alpha,
                                        ceiling = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.f.regression() function. \n")

  return(invisible(pwrss.f.reg.obj))

} # pwrss.f.regression
pwrss.f.reg <- pwrss.f.regression



#####################
# linear regression #
#####################

# if the predictor is binary
# provide sd.predictor = sqrt(p * (1 - p))
# p = proportion of subjects in treatment group
# use defaults if beta is standardized
power.t.regression <- function(beta, null.beta = 0, margin = 0,
                               sd.predictor = 1, sd.outcome = 1,
                               r.squared = (beta * sd.predictor / sd.outcome) ^ 2,
                               k.total = 1,
                               n = NULL, power = NULL, alpha = 0.05,
                               alternative = c("two.sided", "one.sided", "two.one.sided"),
                               ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  check.numeric(beta, null.beta)
  check.proportion(alpha)
  check.logical(ceiling)
  check.sample.size(k.total)
  check.positive(sd.predictor, sd.outcome, k.total)
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.proportion(power)

  alternative <- tolower(match.arg(alternative))

  if (is.null(n) && is.null(power)) stop("`n` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n) && !is.null(power)) stop("Exactly one of the `n` or `power` should be `NULL`.", call. = FALSE)
  if (!is.numeric(r.squared) || r.squared > 1 || r.squared < 0) stop("Incorrect value for `r.squared`, specify `r.squared` explicitly or modify `beta`, `sd.predictor`, `sd.outcome`.", call. = FALSE)
  if (r.squared < (beta * sd.predictor / sd.outcome) ^ 2) warning("`r.squared` is possibly larger.", call. = FALSE)

  if (alternative == "two.one.sided") {
    if (length(margin) != 2) stop("Specify `margin` in the form of margin = c(lower, upper).", call. = FALSE)
  } else {
    if (length(margin) != 1) stop("Specify only one value for the `margin`.", call. = FALSE)
    # if (beta == null.beta) stop("`beta` takes a value different from `null.beta` for 'one.sided' or 'two.sided' tests.", call. = FALSE)
  }

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  pwr.t.reg <- function(beta, null.beta, margin,
                        sd.outcome, sd.predictor,
                        n, k.total, r.squared,
                        alpha, alternative) {

    df <- n - k.total - 1
    lambda <- (beta - null.beta)  / ((sd.outcome / sd.predictor) * sqrt((1 - r.squared) / n))
    null.lambda <- margin / ((sd.outcome / sd.predictor) * sqrt((1 - r.squared) / n))

    pwr.obj <- power.t.test(ncp = lambda, null.ncp = null.lambda, df = df,
                            alpha = alpha, alternative = alternative,
                            plot = FALSE, verbose = FALSE)
    power <- pwr.obj$power
    t.alpha <- pwr.obj$t.alpha

    list(power = power, lambda = lambda, null.lambda = null.lambda, df = df, t.alpha = t.alpha)

  } # pwr.t.reg()


  ss.t.reg <- function(beta, null.beta, margin,
                       sd.outcome, sd.predictor,
                       power, k.total, r.squared,
                       alpha, alternative) {

    n <- try(silent = TRUE,
             suppressWarnings({
               uniroot(function(n) {
                 power - pwr.t.reg(beta = beta, null.beta = null.beta, margin = margin,
                                   sd.outcome = sd.outcome, sd.predictor = sd.predictor,
                                   n = n, k.total = k.total, r.squared = r.squared,
                                   alpha = alpha, alternative =  alternative)$power
               }, interval = c(k.total + 2, 1e10))$root
             }) # supressWarnings
    ) # try

    if (inherits(n, "try-error") || n == 1e10) stop("Design is not feasible.", call. = FALSE)

    n

  } # ss.t.reg()


  if (is.null(power)) {
    pwr.obj <- pwr.t.reg(beta = beta, null.beta = null.beta, margin = margin,
                         sd.outcome = sd.outcome, sd.predictor = sd.predictor,
                         n = n, k.total = k.total, r.squared = r.squared,
                         alpha = alpha, alternative =  alternative)
    power <- pwr.obj$power
    lambda <- pwr.obj$lambda
    null.lambda <- pwr.obj$null.lambda
    df <- pwr.obj$df
    t.alpha <- pwr.obj$t.alpha
  } # pwr

  if (is.null(n)) {
    n <- ss.t.reg(beta = beta, null.beta = null.beta, margin = margin,
                  sd.outcome = sd.outcome, sd.predictor = sd.predictor,
                  power = power, k.total = k.total, r.squared = r.squared,
                  alpha = alpha, alternative =  alternative)

    if (ceiling) n <- ceiling(n)
    pwr.obj <- pwr.t.reg(beta = beta, null.beta = null.beta, margin = margin,
                         sd.outcome = sd.outcome, sd.predictor = sd.predictor,
                         n = n, k.total = k.total, r.squared = r.squared,
                         alpha = alpha, alternative =  alternative)
    power <- pwr.obj$power
    lambda <- pwr.obj$lambda
    null.lambda <- pwr.obj$null.lambda
    df <- pwr.obj$df
    t.alpha <- pwr.obj$t.alpha
  } # ss

  std.beta <- beta * (sd.predictor / sd.outcome)
  std.null.beta <- null.beta * (sd.predictor / sd.outcome)
  std.margin <- margin * (sd.predictor / sd.outcome)

  # verbose check
  if (is.logical(verbose)) {
    ifelse(isTRUE(verbose),
           verbose <- 1,
           verbose <- 0)
  } else if (is.numeric(verbose)) {
    if (length(verbose) == 1 && verbose %% 1 == 0) {
      ifelse(verbose %in% c(0, 1, 2),
             verbose <- verbose,
             verbose <- 1)
    }
  } else {
    verbose <- 1
  } # verbose

  if (verbose != 0) {

    print.obj <- list(requested = requested,
                      test = "Linear Regression Coefficient (T-Test)",
                      alt = alternative,
                      std.beta = std.beta,
                      std.null.beta = std.null.beta,
                      std.margin = std.margin,
                      margin = margin,
                      n = n,
                      df = df,
                      ncp.alternative = lambda,
                      ncp.null = null.lambda,
                      t.alpha = t.alpha,
                      alpha = alpha,
                      power = power)

    if (pretty) {
      .print.pwrss.regression(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.regression(print.obj, verbose = verbose)
    }

  } # verbose

  invisible(structure(list(parms = list(beta = beta, null.beta = null.beta,
                                        margin = margin,
                                        sd.predictor = sd.predictor, sd.outcome = sd.outcome,
                                        r.squared = r.squared, k.total = k.total,
                                        alpha = alpha, alternative = alternative,
                                        ceiling = ceiling, verbose = verbose,
                                        pretty = pretty),
                           test = "t",
                           std.beta = std.beta,
                           std.null.beta = std.null.beta,
                           std.margin = std.margin,
                           df = df,
                           t.alpha = t.alpha,
                           ncp = lambda,
                           null.ncp = null.lambda,
                           power = power,
                           n = n),
                      class = c("pwrss", "t", "regression")))
} # power.t.regression()
power.t.reg <- power.t.regression


pwrss.t.regression <- function(beta1 = 0.25, beta0 = 0, margin = 0,
                               sdx = 1, sdy = 1,
                               k = 1, r2 = (beta1 * sdx / sdy) ^ 2,
                               alpha = 0.05, n = NULL, power = NULL,
                               alternative = c("not equal", "less", "greater",
                                               "non-inferior", "superior", "equivalent"),
                               verbose = TRUE) {

  alternative <- tolower(match.arg(alternative))
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
                                        ceiling = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.t.regression() function. \n")

  return(invisible(pwrss.t.reg.obj))

} # pwrss.t.regression
pwrss.t.reg <- pwrss.t.regression


# defunct
pwrss.z.mean   <- function(...) stop("This function is no longer available. Please use `power.t.student()`.", call. = FALSE)
pwrss.z.2means <- function(...) stop("This function is no longer available. Please use `power.t.student()` or `power.t.welch()`.", call. = FALSE)
pwrss.z.regression <- pwrss.z.reg <- function(...) stop("This function is no longer available. Please use `power.t.regression()`.", call. = FALSE)
