power.f.mixed.anova <- function(eta.squared,
                                null.eta.squared = 0,
                                factor.levels = c(2, 2),
                                factor.type = c("between", "within"),
                                rho.within = 0.50,
                                epsilon = 1,
                                n.total = NULL,
                                power = NULL, alpha = 0.05,
                                effect = c("between", "within", "interaction"),
                                ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  check.proportion(alpha)
  check.logical(ceiling)
  check.nonnegative(eta.squared, null.eta.squared)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n.total)) check.sample.size(n.total)

  for (i in 1:2) {
    factor.type.check <- factor.type[i]
    if (!is.character(factor.type.check) || !(factor.type.check %in% c("between", "within")))
      stop("The 'factor.type' argument must be specified as either c('between', 'within') or c('within', 'between'), indicating the order in which the corresponding values in 'factor.levels' are interpreted - specifically, which factor is treated as between-subjects and which as within-subjects.")
  }

  effect <- tolower(match.arg(effect))

  if (length(factor.levels) != 2 || length(factor.type) != 2) stop("Excatly two factors are allowed in this procedure.", call. = FALSE)
  if (all(factor.type == "within") || all(factor.type == "between")) stop("Two 'within' or two 'between' factors are not allowed in this procedure.", call. = FALSE)

  n.levels.between <- factor.levels[which(tolower(factor.type) == "between")]
  n.levels.within <- factor.levels[which(tolower(factor.type) == "within")]
  if (n.levels.within > 1 && epsilon <  1 / (n.levels.within - 1)) stop("Incorrect value for the non-sphericity correction factor (epsilon).", call. = FALSE)

  if (is.null(n.total) && is.null(power)) stop("`n.total` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n.total) && !is.null(power)) stop("Exactly one of the `n.total` or `power` should be `NULL`.", call. = FALSE)

  f.squared <- eta.squared / (1 - eta.squared)
  null.f.squared <- null.eta.squared / (1 - null.eta.squared)

  if (!is.na(rho.within)) {

    if (effect == "between") {
      f.squared <- f.squared  * (n.levels.within / (1 + (n.levels.within - 1) * rho.within))
      null.f.squared <- null.f.squared  * (n.levels.within / (1 + (n.levels.within - 1) * rho.within))
    } else {
      f.squared <- f.squared  * (n.levels.within / (1 - rho.within))
      null.f.squared <- null.f.squared  * (n.levels.within / (1 - rho.within))
    }

  } else {

    warning("Assuming that 'eta.squared' and 'null.eta.squared' are already adjusted for within-subject correlation.", call. = FALSE)

  }

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n.total")

  ss.mixed <- function(f.squared, null.f.squared, n.levels.between, n.levels.within, epsilon, alpha, power, effect) {

    n.min <- n.levels.between + 1

    n.total <- try(silent = TRUE,
                   uniroot(function(n.total) {
                     if (effect == "between") {
                       df1 <- n.levels.between - 1
                       df2 <- n.total - n.levels.between
                     } else if (effect == "within") {
                       df1 <- (n.levels.within - 1) * epsilon
                       df2 <- (n.total - n.levels.between) * (n.levels.within - 1) * epsilon
                     } else if (effect == "interaction") {
                       df1 <- (n.levels.between - 1) * (n.levels.within - 1) * epsilon
                       df2 <- (n.total - n.levels.between) * (n.levels.within - 1) * epsilon
                     } else {
                       stop("Unknown effect", call. = FALSE)
                     }
                     u <- df1
                     v <- df2
                     lambda <- f.squared * n.total * epsilon
                     null.lambda <- null.f.squared * n.total * epsilon
                     f.alpha <- qf(alpha, df1 = u, df2 = v, ncp = null.lambda, lower.tail = FALSE)

                     power - pf(f.alpha, df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)

                   }, interval = c(n.min, 1e10))$root
                   ) # try

    if (inherits(n.total, "try-error") || n.total == 1e+10) stop("Design is not feasible.", call. = FALSE)

    n.total

  }

  pwr.mixed <- function(f.squared, null.f.squared, n.total,
                        n.levels.between, n.levels.within, epsilon, alpha, effect) {

    if (effect == "between") {
      df1 <- n.levels.between - 1
      df2 <- n.total - n.levels.between
    } else if (effect == "within") {
      df1 <- (n.levels.within - 1) * epsilon
      df2 <- (n.total - n.levels.between) * df1
    } else if (effect == "interaction") {
      df1 <- (n.levels.between - 1) * (n.levels.within - 1) * epsilon
      df2 <- (n.total - n.levels.between) * (n.levels.within - 1) * epsilon
    } else {
      stop("Unknown effect", call. = FALSE)
    }

    u <- df1
    v <- df2
    if (u < 1 || v < 1) stop("Design is not feasible", call. = FALSE)
    lambda <- f.squared * n.total * epsilon
    null.lambda <- null.f.squared * n.total * epsilon
    f.alpha <- qf(alpha, df1 = u, df2 = v, ncp = null.lambda, lower.tail = FALSE)
    power <- pf(f.alpha, df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)

    list(power = power, lambda = lambda, null.lambda = null.lambda,
         df1 = u, df2 = v, f.alpha = f.alpha)

  } # pwr.mixed()


  if (is.null(n.total)) {

    n.total <- ss.mixed(f.squared = f.squared, null.f.squared = null.f.squared,
                        n.levels.between = n.levels.between,
                        n.levels.within = n.levels.within,
                        epsilon = epsilon, alpha = alpha,
                        power = power, effect = effect)

    if (ceiling) n.total <- ceiling(n.total / n.levels.between) * n.levels.between

    pwr.obj <- pwr.mixed(f.squared = f.squared, null.f.squared = null.f.squared,
                         n.total = n.total, n.levels.between = n.levels.between,
                         n.levels.within = n.levels.within,
                         epsilon = epsilon, alpha = alpha, effect = effect)

    power <- ncp <- pwr.obj$power
    df1 <-  pwr.obj$df1
    df2 <-  pwr.obj$df2
    ncp <-  pwr.obj$lambda
    null.ncp <-  pwr.obj$null.lambda
    f.alpha <- pwr.obj$f.alpha

  } # n.total is null

  if (is.null(power)) {

    pwr.obj <- pwr.mixed(f.squared = f.squared, null.f.squared = null.f.squared,
                         n.total = n.total, n.levels.between = n.levels.between,
                         n.levels.within = n.levels.within,
                         epsilon = epsilon, alpha = alpha, effect = effect)

    power <- ncp <- pwr.obj$power
    df1 <-  pwr.obj$df1
    df2 <-  pwr.obj$df2
    ncp <-  pwr.obj$lambda
    null.ncp <-  pwr.obj$null.lambda
    f.alpha <- pwr.obj$f.alpha

  } # power is null

  if (effect == "between")     effect_bw <- paste0(c("B", "W"), "(", c(n.levels.between, n.levels.within), ")", collapse = "|")
  if (effect == "within")      effect_bw <- paste0(c("W", "B"), "(", c(n.levels.within, n.levels.between), ")", collapse = "|")
  if (effect == "interaction") effect_bw <- paste0(c("B", "W"), "(", c(n.levels.between, n.levels.within), ")", collapse = ":")

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

    if (n.levels.between == 1) test <- "Repeated Measures Analysis of Variance (F-Test)"
    if (n.levels.within == 1) test <- "Analysis of Variance (F-Test)"
    if (n.levels.within > 1 && n.levels.between > 1) test <- "Mixed-Effects Analysis of Variance (F-Test)"

    print.obj <- list(test = test, effect = effect_bw, n.total = n.total,
                      requested = requested, factor.levels = factor.levels,
                      power = power, ncp = ncp, null.ncp = null.ncp,
                      alpha = alpha, f.alpha = f.alpha, df1 = df1, df2 = df2)

    if (pretty) {
      .print.pwrss.ancova(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.ancova(print.obj, verbose = verbose)
    }

  } # verbose

  invisible(structure(list(parms = list(eta.squared = eta.squared,
                                        null.eta.squared = null.eta.squared,
                                        factor.levels = factor.levels,
                                        factor.type = factor.type,
                                        rho.within = rho.within,
                                        epsilon = epsilon,
                                        alpha = alpha, effect = effect,
                                        ceiling = ceiling, verbose = verbose,
                                        pretty = pretty),
                           effect = effect_bw,
                           test = "F",
                           df1 = df1,
                           df2 = df2,
                           ncp = ncp,
                           null.ncp = null.ncp,
                           f.alpha = f.alpha,
                           power = power,
                           n.total = n.total),
                      class = c("pwrss", "f", "anova_mixed")))

} # power.f.anova.mixed

pwrss.f.mixed.anova <- power.f.mixed.anova



pwrss.f.rmanova <- function(eta2 = 0.10, f2 = eta2 / (1 - eta2),
                             corr.rm = 0.50, n.levels = 2, n.rm = 2,
                             epsilon = 1, alpha = 0.05,
                             type = c("between", "within", "interaction"),
                             n = NULL, power = NULL, verbose = TRUE) {

  user.parms.names <- names(as.list(match.call()))
  if ("repmeasures.r" %in% user.parms.names) stop("'repmeasures.r' argument is obsolete, use 'corr.rm' instead", call. = FALSE)
  if ("n.measurements" %in% user.parms.names) stop("'n.measurements' argument is obsolete, use 'n.rm' instead", call. = FALSE)
  if (all(c("eta2", "f2") %in% user.parms.names)) warning("'eta2' and 'f2' cannot be specified at the same time \n Using 'f2' as the effect size")
  if ("f2" %in% user.parms.names) eta.squared <- f.to.etasq(f = sqrt(f2), verbose = FALSE)$eta.squared
  if ("eta2" %in% user.parms.names) eta.squared <- eta2

  type <- tolower(match.arg(type))

  mixed.anova.obj <- power.f.mixed.anova(eta.squared = eta.squared,
                                         n.total = n,
                                         factor.levels = c(n.levels, n.rm),
                                         factor.type = c("between", "within"),
                                         rho.within = corr.rm, epsilon = epsilon,
                                         power = power, alpha = alpha,
                                         effect = type,
                                         ceiling = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.f.mixed.anova() function. \n")

  return(invisible(mixed.anova.obj))

} # pwrss.f.rmanova()
