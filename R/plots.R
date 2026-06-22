#################
# plot function #
#################

#' @exportS3Method
plot.pwrss <- function(x, ...) {

  if ("defunct" %in% class(x))
    stop("Plotting is no longer available for this type of object.", call. = FALSE)

  if (all(c("pwrss", "generic") %in% class(x)))
    stop("Use plot = TRUE argument for generic tests.", call. = FALSE)

  if (all(c("pwrss", "t") %in% class(x))) {

    # student, welch, wilcoxon, regression
    power.t.test(ncp = x$ncp,
                 null.ncp = x$null.ncp,
                 df = x$df,
                 alpha = x$parms$alpha,
                 alternative = x$parms$alternative,
                 verbose = 0)

  } else if (all(c("pwrss", "z") %in% class(x))) {

    # proportions, correlations, logistic, poisson, mediation
    power.z.test(mean = x$mean,
                 sd = x$sd,
                 null.mean = x$null.mean,
                 null.sd = x$null.sd,
                 alpha = x$parms$alpha,
                 alternative = x$parms$alternative,
                 verbose = 0)

  } else if (all(c("pwrss", "exact") %in% class(x))) {

    if (any(c("mcnemar", "fisher", "onecor") %in% class(x)))
      stop("Plotting is not available for exact tests.", call. = FALSE)

    # proportions.onetwo (only exact.oneprop)
    power.binom.test(size = ceiling(x$size),
                     prob = x$prob,
                     null.prob = x$null.prob,
                     alpha = x$parms$alpha,
                     alternative = x$parms$alternative,
                     verbose = 0)

  } else if (all(c("pwrss", "f") %in% class(x))) {

    # ancova, keppel, shieh, mixed.anova, regression
    power.f.test(ncp = x$ncp,
                 null.ncp = x$null.ncp,
                 df1 = x$df1,
                 df2 = x$df2,
                 alpha = x$parms$alpha,
                 verbose = 0)

  } else if (all(c("pwrss", "chisq") %in% class(x))) {

    # goodness-of-fit
    power.chisq.test(ncp = x$ncp,
                     null.ncp = x$null.ncp,
                     df = x$df,
                     alpha = x$parms$alpha,
                     verbose = 0)

  } else {

    stop("Not an object of the type 'pwrss'.", call. = FALSE)

  }

} # plot.pwrss
