#################
# plot function #
#################

#' @exportS3Method
plot.pwrss <- function(x, ...) {

  if (inherits(x, "pwrss")) {

    if (inherits(x, "defunct"))
      stop("Plotting is no longer available for this type of object.", call. = FALSE)

    if (inherits(x, "generic"))
      stop("Use plot = TRUE argument for generic tests.", call. = FALSE)

    if (inherits(x, "t")) {

      # student, welch, wilcoxon, regression
      power.t.test(ncp = x$ncp,
                   null.ncp = x$null.ncp,
                   df = x$df,
                   alpha = x$parms$alpha,
                   alternative = x$parms$alternative,
                   verbose = 0)

    } else if (inherits(x, "z")) {

      # proportions, correlations, logistic, poisson, mediation
      power.z.test(mean = x$mean,
                   sd = x$sd,
                   null.mean = x$null.mean,
                   null.sd = x$null.sd,
                   alpha = x$parms$alpha,
                   alternative = x$parms$alternative,
                   verbose = 0)

    } else if (inherits(x, "exact")) {

      if (inherits(x, c("mcnemar", "fisher", "onecor")))
        stop("Plotting is not available for exact tests.", call. = FALSE)

      # proportions.onetwo (only exact.oneprop)
      power.binom.test(size = ceiling(x$size),
                       prob = x$prob,
                       null.prob = x$null.prob,
                       alpha = x$parms$alpha,
                       alternative = x$parms$alternative,
                       verbose = 0)

    } else if (inherits(x, "f")) {

      # ancova, keppel, shieh, mixed.anova, regression
      power.f.test(ncp = x$ncp,
                   null.ncp = x$null.ncp,
                   df1 = x$df1,
                   df2 = x$df2,
                   alpha = x$parms$alpha,
                   verbose = 0)

    } else if (inherits(x, "chisq")) {

      # goodness-of-fit
      power.chisq.test(ncp = x$ncp,
                       null.ncp = x$null.ncp,
                       df = x$df,
                       alpha = x$parms$alpha,
                       verbose = 0)

    }

  } else {

    stop("Not an object of the type 'pwrss'.", call. = FALSE)

  }

} # plot.pwrss
