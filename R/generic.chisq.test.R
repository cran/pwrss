# power for the generic Chi-square test with (optional) type I and type II error plots
power.chisq.test <- function(ncp, null.ncp = 0, df, alpha = 0.05,
                             plot = TRUE, verbose = TRUE, pretty = FALSE) {

  check.sample.size(df)
  check.proportion(alpha)

  ifelse(is.numeric(ncp) || length(ncp) == 1 || ncp >= 0,
         valid.ncp <- TRUE,
         valid.ncp <- FALSE)

  ifelse(is.numeric(null.ncp) || length(null.ncp) == 1 || all(null.ncp >= 0),
         valid.null.ncp <- TRUE,
         valid.null.ncp <- FALSE)

  if (isFALSE(valid.ncp) || isFALSE(valid.null.ncp)) stop("'ncp' or 'null.ncp' must be numeric, positive, and of length one.", call. = FALSE)
  if (ncp < null.ncp) stop("'ncp' should be greater than or equal to 'null.ncp'.")

  chisq.alpha <- qchisq(alpha, df = df, ncp = null.ncp, lower.tail = FALSE)

  power <- pchisq(chisq.alpha, df = df, ncp = ncp, lower.tail = FALSE)

  if (plot) {

    .plot.chisq.t1t2(ncp = ncp, null.ncp = null.ncp, df = df, alpha = alpha)

  }

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

    print.obj <- list(test = "Generic Chi-square Test",
                      requested = "power",
                      power = power, ncp.alternative = ncp, ncp.null = null.ncp,
                      alpha = alpha, chisq.alpha = chisq.alpha, df = df)

    if (pretty) {
      .print.pwrss.chisq(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.chisq(print.obj, verbose = verbose)
    }

  } # end of verbose

  return(invisible(list(power = power, ncp = ncp, null.ncp = null.ncp,
                        alpha = alpha, df = df, chisq.crit = chisq.alpha)))

} # end of power.chisq.test()

power.chisq <- power.chisq.test
