# power for the generic t test with (optional) type I and type II error plots
power.t.test <- function(ncp, null.ncp = 0,
                         df, alpha = 0.05,
                         alternative = c("two.sided", "one.sided", "two.one.sided"),
                         plot = TRUE, verbose = TRUE, pretty = FALSE) {

  check.positive(df)
  check.proportion(alpha)

  alternative <- tolower(match.arg(alternative))

  # calculate statistical power
  if (alternative == "two.sided") {

    ifelse(is.numeric(ncp) && length(ncp) == 1,
           valid.ncp <- TRUE,
           valid.ncp <- FALSE)

    ifelse(is.numeric(null.ncp) && length(null.ncp) == 1,
           valid.null.ncp <- TRUE,
           valid.null.ncp <- FALSE)

    if (isFALSE(valid.ncp) || isFALSE(valid.null.ncp)) stop("'ncp' or 'null.ncp' must be numeric and of length one for the two-sided test.", call. = FALSE)
    # if (ncp < null.ncp) stop("'ncp' must be equal or greater than 'null.ncp' for the two-sided test.", .call = FALSE)

    t.alpha.upper <- qt(alpha / 2, df = df, ncp = null.ncp, lower.tail = FALSE)
    t.alpha.lower <- qt(alpha / 2, df = df, ncp = null.ncp, lower.tail = TRUE)
    t.alpha <- c(t.alpha.lower, t.alpha.upper)
    power <-  pt(t.alpha.lower, df = df, ncp = ncp, lower.tail = TRUE) +
      pt(t.alpha.upper, df = df, ncp = ncp, lower.tail = FALSE)

  } else if (alternative == "one.sided") {

    ifelse(is.numeric(ncp) || length(ncp) == 1,
           valid.ncp <- TRUE,
           valid.ncp <- FALSE)

    ifelse(is.numeric(null.ncp) || length(null.ncp) == 1,
           valid.null.ncp <- TRUE,
           valid.null.ncp <- FALSE)

    if (isFALSE(valid.ncp) || isFALSE(valid.null.ncp)) stop("'ncp' or 'null.ncp' must be numeric and of length one for the one-sided test.", call. = FALSE)
    # if (any(ncp < null.ncp) && alternative == "greater") stop("alternative = 'greater' but ncp < null.ncp.", call. = FALSE)
    # if (any(ncp > null.ncp) && alternative == "less") stop("alternative = 'less' but ncp > null.ncp.", call. = FALSE)

    ifelse(ncp > null.ncp,
           lower.tail <- FALSE,
           lower.tail <- TRUE)
    t.alpha <- qt(alpha, df = df, ncp = null.ncp, lower.tail = lower.tail) # if ncp > null.ncp
    power <- pt(t.alpha, df = df, ncp = ncp, lower.tail = lower.tail) # if ncp > null.ncp

  } else if (alternative == "two.one.sided") {

    ifelse(is.numeric(ncp) && length(ncp) == 1,
           valid.ncp <- TRUE,
           valid.ncp <- FALSE)

    ifelse(is.numeric(null.ncp) && length(null.ncp) %in% c(1, 2),
           valid.null.ncp <- TRUE,
           valid.null.ncp <- FALSE)

    if (isFALSE(valid.ncp)) stop("'ncp' must be numeric and of length one for equivalence tests.", call. = FALSE)
    if (isFALSE(valid.null.ncp)) stop("'null.ncp' must be numeric and of length one (absolute value) or length two (with lower and upper bounds) for the equivalence test.", call. = FALSE)

    if (length(null.ncp) == 1) null.ncp <- c(min(c(-null.ncp, null.ncp)), max(-null.ncp, null.ncp))

    # equivalence test
    if (ncp > min(null.ncp) && ncp < max(null.ncp)) {

      t.alpha.upper <- qt(alpha, df = df, ncp = min(null.ncp), lower.tail = FALSE)
      t.alpha.lower <- qt(alpha, df = df, ncp = max(null.ncp), lower.tail = TRUE)
      t.alpha <- c(t.alpha.upper, t.alpha.lower)

      power <- pt(t.alpha.lower, df = df, ncp = ncp, lower.tail = TRUE) +
        pt(t.alpha.upper, df = df, ncp = ncp, lower.tail = FALSE) - 1

      power[power < 0] <- 0

    }

    # minimum effect test
    if (ncp < min(null.ncp) || ncp > max(null.ncp)) {

      t.alpha.lower <- qt(alpha / 2, df = df, ncp = min(null.ncp), lower.tail = TRUE)
      t.alpha.upper <- qt(alpha / 2, df = df, ncp = max(null.ncp), lower.tail = FALSE)
      t.alpha <- c(t.alpha.lower, t.alpha.upper)

      power <- pt(t.alpha.lower, df = df, ncp = ncp, lower.tail = TRUE) +
        pt(t.alpha.upper, df = df, ncp = ncp, lower.tail = FALSE)

    }

  } else {

    stop("Not a valid hypothesis type.", call. = FALSE)

  }

  if (plot) {

    suppressWarnings({
      .plot.t.t1t2(ncp = ncp, null.ncp = null.ncp, df = df, alpha = alpha, alternative = alternative)
    }) # supressWarnings

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

    print.obj <- list(test = "Generic T-Test",
                      requested = "power",
                      alt = alternative,
                      ncp.alternative = ncp,
                      ncp.null = null.ncp,
                      t.alpha = t.alpha,
                      df = df, alpha = alpha,
                      power = power)

    if (pretty) {
      .print.pwrss.t(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.t(print.obj, verbose = verbose)
    }

  } # verbose


  return(invisible(list(alternative = alternative, ncp = ncp, null.ncp = null.ncp,
                        df = df, alpha = alpha, t.alpha = t.alpha, power = power)))

} # end of power.t.test()

power.t <- power.t.test
