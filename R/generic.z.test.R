# power for the generic t test with (optional) type I and type II error plots
power.z.test <- function(mean = NULL, sd = 1, null.mean = 0, null.sd = 1,
                         alpha = 0.05, alternative = c("two.sided", "one.sided", "two.one.sided"),
                         plot = TRUE, verbose = TRUE, pretty = FALSE, ...) {

  old.args <- list(...) # just arguments in ...
  names.old.args <- names(old.args)

  if ("ncp" %in% names.old.args) {
    if (is.null(mean)) {
      check.numeric(old.args$ncp)
      mean <- old.args$ncp
    } else {
      stop("Both the new argument `mean` and the deprecated argument `ncp` were provided. Please specify only one.", call. = FALSE)
    }
  }

  check.proportion(alpha)
  check.positive(sd)

  alternative <- tolower(match.arg(alternative))

  # calculate statistical power
  if (alternative == "two.sided") {

    ifelse(is.numeric(mean) && length(mean) == 1,
           valid.mean <- TRUE,
           valid.mean <- FALSE)

    ifelse(is.numeric(null.mean) && length(null.mean) == 1,
           valid.null.mean <- TRUE,
           valid.null.mean <- FALSE)

    if (isFALSE(valid.mean) || isFALSE(valid.null.mean)) stop("'mean' or 'null.mean' must be numeric and of length one for the two-sided test.", call. = FALSE)

    # if (mean < null.mean) stop("'mean' must be equal or greater than 'null.mean' for the two-sided test.", .call = FALSE)

    z.alpha.upper <- qnorm(alpha / 2, mean = null.mean, sd = null.sd, lower.tail = FALSE)
    z.alpha.lower <- qnorm(alpha / 2, mean = null.mean, sd = null.sd, lower.tail = TRUE)
    z.alpha <- c(z.alpha.lower, z.alpha.upper)

    power <-  pnorm(z.alpha.lower, mean = mean, sd = sd, lower.tail = TRUE) +
      pnorm(z.alpha.upper, mean = mean, sd = sd, lower.tail = FALSE)


  } else if (alternative == "one.sided") {

    ifelse(is.numeric(mean) && length(mean) == 1,
           valid.mean <- TRUE,
           valid.mean <- FALSE)

    ifelse(is.numeric(null.mean) && length(null.mean) == 1,
           valid.null.mean <- TRUE,
           valid.null.mean <- FALSE)

    if (isFALSE(valid.mean) || isFALSE(valid.null.mean)) stop("'mean' or 'null.mean' must be numeric and of length one for the one-sided test.", call. = FALSE)
    # if (any(mean < null.mean) && alternative == "greater") stop("alternative = 'greater' but mean < null.mean.", call. = FALSE)
    # if (any(mean > null.mean) && alternative == "less") stop("alternative = 'less' but mean > null.mean.", call. = FALSE)

    ifelse(mean > null.mean,
           lower.tail <- FALSE,
           lower.tail <- TRUE)
    z.alpha <- qnorm(alpha, mean = null.mean, sd = null.sd, lower.tail = lower.tail) # if mean > null.mean
    power <- pnorm(z.alpha, mean = mean, sd = sd, lower.tail = lower.tail) # if mean > null.mean

  } else if (alternative == "two.one.sided") {

    ifelse(is.numeric(mean) && length(mean) == 1,
           valid.mean <- TRUE,
           valid.mean <- FALSE)

    ifelse(is.numeric(null.mean) && length(null.mean) %in% c(1, 2),
           valid.null.mean <- TRUE,
           valid.null.mean <- FALSE)

    if (isFALSE(valid.mean)) stop("'mean' must be numeric and of length one for equivalence tests.", call. = FALSE)
    if (isFALSE(valid.null.mean)) stop("'null.mean' must be numeric and of length one (absolute value) or length two (with lower and upper bounds) for the equivalence test.", call. = FALSE)

    if (length(null.mean) == 1) null.mean <- c(min(c(-null.mean, null.mean)), max(-null.mean, null.mean))

    # equivalence test
    if (mean > min(null.mean) && mean < max(null.mean)) {

      z.alpha.upper <- qnorm(alpha, mean = min(null.mean), sd = null.sd, lower.tail = FALSE)
      z.alpha.lower <- qnorm(alpha, mean = max(null.mean), sd = null.sd, lower.tail = TRUE)
      z.alpha <- c(z.alpha.upper, z.alpha.lower)

      power <- pnorm(z.alpha.lower, mean = mean, sd = sd, lower.tail = TRUE) +
        pnorm(z.alpha.upper, mean = mean, sd = sd, lower.tail = FALSE) - 1

      power[power < 0] <- 0

    }

    # minimum effect test
    if (mean < min(null.mean) || mean > max(null.mean)) {

      z.alpha.lower <- qnorm(alpha / 2, mean = min(null.mean), sd = null.sd, lower.tail = TRUE)
      z.alpha.upper <- qnorm(alpha / 2, mean = max(null.mean), sd = null.sd, lower.tail = FALSE)
      z.alpha <- c(z.alpha.lower, z.alpha.upper)

      power <- pnorm(z.alpha.lower, mean = mean, sd = sd, lower.tail = TRUE) +
        pnorm(z.alpha.upper, mean = mean, sd = sd, lower.tail = FALSE)

    }

  } else {

    stop("Not a valid hypothesis type.", call. = FALSE)

  }

  if (plot) {
    # if (sd != 1 || null.sd != 1) stop("Plotting is currently not available when standard deviation of the standard normal distribution deviates from one.", call. = FALSE)

    try(silent = TRUE,
        suppressWarnings({
          .plot.t.t1t2(ncp = mean, null.ncp = null.mean, df = Inf, alpha = alpha, alternative = alternative)
        }) # supressWarnings
    ) # try

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

    print.obj <- list(test = "Generic Z-Test",
                      requested = "power",
                      alt = alternative,
                      mean.alternative = mean,
                      sd.alternative = sd,
                      mean.null = null.mean,
                      sd.null = null.sd,
                      alpha = alpha,
                      z.alpha = z.alpha,
                      power = power)

    if (pretty) {
      .print.pwrss.z(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.z(print.obj, verbose = verbose)
    }

  } # verbose


  return(invisible(list(alternative = alternative,
                        mean = mean,
                        sd = sd,
                        null.mean = null.mean,
                        null.sd = null.sd,
                        alpha = alpha,
                        z.alpha = z.alpha,
                        power = power)))

} # end of power.z.test()

power.z <- power.z.test
