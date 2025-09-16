# one can specify Cohen's d or Hedges' g
# because small sample correction applies to the effect size not its standard error
power.t.student <- function(d, null.d = 0, margin = 0,
                            n2 = NULL, n.ratio = 1, power = NULL, alpha = 0.05,
                            alternative = c("two.sided", "one.sided", "two.one.sided"),
                            design = c("independent", "paired", "one.sample"),
                            claim.basis = c("md.pval", "smd.ci"),
                            ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  check.numeric(d, null.d)
  check.logical(ceiling)
  check.proportion(alpha)
  check.positive(n.ratio)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n2)) check.sample.size(n2)

  alternative <- tolower(match.arg(alternative))
  design <- tolower(match.arg(design))
  claim.basis <- tolower(match.arg(claim.basis))

  if (is.null(n2) && is.null(power)) stop("`n2` and `power` cannot be `NULL` at the same time", call. = FALSE)
  if (!is.null(n2) && !is.null(power)) stop("Exactly one of the `n2` or `power` should be `NULL`", call. = FALSE)

  if (alternative == "two.one.sided") {
    if (length(margin) == 1) margin <- c(min(c(-margin, margin)), max(c(-margin, margin)))
    if (length(margin) > 2) stop("Provide margins in the form of margin = c(lower, upper)", call. = FALSE)
  } else {
    if (isFALSE(all(is.numeric(margin))) || length(margin) != 1) stop("Incorrect value for `margin`", call. = FALSE)
    # if (alternative %in% c("two.sided", "one.sided") && d == null.d) stop("alternative = 'two.sided' or 'one.sided' but d = null.d", call. = FALSE)
  }

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  pwr.student <- function(d, null.d, margin, n2, n.ratio,
                          alpha, alternative,
                          design, claim.basis) {
    n1 <- n2 * n.ratio

    if (design == "independent") {
      df <- n1 + n2 - 2
      if (claim.basis == "md.pval") {
        var.d <- (n1 + n2) / (n1 * n2)
      } else {
        var.d <- (n1 + n2) / (n1 * n2) + d ^ 2 / (2 * (df))
      }
    } # if independent

    if (design %in% c("paired", "one.sample")) {
      df <- n2 - 1
      if (claim.basis == "md.pval") {
        var.d <- 1 / n2
      } else {
        var.d <- (1 / n2 + d ^ 2 / (2 * n))
      }
    } # if paired or one.sample

    se.d <- sqrt(var.d)
    lambda <- (d - null.d) / sqrt(var.d)
    null.lambda <- margin / sqrt(var.d)
    # lambda <- (d - null.d - margin) / sqrt(var.d)
    # null.lambda <- 0

    pwr.obj <- power.t.test(ncp = lambda, null.ncp = null.lambda, df = df,
                            alpha = alpha, alternative = alternative,
                            plot = FALSE, verbose = FALSE)

    list(power = pwr.obj$power,
         t.alpha = pwr.obj$t.alpha,
         ncp = pwr.obj$ncp,
         null.ncp = pwr.obj$null.ncp,
         se.d = se.d,
         df = pwr.obj$df)

  } # pwr.student()

  ss.student <- function(d, null.d, margin, power, n.ratio,
                         alpha, alternative,
                         design, claim.basis) {

    n2 <- try(silent = TRUE,
              suppressWarnings({
                uniroot(function(n2) {
                  power - pwr.student(d = d, null.d = null.d, margin = margin,
                                      n2 = n2, n.ratio = n.ratio,
                                      alpha = alpha, alternative = alternative,
                                      design = design, claim.basis = claim.basis)$power
                }, interval = c(3, 1e10))$root
              }) # supressWarnings
    ) # try

    if (inherits(n2, "try-error") || n2 == 1e10) stop("Design is not feasible.", call. = FALSE)

    n2

  } # ss.student()

  if (is.null(power)) {
    pwr.obj <- pwr.student(d = d, null.d = null.d, margin = margin,
                           n2 = n2, n.ratio = n.ratio,
                           alpha = alpha, alternative = alternative,
                           design = design, claim.basis = claim.basis)
    n1 <- n.ratio * n2
    power <- pwr.obj$power
    ncp <- pwr.obj$ncp
    null.ncp <- pwr.obj$null.ncp
    t.alpha <- pwr.obj$t.alpha
    df <- pwr.obj$df

  } # power is null

  if (is.null(n2)) {
    n2 <- ss.student(d = d, null.d = null.d, margin = margin,
                     power = power, n.ratio = n.ratio,
                     alpha = alpha, alternative = alternative,
                     design = design, claim.basis = claim.basis)
    n1 <- n.ratio * n2

    if (ceiling) {
      n1 <- ceiling(n1)
      n2 <- ceiling(n2)
    }

    pwr.obj <- pwr.student(d = d, null.d = null.d, margin = margin,
                           n2 = n2, n.ratio = n.ratio,
                           alpha = alpha, alternative = alternative,
                           design = design, claim.basis = claim.basis)
    power <- pwr.obj$power
    t.alpha <- pwr.obj$t.alpha
    ncp <- pwr.obj$ncp
    null.ncp <- pwr.obj$null.ncp
    df <- pwr.obj$df

  } # n2 is null

  ifelse(design %in% c("paired", "one.sample"), n <- n2, n <- c(n1 = n1, n2 = n2))

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

    ifelse(design == "independent",
           test <- "Student's T-Test (Independent Samples)",
           ifelse(design == "paired",
                  test <- "Student's T-Test (Paired Samples)",
                  test <- "Student's T-Test (One Sample)"))

    print.obj <- list(requested = requested, test = test,
                      d = d, null.d = null.d, margin = margin,
                      alpha = alpha, t.alpha = t.alpha,
                      alt = alternative, n = n, df = df,
                      ncp.alternative = ncp,
                      ncp.null = null.ncp,
                      power = power)

    if (pretty) {
      .print.pwrss.student(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.student(print.obj, verbose = verbose)
    }

  }

  invisible(structure(list(parms = list(d = d, null.d = null.d, margin = margin,
                                        n2 = n2, n.ratio = n.ratio,
                                        alpha = alpha, alternative = alternative,
                                        design = design, claim.basis = claim.basis,
                                        ceiling = ceiling, verbose = verbose),
                           test = "t",
                           df = df,
                           ncp = ncp,
                           null.ncp = null.ncp,
                           t.alpha = t.alpha,
                           power = power,
                           n = n),
                      class = c("pwrss", "t", "student")))

} # power.t.student()

# provide a table of conversions from d to d.hc2
# mention that it cannot be interpreted same as d
# when both var.ratio and n.ratio deviate from 1
# minor deviations are OK
power.t.welch <- function(d, null.d = 0, margin = 0,
                          var.ratio = 1, n.ratio = 1, n2 = NULL,
                          power = NULL, alpha = 0.05,
                          alternative = c("two.sided", "one.sided", "two.one.sided"),
                          claim.basis = c("md.pval", "smd.ci"),
                          ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  # variance ratio constraint
  vrc <- function(var.ratio, n2, n.ratio) {
    sd2  <- sqrt((n2 * (n.ratio + 1) - 2) / (n2 - 1 + var.ratio * (n2 * n.ratio - 1)))
    sd1 <- sqrt(var.ratio * sd2 ^ 2)
    list(sd1 = sd1, sd2 = sd2)
  } # vrc()

  welch.df <- function(sd1, sd2, n1, n2) {
    (sd1 ^ 2 / n1 + sd2 ^ 2 / n2) ^ 2 /
      (sd1^4 / (n1 ^ 2 * (n1 - 1)) + sd2^4 / (n2 ^ 2 * (n2 - 1)))
  } # welch.df()

  check.numeric(d, null.d)
  check.logical(ceiling)
  check.proportion(alpha)
  check.positive(n.ratio)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n2)) check.sample.size(n2)

  alternative <- tolower(match.arg(alternative))
  claim.basis <- tolower(match.arg(claim.basis))

  if (is.null(n2) && is.null(power)) stop("`n2` and `power` cannot be `NULL` at the same time", call. = FALSE)
  if (!is.null(n2) && !is.null(power)) stop("Exactly one of the `n2` or `power` should be `NULL`", call. = FALSE)

  if (alternative == "two.one.sided") {
    if (length(margin) == 1) margin <- c(min(c(-margin, margin)), max(c(-margin, margin)))
    if (length(margin) > 2) stop("Provide margins in the form of margin = c(lower, upper)", call. = FALSE)
  } else {
    if (isFALSE(all(is.numeric(margin))) || length(margin) != 1) stop("Incorrect value for `null.d`", call. = FALSE)
    # if (alternative %in% c("two.sided", "one.sided") && d == null.d) stop("alternative = 'two.sided' or 'one.sided' but d = null.d", call. = FALSE)
  }

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  pwr.welch <- function(d, null.d, margin, var.ratio,
                        n2, n.ratio, alpha, alternative,
                        claim.basis) {
    n1 <- n.ratio * n2

    sds <- vrc(var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio)
    sd1 <- sds$sd1
    sd2 <- sds$sd2

    gamma1 <- 1 / n1
    gamma2 <- 1 / n2

    df <- welch.df(sd1 = sd1, sd2 = sd2, n1 = n1, n2 = n2)

    var.delta <- gamma1 * sd1 ^ 2 + gamma2 * sd2 ^ 2
    std <- sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
    var.std <- ((n1 - 1) * sd1^4 + (n2 - 1) * sd2^4) / ((2 * (n1 - 1) * sd1 ^ 2 + 2 * (n2 - 1) * sd2 ^ 2) * (n1 + n2 - 2))

    if (claim.basis == "md.pval") {
      var.d.hc2 <- var.delta
    } else {
      var.d.hc2 <- var.delta / std ^ 2 + d ^ 2 * var.std / std ^ 2
    }

    se.d <- sqrt(var.d.hc2)
    lambda <- (d - null.d) / sqrt(var.d.hc2)
    null.lambda <- margin / sqrt(var.d.hc2)

    pwr.obj <- power.t.test(ncp = lambda, null.ncp = null.lambda, df = df,
                            alpha = alpha, alternative = alternative,
                            plot = FALSE, verbose = FALSE)

    list(power = pwr.obj$power,
         t.alpha = pwr.obj$t.alpha,
         ncp = pwr.obj$ncp,
         null.ncp = pwr.obj$null.ncp,
         se.d = se.d,
         df = pwr.obj$df)

  } # pwr.welch()

  ss.welch <- function(d, null.d, margin, var.ratio,
                       power, n.ratio, alpha, alternative,
                       claim.basis) {

    n2 <- try(silent = TRUE,
             suppressWarnings({
               uniroot(function(n2) {
                 power - pwr.welch(d = d, null.d = null.d, margin = margin,
                                   var.ratio = var.ratio,
                                   n2 = n2, n.ratio = n.ratio,
                                   alpha = alpha, alternative = alternative,
                                   claim.basis = claim.basis)$power
               }, interval = c(3, 1e10))$root
             }) # supressWarnings
    ) # try

    if (inherits(n2, "try-error") || n2 == 1e10) stop("Design is not feasible.", call. = FALSE)

    n2

  } # ss.student()


  if (is.null(power)) {
    pwr.obj <- pwr.welch(d = d, null.d = null.d, margin = margin,
                         var.ratio = var.ratio,
                         n2 = n2, n.ratio = n.ratio,
                         alpha = alpha, alternative = alternative,
                         claim.basis = claim.basis)
    n1 <- n.ratio * n2
    power <- pwr.obj$power
    ncp <- pwr.obj$ncp
    null.ncp <- pwr.obj$null.ncp
    t.alpha <- pwr.obj$t.alpha
    se.d <- pwr.obj$se.d
    df <- pwr.obj$df
  } # power is null

  if (is.null(n2)) {
    n2 <- ss.welch(d = d, null.d = null.d, margin = margin,
                   var.ratio = var.ratio,
                   power = power, n.ratio = n.ratio,
                   alpha = alpha, alternative = alternative,
                   claim.basis = claim.basis)
    n1 <- n.ratio * n2

    if (ceiling) {
      n1 <- ceiling(n1)
      n2 <- ceiling(n2)
    }

    pwr.obj <- pwr.welch(d = d, null.d = null.d, margin = margin,
                         var.ratio = var.ratio,
                         n2 = n2, n.ratio = n.ratio,
                         alpha = alpha, alternative = alternative,
                         claim.basis = claim.basis)
    power <- pwr.obj$power
    t.alpha <- pwr.obj$t.alpha
    ncp <- pwr.obj$ncp
    null.ncp <- pwr.obj$null.ncp
    se.d <- pwr.obj$se.d
    df <- pwr.obj$df
  } # n2 is null

  n <- c(n1 = n1, n2 = n2)

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

    test <- "Welch's T-Test (Independent Samples)"

    print.obj <- list(requested = requested, test = test,
                      d = d, se.d = se.d, null.d = null.d, margin = margin,
                      alpha = alpha, t.alpha = t.alpha,
                      alt = alternative, n = n, df = df,
                      ncp.alternative = ncp,
                      ncp.null = null.ncp,
                      power = power)

    if (pretty) {
      .print.pwrss.student(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.student(print.obj, verbose = verbose)
    }

  }

  invisible(structure(list(parms = list(d = d, null.d = null.d,
                                        var.ratio = var.ratio,
                                        n2 = n2, n.ratio = n.ratio,
                                        alpha = alpha, alternative = alternative,
                                        claim.basis = claim.basis,
                                        ceiling = ceiling, verbose = verbose),
                           test = "t",
                           df = df,
                           ncp = ncp,
                           null.ncp = null.ncp,
                           t.alpha = t.alpha,
                           power = power,
                           n = n),
                      class = c("pwrss", "t", "welch")))

} # power.t.welch()




###################
# one mean t test #
###################

pwrss.t.mean <- function(mu, sd = 1, mu0 = 0, margin = 0, alpha = 0.05,
                         alternative = c("not equal", "greater", "less",
                                         "equivalent", "non-inferior", "superior"),
                         n = NULL, power = NULL, verbose = TRUE) {

  check.positive(sd)
  check.numeric(mu, mu0, margin)
  check.logical(verbose)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n)) check.sample.size(n)

  alternative <- tolower(match.arg(alternative))

  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"

  d <- (mu - mu0) / sd

  if (alternative == "two.one.sided") {

    margin.d.lower <- - margin / sd
    margin.d.upper <- + margin / sd
    margin.d <- c(margin.d.lower, margin.d.upper)

  } else if (alternative == "two.sided") {

    if (margin != 0) warning("Margin is forced to be 0 for the 'two.sided' test.", call. = FALSE)
    margin.d <- 0

  } else {

    margin.d <- margin / sd

  }

  student.obj <- power.t.student(d = d, margin = margin.d,
                                 n2 = n, power = power, alpha = alpha,
                                 alternative = alternative,
                                 design = "one.sample",
                                 claim.basis = "md.pval",
                                 ceiling = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.t.student() function. \n")

  return(invisible(student.obj))

} # pwrss.t.mean()


####################
# two means t test #
####################

pwrss.t.2means <- function(mu1, mu2 = 0, margin = 0,
                            sd1 = ifelse(paired, sqrt(1 / (2 * (1 - paired.r))), 1),
                            sd2 = sd1, kappa = 1, paired = FALSE, paired.r = 0.50,
                            alpha = 0.05, welch.df = TRUE,
                            alternative = c("not equal", "greater", "less",
                                            "equivalent", "non-inferior", "superior"),
                            n2 = NULL, power = NULL, verbose = TRUE) {

  if (isFALSE(welch.df)) warning("Forcing welch.df = TRUE.", call. = FALSE)

  check.positive(sd1, sd2)
  check.correlation(paired.r)
  check.numeric(mu1, mu2, margin)
  check.logical(paired, verbose)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n2)) check.sample.size(n2)

  alternative <- tolower(match.arg(alternative))

  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"

  calculate.d <- function(mu1, mu2, sd1, sd2, n2,
                     n.ratio, paired, rho.paired) {
    if (paired) {
      psd <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * sd1 * sd2 * rho.paired)
    } else {
      n1 <- n.ratio * n2
      psd <- sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
    }

    (mu1 - mu2) / psd
  }

  # population d
  d <- calculate.d(mu1 = mu1, mu2 = mu2,
                   sd1 = sd1, sd2 = sd2,
                   n2 = 1e10, n.ratio = kappa,
                   paired = paired, rho.paired = paired.r)

  if (alternative == "two.one.sided") {

    margin.d.lower <- calculate.d(mu1 = -margin, mu2 = 0,
                             sd1 = sd1, sd2 = sd2,
                             n2 = 1e5, n.ratio = kappa,
                             paired = paired, rho.paired = paired.r)
    margin.d.upper <- calculate.d(mu1 = margin, mu2 = 0,
                             sd1 = sd1, sd2 = sd2,
                             n2 = 1e5, n.ratio = kappa,
                             paired = paired, rho.paired = paired.r)
    margin.d <- c(margin.d.lower, margin.d.upper)

  } else if (alternative == "two.sided") {

    if (margin != 0) warning("Margin is forced to be 0 for the 'two.sided' test.", call. = FALSE)
    margin.d <- 0

  } else {

    margin.d <- calculate.d(mu1 = margin, mu2 = 0,
                            sd1 = sd1, sd2 = sd2,
                            n2 = 1e5, n.ratio = kappa,
                            paired = paired, rho.paired = paired.r)

  }

  if (paired) {
    t.obj <- power.t.student(d = d, margin = margin.d,
                             n2 = n2, n.ratio = kappa,
                             power = power, alpha = alpha,
                             alternative = alternative,
                             design = "paired",
                             claim.basis = "md.pval",
                             ceiling = TRUE, verbose = verbose)
  } else {
    t.obj <- power.t.welch(d = d, margin = margin.d,
                           var.ratio = sd1 ^ 2 / sd2 ^ 2,
                           n2 = n2, n.ratio = kappa,
                           power = power, alpha = alpha,
                           alternative = alternative,
                           claim.basis = "md.pval",
                           ceiling = TRUE, verbose = verbose)
  }

  # cat("This function will be removed in the future. \n Please use power.t.student() or power.t.welch() function. \n")

  return(invisible(t.obj))

} # pwrss.t.2means()
