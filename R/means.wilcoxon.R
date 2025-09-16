#############################################################
#################### non-parametric tests ###################
#### one sample, paired samples, and independent samples ####
#############################################################

power.np.wilcoxon <- function(d, null.d = 0, margin = 0,
                              n2 = NULL, n.ratio = 1, power = NULL, alpha = 0.05,
                              alternative = c("two.sided", "one.sided", "two.one.sided"),
                              design = c("independent", "paired", "one.sample"),
                              distribution = c("normal", "uniform", "double.exponential",
                                               "laplace", "logistic"),
                              method = c("guenther", "noether"),
                              ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  check.logical(ceiling)
  check.proportion(alpha)
  check.positive(n.ratio)
  check.numeric(d, null.d)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n2)) check.sample.size(n2)

  alternative <- tolower(match.arg(alternative))
  distribution <- tolower(match.arg(distribution))
  method <- tolower(match.arg(method))
  design <- tolower(match.arg(design))

  if (is.null(n2) && is.null(power)) stop("`n2` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n2) && !is.null(power)) stop("Exactly one of the `n2` or `power` should be `NULL`.", call. = FALSE)

  if (alternative == "two.one.sided") {
    if (isFALSE(all(is.numeric(margin))) || any(margin < -10) || any(margin > 10)) stop("Possibly incorrect value for `margin`.", call. = FALSE)
    if (length(margin) != 2) stop("Provide null margins in the form of margin = c(lower, upper).", call. = FALSE)
    # if (prob - null.prob < margin[1] ||  prob - null.prob > margin[2]) stop("`prob` should be between lower and upper null margins)", call. = FALSE)
  } else {
    if (isFALSE(all(is.numeric(margin))) || length(margin) != 1 || any(margin < -10) || any(margin > 10)) stop("Possibly incorrect value for `margin`.", call. = FALSE)
  }

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  ifelse(design == "independent",
         paired <- FALSE,
         paired <- TRUE)

  if (method == "noether" && isTRUE(paired)) stop("Specify method = 'guenther' to request Wilcoxon signed-rank test for matched pairs.", call. = FALSE)

  pwr.wilcox <- function(d, null.d, margin,
                         n2, n.ratio, alpha,
                         paired, alternative,
                         method) {

    n1 <- n2 * n.ratio

    if (method == "noether") {

      propss <- n.ratio / (n.ratio + 1)

      prob <- d.to.cles(d = d, design = design, verbose = FALSE)$cles
      null.prob <- d.to.cles(d = null.d, design = design, verbose = FALSE)$cles

      if (alternative == "two.one.sided") {

        ignorable.prob.lower <- d.to.cles(d = min(margin) + null.d, design = design, verbose = FALSE)$cles
        margin.prob.lower <- ignorable.prob.lower - null.prob

        ignorable.prob.upper <- d.to.cles(d = max(margin) + null.d, design = design, verbose = FALSE)$cles
        margin.prob.upper <- ignorable.prob.upper - null.prob

        margin.prob <-  c(margin.prob.lower,  margin.prob.upper)

      } else {

        ignorable.prob <- d.to.cles(d = margin + null.d, design = design, verbose = FALSE)$cles
        margin.prob <- ignorable.prob - null.prob
      }

      lambda <- sqrt(n2 + n2 * n.ratio) * sqrt(12 * propss * (1 - propss)) * (prob - null.prob)
      null.lambda <- sqrt(n2 + n2 * n.ratio) * sqrt(12 * propss * (1 - propss)) * (margin.prob)

      pwr.obj <- power.z.test(mean = lambda, null.mean = null.lambda,
                              alpha = alpha, alternative = alternative,
                              plot = FALSE, verbose = FALSE)

      list(power = pwr.obj$power,
           z.alpha = pwr.obj$z.alpha,
           mean = pwr.obj$mean,
           null.mean = pwr.obj$null.mean)

    } else if (method == "guenther") {

      if (paired) {
        df <- n2 - 1
        lambda <- (d - null.d) / sqrt(1 / n2)
        null.lambda <- (margin) / sqrt(1 / n2)
      } else {
        df <- n1 + n2 - 2
        lambda <- (d - null.d) / sqrt(1 / n1 + 1 / n2)
        null.lambda <- (margin) / sqrt(1 / n1 + 1 / n2)
      }

      pwr.obj <- power.t.test(ncp = lambda, null.ncp = null.lambda, df = df,
                              alpha = alpha, alternative = alternative,
                              plot = FALSE, verbose = FALSE)

      list(power = pwr.obj$power,
           t.alpha = pwr.obj$t.alpha,
           ncp = pwr.obj$ncp,
           null.ncp = pwr.obj$null.ncp,
           df = pwr.obj$df)

    }

  } # pwr

  # wilcoxon adjustment for guenther method
  w <- switch(distribution,
              `uniform` = 1,
              `double.exponential` = 2 / 3,
              `laplace` = 2 / 3,
              `logistic` = 9 / pi ^ 2,
              `normal` = pi / 3)

  # ifelse(method == "noether", w.adj <- 1, w.adj <- w)

  # get power or sample size
  if (is.null(power)) {

    # apply wilcoxon adjustment
    n1 <- n2 * n.ratio
    n1 <- n1 / w
    n2 <- n2 / w

    pwr.obj <- pwr.wilcox(d = d, null.d = null.d, margin = margin,
                          n2 = n2, n.ratio = n.ratio, alpha = alpha,
                          paired = paired, alternative = alternative,
                          method = method)

    if (method == "guenther") {
      list.out <- list(power = pwr.obj$power,
                       t.alpha = pwr.obj$t.alpha,
                       ncp = pwr.obj$ncp,
                       null.ncp = pwr.obj$null.ncp,
                       df = pwr.obj$df)
    }

    if (method == "noether") {
      list.out <- list(power = pwr.obj$power,
                       z.alpha = pwr.obj$z.alpha,
                       mean = pwr.obj$mean,
                       sd = 1,
                       null.mean = pwr.obj$null.mean,
                       null.sd = 1,
                       df = Inf)
    }

    if (list.out$power < 0) stop("Design is not feasible.", call. = FALSE)

    # reverse wilcoxon adjustment
    n1.star <- n1 * w
    n2.star <- n2 * w

    if (ceiling) {
      n1.star <- ceiling(n1.star)
      n2.star <- ceiling(n2.star)
    }

  } else if (is.null(n2)) {

    if (method == "noether") {

      H1_H0.min <- 0.001
      lambda.max <- 4
      propss <- n.ratio / (n.ratio + 1)
      n.tot.max <- (lambda.max / (sqrt(12 * propss * (1 - propss)) * (H1_H0.min))) ^ 2
      n2.max <- n.tot.max / (1 + n.ratio)

    } else {

      H1_H0.min <- 0.001
      lambda.max <- 4
      if (paired) {
        n2.max <- (lambda.max / H1_H0.min) ^ 2
      } else {
        n2.max <- (1 + 1 / n.ratio) / (H1_H0.min / lambda.max) ^ 2
      }

    }

    # n2.max <- 1e+08
    # too big of a number throws warning in uniroot()
    # full precision may not have been achieved in 'pnt{final}'
    # because ncp is too large

    n2 <- try(silent = TRUE,
              suppressWarnings({
                uniroot(function(n2) {
                  power - pwr.wilcox(d = d, null.d = null.d, margin = margin,
                                     n2 = n2, n.ratio = n.ratio, alpha = alpha,
                                     paired = paired, alternative = alternative,
                                     method = method)$power
                }, interval = c(2, n2.max))$root
              }) # supressWarnings
    ) # try

    if (inherits(n2, "try-error") || n2 == 1e10) stop("Design is not feasible.", call. = FALSE)

    n1 <- n2 * n.ratio

    # reverse wilcoxon adjustment
    n1.star <- n1 * w
    n2.star <- n2 * w

    if (ceiling) {
      n1.star <- ceiling(n1.star)
      n2.star <- ceiling(n2.star)
    }

    pwr.obj <- pwr.wilcox(d = d, null.d = null.d, margin = margin,
                          n2 = n2.star / w, n.ratio = n.ratio, alpha = alpha,
                          paired = paired, alternative = alternative,
                          method = method)

    if (method == "guenther") {
      list.out <- list(power = pwr.obj$power,
                       t.alpha = pwr.obj$t.alpha,
                       ncp = pwr.obj$ncp,
                       null.ncp = pwr.obj$null.ncp,
                       df = pwr.obj$df)
    }

    if (method == "noether") {
      list.out <- list(power = pwr.obj$power,
                       z.alpha = pwr.obj$z.alpha,
                       mean = pwr.obj$mean,
                       sd = 1,
                       null.mean = pwr.obj$null.mean,
                       null.sd = 1,
                       df = Inf)
    }

    if (list.out$power < 0) stop("Design is not feasible", call. = FALSE)

  } # get power or sample size


  ifelse(design %in% c("paired", "one.sample"), n <- n2.star, n <- c(n1 = n1.star, n2 = n2.star))
  # ifelse(design %in% c("paired", "one.sample"), n <- n2, n <- c(n1 = n1, n2 = n2))

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
           test <- "Wilcoxon Rank-Sum Test (Independent Samples) \n(Wilcoxon-Mann-Whitney or Mann-Whitney U Test)",
           ifelse(design == "paired",
                  test <- "Wilcoxon Signed-Rank Test (Paired Samples)",
                  test <- "Wilcoxon Signed-Rank Test (One Sample)"))

    print.obj <- c(list(requested = requested, test = test,
                      design = design, method = method, dist = distribution,
                      d = d, null.d = null.d, margin = margin,
                      alpha = alpha, alt = alternative, n = n), list.out)

    if (pretty) {
      .print.pwrss.wilcoxon(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.wilcoxon(print.obj, verbose = verbose)
    }

  }

  invisible(structure(c(list(parms = list(d = d, null.d = null.d, margin = margin,
                                          n.ratio = n.ratio, alpha = alpha,
                                          alternative = alternative, design = design,
                                          distribution = distribution,
                                          method = method, ceiling = ceiling,
                                          verbose = verbose, pretty = pretty),
                             test = ifelse(method == "noether", "z", "t"),
                             n = n), list.out),
                      class = c("pwrss", "np", "wilcoxon", ifelse(method == "noether", "z", "t"))))

} # end of pwrss.np.wilcoxon()
power.np.wilcox <- power.np.wilcoxon


pwrss.np.2groups <- function(mu1 = 0.20, mu2 = 0,
                             sd1 = ifelse(paired, sqrt(1 / (2 * (1 - paired.r))), 1),
                             sd2 = sd1, margin = 0, alpha = 0.05, paired = FALSE,
                             paired.r = 0.50, kappa = 1, n2 = NULL, power = NULL,
                             alternative = c("not equal", "greater", "less",
                                             "non-inferior", "superior", "equivalent"),
                             distribution = c("normal", "uniform", "double.exponential",
                                              "laplace", "logistic"),
                             method = c("guenther", "noether"),
                             verbose = TRUE) {

  method <- tolower(match.arg(method))
  alternative <- tolower(match.arg(alternative))
  distribution <- tolower(match.arg(distribution))

  null.d <- 0
  d <- means.to.d(mu1 = mu1, mu2 = mu2,
                  sd1 = sd1, sd2 = sd2,
                  n2 = 1e4, n.ratio = kappa,
                  paired = paired,
                  rho.paired = paired.r,
                  verbose = FALSE)$d

  if (alternative == "equivalent") {
    margin <- c(min(-margin, margin), max(-margin, margin))
    margin.lower <- means.to.d(mu1 = margin[1], mu2 = 0,
                               sd1 = sd1, sd2 = sd2,
                               n2 = 1e10, n.ratio = kappa,
                               paired = paired,
                               rho.paired = paired.r,
                               verbose = FALSE)$d
    margin.upper <- means.to.d(mu1 = margin[2], mu2 = 0,
                               sd1 = sd1, sd2 = sd2,
                               n2 = 1e10, n.ratio = kappa,
                               paired = paired,
                               rho.paired = paired.r,
                               verbose = FALSE)$d
    margin <- c(margin.lower, margin.upper)
  } else {
    margin <- means.to.d(mu1 = margin, mu2 = 0,
                         sd1 = sd1, sd2 = sd2,
                         n2 = 1e10, n.ratio = kappa,
                         paired = paired,
                         rho.paired = paired.r,
                         verbose = FALSE)$d
  }


  if (alternative %in% c("less", "greater", "non-inferior", "superior")) alternative <- "one.sided"
  if (alternative == "equivalent") alternative <- "two.one.sided"
  if (alternative == "not equal") alternative <- "two.sided"

  ifelse(paired, design <- "paired", design <- "independent")

  wilcox.obj <- power.np.wilcoxon(d = d,
                                  null.d = null.d,
                                  margin = margin,
                                  alpha = alpha,
                                  n2 = n2, n.ratio = kappa,
                                  power = power,
                                  design = design,
                                  alternative = alternative,
                                  distribution = distribution,
                                  method = method,
                                  ceiling = TRUE,
                                  verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.np.wilcox() function. \n")

  return(invisible(wilcox.obj))

} # end of pwrss.np.2groups()

pwrss.np.2means <- function(...) stop("This function is no longer available. Please use `power.np.wilcox()`.", call. = FALSE)
