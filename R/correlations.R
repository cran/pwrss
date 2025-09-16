#######################################################################
# Steiger's z-test for dependent correlations (Steiger, 1980, p. 247) #
#######################################################################

power.z.steiger <- function(rho12, rho13, rho23,
                            rho14 = NULL, rho24 = NULL, rho34 = NULL,
                            n = NULL, power = NULL, alpha = 0.05,
                            alternative = c("two.sided", "one.sided"),
                            pooled = TRUE, common.index = FALSE,
                            ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  user.parms <- as.list(match.call())
  names.user.parms <- names(user.parms)

  check.proportion(alpha)
  check.logical(pooled, common.index, ceiling)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n)) check.sample.size(n)

  alternative <- tolower(match.arg(alternative))

  if (is.null(n) && is.null(power)) stop("`n` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n) && !is.null(power)) stop("Exactly one of the `n` or `power` should be `NULL`.", call. = FALSE)

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  pwr.steiger <- function(rho1, rho2,
                          cov.null, cov.alt,
                          n, alpha, alternative) {

    z1 <- cor.to.z(rho1, FALSE)$z
    z2 <- cor.to.z(rho2, FALSE)$z

    sigma.null <- sqrt((2 - 2 * cov.null) / (n - 3))
    sigma.alt <- sqrt((2 - 2 * cov.alt) / (n - 3))

    lambda <- (z1 - z2) /  sigma.null
    sigma.lambda <- sigma.alt / sigma.null

    if (alternative == "two.sided") lambda <- abs(lambda)

    pwr.obj <- power.z.test(mean = lambda, sd = sigma.lambda,
                            null.mean = 0, null.sd = 1,
                            alpha = alpha,
                            alternative = alternative,
                            plot = FALSE, verbose = FALSE)

    pwr.obj

  } # pwr.steiger()

  ss.steiger <- function(rho1, rho2,
                         cov.null, cov.alt,
                         power, alpha, alternative) {

    n <- try(silent = TRUE,
             suppressWarnings({
               uniroot(function(n) {
                 power - pwr.steiger(rho1 = rho1, rho2 = rho2,
                                     cov.null = cov.null, cov.alt = cov.alt,
                                     n = n, alpha = alpha,
                                     alternative = alternative)$power
               }, interval = c(5, 1e+09))$root
             }) # supressWarnings
    ) # try

    if (inherits(n, "try-error") || n == 1e10) stop("Design is not feasible.", call. = FALSE)

    n

  } # ss.steiger()

  if (common.index) {

    if (any(c("rho14", "rho24", "rho34") %in% names.user.parms))
      warning("Ignoring `rho14` `rho24`, or `rho34` because common.index = TRUE.", call. = FALSE)

    check.correlation(rho12, rho13, rho23)

    if (alternative == "two.sided" && rho12 == rho13) stop("alternative = 'two.sided' but rho12 = rho13.", call. = FALSE)

    cor.mat <- matrix(c(1, rho12, rho13,
                           rho12, 1, rho23,
                           rho13, rho23, 1),
                         nrow = 3, ncol = 3)
    check.correlation.matrix(cor.mat)

    # common index
    if (pooled) {

      rho.bar.ab.ac <- (rho12 + rho13) / 2

      ## under null
      psi.ab.ac.0 <- rho23 * (1 - 2 * rho.bar.ab.ac ^ 2) - 0.50 * (rho.bar.ab.ac ^ 2) * (1 - 2 * rho.bar.ab.ac ^ 2 - rho23 ^ 2)
      cov.ab.ac.0 <- psi.ab.ac.0 / (1 - rho.bar.ab.ac ^ 2) ^ 2 # both = (rho12 + rho13) / 2 when pooled
      # sigma.ab.ac.0 <- sqrt((2 - 2 * cov.ab.ac.0) / (n - 3))

      ## under alt
      psi.ab.ac.1 <- rho23 * (1 - rho12 ^ 2 - rho13 ^ 2) - 0.50 * (rho12 * rho13) * (1 - rho12 ^ 2 - rho13 ^ 2 - rho23 ^ 2)
      cov.ab.ac.1 <- psi.ab.ac.1 / ((1 - rho12 ^ 2) * (1 - rho13 ^ 2))
      # sigma.ab.ac.1 <- sqrt((2 - 2 * cov.ab.ac.1) / (n - 3))

    } else {

      ## under null
      psi.ab.ac.0 <- rho23 * (1 - rho12 ^ 2 - rho12 ^ 2) - 0.50 * (rho12 * rho12) * (1 - rho12 ^ 2 - rho12 ^ 2 - rho23 ^ 2) # rho12 = rho13
      cov.ab.ac.0 <- psi.ab.ac.0 / ((1 - rho12 ^ 2) * (1 - rho12 ^ 2)) # rho12 = rho13
      # sigma.ab.ac.0 <- sqrt((2 - 2 * cov.ab.ac.0) / (n - 3))

      ## under alt
      # psi.ab.ac.1 <- rho23 * (1 - rho12 ^ 2 - rho13 ^ 2) - 0.50 * (rho12 * rho13) * (1 - rho12 ^ 2 - rho13 ^ 2 - rho23 ^ 2)
      cov.ab.ac.1 <- psi.ab.ac.1 / ((1 - rho12 ^ 2) * (1 - rho13 ^ 2))
      # sigma.ab.ac.1 <- sqrt((2 - 2 * cov.ab.ac.1) / (n - 3))

    } # if pooled

    # z.ab <- cor.to.z(rho12)
    # z.ac <- cor.to.z(rho13)
    # sigma.ab.ac.0 <- sqrt((2 - 2 * cov.ab.ac.0) / (n - 3))
    # sigma.ab.ac.1 <- sqrt((2 - 2 * cov.ab.ac.1) / (n - 3))
    # lambda <- (z.ac - z.ab) /  sigma.ab.ac.0
    # sigma.lambda <- sigma.ab.ac.1 / sigma.ab.ac.0
    rho1 <- rho12
    rho2 <- rho13
    cov.null <- cov.ab.ac.0
    cov.alt <- cov.ab.ac.1

  } else {

    check.correlation(rho14, rho24, rho34)

    if (alternative == "two.sided" && rho12 == rho34) stop("alternative = 'two.sided' but rho12 = rho34.", call. = FALSE)

    cor.mat <- matrix(c(1, rho12, rho13, rho14,
                             rho12, 1, rho23, rho24,
                             rho13, rho23, 1, rho34,
                             rho14, rho24, rho34, 1),
                           nrow = 4, ncol = 4)
    check.correlation.matrix(cor.mat)

    # no common index
    if (pooled) {

      rho.bar.ab.cd <- (rho12 + rho34) / 2

      ## under null
      psi.ab.cd.0 <- 0.50 * ((rho13 - rho.bar.ab.cd * rho23) * (rho24 - rho23 * rho.bar.ab.cd) +
                               (rho14 - rho13 * rho.bar.ab.cd) * (rho23 - rho.bar.ab.cd * rho13) +
                               (rho13 - rho14 * rho.bar.ab.cd) * (rho24 - rho.bar.ab.cd * rho14) +
                               (rho14 - rho.bar.ab.cd * rho24) * (rho23 - rho24 * rho.bar.ab.cd)) # rho12 = rho34
      cov.ab.cd.0 <- psi.ab.cd.0 / (1 - rho.bar.ab.cd ^ 2) ^ 2
      # sigma.ab.cd.0 <- sqrt((2 - 2 * cov.ab.cd.0) / (n - 3))

      ## under alt
      psi.ab.cd.1 <- 0.50 * ((rho13 - rho12 * rho23) * (rho24 - rho23 * rho34) +
                               (rho14 - rho13 * rho34) * (rho23 - rho12 * rho13) +
                               (rho13 - rho14 * rho34) * (rho24 - rho12 * rho14) +
                               (rho14 - rho12 * rho24) * (rho23 - rho24 * rho34))
      cov.ab.cd.1 <- psi.ab.cd.1 / ((1 - rho12 ^ 2) * (1 - rho34 ^ 2))
      # sigma.ab.cd.1 <- sqrt((2 - 2 * cov.ab.cd.1) / (n - 3))

    } else {

      ## under null
      psi.ab.cd.0 <- 0.50 * ((rho13 - rho12 * rho23) * (rho24 - rho23 * rho12) +
                               (rho14 - rho13 * rho12) * (rho23 - rho12 * rho13) +
                               (rho13 - rho14 * rho12) * (rho24 - rho12 * rho14) +
                               (rho14 - rho12 * rho24) * (rho23 - rho24 * rho12)) # rho12 = rho34
      cov.ab.cd.0 <- psi.ab.cd.0 / ((1 - rho12 ^ 2) * (1 - rho12 ^ 2))
      # sigma.ab.cd.0 <- sqrt((2 - 2 * cov.ab.cd.0) / (n - 3))

      ## under alt
      psi.ab.cd.1 <- 0.50 * ((rho13 - rho12 * rho23) * (rho24 - rho23 * rho34) +
                               (rho14 - rho13 * rho34) * (rho23 - rho12 * rho13) +
                               (rho13 - rho14 * rho34) * (rho24 - rho12 * rho14) +
                               (rho14 - rho12 * rho24) * (rho23 - rho24 * rho34))
      cov.ab.cd.1 <- psi.ab.cd.1 / ((1 - rho12 ^ 2) * (1 - rho34 ^ 2))
      # sigma.ab.cd.1 <- sqrt((2 - 2 * cov.ab.cd.1) / (n - 3))

    } # if pooled

    # z.ab <- cor.to.z(rho12)
    # z.cd <- cor.to.z(rho34)
    # sigma.ab.cd.0 <- sqrt((2 - 2 * cov.ab.cd.0) / (n - 3))
    # sigma.ab.cd.1 <- sqrt((2 - 2 * cov.ab.cd.1) / (n - 3))
    # lambda <- (z.cd - z.ab) /  sigma.ab.cd.0
    # sigma.lambda <- sigma.ab.cd.1 / sigma.ab.cd.0
    rho1 <- rho12
    rho2 <- rho34
    cov.null <- cov.ab.cd.0
    cov.alt <- cov.ab.cd.1

  } # if common.index

  if (is.null(power)) {

    pwr.obj <- pwr.steiger(rho1 = rho1, rho2 = rho2,
                           cov.null = cov.null, cov.alt = cov.alt,
                           n = n, alpha = alpha,
                           alternative = alternative)

    power <- pwr.obj$power
    mean.alternative <- pwr.obj$mean
    sd.alternative <- pwr.obj$sd
    mean.null <- pwr.obj$null.mean
    sd.null <- pwr.obj$null.sd
    z.alpha <- pwr.obj$z.alpha

    if (alternative == "two.sided" && rho1 - rho2 < 0) mean.alternative <- -mean.alternative

  } # power

  if (is.null(n)) {

    n <- ss.steiger(rho1 = rho1, rho2 = rho2,
                          cov.null = cov.null, cov.alt = cov.alt,
                          power = power, alpha = alpha,
                          alternative = alternative)

    if (ceiling) {
      n <- ceiling(n)
    }

    pwr.obj <- pwr.steiger(rho1 = rho1, rho2 = rho2,
                           cov.null = cov.null, cov.alt = cov.alt,
                           n = n, alpha = alpha,
                           alternative = alternative)

    power <- pwr.obj$power
    mean.alternative <- pwr.obj$mean
    sd.alternative <- pwr.obj$sd
    mean.null <- pwr.obj$null.mean
    sd.null <- pwr.obj$null.sd
    z.alpha <- pwr.obj$z.alpha

    if (alternative == "two.sided" && rho1 - rho2 < 0) mean.alternative <- -mean.alternative

  } # ss

  delta <- rho1 - rho2
  q <- cors.to.q(rho1, rho2, FALSE)$q

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

    print.obj <-  list(requested = requested,
                       test = "Dependent Correlations",
                       design = "paired",
                       alpha = alpha,
                       alt = alternative,
                       common = common.index,
                       delta = delta,
                       q = q,
                       mean.alternative = mean.alternative,
                       sd.alternative = sd.alternative,
                       mean.null = mean.null,
                       sd.null = sd.null,
                       z.alpha = z.alpha,
                       power = power,
                       n = n)

    if (pretty) {
      .print.pwrss.steiger(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.steiger(print.obj, verbose = verbose)
    }

  }

  invisible(structure(list(parms = list(rho12 = rho12, rho13 = rho13, rho23 = rho23,
                                        rho14 = rho14, rho24 = rho24, rho34 = rho34,
                                        alpha = alpha, alternative = alternative,
                                        pooled = pooled, common.index = common.index,
                                        ceiling = ceiling, verbose = verbose,
                                        pretty = pretty),
                           test = "z",
                           design = "paired",
                           delta = delta,
                           q = q,
                           mean = mean.alternative,
                           sd = sd.alternative,
                           null.mean = mean.null,
                           null.sd = sd.null,
                           alternative = alternative,
                           z.alpha = z.alpha,
                           power = power,
                           n = n),
                      class = c("pwrss", "z", "twocors", "paired")))

} # power.z.steiger()
power.z.twocors.steiger <- power.z.steiger


power.z.twocors <- function(rho1, rho2,
                            n2 = NULL, n.ratio = 1,
                            power = NULL, alpha = 0.05,
                            alternative = c("two.sided", "one.sided"),
                            ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  # old <- list(...)
  # user.parms <- as.list(match.call(expand.dots = TRUE))
  # names.user.parms <- names(user.parms)

  check.positive(n.ratio)
  check.proportion(alpha)
  check.correlation(rho1, rho2)
  check.logical(ceiling)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n2)) check.sample.size(n2)

  alternative <- tolower(match.arg(alternative))

  if (is.null(n2) && is.null(power)) stop("`n2` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n2) && !is.null(power)) stop("Exactly one of the `n2` or `power` should be `NULL`.", call. = FALSE)

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  if (alternative == "two.sided") {

    if (is.null(n2)) {
      beta <- 1 - power
      z1 <- cor.to.z(rho1, FALSE)$z
      z2 <- cor.to.z(rho2, FALSE)$z
      M <- qnorm(alpha / 2, mean = 0, sd = 1, lower.tail = FALSE) + qnorm(beta, mean = 0, sd = 1, lower.tail = FALSE)
      n2 <- uniroot(function(n2) M ^ 2 - (z1 - z2) ^ 2 / (1 / (n.ratio * n2 - 3) + 1 / (n2 - 3)), interval = c(-1e10, 1e10))$root
      n1 <- n.ratio * n2

      if (ceiling) {
        n2 <- ceiling(n2)
        n1 <- ceiling(n1)
        lambda <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
        z.alpha <- qnorm(alpha / 2, mean = 0, sd = 1, lower.tail = FALSE)
        power <- 1 - pnorm(z.alpha, mean = abs(lambda), sd = 1) + pnorm(-z.alpha, mean = abs(lambda), sd = 1)
        z.alpha <- c(-z.alpha, z.alpha)
      }
    }

    if (is.null(power)) {
      z1 <- cor.to.z(rho1, FALSE)$z
      z2 <- cor.to.z(rho2, FALSE)$z
      n1 <- n.ratio * n2
      if (ceiling) {
        n2 <- ceiling(n2)
        n1 <- ceiling(n1)
      }
      lambda <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
      z.alpha <- qnorm(alpha / 2, mean = 0, sd = 1, lower.tail = FALSE)
      power <- 1 - pnorm(z.alpha, mean = abs(lambda), sd = 1) + pnorm(-z.alpha, mean = abs(lambda), sd = 1)
      z.alpha <- c(-z.alpha, z.alpha)
    }

  }

  if (alternative == "one.sided") {

    if (is.null(n2)) {
      beta <- 1 - power
      z1 <- cor.to.z(rho1, FALSE)$z
      z2 <- cor.to.z(rho2, FALSE)$z
      M <- qnorm(alpha, mean = 0, sd = 1, lower.tail = FALSE) + qnorm(beta, mean = 0, sd = 1, lower.tail = FALSE)
      n2 <- uniroot(function(n2) M ^ 2 - (z1 - z2) ^ 2 / (1 / (n.ratio * n2 - 3) + 1 / (n2 - 3)), interval = c(0, 1e10))$root
      n1 <- n.ratio * n2

      if (ceiling) {
        n2 <- ceiling(n2)
        n1 <- ceiling(n1)
        lambda <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
        z.alpha <- qnorm(alpha, mean = 0, sd = 1, lower.tail = FALSE)
        power <- 1 - pnorm(z.alpha, mean = abs(lambda), sd = 1)
        if (lambda < 0) z.alpha <- -z.alpha
      }
    }

    if (is.null(power)) {
      z1 <- cor.to.z(rho1, FALSE)$z
      z2 <- cor.to.z(rho2, FALSE)$z
      n1 <- n.ratio * n2
      if (ceiling) {
        n2 <- ceiling(n2)
        n1 <- ceiling(n1)
      }
      lambda <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
      z.alpha <- qnorm(alpha,  mean = 0, sd = 1, lower.tail = FALSE)
      power <- 1 - pnorm(z.alpha, mean = abs(lambda), sd = 1)
      if (lambda < 0) z.alpha <- -z.alpha
    }

  }

  delta <- rho1 - rho2
  q <- cors.to.q(rho1, rho2, FALSE)$q

  mean.alternative <- lambda
  sd.alternative <- 1
  mean.null <- 0
  sd.null <- 1

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

    print.obj <-  list(requested = requested,
                       test = "Independent Correlations",
                       design = "independent",
                       alpha = alpha,
                       alt = alternative,
                       delta = delta,
                       q = q,
                       mean.alternative = mean.alternative,
                       sd.alternative = sd.alternative,
                       mean.null = mean.null,
                       sd.null = sd.null,
                       z.alpha = z.alpha,
                       power = power,
                       n = c(n1 = n1, n2 = n2))

    if (pretty) {
      .print.pwrss.twocors(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.twocors(print.obj, verbose = verbose)
    }

  }

  invisible(structure(list(parms = list(rho1 = rho1, rho2 = rho2, n.ratio = n.ratio,
                                        alpha = alpha, alternative = alternative,
                                        ceiling = ceiling, verbose = verbose,
                                        pretty = pretty),
                           test = "z",
                           design = "independent",
                           delta = delta,
                           q = q,
                           mean = mean.alternative,
                           sd = sd.alternative,
                           null.mean = mean.null,
                           null.sd = sd.null,
                           alternative = alternative,
                           z.alpha = z.alpha,
                           n = c(n1 = n1, n2 = n2),
                           power = power),
                      class = c("pwrss", "z", "twocors", "independent")))

} # power.z.twocors
power.z.twocor <- power.z.twocors


##########################
# one correlation z test #
##########################

power.z.onecor <- function(rho, null.rho = 0,
                           n = NULL, power = NULL, alpha = 0.05,
                           alternative = c("two.sided", "one.sided"),
                           ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  # old <- list(...)
  # user.parms <- as.list(match.call(expand.dots = TRUE))
  # names.user.parms <- names(user.parms)

  check.proportion(alpha)
  check.correlation(rho, null.rho)
  check.logical(ceiling)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n)) check.sample.size(n)

  alternative <- tolower(match.arg(alternative))

  if (is.null(n) && is.null(power)) stop("`n` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n) && !is.null(power)) stop("Exactly one of the `n` or `power` should be `NULL`.", call. = FALSE)

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  z <- cor.to.z(rho, FALSE)$z
  null.z <- cor.to.z(null.rho, FALSE)$z

  if (alternative == "two.sided") {

    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha / 2, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n <- M ^ 2 / (z - null.z) ^ 2 + 3
      if (ceiling) {
        n <- ceiling(n)
        lambda <- (z - null.z) / sqrt(1 / (n - 3))
        z.alpha <- qnorm(alpha / 2, lower.tail = FALSE)
        power <- 1 - pnorm(z.alpha, lambda) + pnorm(-z.alpha, lambda)
        z.alpha <- c(-z.alpha, z.alpha)
      }

    }

    if (is.null(power)) {
      lambda <- (z - null.z) / sqrt(1 / (n - 3))
      z.alpha <- qnorm(alpha / 2, lower.tail = FALSE)
      power <- 1 - pnorm(z.alpha, lambda) + pnorm(-z.alpha, lambda)
      z.alpha <- c(-z.alpha, z.alpha)
    }

  }

  if (alternative == "one.sided") {

    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n <- M ^ 2 / (z - null.z) ^ 2 + 3
      if (ceiling) {
        n <- ceiling(n)
        lambda <- (z - null.z) / sqrt(1 / (n - 3))
        z.alpha <- qnorm(alpha, lower.tail = FALSE)
        power <- 1 - pnorm(z.alpha, abs(lambda))
      }
    }

    if (is.null(power)) {
      lambda <- (z - null.z) / sqrt(1 / (n - 3))
      z.alpha <- qnorm(alpha, lower.tail = FALSE)
      power <- 1 - pnorm(z.alpha, abs(lambda))
    }

  }

  delta <- rho - null.rho
  q <- cors.to.q(rho, null.rho, FALSE)$q

  mean.alternative <- lambda
  sd.alternative <- 1
  mean.null <- 0
  sd.null <- 1

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

    print.obj <-  list(requested = requested,
                       test = "One-Sample Correlation",
                       design = "one.sample",
                       alpha = alpha,
                       alt = alternative,
                       delta = delta,
                       q = q,
                       mean.alternative = mean.alternative,
                       sd.alternative = sd.alternative,
                       mean.null = mean.null,
                       sd.null = sd.null,
                       z.alpha = z.alpha,
                       power = power,
                       n = n)

    if (pretty) {
      .print.pwrss.twocors(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.twocors(print.obj, verbose = verbose)
    }

  }

  invisible(structure(list(parms = list(rho = rho, null.rho = null.rho,
                                        alpha = alpha, alternative = alternative,
                                        ceiling = ceiling, verbose = verbose,
                                        pretty = pretty),
                           test = "z",
                           design = "one.sample",
                           delta = delta,
                           q = q,
                           mean = mean.alternative,
                           sd = sd.alternative,
                           null.mean = mean.null,
                           null.sd = sd.null,
                           alternative = alternative,
                           z.alpha = z.alpha,
                           n = n,
                           power = power),
                      class = c("pwrss", "z", "onecor")))

} # power.z.onecor()


pwrss.z.cor <- function(r = 0.50, r0 = 0, alpha = 0.05,
                         alternative = c("not equal", "greater", "less"),
                         n = NULL, power = NULL, verbose = TRUE) {


  check.correlation(r, r0)
  check.logical(verbose)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n)) check.sample.size(n)

  alternative <- tolower(match.arg(alternative))
  if (alternative %in% c("less", "greater")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"

  onecor.obj <- power.z.onecor(rho = r, null.rho = r0,
                               n = n, power = power, alpha = alpha,
                               alternative = alternative,
                               ceiling = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.z.onecor() function. \n")

  return(invisible(onecor.obj))

} # pwrss.z.corr()

pwrss.z.corr <- pwrss.z.cor

pwrss.z.2cors <- function(r1 = 0.50, r2 = 0.30,
                           alpha = 0.05, kappa = 1,
                           alternative = c("not equal", "greater", "less"),
                           n2 = NULL, power = NULL, verbose = TRUE) {


  check.correlation(r1, r2)
  check.logical(verbose)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n2)) check.sample.size(n2)

  alternative <- tolower(match.arg(alternative))
  if (alternative %in% c("less", "greater")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"

  twocors.obj <- power.z.twocors(rho1 = r1, rho2 = r2,
                             n2 = n2, n.ratio = kappa,
                             power = power, alpha = alpha,
                             alternative = alternative,
                             ceiling = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.z.twocors() function. \n")

  return(invisible(twocors.obj))

} # pwrss.z.2corrs()

pwrss.z.2corrs <- pwrss.z.2cors
