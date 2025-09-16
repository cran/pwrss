######################
# poisson regression #
######################

# dist = c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")
# dist = list(dist = "normal", mean = 0, sd = 1)
# dist = list(dist = "poisson", lambda = 1)
# dist = list(dist = "uniform", min = 0, max = 1)
# dist = list(dist = "exponential", rate = 1)
# dist = list(dist = "binomial", size = 1, prob = 0.50)
# dist = list(dist = "bernoulli", prob = 0.50)
# dist = list(dist = "lognormal", meanlog = 0, sdlog = 1)

# mean.exposure is the mean exposure time (should be > 0)
power.z.poisson <- function(base.rate = NULL, rate.ratio = NULL,
                            beta0 = log(base.rate), beta1 = log(rate.ratio),
                            n = NULL, power = NULL,
                            r.squared.predictor = 0, mean.exposure = 1,
                            alpha = 0.05, alternative = c("two.sided", "one.sided"),
                            method = c("demidenko(vc)", "demidenko", "signorini"),
                            distribution = "normal", ceiling = TRUE,
                            verbose = TRUE, pretty = FALSE) {

  # old.args <- list(...) # just arguments in ...
  # names.old.args <- names(old.args)
  user.parms.names <- names(as.list(match.call()))

  check.positive(mean.exposure)
  check.proportion(alpha, r.squared.predictor)
  check.logical(ceiling)
  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.proportion(power)

  alternative <- tolower(match.arg(alternative))
  method <- tolower(match.arg(method))

  if (all(c("beta0", "beta1") %in% user.parms.names)) {
    check.numeric(beta0, beta1)
    base.rate <- exp(beta0)
    rate.ratio <- exp(beta1)
    if (any(c("base.rate", "rate.ratio") %in% user.parms.names))
      message("Ignoring any specifications to 'base.rate' or 'rate.ratio'.")
  }

  if (all(c("base.rate", "rate.ratio") %in% user.parms.names)) {
    check.nonnegative(base.rate, rate.ratio)
    beta0 <- log(base.rate)
    beta1 <- log(rate.ratio)
    if (any(c("beta0", "beta1") %in% user.parms.names))
      message("Ignoring any specifications to 'beta0' or 'beta1'.")
  }

  if (all(c("base.rate", "rate.ratio", "beta0", "beta1") %in% user.parms.names)) {
    stop("Specify 'base.rate' & 'rate.ratio' or 'beta0' & 'beta1'.", call. = FALSE)
  }

  if (beta0 == beta1) stop("`beta0` = `beta1`?", call. = FALSE)
  if (is.null(n) && is.null(power)) stop("`n` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n) && !is.null(power)) stop("Exactly one of the `n` or `power` should be `NULL`.", call. = FALSE)

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  if (length(distribution) == 1 && is.character(distribution)) {
    distribution <- switch(tolower(distribution),
                            `normal` = list(dist = "normal", mean = 0, sd = 1),
                            `poisson` = list(dist = "poisson", lambda = 1),
                            `uniform` = list(dist = "uniform", min = 0, max = 1),
                            `exponential` = list(dist = "exponential", rate = 1),
                            `binomial` = list(dist = "binomial", size = 1, prob = 0.50),
                            `bernoulli` = list(dist = "bernoulli", size = 1, prob = 0.50),
                            `lognormal` = list(dist = "lognormal", meanlog = 0, sdlog = 1))
  } else if (is.list(distribution)) {
    if (length(distribution) > 3) stop("unknown input type for 'distribution' argument", call. = FALSE)
    dist.list.names <- names(distribution)
    dist.attrib <- c(dist.list.names, tolower(distribution$dist))
    dist.invalid <- c(any(is.na(match(dist.attrib, c("dist", "normal", "mean", "sd")))),
                      any(is.na(match(dist.attrib, c("dist", "lognormal", "meanlog", "sdlog")))),
                      any(is.na(match(dist.attrib, c("dist", "uniform", "min", "max")))),
                      any(is.na(match(dist.attrib, c("dist", "exponential", "rate")))),
                      any(is.na(match(dist.attrib, c("dist", "poisson", "lambda")))),
                      any(is.na(match(dist.attrib, c("dist", "binomial", "bernoulli", "size", "prob")))))
    if (all(dist.invalid == TRUE)) stop("unknown input type for 'distribution' argument", call. = FALSE)
  } else {
    stop("Unknown input type for 'distribution'.", call. = FALSE)
  }

  # asymptotic variances
  var.beta <- function(beta0, beta1, distribution) {

    # asymptotic variances
    if (tolower(distribution$dist) == "normal") {

      mean <- distribution$mean
      sd <- distribution$sd

      min.norm <- qnorm(.0000001, mean = mean, sd = sd)
      max.norm <- qnorm(.9999999, mean = mean, sd = sd)

      # variance under null
      mu <- integrate(function(x)  dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x), min.norm, max.norm)$value
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- integrate(function(x)  dnorm(x, mean = mean, sd = sd) * exp(beta0.star + beta1.star * x), min.norm, max.norm)$value
      i01 <- integrate(function(x) x * dnorm(x, mean = mean, sd = sd) * exp(beta0.star + beta1.star * x), min.norm, max.norm)$value
      i11 <- integrate(function(x) x ^ 2 * dnorm(x, mean = mean, sd = sd) * exp(beta0.star + beta1.star * x), min.norm, max.norm)$value
      var.beta0 <- i00 / (i00 * i11 - i01 ^ 2)

      # variance under alternative
      i00 <- integrate(function(x)  dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x), min.norm, max.norm)$value
      i01 <- integrate(function(x) x * dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x), min.norm, max.norm)$value
      i11 <- integrate(function(x) x ^ 2 * dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x), min.norm, max.norm)$value
      var.beta1 <- i00 / (i00 * i11 - i01 ^ 2)

    }

    if (tolower(distribution$dist) == "poisson") {

      lambda <- distribution$lambda
      # maximum value
      max.pois <- qpois(.999999999, lambda = lambda)

      # variance under null
      mu <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:max.pois, function(x) x * dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:max.pois, function(x) x ^ 2 * dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      var.beta0 <- i00 / (i00 * i11 - i01 ^ 2)

      # variance under alternative
      i00 <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:max.pois, function(x) x * dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:max.pois, function(x) x ^ 2 * dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      var.beta1 <- i00 / (i00 * i11 - i01 ^ 2)

    }

    if (tolower(distribution$dist) == "uniform") {

      min <- distribution$min
      max <- distribution$max

      # variance under null
      mu <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x), min, max)$value
      i01 <- integrate(function(x) x * dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x), min, max)$value
      i11 <- integrate(function(x) x ^ 2 * dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x), min, max)$value
      var.beta0 <- i00 / (i00 * i11 - i01 ^ 2)

      # variance under alternative
      i00 <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      i01 <- integrate(function(x) x * dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      i11 <- integrate(function(x) x ^ 2 * dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      var.beta1 <- i00 / (i00 * i11 - i01 ^ 2)

    }

    if (tolower(distribution$dist) == "exponential") {

      rate <- distribution$rate
      max.exp <- qexp(.9999999, rate = rate)

      # variance under null
      mu <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0.star + beta1.star * x), 0, max.exp)$value
      i01 <- integrate(function(x) x * dexp(x, rate = rate) * exp(beta0.star + beta1.star * x), 0, max.exp)$value
      i11 <- integrate(function(x) x ^ 2 * dexp(x, rate = rate) * exp(beta0.star + beta1.star * x), 0, max.exp)$value
      var.beta0 <- i00 / (i00 * i11 - i01 ^ 2)

      # variance under alternative
      i00 <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      i01 <- integrate(function(x) x * dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      i11 <- integrate(function(x) x ^ 2 * dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      var.beta1 <- i00 / (i00 * i11 - i01 ^ 2)

    }

    if (tolower(distribution$dist) %in% c("binomial", "bernoulli")) {

      ifelse(tolower(distribution$dist) == "bernoulli", size <- 1, size <- distribution$size)
      prob <- distribution$prob

      # variance under null
      mu <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:size, function(x) x * dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:size, function(x) x ^ 2 * dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      var.beta0 <- i00 / (i00 * i11 - i01 ^ 2)

      # variance under alternative
      i00 <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:size, function(x) x * dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:size, function(x) x ^ 2 * dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      var.beta1 <- i00 / (i00 * i11 - i01 ^ 2)

    }

    if (tolower(distribution$dist) == "lognormal") {

      meanlog <- distribution$meanlog
      sdlog <- distribution$sdlog
      min.lnorm <- qlnorm(.0000001, meanlog = meanlog, sdlog = sdlog)
      max.lnorm <- qlnorm(.9999999, meanlog = meanlog, sdlog = sdlog)

      # variance under null
      mu <- integrate(function(x)  dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x), min.lnorm, max.lnorm)$value
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- integrate(function(x)  dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0.star + beta1.star * x), min.lnorm, max.lnorm)$value
      i01 <- integrate(function(x) x * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0.star + beta1.star * x), min.lnorm, max.lnorm)$value
      i11 <- integrate(function(x) x ^ 2 * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0.star + beta1.star * x), min.lnorm, max.lnorm)$value
      var.beta0 <- i00 / (i00 * i11 - i01 ^ 2)

      # variance under alternative
      i00 <- integrate(function(x)  dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x), min.lnorm, max.lnorm)$value
      i01 <- integrate(function(x) x * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x), min.lnorm, max.lnorm)$value
      i11 <- integrate(function(x) x ^ 2 * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x), min.lnorm, max.lnorm)$value
      var.beta1 <- i00 / (i00 * i11 - i01 ^ 2)

    }

    list(var.beta0 = var.beta0, var.beta1 = var.beta1, distribution = tolower(distribution$dist))

  } # var.beta()

  pwr.demidenko <- function(beta0, beta1, n,
                            r.squared.predictor,
                            alpha, alternative,
                            method, distribution) {

    # variance correction factor
    if (tolower(method) == "demidenko(vc)") {
      vcf <- switch(tolower(distribution$dist),
                     `normal` = 1,
                     `poisson` = 1,
                     `uniform` = 1,
                     `exponential` = 1,
                     `binomial` = 1, # 0.85
                     `bernoulli` = 1, # 0.85
                     `lognormal` = 0.75)
    } else if (tolower(method) == "demidenko") {
      vcf <- 0
    }

    var.obj <- var.beta(beta0 = beta0, beta1 = beta1, distribution = distribution)
    var.beta0 <- var.obj$var.beta0
    var.beta1 <- var.obj$var.beta1

    # non-centrality parameter and standard deviation of the non-centrality parameter under alternative
    # Signorini, D. F. (1991). Sample size for poisson regression. Biometrika, 78, 446-450.
    # Demidenko, E. (2007). Sample size determination for logistic regression revisited. Statistics in Medicine, 26, 3385-3397.
    if (method == "signorini") {
      ncp <- beta1 / sqrt(var.beta0 / (n * (1 - r.squared.predictor) * mean.exposure))
      sd.ncp <- sqrt(var.beta1 / var.beta0)
      vcf <- NA
    } else {
      ncp <- beta1 / sqrt(var.beta1 / (n * (1 - r.squared.predictor) * mean.exposure))
      sd.ncp <- sqrt((vcf * var.beta0 + (1 - vcf) * var.beta1) / var.beta1)
    } # method

    pwr.obj <- power.z.test(mean = ncp, sd = sd.ncp, null.mean = 0,
                            alpha = alpha, alternative = alternative,
                            plot = FALSE, verbose = FALSE)
    power <- pwr.obj$power
    z.alpha <- pwr.obj$z.alpha

    list(power = power, ncp = ncp, sd.ncp = sd.ncp, vcf = vcf, z.alpha = z.alpha)

  } # pwr.demidenko()

  ss.demidenko <- function(beta0, beta1, power,
                           r.squared.predictor,
                           alpha, alternative,
                           method, distribution) {

    n <- uniroot(function(n) {
      power - pwr.demidenko(beta0 = beta0, beta1 = beta1, n = n,
                            r.squared.predictor = r.squared.predictor,
                            alpha = alpha, alternative = alternative,
                            method = method, distribution = distribution)$power
    }, interval = c(2, 1e10))$root

    n

  } # ss.demidenko()



  if (is.null(power)) {
    pwr.obj <- pwr.demidenko(beta0 = beta0, beta1 = beta1, n = n,
                             r.squared.predictor = r.squared.predictor,
                             alpha = alpha, alternative = alternative,
                             method = method, distribution = distribution)
    power <- pwr.obj$power
    z.alpha <- pwr.obj$z.alpha
    ncp <- pwr.obj$ncp
    sd.ncp <- pwr.obj$sd.ncp
    vcf <- pwr.obj$vcf

  } # pwr

  if (is.null(n)) {

    n <- ss.demidenko(beta0 = beta0, beta1 = beta1, power = power,
                      r.squared.predictor = r.squared.predictor,
                      alpha = alpha, alternative = alternative,
                      method = method, distribution = distribution)

    if (ceiling) {
      n <- ceiling(n)
    }

    pwr.obj <- pwr.demidenko(beta0 = beta0, beta1 = beta1, n = n,
                             r.squared.predictor = r.squared.predictor,
                             alpha = alpha, alternative = alternative,
                             method = method, distribution = distribution)
    power <- pwr.obj$power
    z.alpha <- pwr.obj$z.alpha
    ncp <- pwr.obj$ncp
    sd.ncp <- pwr.obj$sd.ncp
    vcf <- pwr.obj$vcf

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

    print.obj <- list(requested = requested,
                      test = "Poisson Regression Coefficient (Wald's Z-Test)",
                      method = switch(method,
                                      `demidenko(vc)` = "Demidenko (Variance Corrected)",
                                      `demidenko` = "Demidenko",
                                      `signorini` = "Signorini"),
                      dist = switch(tolower(distribution$dist),
                                    `normal` = "Normal",
                                    `poisson` = "Poisson",
                                    `bernoulli` = "Bernoulli",
                                    `binomial` = "Binomial",
                                    `lognormal` = "Log-normal",
                                    `uniform` = "Uniform",
                                    `exponential` = "Exponential"),
                      alt = alternative,
                      base.rate = base.rate,
                      rate.ratio = rate.ratio,
                      n = n,
                      mean.alternative = ncp,
                      sd.alternative = sd.ncp,
                      vcf = vcf,
                      mean.null = 0,
                      sd.null = 1,
                      z.alpha = z.alpha,
                      alpha = alpha,
                      power = power)

    if (pretty) {
      .print.pwrss.poisson(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.poisson(print.obj, verbose = verbose)
    }

  } # verbose

  invisible(structure(list(parms = list(base.rate = base.rate,
                                        rate.ratio = rate.ratio,
                                        beta0 = beta0, beta1 = beta1,
                                        r.squared.predictor = r.squared.predictor,
                                        mean.exposure = mean.exposure,
                                        alpha = alpha, alternative = alternative,
                                        method = method, distribution =  distribution,
                                        ceiling = ceiling, verbose = verbose,
                                        pretty = pretty),
                           test = "z",
                           base.rate = base.rate,
                           rate.ratio = rate.ratio,
                           mean = ncp,
                           sd = sd.ncp,
                           vcf = vcf,
                           null.mean = 0,
                           null.sd = 1,
                           z.alpha = z.alpha,
                           power = power,
                           n = n),
                      class = c("pwrss", "z", "poisson")))
} # power.z.poisson()
power.z.poisreg <- power.z.poisson

pwrss.z.poisson <- function(exp.beta0 = 1.10, exp.beta1 = 1.16,
                            beta0 = log(exp.beta0), beta1 = log(exp.beta1),
                            mean.exposure = 1, n = NULL, power = NULL, r2.other.x = 0,
                            alpha = 0.05, alternative = c("not equal", "less", "greater"),
                            method = c("demidenko(vc)", "demidenko", "signorini"),
                            distribution = "normal", verbose = TRUE) {

  method <- tolower(match.arg(method))
  alternative <- tolower(match.arg(alternative))

  if (alternative %in% c("less", "greater")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"

  poisreg.obj <- power.z.poisson(beta0 = beta0, beta1 = beta1, n = n, power = power,
                                 r.squared.predictor = r2.other.x, mean.exposure = mean.exposure,
                                 alpha = alpha, alternative = alternative, method = method,
                                 distribution = distribution, ceiling = TRUE, verbose = verbose)

  # cat("This function will be removed in the future. \n Please use power.z.poisson() function. \n")

  return(invisible(poisreg.obj))

} # pwrss.z.poisson
pwrss.z.poisreg <- pwrss.z.poisson
