####################
# mediation z test #
####################

## 'cp = 0' by default, implying complete mediation (it increases explanatory power of the covariate
# use 'r.squared.mediator' and 'r.squared.outcome' to adjust standard error for other predictors in mediation and outcome model
power.z.mediation  <- function(beta.a, beta.b, beta.cp = 0,
                               sd.predictor = 1, sd.mediator = 1, sd.outcome = 1,
                               r.squared.mediator = beta.a ^ 2 * sd.predictor ^ 2 / sd.mediator ^ 2,
                               r.squared.outcome = (beta.b ^ 2 * sd.mediator ^ 2 + beta.cp ^ 2 * sd.predictor ^ 2) / sd.outcome ^ 2,
                               n = NULL, power = NULL, alpha = 0.05,
                               alternative = c("two.sided", "one.sided"),
                               method = c("sobel", "aroian", "goodman",
                                          "joint", "monte.carlo"),
                               n.simulation = 1000, n.draws = 1000,
                               ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  user.parms.names <- names(as.list(match.call()))

  alternative <- tolower(match.arg(alternative))
  method <- tolower(match.arg(method))

  check.logical(ceiling)
  check.numeric(beta.a, beta.b, beta.cp)
  check.sample.size(n.simulation, n.draws)
  check.positive(sd.predictor, sd.mediator, sd.outcome)
  check.proportion(alpha, r.squared.mediator, r.squared.outcome)

  if (!is.null(n)) check.sample.size(n)
  if (!is.null(power)) check.proportion(power)
  if (is.null(n) && is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time.", call. = FALSE)
  if (!is.null(n) && !is.null(power))
    stop("Exactly one of the `n` or `power` should be `NULL`.", call. = FALSE)

  if (r.squared.outcome == 0 && "beta.cp" %in% user.parms.names)
    warning("Ignoring any specification to 'beta.cp'.", call. = FALSE)
  if (r.squared.mediator < beta.a ^ 2 * sd.predictor ^ 2 / sd.mediator ^ 2)
    warning("Specified 'r.squared.mediator' is smaller than the base 'r.squared.mediator'.", call. = FALSE)
  if (r.squared.outcome < (beta.b ^ 2 * sd.mediator ^ 2 + beta.cp ^ 2 * sd.predictor ^ 2) / sd.outcome ^ 2)
    warning("Specified 'r.squared.outcome' is smaller than the base 'r.squared.outcome'.", call. = FALSE)

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  se.a <- function(sd.mediator, sd.predictor, r.squared.mediator, n) {
    var.beta.a <- (1 / n) * (sd.mediator ^ 2) * (1 - r.squared.mediator) / (sd.predictor ^ 2)
    sqrt(var.beta.a)
  }

  se.b <- function(sd.outcome, sd.mediator, r.squared.outcome, r.squared.mediator, n) {
    var.beta.b <- (1 / n) * (sd.outcome ^ 2) * (1 - r.squared.outcome) / ((sd.mediator ^ 2) * (1 - r.squared.mediator))
    sqrt(var.beta.b)
  }

  pwr.med <- function(beta.a, beta.b,
                      sd.predictor, sd.mediator, sd.outcome,
                      r.squared.mediator, r.squared.outcome,
                      n, method = c("sobel", "aroian", "goodman",
                                    "joint", "monte.carlo")) {

    method <- tolower(match.arg(method))

    se.beta.a <- se.a(sd.mediator = sd.mediator, sd.predictor = sd.predictor,
                      r.squared.mediator = r.squared.mediator, n = n)
    se.beta.b <- se.b(sd.outcome = sd.outcome, sd.mediator = sd.mediator,
                      r.squared.outcome = r.squared.outcome,
                      r.squared.mediator = r.squared.mediator, n = n)

    if (method == "sobel") {
      sobel.se <- sqrt(beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2)
      lambda <- (beta.a * beta.b) / sobel.se
      pwr.obj <- power.z.test(mean = lambda, alpha = alpha, alternative = alternative,
                              plot = FALSE, verbose = FALSE)
    }

    if (method == "aroian") {
      aroian.se <- sqrt(beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2 + se.beta.a ^ 2 * se.beta.b ^ 2)
      lambda <- (beta.a * beta.b) / aroian.se
      pwr.obj <- power.z.test(mean = lambda, alpha = alpha, alternative = alternative,
                              plot = FALSE, verbose = FALSE)
    }

    if (method == "goodman") {
      goodman.var <- beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2 - se.beta.a ^ 2 * se.beta.b ^ 2
      if (goodman.var <= 0) stop("Design is not feasible for Goodman's Z-Test", call. = FALSE)
      goodman.se <- sqrt(goodman.var)
      lambda <- (beta.a * beta.b) / goodman.se
      pwr.obj <- power.z.test(mean = lambda, alpha = alpha, alternative = alternative,
                              plot = FALSE, verbose = FALSE)
    }

    if (method == "joint") {
      power.beta.a <- power.z.test(mean = beta.a / se.beta.a, alpha = alpha, alternative = alternative,
                                   plot = FALSE, verbose = FALSE)$power
      power.beta.b <- power.z.test(mean = beta.b / se.beta.b, alpha = alpha, alternative = alternative,
                                   plot = FALSE, verbose = FALSE)$power
      power <- power.beta.a * power.beta.b
      pwr.obj <- list(alternative = alternative, mean = NA, sd = NA, null.mean = NA,
                      alpha = alpha, z.alpha = NA, power = power)
    }

    if (method == "monte.carlo") {
      beta.a <- abs(beta.a)
      beta.b <- abs(beta.b)
      reject <- numeric(0)
      for (i in 1:n.simulation) {
        beta.a.star <- rnorm(1, beta.a, se.beta.a)
        beta.b.star <- rnorm(1, beta.b, se.beta.b)
        reject <- c(reject, quantile(rnorm(n.draws, beta.a.star, se.beta.a) * rnorm(n.draws, beta.b.star, se.beta.b),
                                     probs = ifelse(alternative == "two.sided", alpha / 2, alpha), na.rm = TRUE) > 0)
      }
      power <- mean(reject)
      pwr.obj <- list(alternative = alternative, mean = NA, sd = NA, null.mean = NA,
                      alpha = alpha, z.alpha = NA, power = power)
    }

    return(pwr.obj)

  } # pwr.med()

  ss.med <- function(beta.a, beta.b,
                     sd.predictor, sd.mediator, sd.outcome,
                     r.squared.mediator, r.squared.outcome,
                     power, method = c("sobel", "aroian", "goodman")) {

    method <- tolower(match.arg(method))

    if (method == "sobel") {
      n <- uniroot(function(n) {
        se.beta.a <- se.a(sd.mediator = sd.mediator, sd.predictor = sd.predictor,
                          r.squared.mediator = r.squared.mediator, n = n)
        se.beta.b <- se.b(sd.outcome = sd.outcome, sd.mediator = sd.mediator,
                          r.squared.outcome = r.squared.outcome,
                          r.squared.mediator = r.squared.mediator, n = n)
        sobel.se <- sqrt(beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2)
        lambda.sobel <- (beta.a * beta.b) / sobel.se
        power - power.z.test(mean = lambda.sobel, alpha = alpha, alternative = alternative,
                             plot = FALSE, verbose = FALSE)$power
      }, interval = c(10, 1e10))$root
    }

    if (method == "aroian") {
      n <- uniroot(function(n) {
        se.beta.a <- se.a(sd.mediator = sd.mediator, sd.predictor = sd.predictor,
                          r.squared.mediator = r.squared.mediator, n = n)
        se.beta.b <- se.b(sd.outcome = sd.outcome, sd.mediator = sd.mediator,
                          r.squared.outcome = r.squared.outcome,
                          r.squared.mediator = r.squared.mediator, n = n)
        aroian.se <- sqrt(beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2 + se.beta.a ^ 2 * se.beta.b ^ 2)
        lambda.aroian <- (beta.a * beta.b) / aroian.se
        power - power.z.test(mean = lambda.aroian, alpha = alpha, alternative = alternative,
                             plot = FALSE, verbose = FALSE)$power
      }, interval = c(10, 1e10))$root
    }

    if (method == "goodman") {
      n <- uniroot(function(n) {
        se.beta.a <- se.a(sd.mediator = sd.mediator, sd.predictor = sd.predictor,
                          r.squared.mediator = r.squared.mediator, n = n)
        se.beta.b <- se.b(sd.outcome = sd.outcome, sd.mediator = sd.mediator,
                          r.squared.outcome = r.squared.outcome,
                          r.squared.mediator = r.squared.mediator, n = n)
        if (beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2 < se.beta.a ^ 2 * se.beta.b ^ 2) {
          goodman.se <- 1
        } else {
          goodman.se <- sqrt(beta.a ^ 2 * se.beta.b ^ 2  + beta.b ^ 2 * se.beta.a ^ 2 - se.beta.a ^ 2 * se.beta.b ^ 2)
        }
        lambda.goodman <- (beta.a * beta.b) / goodman.se
        power - power.z.test(mean = lambda.goodman, alpha = alpha, alternative = alternative,
                             plot = FALSE, verbose = FALSE)$power
      }, interval = c(10, 1e10))$root
    }

    return(n)

  } # ss.med()

  if (is.null(power)) {

    pwr.obj <- pwr.med(beta.a = beta.a, beta.b = beta.b,
                       sd.predictor = sd.predictor,
                       sd.mediator = sd.mediator, sd.outcome = sd.outcome,
                       r.squared.mediator = r.squared.mediator,
                       r.squared.outcome = r.squared.outcome,
                       n = n, method = method)

    power <- pwr.obj$power
    mean.alternative <- pwr.obj$mean
    sd.alternative <- 1
    mean.null <- pwr.obj$null.mean
    sd.null <- 1
    z.alpha <- pwr.obj$z.alpha

  }

  if (is.null(n)) {

    if (method %in% c("joint", "monte.carlo"))
      stop("Sample size calculation not supported by this method", call. = FALSE)

    n <- ss.med(beta.a = beta.a, beta.b = beta.b,
                sd.predictor = sd.predictor,
                sd.mediator = sd.mediator, sd.outcome = sd.outcome,
                r.squared.mediator = r.squared.mediator,
                r.squared.outcome = r.squared.outcome,
                power = power, method = method)

    if (ceiling) n <- ceiling(n)

    pwr.obj <- pwr.med(beta.a = beta.a, beta.b = beta.b,
                       sd.predictor = sd.predictor,
                       sd.mediator = sd.mediator, sd.outcome = sd.outcome,
                       r.squared.mediator = r.squared.mediator,
                       r.squared.outcome = r.squared.outcome,
                       n = n, method = method)

    power <- pwr.obj$power
    mean.alternative <- pwr.obj$mean
    sd.alternative <- 1
    mean.null <- pwr.obj$null.mean
    sd.null <- 1
    z.alpha <- pwr.obj$z.alpha

  }

  std.beta.a <- beta.a * sd.predictor / sd.mediator
  std.beta.b <- beta.b * sd.mediator / sd.outcome
  std.beta.cp <- beta.cp * sd.predictor / sd.outcome
  std.beta.indirect <-  std.beta.a * std.beta.b

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

    print.obj <- list(test = "Indirect Effect in a Mediation Model",
                      method = method,
                      requested = requested,
                      std.beta.a = std.beta.a,
                      std.beta.b = std.beta.b,
                      std.beta.cp = std.beta.cp,
                      std.beta.indirect = std.beta.indirect,
                      mean.alternative = mean.alternative,
                      sd.alternative = sd.alternative,
                      mean.null = mean.null,
                      sd.null = sd.null,
                      alpha = alpha,
                      alt = alternative,
                      z.alpha = z.alpha,
                      n = n,
                      power = power)

    if (pretty) {
      .print.pwrss.med(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.med(print.obj, verbose = verbose)
    }

  } # verbose

  invisible(structure(list(parms = list(beta.a = beta.a, beta.b = beta.b, beta.cp = beta.cp,
                                        sd.predictor = sd.predictor, sd.mediator = sd.mediator, sd.outcome = sd.outcome,
                                        r.squared.mediator = r.squared.mediator, r.squared.outcome = r.squared.outcome,
                                        alpha = alpha, alternative = alternative, method = method,
                                        ceiling = ceiling, verbose = verbose),
                           test = switch(method, `joint` = "joint", `monte.carlo` = "monte.carlo", "z"),
                           std.beta.a = std.beta.a,
                           std.beta.b = std.beta.b,
                           std.beta.cp = std.beta.cp,
                           std.beta.indirect = std.beta.indirect,
                           mean = mean.alternative,
                           sd = sd.alternative,
                           null.mean = mean.null,
                           null.sd = sd.null,
                           z.alpha = z.alpha,
                           power = power,
                           n = n),
                      class = c("pwrss", "z", "mediation")))

} # end of power.z.mediation()

power.z.med <- power.z.mediation


pwrss.z.mediation  <- function(a, b, cp = 0,
                               sdx = 1, sdm = 1, sdy = 1,
                               r2m.x = a ^ 2 * sdx ^ 2 / sdm ^ 2,
                               r2y.mx = (b ^ 2 * sdm ^ 2 + cp ^ 2 * sdx ^ 2) / sdy ^ 2,
                               n = NULL, power = NULL,
                               alpha = 0.05, alternative = c("not equal", "less", "greater"),
                               mc = TRUE, nsims = 1000, ndraws = 1000,
                               verbose = TRUE) {

  alternative <- tolower(match.arg(alternative))
  if (alternative %in% c("less", "greater")) alternative <- "one.sided"
  if (alternative == "not equal") alternative <- "two.sided"

  if (is.null(power)) {

    sobel.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                   sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                   r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                   n = n, power = power, alpha = alpha,
                                   alternative = alternative,
                                   method = "sobel",
                                   n.simulation = nsims,
                                   n.draws = ndraws,
                                   ceiling = TRUE,
                                   verbose = FALSE)

    aroian.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                    sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                    r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                    n = n, power = power, alpha = alpha,
                                    alternative = alternative,
                                    method = "aroian",
                                    n.simulation = nsims,
                                    n.draws = ndraws,
                                    ceiling = TRUE,
                                    verbose = FALSE)

    goodman.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                     sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                     r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                     n = n, power = power, alpha = alpha,
                                     alternative = alternative,
                                     method = "goodman",
                                     n.simulation = nsims,
                                     n.draws = ndraws,
                                     ceiling = TRUE,
                                     verbose = FALSE)

    joint.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                   sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                   r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                   n = n, power = power, alpha = alpha,
                                   alternative = alternative,
                                   method = "joint",
                                   n.simulation = nsims,
                                   n.draws = ndraws,
                                   ceiling = TRUE,
                                   verbose = FALSE)

    monte.carlo.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                         sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                         r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                         n = n, power = power, alpha = alpha,
                                         alternative = alternative,
                                         method = "monte.carlo",
                                         n.simulation = nsims,
                                         n.draws = ndraws,
                                         ceiling = TRUE,
                                         verbose = FALSE)

    power.out <- data.frame(rbind(
      c(sobel.obj$mean, sobel.obj$n, sobel.obj$power),
      c(aroian.obj$mean, aroian.obj$n, aroian.obj$power),
      c(goodman.obj$mean, goodman.obj$n, goodman.obj$power),
      c(joint.obj$mean, joint.obj$n, joint.obj$power),
      c(monte.carlo.obj$mean, monte.carlo.obj$n, monte.carlo.obj$power)
    ))

    colnames(power.out) <- c("non-centrality", "n", "power")
    power.out$method <- c("Sobel", "Aroian", "Goodman", "Joint", "Monte Carlo")

  } # power

  if (is.null(n)) {

    sobel.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                   sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                   r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                   n = n, power = power, alpha = alpha,
                                   alternative = alternative,
                                   method = "sobel",
                                   n.simulation = nsims,
                                   n.draws = ndraws,
                                   ceiling = TRUE,
                                   verbose = FALSE)

    aroian.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                    sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                    r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                    n = n, power = power, alpha = alpha,
                                    alternative = alternative,
                                    method = "aroian",
                                    n.simulation = nsims,
                                    n.draws = ndraws,
                                    ceiling = TRUE,
                                    verbose = FALSE)

    goodman.obj <- power.z.mediation(beta.a = a, beta.b = b, beta.cp = cp,
                                     sd.predictor = sdx, sd.mediator = sdm, sd.outcome = sdy,
                                     r.squared.mediator = r2m.x, r.squared.outcome = r2y.mx,
                                     n = n, power = power, alpha = alpha,
                                     alternative = alternative,
                                     method = "goodman",
                                     n.simulation = nsims,
                                     n.draws = ndraws,
                                     ceiling = TRUE,
                                     verbose = FALSE)

    power.out <- data.frame(rbind(
      c(sobel.obj$mean, sobel.obj$n, sobel.obj$power),
      c(aroian.obj$mean, aroian.obj$n, aroian.obj$power),
      c(goodman.obj$mean, goodman.obj$n, goodman.obj$power),
      c(NA, NA, NA),
      c(NA, NA, NA)
    ))

    colnames(power.out) <- c("non-centrality", "n", "power")
    power.out$method <- c("Sobel", "Aroian", "Goodman", "Joint", "Monte Carlo")

  } # n


  if (is.null(power)) {
    ncp.vec <- c(sobel = power.out$mean[1], aroian = power.out$mean[2], goodman = power.out$mean[3], joint = power.out$mean[4], mc = power.out$mean[5])
    n.vec <- c(sobel = power.out$n[1], aroian = power.out$n[2], goodman = power.out$n[3], joint = power.out$n[4], mc = power.out$n[5])
    power.vec <- c(sobel = power.out$power[1], aroian = power.out$power[2], goodman = power.out$power[3], joint = power.out$power[4], mc = power.out$power[5])
  } else {
    ncp.vec <- c(sobel = power.out$mean[1], aroian = power.out$mean[2], goodman = power.out$mean[3], joint = power.out$mean[4], mc = power.out$mean[5])
    n.vec <- c(sobel = power.out$n[1], aroian = power.out$n[2], goodman = power.out$n[3], joint = power.out$n[4], mc = power.out$n[5])
    power.vec <- c(sobel = power.out$power[1], aroian = power.out$power[2], goodman = power.out$power[3], joint = power.out$power[4], mc = power.out$power[5])
  }

  if (verbose) {
    cat(" Indirect Effect in a Mediation Model",
        "\n ====================================\n",
        sep = "")

    ifelse(is.null(power),
           print(power.out[, c("method", "non-centrality", "n", "power")], row.names = FALSE),
           print(power.out[1:3, c("method", "non-centrality", "n", "power")], row.names = FALSE))

    cat(" ------------------------------------\n",
        " Type 1 error rate: ", round(alpha, 3), "\n\n",
        sep = "")

    # cat2("beta.indirect := beta.a * beta.b \n", "green")

  } # verbose

# cat("This function will be removed in the future. \n Please use power.z.mediation() function. \n")

invisible(structure(list(parms = list(a = a, b = b, cp = cp,
                                      sdx = sdx, sdm = sdm, sdy = sdy,
                                      r2m.x = r2m.x, r2y.mx = r2y.mx,
                                      alpha = alpha, alternative = alternative, verbose = verbose),
                         test = c("z", "joint", "monte.carlo"),
                         ncp = ncp.vec,
                         power = power.vec,
                         n = n.vec),
                    class = c("pwrss", "z", "med", "defunct")))

} # pwrss.z.mediation()
pwrss.z.med  <- pwrss.z.mediation
