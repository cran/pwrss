power.exact.mcnemar <- function(prob10, prob01, n.paired = NULL,
                                power = NULL,  alpha = 0.05,
                                alternative = c("two.sided", "one.sided"),
                                method = c("exact", "approximate"),
                                ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  check.proportion(prob10, prob01, alpha)
  check.logical(ceiling)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n.paired)) check.sample.size(n.paired)

  method <- tolower(match.arg(method))
  alternative <- tolower(match.arg(alternative))

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  pwr.exact <- function(prob10, prob01, n.paired, alpha, alternative) {

    OR <- prob10 / prob01

    if (alternative == "two.sided") prob <- 1 / (1 + OR)
    if (alternative == "one.sided" && prob10 < prob01) prob <- min(1 / (1 + OR), OR / (1 + OR))
    if (alternative == "one.sided" && prob10 > prob01) prob <- max(1 / (1 + OR), OR / (1 + OR))

    prod1 <- dbinom(x = seq(0, ceiling(n.paired)), size = ceiling(n.paired), prob = prob01 + prob10)
    prod2 <- power.binom.test(prob = prob, null.prob = 0.50,
                              size = seq(0, ceiling(n.paired)), alpha = alpha,
                              alternative = alternative, plot = FALSE, verbose = FALSE)$power
    power <- sum(prod1 * prod2)

    power

  } # pwr.exact()

  ss.exact <- function(prob10, prob01, power, alpha, alternative) {

    achieved.power <- 0
    n.paired <-  ss.approx(prob10, prob01, power, alpha, alternative)

    while (achieved.power < power && n.paired < 1e5) {

      achieved.power <- pwr.exact(prob10 = prob10, prob01 = prob01,
                                    n.paired = n.paired, alpha = alpha,
                                    alternative = alternative)

      if (achieved.power < power) n.paired <- n.paired + 1

    }

    if (n.paired > 1e5) stop("Sample size exceeds 100,000. Check assumptions.", call. = FALSE)

    n.paired

  } #  ss.exact()

  pwr.approx <- function(prob10, prob01, n.paired, alpha, alternative) {

    OR <- prob10 / prob01
    PD <- prob10 + prob01

    # Machin, Campbell, Fayers, and Pinol (1997)
    if (alternative == "two.sided") alpha <- alpha / 2
    z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1, lower.tail = TRUE)
    z.beta <- (sqrt((OR - 1) ^ 2 * PD * n.paired) - z.alpha * (1 + OR)) / sqrt((OR + 1) ^ 2 - (OR - 1) ^ 2 * PD)
    power <- pnorm(z.beta, mean = 0, sd = 1, lower.tail = TRUE)

    mean.alternative <- z.alpha + z.beta
    sd.alternative <- 1
    mean.null <- 0
    sd.null <- 1
    z.alpha <- z.alpha

    if (OR < 1) {
      mean.alternative <- -mean.alternative
      if (alternative == "one.sided") z.alpha <- -z.alpha
    }
    if (alternative == "two.sided") z.alpha <- c(-z.alpha, z.alpha)

    list(power = power,
         mean.alternative =  mean.alternative,
         sd.alternative = sd.alternative,
         mean.null = mean.null,
         sd.null = sd.null,
         z.alpha = z.alpha)

  } # pwr.approx()

  ss.approx <- function(prob10, prob01, power, alpha, alternative) {

    OR <- prob10 / prob01
    PD <- prob10 + prob01
    beta <- 1 - power

    # Machin, Campbell, Fayers, and Pinol (1997)
    if (alternative == "two.sided") alpha <- alpha / 2
    z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1, lower.tail = TRUE)
    z.beta <- qnorm(1 - beta, mean = 0, sd = 1, lower.tail = TRUE)
    n.paired <- (z.alpha * (1 + OR) + z.beta * sqrt((OR + 1) ^ 2 - (OR - 1) ^ 2 * PD)) ^ 2 / ((OR - 1) ^ 2 * PD)

    n.paired

  } # ss.approx()

  # for reasonable round-up
  prob11 <- (1 - (prob10 + prob01)) / 2
  prob00 <- (1 - (prob10 + prob01)) / 2

  # method
  if (method == "exact") {

    if (is.null(power)) {

      power <- pwr.exact(prob10 = prob10, prob01 = prob01,
                         n.paired = n.paired, alpha = alpha,
                         alternative = alternative)

      n.discordant <- sum(n.paired * c(prob10, prob01))

      if (ceiling) {
        n.discordant <- sum(ceiling(n.paired * c(prob10, prob01)))
      }

      mean.alternative <- NA
      sd.alternative <- NA
      mean.null <- NA
      sd.null <- NA
      z.alpha <- NA

    }

    if (is.null(n.paired)) {

      n.paired <- ss.exact(prob10 = prob10, prob01 = prob01,
                           power = power, alpha = alpha,
                           alternative = alternative)
      n.discordant <- sum(n.paired * c(prob10, prob01))

      if (ceiling) {
        n.paired <- sum(ceiling(n.paired * c(prob11, prob10, prob01, prob00)))
        n.discordant <- sum(ceiling(n.paired * c(prob10, prob01)))
      }

      power <- pwr.exact(prob10 = prob10, prob01 = prob01,
                         n.paired = n.paired, alpha = alpha,
                         alternative = alternative)

      mean.alternative <- NA
      sd.alternative <- NA
      mean.null <- NA
      sd.null <- NA
      z.alpha <- NA

    }

  }  else if (method == "approximate") {

    if (is.null(power)) {

      pwr.obj <- pwr.approx(prob10 = prob10, prob01 = prob01,
                            n.paired = n.paired, alpha = alpha,
                            alternative = alternative)

      power <- pwr.obj$power
      mean.alternative <- pwr.obj$mean.alternative
      sd.alternative <- pwr.obj$sd.alternative
      mean.null <- pwr.obj$mean.null
      sd.null <- pwr.obj$sd.null
      z.alpha <- pwr.obj$z.alpha

      n.discordant <- sum(n.paired * c(prob10, prob01))

      if (ceiling) {
        n.discordant <- sum(ceiling(n.paired * c(prob10, prob01)))
      }

    }

    if (is.null(n.paired)) {

      n.paired <- ss.approx(prob10 = prob10, prob01 = prob01,
                            power = power, alpha = alpha,
                            alternative = alternative)
      n.discordant <- sum(n.paired * c(prob10, prob01))

      if (ceiling) {
        n.paired <- sum(ceiling(n.paired * c(prob11, prob10, prob01, prob00)))
        n.discordant <- sum(ceiling(n.paired * c(prob10, prob01)))
      }

      pwr.obj <- pwr.approx(prob10 = prob10, prob01 = prob01,
                            n.paired = n.paired, alpha = alpha,
                            alternative = alternative)

      power <- pwr.obj$power
      mean.alternative <- pwr.obj$mean.alternative
      sd.alternative <- pwr.obj$sd.alternative
      mean.null <- pwr.obj$mean.null
      sd.null <- pwr.obj$sd.null
      z.alpha <- pwr.obj$z.alpha

    }

  }  # method

  # efficient non-centrality parameter for the chi-squared dist
  # f10 <- prob10 * n.paired
  # f01 <- prob01 * n.paired
  # ncp.alternative <- (f10 - f01) ^ 2 / (f10 + f01)
  # df <- 1

  # critical values for the binomial approach
  prob <- prob10 / (prob10 + prob01)
  null.prob <- 0.50
  size <- ceiling(n.discordant)
  if (alternative == "one.sided") {
    if (prob < null.prob) {
      # less
      q.binom.alpha <- qbinom(alpha, size = size, prob = null.prob, lower.tail = TRUE)
      p.binom.alpha <- pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = TRUE)
      if (p.binom.alpha > alpha) q.binom.alpha <- q.binom.alpha - 1
      if (method == "exact") {
        approx.alpha <- pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = TRUE)
      } else {
        approx.alpha <- alpha
      }
    } else {
      # greater
      q.binom.alpha <- qbinom(alpha, size = size, prob = null.prob, lower.tail = FALSE)
      p.binom.alpha <- pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = FALSE)
      if (p.binom.alpha > alpha) q.binom.alpha <- q.binom.alpha + 1
      if (method == "exact") {
        approx.alpha <- pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = FALSE)
      } else {
        approx.alpha <- alpha
      }
    }
  }

  if (alternative == "two.sided") {
    q.binom.lower <- qbinom(alpha / 2, size = size, prob = null.prob, lower.tail = TRUE)
    p.binom.lower <- pbinom(q.binom.lower, size = size, prob = null.prob, lower.tail = TRUE)
    if (p.binom.lower > alpha / 2) q.binom.lower <- q.binom.lower - 1
    p.binom.lower <- pbinom(q.binom.lower, size = size, prob = null.prob, lower.tail = TRUE)

    q.binom.upper <- qbinom(alpha / 2, size = size, prob = null.prob, lower.tail = FALSE)
    p.binom.upper <- pbinom(q.binom.upper, size = size, prob = null.prob, lower.tail = FALSE)
    if (p.binom.upper > alpha / 2) q.binom.upper <- q.binom.upper + 1
    p.binom.upper <- pbinom(q.binom.upper, size = size, prob = null.prob, lower.tail = FALSE)

    q.binom.alpha <- c(q.binom.lower, q.binom.upper)
    if (method == "exact") {
      approx.alpha <- p.binom.upper + p.binom.lower
    } else {
      approx.alpha <- alpha
    }
  }


  ifelse(method == "exact",
         class <- c("pwrss", "exact", "mcnemar"),
         class <- c("pwrss", "z", "twoprops"))

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
                      test = "Paired Proportions",
                      alpha = approx.alpha,
                      alt = alternative,
                      method = ifelse(method == "exact", "exact", "z"),
                      delta = prob10 - prob01,
                      odds.ratio = prob10 / prob01,
                      size = n.discordant,
                      prob.alternative = prob10 / (prob10 + prob01),
                      prob.null = 0.50,
                      binom.alpha = q.binom.alpha,
                      mean.alternative = mean.alternative,
                      sd.alternative = sd.alternative,
                      mean.null = mean.null,
                      sd.null = sd.null,
                      z.alpha = z.alpha,
                      power = power,
                      n.paired = n.paired)

    if (pretty) {
      .print.pwrss.mcnemar(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.mcnemar(print.obj, verbose = verbose)
    }

  } # verbose

  invisible(structure(list(parms = list(prob10 = prob10, prob01 = prob01,
                                        alpha = approx.alpha, alternative = alternative,
                                        method = method, ceiling = ceiling,
                                        verbose = verbose, pretty = pretty),
                           test = ifelse(method == "exact", "exact", "z"),
                           delta = prob10 - prob01,
                           odds.ratio = prob10 / prob01,
                           size = n.discordant,
                           prob = prob10 / (prob10 + prob01),
                           null.prob = 0.50,
                           binom.alpha = q.binom.alpha,
                           mean = mean.alternative,
                           sd = sd.alternative,
                           null.mean = mean.null,
                           null.sd = sd.null,
                           z.alpha = z.alpha,
                           alpha = approx.alpha,
                           power = power,
                           n.paired = n.paired),
                      class = class))

} # power.exact.mcnemar()

power.exact.twoprops.mcnemar <- power.exact.mcnemar
