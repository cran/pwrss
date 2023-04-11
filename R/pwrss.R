#########################
# one proportion z test #
#########################

pwrss.z.prop <- function (p, p0 = 0, margin = 0, arcsin.trans = FALSE,
                          alpha = 0.05,
                          alternative = c("not equal", "greater", "less",
                                          "equivalent", "non-inferior", "superior"),
                          n = NULL, power = NULL, verbose = TRUE)
{

  alternative <- match.arg(alternative)

  if (is.null(n) & is.null(power)) stop("`n` and `power` cannot be `NULL` at the same time", call. = FALSE)
  if (!is.null(n) & !is.null(power)) stop("one of the `n` or `power` should be `NULL`", call. = FALSE)
  if ((alternative == "not equal" | alternative == "not equal" | alternative == "not equal") & margin != 0) warning("`margin` argument is ignored", call. = FALSE)

  if(arcsin.trans) {
    var.num <- 1
    h <- switch(alternative,
                `not equal` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)),
                `greater` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)),
                `less` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)),
                `non-inferior` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0 + margin)),
                `superior` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0 + margin)),
                `equivalent` = c(2*asin(sqrt(p)) - 2*asin(sqrt(p0 + margin)), 2*asin(sqrt(p)) - 2*asin(sqrt(p0 - margin))))
    if(verbose) cat(" Approach: Arcsine Transformation \n")
  } else {
    var.num <- p * (1 - p)
    h <- switch(alternative,
                `not equal` = p - p0,
                `greater` = p - p0,
                `less` = p - p0,
                `non-inferior` = p - p0 - margin,
                `superior` = p - p0 - margin,
                `equivalent` = c(p - p0 - margin, p - p0 + margin))
    if(verbose) cat(" Approach: Normal Approximation \n")
  }

  if (alternative == "not equal") {

    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha / 2, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n <- M^2 * var.num / h^2
      lambda <- h / sqrt(var.num / n)
    }
    if (is.null(power)) {
      lambda <- h / sqrt(var.num / n)
      power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), lambda) +
        pnorm(-qnorm(alpha / 2, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "greater" | alternative ==
           "less") {


    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n <- M^2 * var.num / h^2
      lambda <- h / sqrt(var.num / n)

    }
    if (is.null(power)) {
      lambda <- h / sqrt(var.num / n)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)

    }

  }
  else if (alternative == "non-inferior" | alternative ==
           "superior") {

    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n <- M^2 * var.num / h^2
      lambda <- h / sqrt(var.num / n)

    }
    if (is.null(power)) {
      lambda <- h / sqrt(var.num / n)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  } else if (alternative == "equivalent") {

    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta / 2, lower.tail = FALSE)
      # h <- min(abs(h))
      n <- M^2 * var.num / h^2
      n <- max(n)
      lambda <- h / sqrt(var.num / n)
    }

    if (is.null(power)) {
      lambda <- h / sqrt(var.num / n)
      power <- 1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = abs(lambda[1])) +
        1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = abs(lambda[2])) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }

  } else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }

  if(verbose) {
    cat(" A Proportion against a Constant (z Test) \n",
        switch(alternative,
               `not equal` = "H0: p = p0 \n HA: p != p0 \n",
               `greater` = "H0: p = p0 \n HA: p > p0 \n",
               `less` = "H0: p = p0 \n HA: p < p0 \n",
               `non-inferior` = "H0: p - p0 <= margin \n HA: p - p0 > margin \n",
               `superior` = "H0: p - p0 <= margin \n HA: p - p0 > margin \n",
               `equivalent` = "H0: |p - p0| >= margin \n HA: |p - p0| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(lambda, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(call = call,
                           parms = list(p = p, p0 = p0, arcsin.trans = arcsin.trans,
                                        alpha = alpha, margin = margin, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = lambda,
                           power = power,
                           n = n),
                      class = c("pwrss", "z", "prop")))
}

##########################
# two proportions z test #
##########################

pwrss.z.2props <- function (p1, p2, margin = 0, arcsin.trans = FALSE,
                            kappa = 1, alpha = 0.05,
                            alternative = c("not equal", "greater", "less",
                                            "equivalent", "non-inferior", "superior"),
                            n2 = NULL, power = NULL, verbose = TRUE)
{

  alternative <- match.arg(alternative)

  if (is.null(n2) & is.null(power)) stop("`n2` and `power` cannot be `NULL` at the same time", call. = FALSE)
  if (!is.null(n2) & !is.null(power)) stop("one of the `n2` or `power` should be `NULL`", call. = FALSE)
  if ((alternative == "not equal" | alternative == "not equal" | alternative == "not equal") & margin != 0)
    warning("`margin` argument is ignored", call. = FALSE)
  # if (alternative == "superior" & margin < 0) warning("expecting `margin > 0`", call. = FALSE)
  # if (alternative == "non-inferior" & margin > 0) warning("expecting `margin < 0`", call. = FALSE)
  # if (alternative == "greater" & (p1 < p2)) stop("alternative = 'greater' but p1 < p2", call. = FALSE)
  # if (alternative == "less" & (p1 > p2)) stop("alternative = 'less' but p1 > p2", call. = FALSE)

  if(arcsin.trans) {
    h <- switch(alternative,
                `not equal` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)),
                `greater` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)),
                `less` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)),
                `non-inferior` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2 + margin)),
                `superior` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2 + margin)),
                `equivalent` = c(2*asin(sqrt(p1)) - 2*asin(sqrt(p2 + margin)),
                                 2*asin(sqrt(p1)) - 2*asin(sqrt(p2 - margin))))

    # 2*asin(sqrt(p1)) - 2*asin(sqrt(p2 + margin)) is same as
    # 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)) - (2*asin(sqrt(p2 + margin)) - 2*asin(sqrt(p2)))
    if(verbose) cat(" Approach: Arcsine Transformation \n", "Cohen's h =", round(h,3), "\n")
  } else {
    h <- switch(alternative,
                `not equal` = p1 - p2,
                `greater` = p1 - p2,
                `less` = p1 - p2,
                `non-inferior` = p1 - p2 - margin,
                `superior` = p1 - p2 - margin,
                `equivalent` = c(p1 - p2 - margin,
                                 p1 - p2 + margin))
    if(verbose) cat(" Approach: Normal Approximation \n")
  }

  if (alternative == "not equal") {

    if (is.null(n2)) {
      beta <- 1 - power
      # kappa = n1/n2
      M <- qnorm(alpha / 2, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      if(arcsin.trans) {
        n2 <- (M^2 / h^2) * (kappa + 1) / kappa
        n1 <- kappa * n2
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        n1 = M^2 * (p1 * (1 - p1) + p2 * (1 - p2) * kappa) / h^2
        n2 <- n1 / kappa
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }

    }
    if (is.null(power)) {
      # kappa = n1/n2
      n1 <- kappa*n2

      if(arcsin.trans) {
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }

      power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), abs(lambda)) +
        pnorm(-qnorm(alpha / 2, lower.tail = FALSE), abs(lambda))
    }

  }
  else if (alternative == "greater" | alternative == "less") {

    if (is.null(n2)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      if(arcsin.trans) {
        n2 <- (M^2 / h^2) * (kappa + 1) / kappa
        n1 <- kappa * n2
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        n1 = M^2 * (p1 * (1 - p1) + p2 * (1 - p2) * kappa) / h^2
        n2 <- n1 / kappa
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }
    }
    if (is.null(power)) {
      # kappa = n1/n2
      n1 <- kappa*n2

      if(arcsin.trans) {
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }

      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), abs(lambda))
    }

  }
  else if (alternative == "non-inferior" | alternative == "superior") {

    if (is.null(n2)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      if(arcsin.trans) {
        n2 <- (M^2 / h^2) * (kappa + 1) / kappa
        n1 <- kappa * n2
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        n1 = M^2 * (p1 * (1 - p1) + p2 * (1 - p2) * kappa) / h^2
        n2 <- n1 / kappa
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }
    }
    if (is.null(power)) {
      # kappa = n1/n2
      n1 <- kappa*n2

      if(arcsin.trans) {
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }

      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), abs(lambda))
    }

  }
  else if (alternative == "equivalent") {

    if (is.null(n2)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta / 2, lower.tail = FALSE)
      # h <- min(abs(h))
      if(arcsin.trans) {
        n2 <- (M^2 / h^2) * (kappa + 1) / kappa
        n2 <- max(n2)
        n1 <- kappa * n2
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        n1 = M^2 * (p1 * (1 - p1) + p2 * (1 - p2) * kappa) / h^2
        n1 <- max(n1)
        n2 <- n1 / kappa
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }
    }
    if (is.null(power)) {
      # kappa = n1/n2
      n1 <- kappa*n2

      if(arcsin.trans) {
        lambda = h / sqrt(1 / n1 + 1 / n2)
      } else {
        lambda = h / sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
      }

      power <- 1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = abs(lambda[1])) +
        1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = abs(lambda[2])) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }

  }
  else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }

  if(verbose) {
    cat(" Difference between Two Proportions \n (Independent Samples z Test) \n",
        switch(alternative,
               `not equal` = "H0: p1 = p2 \n HA: p1 != p2 \n",
               `greater` = "H0: p1 = p2 \n HA: p1 > p2 \n",
               `less` = "H0: p1 = p2 \n HA: p1 < p2 \n",
               `non-inferior` = "H0: p1 - p2 <= margin \n HA: p1 - p2 > margin \n",
               `superior` = "H0: p1 - p2 <= margin \n HA: p1 - p2 > margin \n",
               `equivalent` = "H0: |p1 - p2| >= margin \n HA: |p1 - p2| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n1 =", ceiling(n1), "\n",
        " n2 =", ceiling(n2), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(lambda, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(p1 = p1, p2 = p2, kappa = kappa, arcsin.trans = arcsin.trans,
                                        alpha = alpha, margin = margin, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = lambda,
                           power = power,
                           n = c(n1 = n1, n2 = n2)),
                      class = c("pwrss", "z", "2props")))
}

###################
# one mean z test #
###################

pwrss.z.mean <- function (mu, sd = 1, mu0 = 0, margin = 0, alpha = 0.05,
                          alternative = c("not equal", "greater", "less",
                                          "equivalent", "non-inferior", "superior"),
                          n = NULL, power = NULL, verbose = TRUE)
{

  if (length(alternative) > 1)
    alternative <- alternative[1]
  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)
  if (!is.null(n) & is.null(power))
    requested <- "power"
  if (is.null(n) & !is.null(power))
    requested <- "n"
  if (alternative == "not equal") {
    if (margin != 0)
      warning("`margin` argument is ignored")

    HA_H0 <- mu - mu0
    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha / 2, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n = M^2 * sd^2 / (HA_H0)^2
    }
    if (is.null(power)) {
      lambda = (HA_H0) / sqrt(sd^2 / n)
      power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), lambda) +
        pnorm(-qnorm(alpha / 2, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "greater" | alternative ==
           "less") {
    if (margin != 0)
      warning("`margin` argument is ignored")

    if (alternative == "greater" & (mu < mu0))
      stop("alternative = 'greater' but mu < mu0",
           call. = FALSE)
    if (alternative == "less" & (mu > mu0))
      stop("alternative = 'less' but mu > mu0", call. = FALSE)

    HA_H0 <- mu - mu0
    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n = M^2 * sd^2 / (HA_H0)^2
    }
    if (is.null(power)) {
      lambda = (HA_H0) / sqrt(sd^2 / n)
      if(alternative == "less") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "non-inferior" | alternative ==
           "superior") {

    HA_H0 <- mu - mu0 - margin
    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n = M^2 * sd^2 / (HA_H0)^2
    }
    if (is.null(power)) {
      lambda = (HA_H0) / sqrt(sd^2 / n)
      if(alternative == "non-inferior") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "equivalent") {

    HA_H0 <- abs(mu - mu0) - margin
    if (is.null(n)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta / 2, lower.tail = FALSE)
      n = M^2 * sd^2 / (HA_H0)^2
    }
    if (is.null(power)) {
      lambda = (HA_H0) / sqrt(sd^2 / n)
      power <- 2*(1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }

  }
  else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }

  ncp <- (HA_H0) / sqrt(sd^2 / n)
  hypothesis <- alternative

  if(verbose) {
    cat(" A Mean against a Constant (z Test) \n",
        switch(hypothesis,
               `not equal` = "H0: mu = mu0 \n HA: mu != mu0 \n",
               `greater` = "H0: mu = mu0 \n HA: mu > mu0 \n",
               `less` = "H0: mu = mu0 \n HA: mu < mu0 \n",
               `non-inferior` = "H0: mu - mu0 <= margin \n HA: mu - mu0 > margin \n",
               `superior` = "H0: mu - mu0 <= margin \n HA: mu - mu0 > margin \n",
               `equivalent` = "H0: |mu - mu0| >= margin \n HA: |mu - mu0| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(mu = mu, sd = sd, mu0 = mu0,
                                        alpha = alpha, margin = margin, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = ncp,
                           power = power,
                           n = n),
                      class = c("pwrss", "z", "mean")))
}

###################
# one mean t test #
###################

pwrss.t.mean <- function (mu, sd = 1, mu0 = 0, margin = 0, alpha = 0.05,
                          alternative = c("not equal", "greater", "less",
                                          "equivalent", "non-inferior", "superior"),
                          n = NULL, power = NULL, verbose = TRUE)
{

  if (length(alternative) > 1)
    alternative <- alternative[1]
  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)
  if (!is.null(n) & is.null(power))
    requested <- "power"
  if (is.null(n) & !is.null(power))
    requested <- "n"

  if (alternative == "not equal") {
    if (margin != 0)
      warning("`margin` argument is ignored")

    HA_H0 <- mu - mu0

    pwr.fun.body <- quote({
      lambda = (HA_H0) / sqrt(sd^2 / n)
      power <- 1 - pt(qt(alpha / 2, df = n - 1, lower.tail = FALSE), df = n - 1, ncp = abs(lambda)) +
        pt(-qt(alpha / 2, df = n - 1, lower.tail = FALSE), df = n - 1, ncp = abs(lambda))
    })

  }
  else if (alternative == "greater" | alternative ==
           "less") {
    if (margin != 0)
      warning("`margin` argument is ignored")
    if (alternative == "greater" & (mu < mu0))
      stop("alternative = 'greater' but mu < mu0",
           call. = FALSE)
    if (alternative == "less" & (mu > mu0))
      stop("alternative = 'less' but mu > mu0", call. = FALSE)

    HA_H0 <- mu - mu0

    pwr.fun.body <- quote({
      lambda = (HA_H0) / sqrt(sd^2 / n)
      power <- 1 - pt(qt(alpha, df = n - 1, lower.tail = FALSE), df = n - 1, ncp = abs(lambda))
    })

  }
  else if (alternative == "non-inferior" | alternative ==
           "superior") {

    HA_H0 <- mu - mu0 - margin

    pwr.fun.body <- quote({
      lambda = (HA_H0) / sqrt(sd^2 / n)
      power <- 1 - pt(qt(alpha, df = n - 1, lower.tail = FALSE), df = n - 1, ncp = abs(lambda))
    })

  }
  else if (alternative == "equivalent") {

    HA_H0 <- abs(mu - mu0) - margin

    pwr.fun.body <- quote({
      lambda = (HA_H0) / sqrt(sd^2 / n)
      power <- 2*(1 - pt(qt(alpha, df = n - 1, lower.tail = FALSE), df = n - 1, ncp = abs(lambda))) - 1
    })

  } else {

    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)

  }

  if (is.null(power)) {

    power <- eval(pwr.fun.body)
    if(power < 0) power <- 0

  } else if(is.null(n)) {

    n <- uniroot(function(n){
      power - eval(pwr.fun.body)
    }, interval = c(2, 1e+09))$root

    lambda <- (HA_H0) / sqrt(sd^2 / n)

  }

  n <- ceiling(n)
  ncp <- (HA_H0) / sqrt(sd^2 / n)
  hypothesis <- alternative

  if(verbose) {
    cat(" A Mean against a Constant (t Test) \n",
        switch(hypothesis,
               `not equal` = "H0: mu = mu0 \n HA: mu != mu0 \n",
               `greater` = "H0: mu = mu0 \n HA: mu > mu0 \n",
               `less` = "H0: mu = mu0 \n HA: mu < mu0 \n",
               `non-inferior` = "H0: mu - mu0 <= margin \n HA: mu - mu0 > margin \n",
               `superior` = "H0: mu - mu0 <= margin \n HA: mu - mu0 > margin \n",
               `equivalent` = "H0: |mu - mu0| >= margin \n HA: |mu - mu0| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Degrees of freedom =", round(ceiling(n) - 1, 3),"\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(mu = mu, sd = sd, mu0 = mu0,
                                        alpha = alpha, margin = margin, alternative = alternative, verbose = verbose),
                           test = "t",
                           ncp = ncp,
                           df = n - 1,
                           power = power,
                           n = n),
                      class = c("pwrss", "t", "mean")))
}

####################
# two means z test #
####################

pwrss.z.2means <- function (mu1, mu2 = 0, sd1 = 1, sd2 = sd1, margin = 0,
                            kappa = 1, alpha = 0.05,
                            alternative = c("not equal", "greater", "less",
                                            "equivalent", "non-inferior", "superior"),
                            n2 = NULL, power = NULL, verbose = TRUE)
{

  # sd1 <- sd2 <- sd # pooled standard devation
  if (length(alternative) > 1)
    alternative <- alternative[1]
  if (is.null(n2) & is.null(power))
    stop("`n2` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n2) & !is.null(power))
    stop("one of the `n2` or `power` should be `NULL`",
         call. = FALSE)
  if (!is.null(n2) & is.null(power))
    requested <- "power"
  if (is.null(n2) & !is.null(power))
    requested <- "n2"
  if (alternative == "not equal") {
    if (margin != 0)
      warning("`margin` argument is ignored")

    HA_H0 <- mu1 - mu2
    if (is.null(n2)) {
      beta <- 1 - power
      # kappa = n1/n2
      M <- qnorm(alpha / 2, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
      n1 <- n2 * kappa
    }
    if (is.null(power)) {
      n1 <- n2 * kappa
      lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), lambda) +
        pnorm(-qnorm(alpha / 2, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "greater" | alternative ==
           "less") {
    if (margin != 0)
      warning("`margin` argument is ignored")
    if (alternative == "greater" & (mu1 < mu2))
      stop("alternative = 'greater' but mu1 < mu2",
           call. = FALSE)
    if (alternative == "less" & (mu1 > mu2))
      stop("alternative = 'less' but mu1 > mu2",
           call. = FALSE)

    HA_H0 <- mu1 - mu2
    if (is.null(n2)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
      n1 <- n2 * kappa
    }
    if (is.null(power)) {
      n1 <- n2 * kappa
      lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      if(alternative == "less") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "non-inferior" | alternative ==
           "superior") {
    if (alternative == "superior" & margin < 0)
      warning("expecting 'margin > 0' when mu1 - mu2 > 0", call. = FALSE)
    if (alternative == "non-inferior" & margin > 0)
      warning("expecting 'margin < 0' when mu1 - mu2 > 0", call. = FALSE)

    HA_H0 <- mu1 - mu2 - margin
    if (is.null(n2)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
      n1 <- n2 * kappa
    }
    if (is.null(power)) {
      n1 <- n2 * kappa
      lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      if(alternative == "non-inferior") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "equivalent") {

    HA_H0 <- abs(mu1 - mu2) - margin
    if (is.null(n2)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta / 2, lower.tail = FALSE)
      n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
      n1 <- n2 * kappa
    }
    if (is.null(power)) {
      n1 <- n2 * kappa
      lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      power <- 2*(1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }

  }
  else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }
  n1 <- kappa * n2
  ncp <- (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
  hypothesis <- alternative

  if(verbose) {
    cat(" Difference between Two Means (Independent Samples z Test) \n",
        switch(hypothesis,
               `not equal` = "H0: mu1 = mu2 \n HA: mu1 != mu2 \n",
               `greater` = "H0: mu1 = mu2 \n HA: mu1 > mu2 \n",
               `less` = "H0: mu1 = mu2 \n HA: mu1 < mu2 \n",
               `non-inferior` = "H0: mu1 - mu2 <= margin \n HA: mu1 - mu2 > margin \n",
               `superior` = "H0: mu1 - mu2 <= margin \n HA: mu1 - mu2 > margin \n",
               `equivalent` = "H0: |mu1 - mu2| >= margin \n HA: |mu1 - mu2| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n1 =", ceiling(n1), "\n",
        " n2 =", ceiling(n2), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, kappa = kappa, alpha = alpha,
                                        margin = margin, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = ncp,
                           power = power,
                           n = c(n1 = n1, n2 = n2)),
                      class = c("pwrss", "z", "2means")))
}

####################
# two means t test #
####################

# use welch.df = TRUE for unequal variances and unequal n in independent samples t test
# allows r between pre and post to be different from 0, pwr.t.test(...,paired = TRUE) assumes paired.r = 0
# allows diferent sd for pre and post
# default specification assumes standardized means (or mean difference)
# margin is defined as the minimum mu1 - mu2 that is practically relevant
# suggest margin = t_critical * SE by default?
pwrss.t.2means <- function (mu1, mu2 = 0, margin = 0,
                            sd1 = ifelse(paired, sqrt(1/(2*(1-paired.r))), 1), sd2 = sd1,
                            kappa = 1, paired = FALSE, paired.r = 0.50,
                            alpha = 0.05, welch.df = FALSE,
                            alternative = c("not equal", "greater", "less",
                                            "equivalent", "non-inferior", "superior"),
                            n2 = NULL, power = NULL, verbose = TRUE)
{

  welch_df <- function(sd1, sd2, n1, n2) {
    (sd1^2 / n1 + sd2^2 / n2)^2 /
      (sd1^4 / (n1^2 * (n1 - 1)) + sd2^4 / (n2^2 * (n2 - 1)))
  }

  if (length(alternative) > 1)
    alternative <- alternative[1]
  if (is.null(n2) & is.null(power))
    stop("`n2` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n2) & !is.null(power))
    stop("one of the `n2` or `power` should be `NULL`",
         call. = FALSE)
  if (paired & isTRUE(welch.df))
    warning("Welch test does not apply to paired samples",
            call. = FALSE)
  if (!is.null(n2) & is.null(power))
    requested <- "power"
  if (is.null(n2) & !is.null(power))
    requested <- "n2"

  if (alternative == "not equal") {

    if (margin != 0)
      warning("`margin` argument is ignored")

    # sample size
    if (is.null(n2)) {
      HA_H0 <- mu1 - mu2
      beta <- 1 - power
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){
        n1.0 <- n2.0 * kappa

        if(paired) {
          df <- n2.0 - 1
        } else {
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }

        if(df<= 0 | is.infinite(df)){break}
        M <- qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta, df = df, ncp = 0, lower.tail = FALSE)

        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }


        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }

    # power
    if (is.null(power)) {

      HA_H0 <- mu1 - mu2
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        n1 <- n2 * kappa
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }

      power <- 1 - pt(qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda) +
        pt(-qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)
    }

  }
  else if (alternative == "greater" | alternative ==
           "less") {

    if (margin != 0)
      warning("`margin` argument is ignored", call. = FALSE)
    if (alternative == "greater" & (mu1 < mu2))
      stop("alternative = 'greater' but mu1 < mu2",
           call. = FALSE)
    if (alternative == "less" & (mu1 > mu2))
      stop("alternative = 'less' but mu1 > mu2",
           call. = FALSE)

    # sample size
    if (is.null(n2)) {
      HA_H0 <- mu1 - mu2
      beta <- 1 - power
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){

        if(paired) {
          df <- n2.0 - 1
        } else {
          n1.0 <- n2.0 * kappa
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }

        if(df <= 0 | is.infinite(df)){break}

        M <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta, df = df, ncp = 0, lower.tail = FALSE)

        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }

        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }

    # power
    if (is.null(power)) {
      HA_H0 <- mu1 - mu2
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        n1 <- n2 * kappa
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }


      if(alternative == "less") lambda <- abs(lambda)
      power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)
    }

  }
  else if (alternative == "non-inferior" | alternative ==
           "superior") {

    # sample size
    if (is.null(n2)) {
      beta <- 1 - power
      HA_H0 <- mu1 - mu2 - margin
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){
        n1.0 <- n2.0 * kappa

        if(paired) {
          df <- n2.0 - 1
        } else {
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }

        if(df <= 0 | is.infinite(df)){break}
        M <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta, df = df, ncp = 0, lower.tail = FALSE)

        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }

        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }

    # power
    if (is.null(power)) {
      HA_H0 <- mu1 - mu2 - margin
      n1 <- n2 * kappa
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }

      if(alternative == "non-inferior") lambda <- abs(lambda)
      power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)
    }

  }
  else if (alternative == "equivalent") {

    # sample size
    if (is.null(n2)) {
      beta <- 1 - power
      HA_H0 <- abs(mu1 - mu2) - margin
      i <- 0
      tol <- 0.01
      conv <- FALSE
      n2.0 <- 30
      while(i <= 100 & conv == FALSE){
        n1.0 <- n2.0 * kappa

        if(paired) {
          df <- n2.0 - 1
        } else {
          ifelse(welch.df,
                 df <- welch_df(sd1, sd2, n1.0, n2.0),
                 df <- n1.0 + n2.0 - 2)
        }

        if(df <= 0 | is.infinite(df)){break}
        M <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE) +
          qt(beta / 2, df = df, ncp = 0, lower.tail = FALSE)

        if(paired) {
          n2 = M^2 * (sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / (HA_H0)^2
        } else {
          n2 = M^2 * (sd1^2 / kappa + sd2^2) / (HA_H0)^2
        }

        if(abs(n2 - n2.0) < tol){conv <- TRUE}
        n2.0 <- (n2 + n2.0)/2
        i <- i + 1
      }
      ifelse(df > 0, n2 <- n2.0, n2 <- NA)
      n1 <- n2 * kappa
    }

    # power
    if (is.null(power)) {
      HA_H0 <- abs(mu1 - mu2) - margin
      if(paired) {
        df <- n2 - 1
        lambda <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
      } else {
        n1 <- n2 * kappa
        ifelse(welch.df,
               df <- welch_df(sd1, sd2, n1, n2),
               df <- n1 + n2 - 2)
        lambda = (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
      }

      power <- 2 * (1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }

  }
  else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }


  if(paired) {
    n2 <- ceiling(n2)
    ncp <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
    df <- n2 - 1
  } else {
    n1 <- ceiling(n1)
    n2 <- ceiling(n2)
    ncp <- (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
    ifelse(welch.df,
           df <- welch_df(sd1, sd2, n1, n2),
           df <- n1 + n2 - 2)
  }
  ifelse(paired, n <- n2, n <- c(n1 = n1, n2 = n2))

  hypothesis <- alternative

  if(verbose) {
    cat(ifelse(paired,
               " Difference between Two means \n (Paired Samples t Test) \n",
               " Difference between Two means \n (Independent Samples t Test) \n"),
        switch(hypothesis,
               `not equal` = "H0: mu1 = mu2 \n HA: mu1 != mu2 \n",
               `greater` = "H0: mu1 = mu2 \n HA: mu1 > mu2 \n",
               `less` = "H0: mu1 = mu2 \n HA: mu1 < mu2 \n",
               `non-inferior` = "H0: mu1 - mu2 <= margin \n HA: mu1 - mu2 > margin \n",
               `superior` = "H0: mu1 - mu2 <= margin \n HA: mu1 - mu2 > margin \n",
               `equivalent` = "H0: |mu1 - mu2| >= margin \n HA: |mu1 - mu2| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        if(paired) {
          c(" n =", ceiling(n))
        } else {
          c(" n1 =", ceiling(n1), "\n  n2 =", ceiling(n2))
        }, "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Degrees of freedom =", round(df, 2), "\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, kappa = kappa, welch.df = welch.df,
                                        paired = paired, paired.r = paired.r,
                                        alpha = alpha, margin = margin, alternative = alternative, verbose = verbose),
                           test = "t",
                           df = df,
                           ncp =  ncp,
                           power = power,
                           n = n),
                      class = c("pwrss", "t", "2means")))
}

##########################
# one correlation z test #
##########################

pwrss.z.cor <- pwrss.z.corr <- function (r = 0.50, r0 = 0, alpha = 0.05,
                          alternative = c("not equal", "greater", "less"),
                          n = NULL, power = NULL, verbose = TRUE)
{

  if (length(alternative) > 1)
    alternative <- alternative[1]

  if (!is.null(n) & is.null(power))
    requested <- "power"
  if (is.null(n) & !is.null(power))
    requested <- "n"

  if (r < 0 & alternative == "greater")
    stop("r < 0 but alternative  = 'greater' than 0",
         call. = FALSE)
  if (r > 0 & alternative == "less")
    stop("r > 0 but alternative  = 'less' than 0",
         call. = FALSE)
  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)

  if (alternative == "not equal") {

    if (is.null(n)) {
      beta <- 1 - power
      z <- 0.5 * log((1 + r) / (1 - r))
      z0 <- 0.5 * log((1 + r0) / (1 - r0))
      M <- qnorm(alpha / 2, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n <- M^2 / (z - z0)^2 + 3
    }
    if (is.null(power)) {
      z <- 0.5 * log((1 + r) / (1 - r))
      z0 <- 0.5 * log((1 + r0) / (1 - r0))
      lambda = (z - z0) / sqrt(1 / (n - 3))
      power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), lambda) +
        pnorm(-qnorm(alpha / 2, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "greater" | alternative ==
           "less") {

    if (alternative == "greater" & (r < r0))
      stop("alternative = 'greater' but r < r0",
           call. = FALSE)
    if (alternative == "less" & (r > r0))
      stop("alternative = 'less' but r > r0", call. = FALSE)

    if (is.null(n)) {
      beta <- 1 - power
      z <- 0.5 * log((1 + r) / (1 - r))
      z0 <- 0.5 * log((1 + r0) / (1 - r0))
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n <- M^2 / (z - z0)^2 + 3
    }
    if (is.null(power)) {
      z <- 0.5 * log((1 + r) / (1 - r))
      z0 <- 0.5 * log((1 + r0) / (1 - r0))
      lambda = (z - z0) / sqrt(1 / (n - 3))
      if(alternative == "less") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  }

  ncp <- (z - z0) / sqrt(1 / (n - 3))
  hypothesis <- alternative

  if(verbose) {
    cat(" A Correlation against a Constant (z Test) \n",
        switch(hypothesis,
               `not equal` = "H0: r = r0 \n HA: r != r0 \n",
               `greater` = "H0: r = r0 \n HA: r > r0 \n",
               `less` = "H0: r = r0 \n HA: r < r0 \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(r = r, r0 = r0, alpha = alpha, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = ncp,
                           power = power,
                           n = n),
                      class = c("pwrss", "z", "corr")))

}

###########################
# two correlations z test #
###########################

pwrss.z.2cors <- pwrss.z.2corrs <- function (r1 = 0.50, r2 = 0.30,
                            alpha = 0.05, kappa = 1,
                            alternative = c("not equal", "greater", "less"),
                            n2 = NULL, power = NULL, verbose = TRUE)
{

  if (length(alternative) > 1)
    alternative <- alternative[1]

  if (!is.null(n2) & is.null(power))
    requested <- "power"
  if (is.null(n2) & !is.null(power))
    requested <- "n2"

  if (r1 - r2 < 0 & alternative == "greater")
    stop("r1 - r2 < 0 but alternative  = 'greater' than 0",
         call. = FALSE)
  if (r1 - r2 > 0 & alternative == "less")
    stop("r1 - r2 > 0 but alternative  = 'less' than 0",
         call. = FALSE)
  if (is.null(n2) & is.null(power))
    stop("`n2` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n2) & !is.null(power))
    stop("one of the `n2` or `power` should be `NULL`",
         call. = FALSE)

  if (alternative == "not equal") {

    if (is.null(n2)) {
      beta <- 1 - power
      z1 <- 0.5 * log((1 + r1) / (1 - r1))
      z2 <- 0.5 * log((1 + r2) / (1 - r2))
      M <- qnorm(alpha / 2, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n2 <- uniroot(function(n2) M^2 - (z1 - z2)^2 / (1 / (kappa*n2 - 3) + 1 / (n2 - 3)), interval = c(-1e10,1e10))$root
      n1 <- kappa*n2
    }
    if (is.null(power)) {
      z1 <- 0.5 * log((1 + r1) / (1 - r1))
      z2 <- 0.5 * log((1 + r2) / (1 - r2))
      n1 <- kappa*n2
      lambda = (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
      power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), lambda) +
        pnorm(-qnorm(alpha / 2, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "greater" | alternative ==
           "less") {

    if (alternative == "greater" & (r1 < r2))
      stop("alternative = 'greater' but r1 < r2",
           call. = FALSE)
    if (alternative == "less" & (r1 > r2))
      stop("alternative = 'less' but r1 > r2", call. = FALSE)

    if (is.null(n2)) {
      beta <- 1 - power
      z1 <- 0.5 * log((1 + r1) / (1 - r1))
      z2 <- 0.5 * log((1 + r2) / (1 - r2))
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta, lower.tail = FALSE)
      n2 <- uniroot(function(n2) M^2 - (z1 - z2)^2 / (1 / (kappa*n2 - 3) + 1 / (n2 - 3)), interval = c(-1e10,1e10))$root
      n1 <- kappa*n2
    }
    if (is.null(power)) {
      z1 <- 0.5 * log((1 + r1) / (1 - r1))
      z2 <- 0.5 * log((1 + r2) / (1 - r2))
      n1 <- kappa*n2
      lambda = (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
      if(alternative == "less") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  }

  ncp <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  hypothesis <- alternative

  if(verbose) {
    cat(" Difference between Two Correlations \n (Independent Samples z Test) \n",
        switch(hypothesis,
               `not equal` = "H0: r1 = r2 \n HA: r1 != r2 \n",
               `greater` = "H0: r1 = r2 \n HA: r1 > r2 \n",
               `less` = "H0: r1 = r2 \n HA: r1 < r2 \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n1 =", ceiling(n1), "\n",
        " n2 =", ceiling(n2), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(r1 = r1, r2 = r2, kappa = kappa, alpha = alpha, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = ncp,
                           power = power,
                           n = c(n1 = n1, n2 = n2)),
                      class = c("pwrss", "z", "2corrs")))

}

############################
# linear regression z test #
############################

# when k = 1 and predictor is binary
# d <- 0.20
# r2 <- d^2 / (d^2 + 4)
# f2 <- r2 /(1 - r2)
# results will be same as pwrss.t.2means(mu1 = d,...)
# specify k = m (the default) to test r2 difference from zero
# specify k > m to test r2 change from zero
pwrss.f.regression <- pwrss.f.reg <- function (r2 = 0.10, f2 = r2 /(1 - r2),
                         k = 1, m = k, alpha = 0.05,
                         n = NULL, power = NULL, verbose = TRUE)
{

  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)
  if(m > k) stop("'m' cannot be greater than 'k'", call. = FALSE)

  if (is.null(n)) {

    n <- uniroot(function(n) {
      u <- m
      v <- n - k - 1
      lambda <- f2 * n
      power - pf(qf(alpha, df1 = u, df2 = v, lower.tail = FALSE),
                 df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)
    }, interval = c(k + 2, 1e10))$root
  }
  if (is.null(power)) {
    u <- m
    v <- n - k - 1
    lambda <- f2 * n
    power <- pf(qf(alpha, df1 = u, df2 = v, lower.tail = FALSE),
                df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)

  }

  ncp <- f2 * n
  df1 <- m
  df2 <- n - k - 1

  if(verbose) {
    cat(ifelse(m == k,
               " Linear Regression (F test) \n R-squared Deviation from 0 (zero) \n",
               " Hierarchical Linear Regression (F test) \n R-squared Change \n"),
        "H0: r2 = 0 \n HA: r2 > 0 \n",
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Numerator degrees of freedom =", round(df1, 3), "\n",
        "Denominator degrees of freedom =", round(df2, 3), "\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(r2 = r2, k = k, m = m, f2 = f2,
                                        n = ceiling(n), power = power, alpha = alpha, verbose = verbose),
                           test = "F",
                           df1 = df1,
                           df2 = df2,
                           ncp = ncp,
                           power = power,
                           n = ceiling(n)),
                      class = c("pwrss", "f", "reg")))


}

#######################
# ANOVA/ANCOVA F test #
#######################

pwrss.f.ancova <- function(eta2 = 0.01, f2 = eta2 / (1 - eta2),
                           n.way = length(n.levels),
                           n.levels = 2, n.covariates = 0, alpha = 0.05,
                           n = NULL, power = NULL, verbose = TRUE)
{

  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)

  if (!is.null(n) & is.null(power))
    requested <- "power"
  if (is.null(n) & !is.null(power))
    requested <- "n"

  ss <- function(df1, n.groups, n.covariates, f2, alpha, power) {
    n <- uniroot(function(n) {
      u <- df1
      v <- n - n.groups - n.covariates
      lambda <- f2 * n
      power - pf(qf(alpha, df1 = u, df2 = v, lower.tail = FALSE),
                 df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)
    }, interval = c(n.groups + n.covariates + 2, 1e10))$root
    n
  }

  pwr <- function(df1, n, n.groups, n.covariates, f2, alpha) {
    u <- df1
    v <- n - n.groups - n.covariates
    lambda <- f2 * n
    power <- pf(qf(alpha, df1 = u, df2 = v, lower.tail = FALSE),
                df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)
    power
  }

  if("f2" %in% names(as.list(match.call()))) eta2 <- f2 / (1 + f2)

  nway <- length(n.levels)
  n.groups <- prod(n.levels)

  if(nway != n.way) {
    warning("'n.way' does not match the length of 'n.levels' \n using 'n.way = length(n.levels)'", call. = FALSE)
    n.way <- nway
  }

  if(verbose) {
    cat(" ",
        switch(n.way,
               `1` = "One",
               `2` = "Two",
               `3` = "Three"),
        ifelse(n.covariates > 0,
               "-way Analysis of Covariance (ANCOVA) \n ",
               "-way Analysis of Variance (ANOVA) \n "),
        " H0: 'eta2' or 'f2' = 0 \n  HA: 'eta2' or 'f2' > 0 \n --------------------------------------\n",
        switch(n.way,
               `1` = c(" Factor A: ", n.levels, " levels \n"),
               `2` = c(" Factor A: ", n.levels[1], " levels \n", " Factor B: ", n.levels[2], " levels \n"),
               `3` = c(" Factor A: ", n.levels[1], " levels \n", " Factor B: ", n.levels[2], " levels \n", " Factor C: ", n.levels[3], " levels \n")),
        " --------------------------------------\n",
        sep = "")

  }

  if(n.way == 1) {

    df1 <- n.levels[1] - 1
    if(requested == "n") {
      n <- ss(df1 = df1, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      df2 <- n - n.groups - n.covariates
      ncp <- n*f2

      if(verbose) {
        effect <- "A"
        print(data.frame(effect = effect, power = round(power,3), n.total = ceiling(n),
                         ncp = round(ncp,3), df1 = df1, df2 = round(df2,3)),
              row.names = FALSE)
      }

    } else if (requested == "power") {
      power <- pwr(df1 = df1, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      df2 <- n - n.groups - n.covariates
      ncp <- n*f2

      if(verbose) {
        effect <- "A"
        print(data.frame(effect = effect, power = round(power,3), n.total = n,
                         ncp = round(ncp,3), df1 = df1, df2 = round(df2,3)),
              row.names = FALSE)
      }

    } else {
      stop("Invalid 'n' or 'power'", call. = FALSE)
    }

  } else if(n.way == 2) {

    df1.f1 <- n.levels[1] - 1
    df1.f2 <- n.levels[2] - 1
    df1.f1f2 <- prod(n.levels - 1)

    df1 <- c(A = df1.f1,
             B  = df1.f2 ,
             AxB = df1.f1f2)

    if(requested == "n") {
      ss.f1 <- ss(df1 = df1.f1, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f2 <- ss(df1 = df1.f2, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f1f2 <- ss(df1 = df1.f1f2, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      n <- c(A = ss.f1,
             B = ss.f2,
             AxB = ss.f1f2)
      df2 <- c(A = ss.f1 - n.groups - n.covariates,
               B  = ss.f2 - n.groups - n.covariates,
               AxB = ss.f1f2 - n.groups - n.covariates)
      ncp <- c(A = ss.f1*f2,
               B = ss.f2*f2,
               AxB = ss.f1f2*f2)
      power <- c(A = power,
                 B = power,
                 AxB = power)

      if(verbose) {
        effect <- c("A", "B", "A x B")
        n.tot <- c(ss.f1, ss.f2, ss.f1f2)
        print(data.frame(effect = effect, power = round(power,3), n.total = ceiling(n.tot),
                         ncp = round(ncp,3), df1 = df1, df2 = round(df2,3)),
              row.names = FALSE)
      }

    } else if (requested == "power") {
      pwr.f1 <- pwr(df1 = df1.f1, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f2 <- pwr(df1 = df1.f2, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f1f2 <- pwr(df1 = df1.f1f2, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      df2 <- c(A = n - n.groups - n.covariates,
               B  = n - n.groups - n.covariates,
               AxB = n - n.groups - n.covariates)
      ncp <- c(A = n*f2,
               B = n*f2,
               AxB = n*f2)
      n <- c(A  = n,
             B  = n,
             AxB = n)
      power <- c(A = pwr.f1,
                 B = pwr.f2,
                 AxB = pwr.f1f2)

      if(verbose) {
        effect <- c("A", "B", "A x B")
        power <- c(pwr.f1, pwr.f2, pwr.f1f2)
        print(data.frame(effect = effect, power = round(power,3), n.total = n,
                         ncp = round(ncp,3), df1 = df1, df2 = round(df2,3)),
              row.names = FALSE)
      }

    } else {

      stop("Invalid 'n' or 'power'", call. = FALSE)

    }

  } else if(n.way == 3) {

    df1.f1 <- n.levels[1] - 1
    df1.f2 <- n.levels[2] - 1
    df1.f3 <- n.levels[3] - 1
    df1.f1f2 <- prod(n.levels[c(1,2)] - 1)
    df1.f1f3 <- prod(n.levels[c(1,3)] - 1)
    df1.f2f3 <- prod(n.levels[c(2,3)] - 1)
    df1.f1f2f3 <- prod(n.levels - 1)

    df1 <- c(A = df1.f1,
             B = df1.f2,
             C = df1.f3,
             AxB = df1.f1f2,
             AxC = df1.f1f3,
             BxC = df1.f2f3,
             AxBxC = df1.f1f2f3)

    if(requested == "n") {
      ss.f1 <- ss(df1 = df1.f1, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f2 <- ss(df1 = df1.f2, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f3 <- ss(df1 = df1.f3, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f1f2 <- ss(df1 = df1.f1f2, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f1f3 <- ss(df1 = df1.f1f3, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f2f3 <- ss(df1 = df1.f2f3, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      ss.f1f2f3 <- ss(df1 = df1.f1f2f3, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      n <- c(A = ss.f1,
             B = ss.f2,
             C = ss.f3,
             AxB = ss.f1f2,
             AxC = ss.f1f3,
             BxC = ss.f2f3,
             AxBxC = ss.f1f2f3)
      df2 <- c(A = ss.f1 - n.groups - n.covariates,
               B = ss.f2 - n.groups - n.covariates,
               C = ss.f3 - n.groups - n.covariates,
               AxB = ss.f1f2 - n.groups - n.covariates,
               AxC = ss.f1f3 - n.groups - n.covariates,
               BxC = ss.f2f3 - n.groups - n.covariates,
               AxBxC = ss.f1f2f3 - n.groups - n.covariates)
      ncp <- c(A = ss.f1*f2,
               B = ss.f2*f2,
               C = ss.f3*f2,
               AxB = ss.f1f2*f2,
               AxC = ss.f1f3*f2,
               BxC = ss.f2f3*f2,
               AxBxC = ss.f1f2f3*f2)
      power <- c(A = power,
                 B = power,
                 C = power,
                 AxB = power,
                 AxC = power,
                 BxC = power,
                 AxBxC = power)

      if(verbose) {
        effect <- c("A", "B", "C", "A x B", "A x C", "B x C", "A x B x C")
        n.tot <- c(ss.f1, ss.f2, ss.f3, ss.f1f2, ss.f1f3, ss.f2f3, ss.f1f2f3)
        print(data.frame(effect = effect, power = round(power,3), n.total = ceiling(n.tot),
                         ncp = round(ncp,3), df1 = df1, df2 = round(df2,3)),
              row.names = FALSE)
      }

    } else if (requested == "power") {
      pwr.f1 <- pwr(df1 = df1.f1, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f2 <- pwr(df1 = df1.f2, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f3 <- pwr(df1 = df1.f3, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f1f2 <- pwr(df1 = df1.f1f2, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f1f3 <- pwr(df1 = df1.f1f3, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f2f3 <- pwr(df1 = df1.f2f3, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      pwr.f1f2f3 <- pwr(df1 = df1.f1f2f3, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      power <-c(A = pwr.f1,
                B = pwr.f2,
                C = pwr.f3,
                AxB = pwr.f1f2,
                AxC = pwr.f1f3,
                BxC = pwr.f2f3,
                AxBxC = pwr.f1f2f3)
      ncp <-c(A = n*f2,
              B = n*f2,
              C = n*f2,
              AxB = n*f2,
              AxC = n*f2,
              BxC = n*f2,
              AxBxC = n*f2)
      df2 <- c(A = n - n.groups - n.covariates,
               B = n - n.groups - n.covariates,
               C = n - n.groups - n.covariates,
               AxB = n - n.groups - n.covariates,
               AxC = n - n.groups - n.covariates,
               BxC = n - n.groups - n.covariates,
               AxBxC = n - n.groups - n.covariates)
      n <- c( A = n,
              B = n,
              C = n,
              AxB = n,
              AxC = n,
              BxC = n,
              AxBxC = n)

      if(verbose) {
        effect <- c("A", "B", "C", "A x B", "A x C", "B x C", "A x B x C")
        power <- c(pwr.f1, pwr.f2, pwr.f3, pwr.f1f2, pwr.f1f3, pwr.f2f3, pwr.f1f2f3)
        print(data.frame(effect = effect, power = round(power,3), n.total = n,
                         ncp = round(ncp,3), df1 = df1, df2 = round(df2,3)),
              row.names = FALSE)
      }

    } else {

      stop("Invalid 'n' or 'power'", call. = FALSE)

    }

  } else {

    stop("More than three-way ANOVA or ANCOVA is not allowed")

  }

  if(verbose) {
    cat(" --------------------------------------\n",
        "Type I error rate:", round(alpha, 3))
  }

  invisible(structure(list(parms = list(eta2 = eta2, f2 = f2, n.way = n.way, n.levels = n.levels,
                                        n.covariates = n.covariates, alpha = alpha, verbose = verbose),
                           test = "F",
                           df1 = df1,
                           df2 = df2,
                           ncp = ncp,
                           power = power,
                           n = n),
                      class = c("pwrss", "f", "ancova")))
} # end of pwrss.f.ancova()

##################################
# Repeated measures ANOVA F test #
##################################

# 1 / (n.rm - 1) < epsilon < 1 when spechiricty assumptions does not hold (or when compound symmetry is violated)
# Greenhouse and Geisser (1959), the other is by Huynh and Feldt (1976).
pwrss.f.rmanova <- function (eta2 = 0.10, f2 = eta2/(1 - eta2),
                             corr.rm = 0.50, n.levels = 2, n.rm = 2,
                             epsilon = 1, alpha = 0.05,
                             type = c("between", "within", "interaction"),
                             n = NULL, power = NULL, verbose = TRUE)
{

  user.parms.names <- names(as.list(match.call()))
  if("repmeasures.r" %in% user.parms.names) stop("'repmeasures.r' argument is obsolete, use 'corr.rm' instead", call. = FALSE)
  if("n.measurements" %in% user.parms.names) stop("'n.measurements' argument is obsolete, use 'n.rm' instead", call. = FALSE)

  if (length(type > 1)) type <- type[1]
  effect <- type
  type <- switch(type, between = 0, within = 1, interaction = 2)

  if(type == 0 & "epsilon" %in% names(as.list(match.call())))
    warning("non-spehericity correction does not apply to 'between' effect", call. = FALSE)

  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)

  ss <- function(f2, n.levels, n.rm, epsilon, alpha, power, type) {
    n.min <- n.levels + 1
    n <- try(silent = TRUE,
             uniroot(function(n) {
               if(type == 0) {
                 df1 <- n.levels - 1
                 df2 <- n - n.levels
               } else if(type == 1) {
                 df1 <- (n.rm - 1) * epsilon
                 df2 <- (n - n.levels) * (n.rm - 1) * epsilon
               } else if(type == 2) {
                 df1 <- (n.levels - 1) * (n.rm - 1) * epsilon
                 df2 <- (n - n.levels) * (n.rm - 1) * epsilon
               } else {
                 stop("Unknown type of effect", call. = FALSE)
               }
               u <- df1
               v <- df2
               lambda <- f2 * n * epsilon
               power - pf(qf(alpha, df1 = u, df2 = v, lower.tail = FALSE),
                          df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)
             }, interval = c(n.min, 1e10))$root
    ) # try
    if(inherits(n, "try-error") | n == 1e+10) stop("design is not feasible", call. = FALSE)
    n
  }

  pwr <- function(f2, n, n.levels, n.rm, epsilon, alpha, type) {
    if(type == 0) {
      df1 <- n.levels - 1
      df2 <- n - n.levels
    } else if(type == 1) {
      df1 <- (n.rm - 1) * epsilon
      df2 <- (n - n.levels) * df1
    } else if(type == 2) {
      df1 <- (n.levels - 1) * (n.rm - 1) * epsilon
      df2 <- (n - n.levels) * (n.rm - 1) * epsilon
    } else {
      stop("Unknown type of effect", call. = FALSE)
    }
    u <- df1
    v <- df2
    if(u < 1 | v < 1) stop("design is not feasible", call. = FALSE)
    lambda <- f2 * n * epsilon
    power <- pf(qf(alpha, df1 = u, df2 = v, lower.tail = FALSE),
                df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)
    power
  }

  if (type == 0) {
    f2 <- f2  * (n.rm / (1 + (n.rm - 1) * corr.rm))
  } else {
    f2 <- f2  * (n.rm / (1 - corr.rm))
  }

  if (is.null(n)) {
    n <- ss(f2 = f2, n.levels = n.levels, n.rm = n.rm,
            epsilon = epsilon, alpha = alpha, power = power, type = type)
  }

  if (is.null(power)) {
    power <- pwr(f2 = f2, n = n, n.levels = n.levels, n.rm = n.rm,
                 epsilon = epsilon, alpha = alpha, type = type)
  }

  if(type == 0) {
    df1 <- n.levels - 1
    df2 <- n - n.levels
  } else if(type == 1) {
    df1 <- (n.rm - 1) * epsilon
    df2 <- (n - n.levels) * df1
  } else if(type == 2) {
    df1 <- (n.levels - 1) * (n.rm - 1) * epsilon
    df2 <- (n - n.levels) * (n.rm - 1) * epsilon
  }

  ncp <- f2 * n * epsilon

  if(verbose) {
    cat(" One-way Repeated Measures \n Analysis of Variance (F test) \n",
        "H0: eta2 = 0 (or f2 = 0) \n HA: eta2 > 0 (or f2 > 0) \n",
        "------------------------------ \n",
        "Number of levels (groups) =", n.levels, "\n",
        "Number of repeated measurements =", n.rm, "\n",
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " Total n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Type of the effect =", dQuote(effect), "\n",
        "Numerator degrees of freedom =", round(df1, 3), "\n",
        "Denominator degrees of freedom =", round(df2, 3), "\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(f2 = f2, n.levels = n.levels,
                                        n.rm = n.rm, corr.rm = corr.rm,
                                        epsilon = epsilon, alpha = alpha, verbose = verbose),
                           test = "F",
                           df1 = df1,
                           df2 = df2,
                           ncp = f2 * n * epsilon,
                           power = power,
                           n = n),
                      class = c("pwrss", "f", "rmanova")))

}

#############################
# linear regresstion t test #
#############################

# if the predictor is binary
# provide sdx = sqrt(p * (1 - p))
# p = proportion of subjects in treatment group
# use defaults if beta1 is standardized
pwrss.t.regression <- pwrss.t.reg <- function (beta1 = 0.25, beta0 = 0, margin = 0,
                                               sdx = 1, sdy = 1,
                                               k = 1, r2 = (beta1 * sdx / sdy)^2,
                                               alpha = 0.05, n = NULL, power = NULL,
                                               alternative = c("not equal", "less", "greater",
                                                               "non-inferior", "superior", "equivalent"),
                                               verbose = TRUE) {


  alternative <- tolower(match.arg(alternative))
  user.parms.names <- names(as.list(match.call()))

  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)
  if(r2 > 0.99) stop("unreasonable value for r2, specify r2 explicitly or modify 'beta', 'sdx', 'sdy'", call. = FALSE)
  if(r2 < (beta1 * sdx / sdy)^2)
    warning("possibly incongruent arguments, need a larger 'r2'", call. = FALSE)


  if(alternative == "equivalent") {

    if(margin < 0) stop("'margin' should be positive in equivalence design", call. = FALSE)

  } else if(alternative == "non-inferior") {

    if(margin > 0 & (beta1 > beta0)) stop("expecting 'beta1 > beta0' and negative 'margin' if higher values of the outcome are better", call. = FALSE)
    if(margin < 0 & (beta1 < beta0)) stop("expecting 'beta1 < beta0' and positive 'margin' if lower values of the outcome are better", call. = FALSE)

  } else if(alternative == "superior") {

    if(margin < 0 & (beta1 > beta0)) stop("expecting 'beta1 < beta0' and positive 'margin' if lower values of the outcome are better", call. = FALSE)
    if(margin > 0 & (beta1 < beta0)) stop("expecting 'beta1 > beta0' and negative 'margin' if higher values of the outcome are better", call. = FALSE)

  } else if(alternative == "greater") {

    if("margin" %in% user.parms.names) message("ignoring any specifications to 'margin'")
    if(beta1 > beta0) stop("alternative = 'less' but beta1 > beta0", call. = FALSE)

  } else if(alternative == "less") {

    if("margin" %in% user.parms.names) message("ignoring any specifications to 'margin'")
    if(beta1 > beta0) stop("alternative = 'less' but beta1 > beta0", call. = FALSE)

  }

  HA_H0 <- switch(alternative,
                  `greater` = beta1 - beta0,
                  `less` = beta1 - beta0,
                  `not equal` = beta1 - beta0,
                  `non-inferior` = beta1 - beta0 - margin,
                  `superior` = beta1 - beta0 - margin,
                  `equivalent` = c(abs(beta1 - beta0) - margin, abs(beta1 - beta0) + margin))

  pwr.fun.body <- quote({
    if(alternative == "not equal") {
      power <- 1 - pt(qt(alpha / 2, df = v, ncp = 0, lower.tail = FALSE), df = v, ncp = abs(lambda)) +
        pt(-qt(alpha / 2, df = v, ncp = 0, lower.tail = FALSE), df = v, ncp = abs(lambda))
    } else if(alternative %in% c("greater", "less", "superior", "non-inferior")) {
      power <- 1 - pt(qt(alpha, df = v, ncp = 0, lower.tail = FALSE), df = v, ncp = abs(lambda))
    } else if(alternative == "equivalent") {
      power <- 1 - pt(qt(alpha, df = v, ncp = 0, lower.tail = FALSE), df = v, ncp = abs(lambda[1])) +
        1 - pt(qt(alpha, df = v, ncp = 0, lower.tail = FALSE), df = v, ncp = abs(lambda[2])) - 1
    } else {
      stop("not a valid test type", call. = FALSE)
    }
    power
  })


  if(is.null(power)) {
    ifelse(alternative == "equivalent",
           lambda <- cbind(left = HA_H0[1] / ((sdy / sdx) * sqrt((1 - r2) / n)),
                           right = HA_H0[2] / ((sdy / sdx) * sqrt((1 - r2) / n))),
           lambda <- HA_H0 / ((sdy / sdx) * sqrt((1 - r2) / n)))
    v <- n - k - 1
    power <-  eval(pwr.fun.body)
  }

  if(is.null(n)) {
    n <- uniroot(function(n) {
      ifelse(alternative == "equivalent",
             lambda <- cbind(left = HA_H0[1] / ((sdy / sdx) * sqrt((1 - r2) / n)),
                             right = HA_H0[2] / ((sdy / sdx) * sqrt((1 - r2) / n))),
             lambda <- HA_H0 / ((sdy / sdx) * sqrt((1 - r2) / n)))
      v <- n - k - 1
      power - eval(pwr.fun.body)
    }, interval = c(k + 2, 1e10))$root

    ifelse(alternative == "equivalent",
           lambda <- cbind(left = HA_H0[1] / ((sdy / sdx) * sqrt((1 - r2) / n)),
                           right = HA_H0[2] / ((sdy / sdx) * sqrt((1 - r2) / n))),
           lambda <- HA_H0 / ((sdy / sdx) * sqrt((1 - r2) / n)))
    v <- n - k - 1
  }

  if(verbose) {
    cat(" Linear Regression Coefficient (t Test) \n",
        switch(alternative,
               `not equal` = "H0: beta1 = beta0 \n HA: beta1 != beta0 \n",
               `greater` = "H0: beta1 = beta0 \n HA: beta1 > beta0 \n",
               `less` = "H0: beta1 = beta0 \n HA: beta1 < beta0 \n",
               `non-inferior` = ifelse(beta1 > beta0,
                                       "H0: beta1 - beta0 <= margin \n HA: beta1 - beta0 > margin \n",
                                       "H0: beta1 - beta0 >= margin \n HA: beta1 - beta0 < margin \n"),
               `superior` = ifelse(beta1 > beta0,
                                   "H0: beta1 - beta0 <= margin \n HA: beta1 - beta0 > margin \n",
                                   "H0: beta1 - beta0 >= margin \n HA: beta1 - beta0 < margin \n"),
               `equivalent` = "H0: |beta1 - beta0| >= margin \n HA: |beta1 - beta0| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Degrees of freedom =", round(v, 3), "\n",
        "Non-centrality parameter =", round(lambda, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(beta1 = beta1, beta0, beta0, margin = margin,
                                        sdx = sdx, sdy = sdy, k = k, r2 = r2,
                                        alpha = alpha, alternative = alternative, verbose = verbose),
                           test = "t",
                           df = v,
                           ncp = lambda,
                           power = power,
                           n = n),
                      class = c("pwrss", "t", "reg")))
}

############################
# linear regression z test #
############################

# if the predictor is binary
# provide sdx = sqrt(p * (1 - p))
# p = proportion of subjects in treatment group
# use defaults if beta1 is standardized
pwrss.z.regression <- pwrss.z.reg <- function (beta1 = 0.25, beta0 = 0, margin = 0,
                                               sdx = 1, sdy = 1,
                                               k = 1, r2 = (beta1 * sdx / sdy)^2,
                                               alpha = 0.05, n = NULL, power = NULL,
                                               alternative = c("not equal", "less", "greater",
                                                               "non-inferior", "superior", "equivalent"),
                                               verbose = TRUE) {


  alternative <- tolower(match.arg(alternative))
  user.parms.names <- names(as.list(match.call()))

  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)
  if(r2 > 0.99) stop("unreasonable value for r2, specify r2 explicitly or modify 'beta', 'sdx', 'sdy'", call. = FALSE)
  if(r2 < (beta1 * sdx / sdy)^2)
    warning("possibly incongruent arguments, need a larger 'r2'", call. = FALSE)


  if(alternative == "equivalent") {

    if(margin < 0) stop("'margin' should be positive in equivalence design", call. = FALSE)

  } else if(alternative == "non-inferior") {

    if(margin > 0 & (beta1 > beta0)) stop("expecting 'beta1 > beta0' and negative 'margin' if higher values of the outcome are better", call. = FALSE)
    if(margin < 0 & (beta1 < beta0)) stop("expecting 'beta1 < beta0' and positive 'margin' if lower values of the outcome are better", call. = FALSE)

  } else if(alternative == "superior") {

    if(margin < 0 & (beta1 > beta0)) stop("expecting 'beta1 < beta0' and positive 'margin' if lower values of the outcome are better", call. = FALSE)
    if(margin > 0 & (beta1 < beta0)) stop("expecting 'beta1 > beta0' and negative 'margin' if higher values of the outcome are better", call. = FALSE)

  } else if(alternative == "greater") {

    if("margin" %in% user.parms.names) message("ignoring any specifications to 'margin'")
    if(beta1 > beta0) stop("alternative = 'less' but beta1 > beta0", call. = FALSE)

  } else if(alternative == "less") {

    if("margin" %in% user.parms.names) message("ignoring any specifications to 'margin'")
    if(beta1 > beta0) stop("alternative = 'less' but beta1 > beta0", call. = FALSE)

  }

  HA_H0 <- switch(alternative,
                  `greater` = beta1 - beta0,
                  `less` = beta1 - beta0,
                  `not equal` = beta1 - beta0,
                  `non-inferior` = beta1 - beta0 - margin,
                  `superior` = beta1 - beta0 - margin,
                  `equivalent` = c(abs(beta1 - beta0) - margin, abs(beta1 - beta0) + margin))

  pwr.fun.body <- quote({
    if(alternative == "not equal") {
      power <- 1 - pnorm(qnorm(alpha / 2, mean = 0, lower.tail = FALSE), mean = abs(lambda)) +
        pnorm(-qnorm(alpha / 2, mean = 0, lower.tail = FALSE), mean = abs(lambda))
    } else if(alternative %in% c("greater", "less", "superior", "non-inferior")) {
      power <- 1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = abs(lambda))
    } else if(alternative == "equivalent") {
      power <- 1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = abs(lambda[1])) +
        1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = abs(lambda[2])) - 1
    } else {
      stop("not a valid test type", call. = FALSE)
    }
    power
  })


  if(is.null(power)) {
    ifelse(alternative == "equivalent",
           lambda <- cbind(left = HA_H0[1] / ((sdy / sdx) * sqrt((1 - r2) / n)),
                           right = HA_H0[2] / ((sdy / sdx) * sqrt((1 - r2) / n))),
           lambda <- HA_H0 / ((sdy / sdx) * sqrt((1 - r2) / n)))
    v <- n - k - 1
    power <-  eval(pwr.fun.body)
  }

  if(is.null(n)) {
    n <- uniroot(function(n) {
      ifelse(alternative == "equivalent",
             lambda <- cbind(left = HA_H0[1] / ((sdy / sdx) * sqrt((1 - r2) / n)),
                             right = HA_H0[2] / ((sdy / sdx) * sqrt((1 - r2) / n))),
             lambda <- HA_H0 / ((sdy / sdx) * sqrt((1 - r2) / n)))
      v <- n - k - 1
      power - eval(pwr.fun.body)
    }, interval = c(k + 2, 1e10))$root

    ifelse(alternative == "equivalent",
           lambda <- cbind(left = HA_H0[1] / ((sdy / sdx) * sqrt((1 - r2) / n)),
                           right = HA_H0[2] / ((sdy / sdx) * sqrt((1 - r2) / n))),
           lambda <- HA_H0 / ((sdy / sdx) * sqrt((1 - r2) / n)))
    v <- n - k - 1
  }

  if(verbose) {
    cat(" Linear Regression Coefficient (z Test) \n",
        switch(alternative,
               `not equal` = "H0: beta1 = beta0 \n HA: beta1 != beta0 \n",
               `greater` = "H0: beta1 = beta0 \n HA: beta1 > beta0 \n",
               `less` = "H0: beta1 = beta0 \n HA: beta1 < beta0 \n",
               `non-inferior` = ifelse(beta1 > beta0,
                                       "H0: beta1 - beta0 <= margin \n HA: beta1 - beta0 > margin \n",
                                       "H0: beta1 - beta0 >= margin \n HA: beta1 - beta0 < margin \n"),
               `superior` = ifelse(beta1 > beta0,
                                   "H0: beta1 - beta0 <= margin \n HA: beta1 - beta0 > margin \n",
                                   "H0: beta1 - beta0 >= margin \n HA: beta1 - beta0 < margin \n"),
               `equivalent` = "H0: |beta1 - beta0| >= margin \n HA: |beta1 - beta0| < margin \n"),
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Degrees of freedom =", round(v, 3), "\n",
        "Non-centrality parameter =", round(lambda, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(beta1 = beta1, beta0, beta0, margin = margin,
                                        sdx = sdx, sdy = sdy, k = k, r2 = r2,
                                        alpha = alpha, alternative = alternative, verbose = verbose),
                           test = "z",
                           df = v,
                           ncp = lambda,
                           power = power,
                           n = n),
                      class = c("pwrss", "z", "reg")))
}

####################
# mediation z test #
####################

## 'cp = 0' by default, implying complete mediation (it increases explanatory power of the covariate
# use 'r2m.x' and 'r2y.mx' to adjust standard error for other predictors in mediation and outcome model
pwrss.z.mediation  <- pwrss.z.med  <- function(a, b, cp = 0,
                                               sdx = 1, sdm = 1, sdy = 1,
                                               r2m.x = a^2 * sdx^2 / sdm^2,
                                               r2y.mx = (b^2 * sdm^2 + cp^2 * sdx^2) / sdy^2,
                                               n = NULL, power = NULL,
                                               alpha = 0.05, alternative = c("not equal", "less", "greater"),
                                               mc = TRUE, nsims = 1000, ndraws = 1000,
                                               verbose = TRUE) {


  user.parms.names <- names(as.list(match.call()))

  if (is.null(n) & is.null(power))
    stop("`n` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n) & !is.null(power))
    stop("one of the `n` or `power` should be `NULL`",
         call. = FALSE)

  alternative <- tolower(match.arg(alternative))
  if (alternative == "greater" & (a * b < 0)) stop("alternative = 'greater' but a * b < 0", call. = FALSE)
  if (alternative == "less" & (a * b > 0)) stop("alternative = 'less' but a * b > 0", call. = FALSE)

  if(r2y.mx == 0 & "cp" %in% user.parms.names)
    warning("ignoring any specification to 'cp'", call. = FALSE)

  if(r2m.x < a^2 * sdx^2 / sdm^2)
    warning("specified 'r2m.x' is smaller than the base 'r2m.x'", call. = FALSE)

  if(r2y.mx < (b^2 * sdm^2 + cp^2 * sdx^2) / sdy^2)
    warning("specified 'r2y.mx' is smaller than the base 'r2y.mx'", call. = FALSE)

  .se.a <- function(sdm, sdx, r2m.x, n) {
    var.a <- (1 / n) * (sdm^2) * (1 - r2m.x) / (sdx^2)
    sqrt(var.a)
  }

  .se.b <- function(sdy, sdm, r2y.mx, r2m.x, n) {
    var.b <- (1 / n) * (sdy^2) * (1 - r2y.mx) / ((sdm^2) * (1 - r2m.x))
    sqrt(var.b)
  }

  if(is.null(power)) {

    se.a <- .se.a(sdm, sdx, r2m.x, n)
    se.b <- .se.b(sdy, sdm, r2y.mx, r2m.x, n)

    sobel.se <- sqrt(a^2 * se.b^2  + b^2 * se.a^2)
    aroian.se <- sqrt(a^2 * se.b^2  + b^2 * se.a^2 + se.a^2 * se.b^2)
    goodman.se <- sqrt(a^2 * se.b^2  + b^2 * se.a^2 - se.a^2 * se.b^2)
    lambda.sobel <- (a * b) / sobel.se
    lambda.aroian <- (a * b) / aroian.se
    lambda.goodman <- (a * b) / goodman.se
    power.sobel <- power.z.test(ncp = lambda.sobel, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)
    power.aroian <- power.z.test(ncp = lambda.aroian, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)
    power.goodman <- power.z.test(ncp = lambda.goodman, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)

    # joint test
    power.a <- power.z.test(ncp = a / se.a, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)
    power.b <- power.z.test(ncp = b / se.b, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)
    power.joint <- power.a * power.b

    # monte carlo interval test
    if(mc) {
      a <- abs(a)
      b <- abs(b)
      rejmc <- NULL
      for (i in 1:nsims){
        a.star <- rnorm(1, a, se.a)
        b.star <- rnorm(1, b, se.b)
        rejmc <- c(rejmc, quantile(rnorm(ndraws, a.star, se.a) * rnorm(ndraws, b.star, se.b),
                                   probs = ifelse(alternative == "not equal", alpha / 2, alpha), na.rm = TRUE) > 0)
      }
      power.mc <- mean(rejmc)
    } else {
      power.mc <- NA
    }

    n.sobel <- n.aroian <- n.goodman <- n.joint <- n.mc <- n

  }


  if(is.null(n)) {

    n.sobel <- uniroot(function(n) {
      se.a <- .se.a(sdm, sdx, r2m.x, n)
      se.b <- .se.b(sdy, sdm, r2y.mx, r2m.x, n)
      sobel.se <- sqrt(a^2 * se.b^2  + b^2 * se.a^2)
      lambda.sobel <- (a * b) / sobel.se
      power - power.z.test(ncp = lambda.sobel, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)
    }, interval = c(10, 1e10))$root

    se.a <- .se.a(sdm, sdx, r2m.x, n.sobel)
    se.b <- .se.b(sdy, sdm, r2y.mx, r2m.x, n.sobel)
    sobel.se <- sqrt(a^2 * se.b^2  + b^2 * se.a^2)
    lambda.sobel <- (a * b) / sobel.se

    n.aroian <- uniroot(function(n) {
      se.a <- .se.a(sdm, sdx, r2m.x, n)
      se.b <- .se.b(sdy, sdm, r2y.mx, r2m.x, n)
      aroian.se <- sqrt(a^2 * se.b^2  + b^2 * se.a^2 + se.a^2 * se.b^2)
      lambda.aroian <- (a * b) / aroian.se
      power - power.z.test(ncp = lambda.aroian, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)

    }, interval = c(10, 1e10))$root

    n.goodman <- uniroot(function(n) {
      se.a <- .se.a(sdm, sdx, r2m.x, n)
      se.b <- .se.b(sdy, sdm, r2y.mx, r2m.x, n)
      if(a^2 * se.b^2  + b^2 * se.a^2 < se.a^2 * se.b^2) {
        goodman.se <- 1
      } else {
        goodman.se <- sqrt(a^2 * se.b^2  + b^2 * se.a^2 - se.a^2 * se.b^2)
      }
      lambda.goodman <- (a * b) / goodman.se
      power - power.z.test(ncp = lambda.goodman, alpha = alpha, alternative = alternative, plot = FALSE, verbose = FALSE)

    }, interval = c(10, 1e10))$root

    n.joint <- n.mc <- NA
    lambda.goodman <- lambda.aroian <- lambda.sobel
    power.sobel <- power.aroian <- power.goodman <- power
    power.joint <- power.mc <- NA
  }

  if(verbose) {

    if(alternative == "not equal") H0HA.test <- "H0: a * b = 0 \n HA: a * b != 0 \n"
    if(alternative == "greater") H0HA.test <- "H0: a * b <= 0 \n HA: a * b > 0 \n"
    if(alternative == "less") H0HA.test <- "H0: a * b >= 0 \n HA: a * b < 0 \n"
    cat(" Indirect Effect in Mediation Model\n", H0HA.test, "--------------------------- \n")

    test <- c("Sobel", "Aroian", "Goodman", "Joint", "Monte Carlo")
    power <- c(round(power.sobel, 3),
               round(power.aroian, 3),
               round(power.goodman, 3),
               round(power.joint, 3),
               round(power.mc, 3))

    n <- c(ceiling(n.sobel),
           ceiling(n.aroian),
           ceiling(n.goodman),
           ceiling(n.joint),
           ceiling(n.mc))

    ncp <- c(round(lambda.sobel, 3),
             round(lambda.aroian, 3),
             round(lambda.goodman, 3),
             NA,
             NA)

    print(data.frame(test = test, power = power, n = n,
                     ncp = ncp),
          row.names = FALSE)

    cat("---------------------------- \n", "Type I error rate =", round(alpha, 3), "\n")

  }

  invisible(structure(list(parms = list(a = a, b = b, cp = cp,
                                        sdx = sdx, sdm = sdm, sdy = sdy,
                                        r2m.x = r2m.x, r2y.mx = r2y.mx,
                                        alpha = alpha, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = c(sobel = lambda.sobel, aroian = lambda.aroian, goodman = lambda.goodman),
                           power = c(sobel = power.sobel, aroian = power.aroian,
                                     goodman = power.goodman, joint = power.joint, mc = power.mc),
                           n = c(sobel = n.sobel, aroian = n.aroian, goodman = n.goodman)),
                      class = c("pwrss", "z", "med")))

} # end of pwrss.z.med()

#################
# plot function #
#################

plot.pwrss <- function(x, ...) {

   if(all(c("pwrss","t") %in% class(x))) {
    power.t.test(ncp = abs(max(x$ncp)),
                 df = x$df,
                 alpha = x$parms$alpha,
                 alternative = x$parms$alternative,
                 verbose = FALSE)
  } else if(all(c("pwrss","z") %in% class(x))) {

    if("med" %in% class(x)) {
      layout(matrix(c(1,3,
                      2,3), nrow = 2, ncol = 2))
      default.pars <- par(no.readonly = TRUE)
      on.exit(par(default.pars))

      power.z.test(ncp = abs(x$ncp["sobel"]), plot.main = "Sobel",
                   alpha = x$parms$alpha,
                   alternative = x$parms$alternative,
                   verbose = FALSE)
      power.z.test(ncp = abs(x$ncp["aroian"]), plot.main = "Aroian",
                   alpha = x$parms$alpha,
                   alternative = x$parms$alternative,
                   verbose = FALSE)
      power.z.test(ncp = abs(x$ncp["goodman"]), plot.main = "Goodman",
                   alpha = x$parms$alpha,
                   alternative = x$parms$alternative,
                   verbose = FALSE)

      par(mfrow = c(1,1))
    } else {
      power.z.test(ncp = abs(max(x$ncp)),
                   alpha = x$parms$alpha,
                   alternative = x$parms$alternative,
                   verbose = FALSE)
    }

  } else if(all(c("pwrss","f") %in% class(x))) {

    if("reg" %in% class(x)) {
      power.f.test(ncp = abs(x$ncp),
                   df1 = x$df1,
                   df2 = x$df2,
                   alpha = x$parms$alpha,
                   verbose = FALSE)
    } else if("ancova" %in% class(x)) {
      if(x$parms$n.way == 1) {
        power.f.test(ncp = abs(x$ncp),
                     df1 = x$df1,
                     df2 = x$df2,
                     alpha = x$parms$alpha,
                     verbose = FALSE)
      } else if(x$parms$n.way == 2) {

        layout(matrix(c(1,3,
                        2,3), nrow = 2, ncol = 2))
        default.pars <- par(no.readonly = TRUE)
        on.exit(par(default.pars))

        power.f.test(ncp = abs(x$ncp["A"]), plot.main = "A",
                     df1 = x$df1["A"],
                     df2 = x$df2["A"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["B"]), plot.main = "B",
                     df1 = x$df1["B"],
                     df2 = x$df2["B"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["AxB"]), plot.main = "A x B",
                     df1 = x$df1["AxB"],
                     df2 = x$df2["AxB"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)

        par(mfrow = c(1,1))

      } else if(x$parms$n.way == 3) {

        layout(matrix(c(1,4,7,
                        2,5,7,
                        3,6,7), nrow = 3, ncol = 3))
        default.pars <- par(no.readonly = TRUE)
        on.exit(par(default.pars))

        power.f.test(ncp = abs(x$ncp["A"]), plot.main = "A",
                     df1 = x$df1["A"],
                     df2 = x$df2["A"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["B"]), plot.main = "B",
                     df1 = x$df1["B"],
                     df2 = x$df2["B"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["C"]), plot.main = "C",
                     df1 = x$df1["C"],
                     df2 = x$df2["C"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["AxB"]), plot.main = "A x B",
                     df1 = x$df1["AxB"],
                     df2 = x$df2["AxB"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["AxC"]), plot.main = "A x C",
                     df1 = x$df1["AxC"],
                     df2 = x$df2["AxC"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["BxC"]), plot.main = "B x C",
                     df1 = x$df1["BxC"],
                     df2 = x$df2["BxC"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)
        power.f.test(ncp = abs(x$ncp["AxBxC"]), plot.main = "A x B x C",
                     df1 = x$df1["AxBxC"],
                     df2 = x$df2["AxBxC"],
                     alpha = x$parms$alpha,
                     verbose = FALSE)

        par(mfrow = c(1,1))

      }

    } else if("rmanova" %in% class(x)) {
      power.f.test(ncp = abs(x$ncp),
                   df1 = x$df1,
                   df2 = x$df2,
                   alpha = x$parms$alpha,
                   verbose = FALSE)
    }

  } else  if(all(c("pwrss","chisq") %in% class(x))) {

    power.chisq.test(ncp = abs(x$ncp),
                     df = x$df,
                     alpha = x$parms$alpha,
                     verbose = FALSE)

  } else {

    stop("not an object of the type 'pwrss'", call. = FALSE)

  }

}

#################################
# two means non-parametric test #
#################################

# margin should be on the same scale as mu1 and mu2
# just specify mu1 for cohen's d
# just specify sd1 for pooled standard deviation
pwrss.np.2groups <- pwrss.np.2means <- function(mu1 = 0.20, mu2 = 0,
                            sd1 = ifelse(paired, sqrt(1/(2*(1-paired.r))), 1), sd2 = sd1,
                            margin = 0, alpha = 0.05, paired = FALSE, paired.r = 0.50,
                            kappa = 1, n2 = NULL, power = NULL,
                            alternative = c("not equal", "greater", "less",
                                            "non-inferior", "superior", "equivalent"),
                            distribution = c("normal", "uniform", "double exponential",
                                             "laplace", "logistic"),
                            method = c("guenther", "noether"),
                            verbose = TRUE) {


  alternative <- tolower(match.arg(alternative))
  distribution <- tolower(match.arg(distribution))
  method <- tolower(match.arg(method))

  if (is.null(n2) & is.null(power))
    stop("`n2` and `power` cannot be `NULL` at the same time",
         call. = FALSE)
  if (!is.null(n2) & !is.null(power))
    stop("one of the `n2` or `power` should be `NULL`",
         call. = FALSE)
  if (!is.null(n2) & is.null(power))
    requested <- "power"
  if (is.null(n2) & !is.null(power))
    requested <- "n2"

  # wilcoxon adjustment for guenther method
  w <- switch(distribution,
              `uniform` = 1,
              `double exponential`= 2/3,
              `laplace`= 2/3,
              `logistic` = 9 / pi^2,
              `normal` = pi / 3)
  ifelse(method == "noether", w.adj <- 1, w.adj <- w)

  pwr.fun.body <- quote({

    n1 <- n2 * kappa
    if(paired) {
      psd <- sqrt(sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r)
    } else {
      psd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
    }

    d <- (mu1 - mu2) / psd
    dmargin <- margin / psd # standardize margin

    if(method == "noether") {

      if(paired) stop("specify method = 'guenther' to request Wilcoxon signed-rank test for matched pairs", call. = FALSE)

      prob0 <- 0.50 # null value
      prob1 <- switch(distribution,
                      `uniform` = d,
                      `double exponential`= 1 - exp(-(d*sqrt(2)))*(1 + (d*sqrt(2))/2) / 2,
                      `laplace`= 1 - exp(-(d*sqrt(2)))*(1 + (d*sqrt(2))/2) / 2,
                      `logistic` = (1 - (1 + (d*pi/sqrt(3))) * exp(-(d*pi/sqrt(3)))) / (1 - exp(-(d*pi/sqrt(3))))^2,
                      `normal` = pnorm(d/sqrt(2)))

      prob1margin <- switch(distribution,
                            `uniform` = dmargin,
                            `double exponential`= 1 - exp(-(dmargin*sqrt(2)))*(1 + (dmargin*sqrt(2))/2) / 2,
                            `laplace`= 1 - exp(-(dmargin*sqrt(2)))*(1 + (dmargin*sqrt(2))/2) / 2,
                            `logistic` = (1 - (1 + (dmargin*pi/sqrt(3))) * exp(-(dmargin*pi/sqrt(3)))) / (1 - exp(-(dmargin*pi/sqrt(3))))^2,
                            `normal` = pnorm(dmargin/sqrt(2)))

      propss <- kappa / (kappa + 1)

      if(alternative == "not equal") {

        HA_H0 <- prob1 - prob0
        lambda <- sqrt(n2 + n2 * kappa) * sqrt(12 * propss * (1 - propss)) * (HA_H0)
        lambda <- abs(lambda)
        power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), mean = lambda) +
          pnorm(-qnorm(alpha / 2, lower.tail = FALSE), mean = lambda)

      } else if(alternative == "greater" | alternative == "less") {

        HA_H0 <- prob1 - prob0
        lambda <- sqrt(n2 + n2 * kappa) * sqrt(12 * propss * (1 - propss)) * (HA_H0)
        lambda <- abs(lambda)
        if(alternative == "less") lambda <- abs(lambda)
        power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), mean = lambda)

      } else if(alternative == "non-inferior" | alternative == "superior") {

        HA_H0 <- prob1 - prob0 - prob1margin
        lambda <- sqrt(n2 + n2 * kappa) * sqrt(12 * propss * (1 - propss)) * (HA_H0)
        lambda <- abs(lambda)
        if(alternative == "non-inferior") lambda <- abs(lambda)
        power <- 1 - pnorm(qnorm(alpha, ncp = 0, lower.tail = FALSE), mean = lambda)

      } else if(alternative == "equivalent") {

        HA_H0 <- abs(prob1 - prob0) - prob1margin
        lambda <- sqrt(n2 + n2 * kappa) * sqrt(12 * propss * (1 - propss)) * (HA_H0)
        lambda <- abs(lambda)
        power <- 2 * (1 - pnorm(qnorm(alpha, ncp = 0, lower.tail = FALSE), mean = lambda)) - 1

      }

    } else if(method == "guenther") {

      if(paired){
        df <- n2 - 1
      } else {
        df <- n1 + n2 - 2
      }

      if(alternative == "not equal") {

        HA_H0 <- d
        if(paired){
          lambda <- HA_H0 / sqrt(1 / n2)
        } else {
          lambda <- HA_H0 / sqrt(1 / n1 + 1 / n2)
        }
        lambda <- abs(lambda)
        power <- 1 - pt(qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda) +
          pt(-qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)

      } else if(alternative == "greater" | alternative == "less") {

        HA_H0 <- d
        if(paired){
          lambda <- HA_H0 / sqrt(1 / n2)
        } else {
          lambda <- HA_H0 / sqrt(1 / n1 + 1 / n2)
        }
        lambda <- abs(lambda)
        power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)

      } else if(alternative == "non-inferior" | alternative == "superior") {

        HA_H0 <- d - dmargin
        if(paired){
          lambda <- HA_H0 / sqrt(1 / n2)
        } else {
          lambda <- HA_H0 / sqrt(1 / n1 + 1 / n2)
        }
        lambda <- abs(lambda)
        power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)

      } else if(alternative == "equivalent") {

        HA_H0 <- abs(d) - dmargin
        if(paired){
          lambda <- HA_H0 / sqrt(1 / n2)
        } else {
          lambda <- HA_H0 / sqrt(1 / n1 + 1 / n2)
        }
        lambda <- abs(lambda)
        power <- 2 * (1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = lambda)) - 1

      }

    }

    return(c(lambda, power))

  }) # end of pwr.fun.body

  # get power or sample size
  if (is.null(power)) {

    # apply wilcoxon adjustment
    n1 <- n2 * kappa
    n1 <- n1 / w
    n2 <- n2 / w

    power <- eval(pwr.fun.body)[2]
    if(power < 0) stop("design is not feasible", call. = FALSE)

    # reverse wilcoxon adjustment
    n1.star <- n1 * w
    n2.star <- n2 * w

  } else if(is.null(n2)) {

    if(tolower(method) == "noether") {
      HA_H0.min <- 0.01
      lambda.max <- 4
      propss <- kappa / (kappa + 1)
      n.tot.max <- (lambda.max / (sqrt(12 * propss * (1 - propss)) * (HA_H0.min )))^2
      n2.max <- n.tot.max / (1 + kappa)
    } else {
      HA_H0.min <- 0.01
      lambda.max <- 4
      if(paired) {
        n2.max <- (lambda.max / HA_H0.min)^2
      } else {
        n2.max <- (1 + 1 / kappa) / (HA_H0.min / lambda.max)^2
      }
    }
    # n2.max <- 1e+08
    # too big of a number throw warning in uniroot()
    # full precision may not have been achieved in 'pnt{final}'
    # because ncp is too large

    n2 <- uniroot(function(n2){
      power - eval(pwr.fun.body)[2]
    }, interval = c(2, n2.max))$root
    n1 <- n2 * kappa

    # reverse wilcoxon adjustment
    n1.star <- n1 * w
    n2.star <- n2 * w

  } # get power or sample size

  # get non-centrality parameter
  ncp <-  eval(pwr.fun.body)[1]

  # get the common language effect size (probability of superiority)
  # n1 <- n2 * kappa
  if(paired) {
    psd <- sqrt(sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r)
  } else {
    psd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  }
  d <- (mu1 - mu2) / psd
  dmargin <- margin / psd # standardize margin
  if(d - dmargin < 0) ncp <- -ncp
  prob1 <- switch(distribution,
                  `uniform` = d,
                  `double exponential`= 1 - exp(-(d*sqrt(2)))*(1 + (d*sqrt(2))/2) / 2,
                  `laplace`= 1 - exp(-(d*sqrt(2)))*(1 + (d*sqrt(2))/2) / 2,
                  `logistic` = (1 - (1 + (d*pi/sqrt(3))) * exp(-(d*pi/sqrt(3)))) / (1 - exp(-(d*pi/sqrt(3))))^2,
                  `normal` = pnorm(d/sqrt(2)))
  prob1margin <- switch(distribution,
                        `uniform` = dmargin,
                        `double exponential`= 1 - exp(-(dmargin*sqrt(2)))*(1 + (dmargin*sqrt(2))/2) / 2,
                        `laplace`= 1 - exp(-(dmargin*sqrt(2)))*(1 + (dmargin*sqrt(2))/2) / 2,
                        `logistic` = (1 - (1 + (dmargin*pi/sqrt(3))) * exp(-(dmargin*pi/sqrt(3)))) / (1 - exp(-(dmargin*pi/sqrt(3))))^2,
                        `normal` = pnorm(dmargin/sqrt(2)))
  ifelse(paired, n <- n2.star, n <- c(n1 = n1.star, n2 = n2.star))

  if(verbose) {
    cat(ifelse(paired,
               " Non-parametric Difference between Two Groups (Dependent samples) \n Wilcoxon signed-rank Test for Matched Pairs \n",
               " Non-parametric Difference between Two Groups (Independent samples) \n Mann-Whitney U or Wilcoxon Rank-sum Test \n (a.k.a Wilcoxon-Mann-Whitney Test) \n"),
        "Method:", toupper(method), "\n",
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        ifelse(paired,
               paste(" n =", ceiling(n2.star)),
               paste(" n1 =", ceiling(n1.star), "\n  n2 =", ceiling(n2.star))), "\n",
        "------------------------------ \n",
        "Alternative =", dQuote(alternative),"\n",
        "Non-centrality parameter =", round(ncp, 3), "\n",
        "Degrees of freedom =", ifelse(method == "noether", NA, round(df, 2)), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(prob1 = prob1, d = d, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2,
                                        prob1margin = prob1margin, margin = margin,
                                        kappa = kappa, alpha = alpha, paired = paired, paired.r = paired.r,
                                        alternative = alternative, verbose = verbose),
                           test = ifelse(method == "noether", "z", "t"),
                           ncp = ncp,
                           power = power,
                           n = n),
                      class = c("pwrss", "np", "2means", ifelse(method == "noether", "z", "t"))))
} # end of pwrss.np.2means()


##############################
# logistic regression z test #
##############################

# dist = c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")
# dist = list(dist = "normal", mean = 0, sd = 1)
# dist = list(dist = "poisson", lambda = 1)
# dist = list(dist = "uniform", min = 0, max = 1)
# dist = list(dist = "exponential", rate = 1)
# dist = list(dist = "binomial", size = 1, prob = 0.50)
# dist = list(dist = "bernoulli", prob = 0.50)
# dist = list(dist = "lognormal", meanlog = 0, sdlog = 1)
pwrss.z.logistic <- pwrss.z.logreg <-
  function(p1 = 0.10, p0 = 0.15,
           odds.ratio  = (p1 / (1 - p1)) / (p0 / (1 - p0)),
           beta0 = log(p0 / (1 - p0)), beta1 = log(odds.ratio),
           n = NULL, power = NULL, r2.other.x = 0,
           alpha = 0.05, alternative = c("not equal", "less", "greater"),
           method = c("demidenko(vc)", "demidenko", "hsieh"),
           distribution = "normal", verbose = TRUE) {


    user.parms.names <- names(as.list(match.call()))

    if(all(c("p0","odds.ratio", "beta0", "beta1") %in% user.parms.names)) {
      stop("specify 'p0' & 'p1' \n  or 'p0' & 'odds.ratio' \n  or 'beta0' & 'beta1'", call. = FALSE)
    }

    if(all(c("p0","p1") %in% user.parms.names)) {
      if(any(c("odds.ratio", "beta0", "beta1") %in% user.parms.names))
        stop("specify 'p0' & 'p1' \n  or 'p0' & 'odds.ratio' \n  or 'beta0' & 'beta1'", call. = FALSE)
    }

    if(all(c("p0","odds.ratio") %in% user.parms.names)) {
      if(any(c("p1","beta0", "beta1") %in% user.parms.names))
        message("ignoring any specifications to 'p1', 'beta0', or 'beta1'")
      beta0 <- log(p0 / (1 - p0))
      beta1 <- log(odds.ratio)
      p1 <- odds.ratio * (p0 / (1 - p0)) / (1 + odds.ratio * (p0 / (1 - p0)))
    }

    if(all(c("beta0","beta1") %in% user.parms.names)) {
      if(any(c("p0", "p1", "odds.ratio") %in% user.parms.names))
        message("ignoring any specifications to 'p0', 'p1', or 'odds.ratio'")
      p0 <- exp(beta0) / (1 + exp(beta0))
      odds.ratio <- exp(beta1)
      p1 <- odds.ratio * (p0 / (1 - p0)) / (1 + odds.ratio * (p0 / (1 - p0)))
    }

    alternative <- match.arg(alternative)
    method <- match.arg(method)

    if (alternative == "less" & p1 > p0)
      warning("expecting 'p1 < p0'", call. = FALSE)
    if (alternative == "greater" & p1 < p0)
      warning("expecting 'p1 > p0'", call. = FALSE)

    if(length(distribution) == 1 & is.character(distribution)) {
      distribution <- switch (tolower(distribution),
                              `normal` = list(dist = "normal", mean = 0, sd = 1),
                              `poisson` = list(dist = "poisson", lambda = 1),
                              `uniform` = list(dist = "uniform", min = 0, max = 1),
                              `exponential` = list(dist = "exponential", rate = 1),
                              `binomial` = list(dist = "binomial", size = 1, prob = 0.50),
                              `bernoulli` = list(dist = "bernoulli", size = 1, prob = 0.50),
                              `lognormal` = list(dist = "lognormal", meanlog = 0, sdlog = 1))
    } else if(is.list(distribution)){
      if(length(distribution) > 3) stop("unknown input type for 'distribution' argument", call. = FALSE)
      dist.list.names <- names(distribution)
      dist.attrib <- c(dist.list.names, tolower(distribution$dist))
      dist.invalid <- c(any(is.na(match(dist.attrib, c("dist", "normal", "mean", "sd")))),
                        any(is.na(match(dist.attrib, c("dist", "lognormal", "meanlog", "sdlog")))),
                        any(is.na(match(dist.attrib, c("dist", "uniform", "min", "max")))),
                        any(is.na(match(dist.attrib, c("dist", "exponential", "rate")))),
                        any(is.na(match(dist.attrib, c("dist", "poisson", "lambda")))),
                        any(is.na(match(dist.attrib, c("dist", "binomial", "bernoulli", "size","prob")))))
      if(all(dist.invalid == TRUE)) stop("unknown input type for 'distribution' argument", call. = FALSE)
    } else {
      stop("unknown input type for 'distribution' argument", call. = FALSE)
    }

    # asymptotic variances
    if(tolower(distribution$dist) == "normal"){

      mean <- distribution$mean
      sd <- distribution$sd

      min.norm <- qnorm(.0000001, mean = mean, sd = sd)
      max.norm <- qnorm(.9999999, mean = mean, sd = sd)

      # variance under null
      mu <- integrate(function(x)  dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x)), min.norm, max.norm)$value
      beta0.star <- log(mu / (1 - mu))
      beta1.star <- 0
      i00 <- integrate(function(x)  dnorm(x, mean = mean, sd = sd) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min.norm, max.norm)$value
      i01 <- integrate(function(x) x * dnorm(x, mean = mean, sd = sd) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min.norm, max.norm)$value
      i11 <- integrate(function(x) x^2 * dnorm(x, mean = mean, sd = sd) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min.norm, max.norm)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min.norm, max.norm)$value
      i01 <- integrate(function(x) x * dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min.norm, max.norm)$value
      i11 <- integrate(function(x) x^2 * dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min.norm, max.norm)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "poisson"){

      lambda <- distribution$lambda
      # maximum value
      max.pois <- qpois(.9999999, lambda = lambda)

      # variance under null
      mu <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))), na.rm = TRUE)
      beta0.star <- log(mu / (1 - mu))
      beta1.star <- 0
      i00 <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0.star+ beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2), na.rm = TRUE)
      i01 <- sum(sapply(0:max.pois, function(x) x * dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2), na.rm = TRUE)
      i11 <- sum(sapply(0:max.pois, function(x) x^2 * dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2), na.rm = TRUE)
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2), na.rm = TRUE)
      i01 <- sum(sapply(0:max.pois, function(x) x * dpois(x, lambda = lambda) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2), na.rm = TRUE)
      i11 <- sum(sapply(0:max.pois, function(x) x^2 * dpois(x, lambda = lambda) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2), na.rm = TRUE)
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "uniform"){

      min <- distribution$min
      max <- distribution$max

      # variance under null
      mu <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x)), min, max)$value
      beta0.star <- log(mu / (1 - mu))
      beta1.star <- 0
      i00 <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min, max)$value
      i01 <- integrate(function(x) x * dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min, max)$value
      i11 <- integrate(function(x) x^2 * dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min, max)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min, max)$value
      i01 <- integrate(function(x) x * dunif(x, min = min, max = max) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min, max)$value
      i11 <- integrate(function(x) x^2 * dunif(x, min = min, max = max) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min, max)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "exponential"){

      rate <- distribution$rate
      max.exp <- qexp(.9999999, rate = rate)

      # variance under null
      mu <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x)), 0, max.exp)$value
      beta0.star <- log(mu / (1 - mu))
      beta1.star <- 0
      i00 <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, 0, max.exp)$value
      i01 <- integrate(function(x) x * dexp(x, rate = rate) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, 0, max.exp)$value
      i11 <- integrate(function(x) x^2 * dexp(x, rate = rate) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, 0, max.exp)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, 0, max.exp)$value
      i01 <- integrate(function(x) x * dexp(x, rate = rate) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, 0, max.exp)$value
      i11 <- integrate(function(x) x^2 * dexp(x, rate = rate) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, 0, max.exp)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) %in% c("binomial","bernoulli")){

      ifelse(tolower(distribution$dist) == "bernoulli", size <- 1, size <- distribution$size)
      prob <- distribution$prob

      # variance under null
      mu <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))), na.rm = TRUE)
      beta0.star <- log(mu / (1 - mu))
      beta1.star <- 0
      i00 <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2), na.rm = TRUE)
      i01 <- sum(sapply(0:size, function(x) x * dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2), na.rm = TRUE)
      i11 <- sum(sapply(0:size, function(x) x^2 * dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2), na.rm = TRUE)
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2), na.rm = TRUE)
      i01 <- sum(sapply(0:size, function(x) x * dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2), na.rm = TRUE)
      i11 <- sum(sapply(0:size, function(x) x^2 * dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2), na.rm = TRUE)
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "lognormal"){

      meanlog <- distribution$meanlog
      sdlog <- distribution$sdlog
      min.lnorm <- qlnorm(.0000001, meanlog = meanlog, sdlog = sdlog)
      max.lnorm <- qlnorm(.9999999, meanlog = meanlog, sdlog = sdlog)

      # variance under null
      mu <- integrate(function(x)  dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x)), min.lnorm, max.lnorm)$value
      beta0.star <- log(mu / (1 - mu))
      beta1.star <- 0
      i00 <- integrate(function(x)  dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min.lnorm, max.lnorm)$value
      i01 <- integrate(function(x) x * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min.lnorm, max.lnorm)$value
      i11 <- integrate(function(x) x^2 * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0.star + beta1.star * x) / (1 + exp(beta0.star + beta1.star * x))^2, min.lnorm, max.lnorm)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min.lnorm, max.lnorm)$value
      i01 <- integrate(function(x) x * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min.lnorm, max.lnorm)$value
      i11 <- integrate(function(x) x^2 * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x))^2, min.lnorm, max.lnorm)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    # Demidenko, E. (2007). Sample size determination for logistic
    # regression revisited. Statistics in Medicine, 26, 3385-3397.
    pwr.fun.body <- quote({

      if(tolower(method) == "demidenko(vc)") {
        # correction factor
        vcf <- switch (tolower(distribution$dist),
                       `normal` = 1,
                       `poisson` = 1,
                       `uniform` = 1,
                       `exponential` = 1,
                       `binomial` = 0.85,
                       `bernoulli` = 0.85,
                       `lognormal` = 0.75)
      } else if (tolower(method) == "demidenko") {
        vcf <- 0
      }

      # non-centrality parameter and standard deviation of the non-centrality parameter under alternative
      ncp <- beta1 / sqrt(var.beta1 / (n * (1 - r2.other.x)))
      sd.ncp <- sqrt((vcf * var.beta0 + (1 - vcf) * var.beta1) / var.beta1)

      # compute power
      if(alternative == "not equal") {
        power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), mean = ncp, sd = sd.ncp) +
          pnorm(-qnorm(alpha / 2, lower.tail = FALSE), mean = ncp, sd = sd.ncp)

      } else if(alternative == "greater" | alternative == "less") {
        if(alternative == "less") ncp <- abs(ncp)
        power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), mean = ncp, sd = sd.ncp)

      }

      power

    })

    # Hsieh, F. Y., Bloch, D. A., & Larsen, M. D. (1998). A simple
    # method of sample size calculation for linear and logistic
    # regression. Statistics in Medicine, 17, 1623-1634.
    n.fun.body <- quote({

      if(tolower(distribution$dist) %in% c("binomial","bernoulli")) {

        if(tolower(distribution$dist) == "binomial" & distribution$size > 1)
          stop("Hsieh et al. (1998) is valid only for a binary covariate or a continuous covariate following normal distribution")
        prob <- distribution$prob
        beta <- 1 - power
        ifelse(tolower(alternative) == "not equal",
               z.alpha <- qnorm(alpha / 2, lower.tail = FALSE),
               z.alpha <- qnorm(alpha, lower.tail = FALSE))
        z.beta <- qnorm(beta, lower.tail = FALSE)
        p.bar <- (1 - prob) * p0 + prob * p1
        n <- (z.alpha * sqrt(p.bar * (1 - p.bar) / prob) + z.beta * sqrt(p0 * (1 - p0) + p1 * (1 - p1) * (1 - prob) / prob))^2 / ((p0 - p1)^2 * (1 - prob))
        n <- n / (1 - r2.other.x)

      } else if (tolower(distribution$dist) == "normal") {

        beta <- 1 - power
        ifelse(tolower(alternative) == "not equal",
               z.alpha <- qnorm(alpha / 2, lower.tail = FALSE),
               z.alpha <- qnorm(alpha, lower.tail = FALSE))
        z.beta <- qnorm(beta, lower.tail = FALSE)
        odds.ratio <- (p1 / (1 - p1)) / (p0 / (1 - p0))
        beta1 <- log(odds.ratio)
        n <- (z.alpha + z.beta)^2 / (p0 * (1 - p0) * beta1^2)
        n <- n / (1 - r2.other.x)

      } else {

        stop("not a valid distribution for Hsieh et al. (1998) procedure", call. = FALSE)

      }

      n

    })


    if(tolower(method) == "demidenko(vc)" | tolower(method) == "demidenko") {

      if (is.null(power)) {
        power <- eval(pwr.fun.body)
      } else if(is.null(n)) {
        n <- uniroot(function(n){
          power - eval(pwr.fun.body)
        }, interval = c(2, 1e+09))$root
      }

      ncp <- beta1 / sqrt(var.beta1 / (n * (1 - r2.other.x)))

    } else if (tolower(method) == "hsieh") {

      if (is.null(n)) {
        n <- eval(n.fun.body)
      } else if(is.null(power)) {
        power <- uniroot(function(power){
          n - eval(n.fun.body)
        }, interval = c(0.01, 0.999))$root
      }

      eval(n.fun.body)
      ifelse(tolower(alternative) == "not equal",
             z.alpha <- qnorm(alpha / 2, lower.tail = FALSE),
             z.alpha <- qnorm(alpha, lower.tail = FALSE))
      z.beta <- qnorm(beta, lower.tail = FALSE)
      ncp <- z.alpha + z.beta
      ifelse(p1 < p0, ncp <- -ncp, ncp <- ncp)

    } else {

      stop("unknown method", call. = FALSE)

    }

    if(verbose) {
      cat(" Logistic Regression Coefficient \n (Large Sample Approx. Wald's z Test) \n",
          switch(alternative,
                 `not equal` = "H0: beta1 = 0 \n HA: beta1 != 0 \n",
                 `greater` = "H0: beta1 = 0 \n HA: beta1 > 0 \n",
                 `less` = "H0: beta1 = 0 \n HA: beta1 < 0 \n"),
          "Distribution of X =", sQuote(tolower(distribution$dist)), "\n",
          "Method =", toupper(method), "\n",
          "------------------------------ \n",
          " Statistical power =", round(power, 3), "\n",
          " n =", ceiling(n), "\n",
          "------------------------------ \n",
          "Alternative =", dQuote(alternative),"\n",
          "Non-centrality parameter =", round(ncp, 3), "\n",
          "Type I error rate =", round(alpha, 3), "\n",
          "Type II error rate =", round(1 - power, 3), "\n")
    }

    invisible(structure(list(parms = list(p0 = p0, p1 = p1, beta0 = beta0, beta1 = beta1,
                                          odds.ratio = odds.ratio, r2.other.x = r2.other.x,
                                          alpha = alpha, alternative = alternative, method = method,
                                          distribution =  distribution, verbose = verbose),
                             test = "z",
                             ncp = ncp,
                             power = power,
                             n = n),
                        class = c("pwrss", "z", "logreg")))

  } # end of pwrss.z.logistic()

#############################
# poisson regression z test #
#############################

# dist = c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")
# dist = list(dist = "normal", mean = 0, sd = 1)
# dist = list(dist = "poisson", lambda = 1)
# dist = list(dist = "uniform", min = 0, max = 1)
# dist = list(dist = "exponential", rate = 1)
# dist = list(dist = "binomial", size = 1, prob = 0.50)
# dist = list(dist = "bernoulli", prob = 0.50)
# dist = list(dist = "lognormal", meanlog = 0, sdlog = 1)

# mean.exposure is the mean exposure time (should be > 0)
# lambda is the mean event rate for the exposure time
pwrss.z.poisson <- pwrss.z.poisreg <-
  function(exp.beta0 = 1.10, exp.beta1 = 1.16,
           beta0 = log(exp.beta0), beta1 = log(exp.beta1),
           mean.exposure = 1, n = NULL, power = NULL, r2.other.x = 0,
           alpha = 0.05, alternative = c("not equal", "less", "greater"),
           method = c("demidenko(vc)", "demidenko", "signorini"),
           distribution = "normal", verbose = TRUE) {


    user.parms.names <- names(as.list(match.call()))

    if(all(c("beta0","betap1") %in% user.parms.names)) {
      if(any(c("exp.beta0", "exp.beta1") %in% user.parms.names))
        message("ignoring any specifications to 'exp.beta0', or 'exp.beta1'")
    }

    alternative <- match.arg(alternative)
    method <- match.arg(method)

    if (alternative == "less" & beta1 > 0)
      warning("expecting 'beta1 < 0'", call. = FALSE)
    if (alternative == "greater" & beta1 < 0)
      warning("expecting 'beta1 > 0'", call. = FALSE)

    if(length(distribution) == 1 & is.character(distribution)) {
      distribution <- switch (tolower(distribution),
                              `normal` = list(dist = "normal", mean = 0, sd = 1),
                              `poisson` = list(dist = "poisson", lambda = 1),
                              `uniform` = list(dist = "uniform", min = 0, max = 1),
                              `exponential` = list(dist = "exponential", rate = 1),
                              `binomial` = list(dist = "binomial", size = 1, prob = 0.50),
                              `bernoulli` = list(dist = "bernoulli", size = 1, prob = 0.50),
                              `lognormal` = list(dist = "lognormal", meanlog = 0, sdlog = 1))
    } else if(is.list(distribution)){
      if(length(distribution) > 3) stop("unknown input type for 'distribution' argument", call. = FALSE)
      dist.list.names <- names(distribution)
      dist.attrib <- c(dist.list.names, tolower(distribution$dist))
      dist.invalid <- c(any(is.na(match(dist.attrib, c("dist", "normal", "mean", "sd")))),
                        any(is.na(match(dist.attrib, c("dist", "lognormal", "meanlog", "sdlog")))),
                        any(is.na(match(dist.attrib, c("dist", "uniform", "min", "max")))),
                        any(is.na(match(dist.attrib, c("dist", "exponential", "rate")))),
                        any(is.na(match(dist.attrib, c("dist", "poisson", "lambda")))),
                        any(is.na(match(dist.attrib, c("dist", "binomial", "bernoulli", "size","prob")))))
      if(all(dist.invalid == TRUE)) stop("unknown input type for 'distribution' argument", call. = FALSE)
    } else {
      stop("unknown input type for 'distribution' argument", call. = FALSE)
    }

    # asymptotic variances
    if(tolower(distribution$dist) == "normal"){

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
      i11 <- integrate(function(x) x^2 * dnorm(x, mean = mean, sd = sd) * exp(beta0.star + beta1.star * x), min.norm, max.norm)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x), min.norm, max.norm)$value
      i01 <- integrate(function(x) x * dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x), min.norm, max.norm)$value
      i11 <- integrate(function(x) x^2 * dnorm(x, mean = mean, sd = sd) * exp(beta0 + beta1 * x), min.norm, max.norm)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "poisson"){

      lambda <- distribution$lambda
      # maximum value
      max.pois <- qpois(.999999999, lambda = lambda)

      # variance under null
      mu <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:max.pois, function(x) x * dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:max.pois, function(x) x^2 * dpois(x, lambda = lambda) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- sum(sapply(0:max.pois, function(x)  dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:max.pois, function(x) x * dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:max.pois, function(x) x^2 * dpois(x, lambda = lambda) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "uniform"){

      min <- distribution$min
      max <- distribution$max

      # variance under null
      mu <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x), min, max)$value
      i01 <- integrate(function(x) x * dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x), min, max)$value
      i11 <- integrate(function(x) x^2 * dunif(x, min = min, max = max) * exp(beta0.star + beta1.star * x), min, max)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      i01 <- integrate(function(x) x * dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      i11 <- integrate(function(x) x^2 * dunif(x, min = min, max = max) * exp(beta0 + beta1 * x), min, max)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "exponential"){

      rate <- distribution$rate
      max.exp <- qexp(.9999999, rate = rate)

      # variance under null
      mu <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0.star + beta1.star * x), 0, max.exp)$value
      i01 <- integrate(function(x) x * dexp(x, rate = rate) * exp(beta0.star + beta1.star * x), 0, max.exp)$value
      i11 <- integrate(function(x) x^2 * dexp(x, rate = rate) * exp(beta0.star + beta1.star * x), 0, max.exp)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      i01 <- integrate(function(x) x * dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      i11 <- integrate(function(x) x^2 * dexp(x, rate = rate) * exp(beta0 + beta1 * x), 0, max.exp)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) %in% c("binomial","bernoulli")){

      ifelse(tolower(distribution$dist) == "bernoulli", size <- 1, size <- distribution$size)
      prob <- distribution$prob

      # variance under null
      mu <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      beta0.star <- log(mu)
      beta1.star <- 0
      i00 <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:size, function(x) x * dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:size, function(x) x^2 * dbinom(x, size = size, prob = prob) * exp(beta0.star + beta1.star * x)), na.rm = TRUE)
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- sum(sapply(0:size, function(x) dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i01 <- sum(sapply(0:size, function(x) x * dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      i11 <- sum(sapply(0:size, function(x) x^2 * dbinom(x, size = size, prob = prob) * exp(beta0 + beta1 * x)), na.rm = TRUE)
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    if(tolower(distribution$dist) == "lognormal"){

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
      i11 <- integrate(function(x) x^2 * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0.star + beta1.star * x), min.lnorm, max.lnorm)$value
      var.beta0 <- i00 / (i00 * i11 - i01^2)

      # variance under alternative
      i00 <- integrate(function(x)  dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x), min.lnorm, max.lnorm)$value
      i01 <- integrate(function(x) x * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x), min.lnorm, max.lnorm)$value
      i11 <- integrate(function(x) x^2 * dlnorm(x, meanlog = meanlog, sdlog = sdlog) * exp(beta0 + beta1 * x), min.lnorm, max.lnorm)$value
      var.beta1 <- i00 / (i00 * i11 - i01^2)

    }

    pwr.fun.body <- quote({

      if(tolower(method) == "demidenko(vc)") {
        # correction factor
        vcf <- switch (tolower(distribution$dist),
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

      # non-centrality parameter and standard deviation of the non-centrality parameter under alternative
      if (tolower(method) == "signorini") {
        # Signorini, D. F. (1991). Sample size for poisson regression.
        # Biometrika, 78, 446-450.
        ncp <- beta1 / sqrt(var.beta0 / (n * (1 - r2.other.x) * mean.exposure))
        sd.ncp <- sqrt(var.beta1 / var.beta0)
      } else {
        # Demidenko, E. (2007). Sample size determination for logistic
        # regression revisited. Statistics in Medicine, 26, 3385-3397.
        ncp <- beta1 / sqrt(var.beta1 / (n * (1 - r2.other.x) * mean.exposure))
        sd.ncp <- sqrt((vcf * var.beta0 + (1 - vcf) * var.beta1) / var.beta1)
      }

      # compute power
      if(alternative == "not equal") {
        power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), mean = ncp, sd = sd.ncp) +
          pnorm(-qnorm(alpha / 2, lower.tail = FALSE), mean = ncp, sd = sd.ncp)

      } else if(alternative == "greater" | alternative == "less") {
        if(alternative == "less") ncp <- abs(ncp)
        power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), mean = ncp, sd = sd.ncp)

      }

      power

    })

    if (is.null(power)) {
      power <- eval(pwr.fun.body)
    } else if(is.null(n)) {
      n <- uniroot(function(n){
        power - eval(pwr.fun.body)
      }, interval = c(2, 1e+09))$root
    }

    ncp <- beta1 / sqrt(var.beta1 / (n * (1 - r2.other.x)))

    if(verbose) {
      cat(" Poisson Regression Coefficient \n (Large Sample Approx. Wald's z Test) \n",
          switch(alternative,
                 `not equal` = "H0: beta1 = 0 \n HA: beta1 != 0 \n",
                 `greater` = "H0: beta1 = 0 \n HA: beta1 > 0 \n",
                 `less` = "H0: beta1 = 0 \n HA: beta1 < 0 \n"),
          "Distribution of X =", sQuote(tolower(distribution$dist)), "\n",
          "Method =", toupper(method), "\n",
          "------------------------------ \n",
          " Statistical power =", round(power, 3), "\n",
          " n =", ceiling(n), "\n",
          "------------------------------ \n",
          "Alternative =", dQuote(alternative),"\n",
          "Non-centrality parameter =", round(ncp, 3), "\n",
          "Type I error rate =", round(alpha, 3), "\n",
          "Type II error rate =", round(1 - power, 3), "\n")
    }

    invisible(structure(list(parms = list(beta0 = beta0, beta1 = beta1, mean.exposure = mean.exposure,
                                          r2.other.x = r2.other.x,
                                          alpha = alpha, alternative = alternative, method = method,
                                          distribution =  distribution, verbose = verbose),
                             test = "z",
                             ncp = ncp,
                             power = power,
                             n = n),
                        class = c("pwrss", "z", "poisreg")))

  } # end of pwrss.z.poisson()

####################################################
# Chi-square test for independence/goodness-of-fit #
####################################################

# internal function to get some chisq stat
.chisq.fun <- function(p1) {
  ifelse(is.vector(p1),
         p0 <- rep(1/length(p1), length(p1)),
         p0 <- outer(rowSums(p1), colSums(p1)) / sum(p1))
  ifelse(is.vector(p1),
         df <- length(p1) - 1,
         df <- (nrow(p1) - 1) * (ncol(p1) - 1))
  chisq <- sum((p1 - p0)^2 / p0)
  ifelse(is.vector(p1),
         mdf <- 1,
         mdf <- min(nrow(p1) - 1, ncol(p1) - 1))
  w <- sqrt(chisq / (sum(p1) * mdf))

  list(p1 = p1, p0 = p0, df = df, w = w)
}

pwrss.chisq.gofit <- function (p1 = c(0.50, 0.50),
                               p0 = .chisq.fun(p1)$p0,
                               w = .chisq.fun(p1)$w,
                               df = .chisq.fun(p1)$df,
                               n = NULL, power = NULL,
                               alpha = 0.05, verbose = TRUE)
{

  user.parms.names <- names(as.list(match.call()))
  if("p1" %in% user.parms.names & "w" %in% user.parms.names) {
    warning("ignoring any specifications to 'p1', or 'p0'")
  }
  if("w" %in% user.parms.names & !("df" %in% user.parms.names)) {
    stop("specify 'df'", call. = FALSE)
  }

  if(is.vector(p1)) {

    if(length(p1) != length(p0)) stop("length of 'p1' and 'p0' should match", call. = FALSE)
    if(sum(p1) != 1 | sum(p0) != 1) stop("cell probabilities should sum to 1", call. = FALSE)

  } else if(is.matrix(p1)) {

    if(any(dim(p1) != dim(p0))) stop("dimensions for 'p1' and 'p0' differ", call. = FALSE)

  } else {

    stop("not a valid 'p1' argument", call. = FALSE)

  }


  pwr.fun.body <- quote({
    lambda <- n * w^2
    pchisq(qchisq(alpha, df = df, lower = FALSE), df = df, ncp = lambda, lower = FALSE)
  })


  if (is.null(power)) {
    power <- eval(pwr.fun.body)
  } else if(is.null(n)) {
    n <- uniroot(function(n){
      power - eval(pwr.fun.body)
    }, interval = c(df + 1, 1e+09))$root
    lambda <- n * w^2
  }

  if(verbose) {
    cat(" Pearson's Chi-square Goodness-of-fit Test \n for Contingency Tables \n",
        "------------------------------ \n",
        " Statistical power =", round(power, 3), "\n",
        " Total n =", ceiling(n), "\n",
        "------------------------------ \n",
        "Degrees of freedom =", round(df, 3), "\n",
        "Non-centrality parameter =", round(lambda, 3), "\n",
        "Type I error rate =", round(alpha, 3), "\n",
        "Type II error rate =", round(1 - power, 3), "\n")
  }

  invisible(structure(list(parms = list(w = w, df = df, alpha = alpha, verbose = verbose),
                           test = "chisq",
                           df = df,
                           ncp = lambda,
                           power = power,
                           n = n),
                      class = c("pwrss", "chisq", "gofit")))

}
# end of pwrss.chisq.gofit()

##################################
# generic t and z test functions #
##################################

# type = 1 for light red, 2 for light blue
.plot.t.dist <- function(ncp = 0, df = Inf, xlim, type = 1) {


  plot.window.dim <- dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  ifelse(type == 1,
         color <- adjustcolor(2, alpha.f = 1),
         color <- adjustcolor(4, alpha.f = 1))

  # non-central t function
  funt <- function(x){
    dt(x, df = df, ncp = ncp)
  }

  # plot central t distribution
  ylim <- c(0, 0.50)
  plot(funt, xlim = xlim, ylim = ylim,
       xaxs = "i", yaxs = "i", bty = "l",
       col = color, lwd = 2, lty = type,
       xlab = "", ylab = "", cex.axis = cex.axis)

}

# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.t.dist<- function(ncp = 0, df = Inf, xlim, type = 1) {

  color <- switch (type,
                   `1` = adjustcolor(2, alpha.f = 0.3),
                   `2` = adjustcolor(4, alpha.f = 0.3),
                   `3` = adjustcolor(1, alpha.f = 0.3))

  # non-central t function
  funt <- function(x){
    dt(x, df = df, ncp = ncp)
  }

  x <- seq(min(xlim), max(xlim), by = .001)
  y <- funt(x)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if(type == 1 | type == 2) {
    polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    polygon(x = xs, y = ys, col = color, density = 20, angle = 45, border = NA)
  }

  prob <- pt(max(xlim), df = df, ncp = ncp, lower.tail = TRUE) -
    pt(min(xlim), df = df, ncp = ncp, lower.tail = TRUE)

  return(invisible(prob))

}

.plot.t.t1t2 <- function(ncp, df, alpha, alternative,
                         plot.main = NULL, plot.sub = NULL) {

  # critical t line segment coordinates
  if(alternative == "equivalent") {

    if(length(ncp) == 1) ncp <- rbind(min(-ncp, ncp), max(-ncp, ncp))
    if(length(ncp) == 2 & is.vector(ncp)) ncp <- rbind(min(ncp), max(ncp))
    if(length(ncp) > 2) stop("not a valid plotting option", call. = FALSE)

    talpha <- c(qt(alpha, df = df, ncp = ncp[1], lower.tail = FALSE),
                qt(alpha, df = df, ncp = ncp[2], lower.tail = TRUE))
    ytalpha <- dt(talpha, df = df, ncp = ncp)

  } else if(alternative == "not equal") {

    if(length(ncp) > 1) stop("not a valid plotting option", call. = FALSE)

    talpha <- qt(c(1 - alpha / 2, alpha / 2), df = df, ncp = 0, lower.tail = FALSE)
    ytalpha <- dt(talpha, df = df, ncp = 0)

  } else {

    if(length(ncp) > 1) stop("not a valid plotting option", call. = FALSE)

    talpha <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE)
    ytalpha <- dt(talpha, df = df, ncp = 0)
    if(ncp < 0) talpha <- -talpha

  }

  # x-axis limits
  ifelse(df < 20, prob.extreme <- 0.001, prob.extreme <- 0.0001)
  lower <- min(min(qt(prob.extreme, df = df, ncp = ncp, lower.tail = TRUE)),
               qt(prob.extreme, df = df, ncp = 0, lower.tail = TRUE))
  upper <- max(max(qt(1-prob.extreme, df = df, ncp = ncp, lower.tail = TRUE)),
               qt(1-prob.extreme, df = df, ncp = 0, lower.tail = TRUE))
  xlim <- c(lower, upper )

  plot.window.dim <- dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  if(alternative == "equivalent") {

    .plot.t.dist(ncp = ncp[1], df = df, xlim = xlim, type = 1)
    par(new = TRUE)
    .plot.t.dist(ncp = ncp[2], df = df, xlim = xlim, type = 1)
    par(new = TRUE)
    .plot.t.dist(ncp = 0, df = df, xlim = xlim, type = 2)

    text(0, dt(0, df = df, ncp = 0) + 0.05,
         labels = "Alternative \n Hypothesis",
         cex = cex.legend, col = adjustcolor(4, alpha.f = 1))

  } else {

    .plot.t.dist(ncp = ncp, df = df, xlim = xlim, type = 2)
    par(new = TRUE)
    .plot.t.dist(ncp = 0, df = df, xlim = xlim, type = 1)

    text(ncp, dt(ncp, df = df, ncp = ncp) + 0.05,
         labels = "Alternative \n Hypothesis",
         cex = cex.legend, col = adjustcolor(4, alpha.f = 1))

  } # end of plots

  # draw vertical lines for critical region
  segments(x0 = talpha, y0 = rep(0, length(talpha)),
           x1 = talpha, y1 = ytalpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  if(alternative == "equivalent") {

    .paint.t.dist(ncp = ncp[1], df = df, xlim = c(talpha[1], max(xlim)), type = 1)
    .paint.t.dist(ncp = ncp[2], df = df, xlim = c(talpha[2], min(xlim)), type = 1)

    .paint.t.dist(ncp = 0, df = df, xlim = c(talpha[1], min(xlim)), type = 2)
    .paint.t.dist(ncp = 0, df = df, xlim = c(talpha[2], max(xlim)), type = 2)

    ifelse(talpha[2] > talpha[1],
           power <- .paint.t.dist(ncp = 0, df = df, xlim = talpha, type = 3),
           power <- 0)

  } else if(alternative == "not equal") {

    .paint.t.dist(ncp = 0, df = df, xlim = c(talpha[1], min(xlim)), type = 1)
    .paint.t.dist(ncp = 0, df = df, xlim = c(talpha[2], max(xlim)), type = 1)

    ifelse(ncp < 0,
           .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha[1], max(xlim)), type = 2),
           .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha[2], min(xlim)), type = 2))

    ifelse(ncp < 0,
           power <- .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha[1], min(xlim)), type = 3),
           power <- .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha[2], max(xlim)), type = 3))

  } else {

    ifelse(ncp < 0,
           .paint.t.dist(ncp = 0, df = df, xlim = c(talpha, min(xlim)), type = 1),
           .paint.t.dist(ncp = 0, df = df, xlim = c(talpha, max(xlim)), type = 1))

    ifelse(ncp < 0,
           .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha, max(xlim)), type = 2),
           .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha, min(xlim)), type = 2))

    ifelse(ncp < 0,
           power <- .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha, min(xlim)), type = 3),
           power <- .paint.t.dist(ncp = ncp, df = df, xlim = c(talpha, max(xlim)), type = 3))

  } # end of paint regions

  # axes labels and subtitle
  title(main = plot.main, line = 2, cex.main = cex.title)
  title(sub = plot.sub, line = 3, cex.sub = cex.title)
  title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))
  if(is.finite(df)) {
    title(xlab = paste0("t Value (df = ", round(df, digits = 2), ")"),
          line = 2.2, cex.lab = cex.label,
          col.lab = adjustcolor(1, alpha.f = 0.8))
  } else {
    title(xlab = "z Value", line = 2.2, cex.lab = cex.label,
          col.lab = adjustcolor(1, alpha.f = 0.8))
  }

  if(power < 0) power <- 0
  alpha <- round(alpha, 2)
  beta <- round(1 - power, 2)
  power <- round(power, 2)

  legend("topright", cex = cex.legend,
         c(as.expression(bquote(Power == .(power))),
           as.expression(bquote(alpha == .(alpha))),
           as.expression(bquote(beta == .(beta)))),
         fill = c(adjustcolor(1, alpha.f = 0.3),
                  adjustcolor(2, alpha.f = 0.3),
                  adjustcolor(4, alpha.f = 0.3)),
         border = c(adjustcolor(1, alpha.f = 0.15),
                    adjustcolor(2, alpha.f = 0.15),
                    adjustcolor(4, alpha.f = 0.15)),
         bg = adjustcolor(1, alpha.f = 0.08),
         box.col = adjustcolor(1, alpha.f = 0),
         density = c(30, NA, NA),
         angle = c(45, NA, NA))

} # end of .plot.t.t1t2()

# power for the generic t test with (optional) type I and type II error plots
power.t.test <- function(ncp, df, alpha = 0.05,
                         alternative = c("not equal", "greater", "less",
                                         "non-inferior", "superior", "equivalent"),
                         plot = TRUE, plot.main = NULL, plot.sub = NULL,
                         verbose = TRUE){

  alternative <- tolower(match.arg(alternative))

  if(sum(c(length(ncp), length(df), length(alpha)) > 1) > 1)
    stop("only one of the 'ncp', 'df', or 'alpha' arguments can take multiple values", call. = FALSE)

  if(any(df < 0))
    stop("'df' sould be positive", call. = FALSE)

  if(any(alpha < 0) | any(alpha > .999))
    stop("'alpha' sould be between 0 and 0.999", call. = FALSE)

  # calculate statistical power
  if(alternative == "not equal") {

    if(any(c(length(ncp), length(df), length(alpha)) > 1) & isTRUE(plot))
      stop("plotting is not available for multiple values", call. = FALSE)

    power <- 1 - pt(qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = abs(ncp)) +
      pt(-qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = abs(ncp))

    talpha <- rbind(qt(1 - alpha / 2, df = df, ncp = 0, lower.tail = FALSE),
                    qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE))

  } else if(alternative %in% c("greater", "less", "superior", "non-inferior")) {

    if(any(ncp < 0) & alternative == "greater")
      warning("alternative = 'greater' but non-centrality parameter is negative", call. = FALSE)

    if(any(ncp > 0) & alternative == "less")
      warning("alternative = 'less' but non-centrality parameter is positive", call. = FALSE)

    if(any(c(length(ncp), length(df), length(alpha)) > 1) & isTRUE(plot))
      stop("plotting is not available for multiple values", call. = FALSE)

    power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = abs(ncp))

    talpha <- qt(alpha, df = df, ncp = 0, lower.tail = FALSE)

  } else if(alternative == "equivalent") {

    if(any(c(length(df), length(alpha)) > 1) & isTRUE(plot))
      stop("plotting is not available for multiple values", call. = FALSE)

    if(length(ncp) == 1) {

      ncp <- rbind(min(-ncp, ncp), max(-ncp, ncp))

    } else {

      if(is.vector(ncp)) {
        if(length(ncp) > 1 & isTRUE(plot))
          stop("plotting is not available for multiple values", call. = FALSE)
        ncp <- rbind(-ncp, ncp)
      }

      if(is.matrix(ncp)) {
        if(dim(ncp)[1] > 2)
          stop("equivalence testing allows multiple 'ncp' parameters in the form of 2 x k matrix \n the first row is the lower bound of 'ncp' and the second row is the upper bound of 'ncp'", call. = FALSE)
        if(dim(ncp)[2] > 1 & isTRUE(plot))
          stop("plotting is not available for multiple values", call. = FALSE)
      }

    }

    power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = abs(ncp[1,])) +
      1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = abs(ncp[2,])) - 1

    # do not rely on lower.tail defaults
    # the power function above is same as the following
    ## beta1 <- 1 - pt(qt(alpha, df = df, ncp = ncp[1,], lower.tail = FALSE),
    ##                 df = df, ncp = 0, lower.tail = FALSE)
    ## beta2 <- 1 - pt(qt(alpha, df = df, ncp = ncp[2,], lower.tail = TRUE),
    ##                 df = df, ncp = 0, lower.tail = TRUE)
    ## power <- 1 - beta1 + 1 - beta2 - 1
    # which is more intuative based on type I and type II error plots

    power[power < 0] <- 0

    talpha <- rbind(qt(alpha, df = df, ncp = ncp[1,], lower.tail = FALSE),
                    qt(alpha, df = df, ncp = ncp[2,], lower.tail = TRUE))

  } else {

    stop("not a valid hypothesis type", call. = FALSE)

  }

  if(plot) {

    .plot.t.t1t2(ncp = ncp, df = df, alpha = alpha, alternative = alternative,
                 plot.main = plot.main, plot.sub = plot.sub)

  }

  if(verbose) {

    if(alternative == "equivalent") {
      if(any(c(is.matrix(df), is.matrix(alpha))))
        stop("'verbose' argument is not available for matrix type input", call. = FALSE)
      ncp.alt <- 0
      ncp.null <- round(ncp, 3)
      print(
        data.frame(power = power, ncp.alt = ncp.alt,
                   ncp.null.1 = ncp.null[1,], ncp.null.2 = ncp.null[2,],
                   alpha = alpha, df = df,
                   t.crit.1 = talpha[1,], t.crit.2 = talpha[2,]),
        row.names = FALSE)
    } else {
      if(any(c(is.matrix(ncp), is.matrix(df), is.matrix(alpha))))
        stop("'verbose' argument is not available for matrix type input", call. = FALSE)
      ncp.alt <- round(ncp, 3)
      ncp.null <- 0
      ifelse(alternative == "not equal",
             print(
               data.frame(power = power, ncp.alt = ncp.alt, ncp.null = ncp.null,
                          alpha = alpha, df = df,
                          t.crit.1 = talpha[1,], t.crit.2 = talpha[2,]),
               row.names = FALSE),
             print(
               data.frame(power = power, ncp.alt = ncp.alt, ncp.null = ncp.null,
                          alpha = alpha, df = df, t.crit = talpha),
               row.names = FALSE))
    }

  } # end of verbose

  return(invisible(power))

} # end of plot.t.test()


# power for the generic z test with (optional) type I and type II error plots
power.z.test <- function(ncp, alpha = 0.05,
                         alternative = c("not equal", "greater", "less",
                                         "non-inferior", "superior", "equivalent"),
                         plot = TRUE, plot.main = NULL, plot.sub = NULL,
                         verbose = TRUE){

  alternative <- tolower(match.arg(alternative))

  if(sum(c(length(ncp), length(alpha)) > 1) > 1)
    stop("only one of the 'ncp' or 'alpha' arguments can take multiple values", call. = FALSE)

  if(any(alpha < 0) | any(alpha > .999))
    stop("'alpha' sould be between 0 and 0.999", call. = FALSE)

  # calculate statistical power
  if(alternative == "not equal") {

    if(any(c(length(ncp), length(alpha)) > 1) & isTRUE(plot))
      stop("plotting is not available for multiple values", call. = FALSE)

    power <- 1 - pnorm(qnorm(alpha / 2, mean = 0, sd = 1, lower.tail = FALSE), sd = 1, mean = abs(ncp)) +
      pnorm(-qnorm(alpha / 2, mean = 0, sd = 1, lower.tail = FALSE), sd = 1, mean = abs(ncp))

    zalpha <- rbind(qnorm(1 - alpha / 2, mean = 0, lower.tail = FALSE),
                    qnorm(alpha / 2, mean = 0, lower.tail = FALSE))

  } else if(alternative %in% c("greater", "less", "superior", "non-inferior")) {

    if(any(ncp < 0) & alternative == "greater")
      warning("alternative = 'greater' but non-centrality parameter is negative", call. = FALSE)

    if(any(ncp > 0) & alternative == "less")
      warning("alternative = 'less' but non-centrality parameter is positive", call. = FALSE)

    if(any(c(length(ncp), length(alpha)) > 1) & isTRUE(plot))
      stop("plotting is not available for multiple values", call. = FALSE)

    power <- 1 - pnorm(qnorm(alpha, mean = 0, sd = 1, lower.tail = FALSE), sd = 1, mean = abs(ncp))

    zalpha <- qnorm(alpha, mean = 0, lower.tail = FALSE)

  } else if(alternative == "equivalent") {

    if(length(alpha) > 1 & isTRUE(plot))
      stop("plotting is not available for multiple values", call. = FALSE)

    if(length(ncp) == 1) {

      ncp <- rbind(min(-ncp, ncp), max(-ncp, ncp))

    } else {

      if(is.vector(ncp)) {
        if(length(ncp) > 1 & isTRUE(plot))
          stop("plotting is not available for multiple values", call. = FALSE)
        ncp <- rbind(-ncp, ncp)
      }

      if(is.matrix(ncp)) {
        if(dim(ncp)[1] > 2)
          stop("equivalence testing allows multiple 'ncp' parameters in the form of 2 x k matrix \n the first row is the lower bound of 'ncp' and the second row is the upper bound of 'ncp'", call. = FALSE)
        if(dim(ncp)[2] > 1 & isTRUE(plot))
          stop("plotting is not available for multiple values", call. = FALSE)
      }

    }

    power <- 1 - pnorm(qnorm(alpha, mean = 0, sd = 1, lower.tail = FALSE), sd = 1, mean = abs(ncp[1,])) +
      1 - pnorm(qnorm(alpha, mean = 0, sd = 1, lower.tail = FALSE), sd = 1, mean = abs(ncp[2,])) - 1

    power[power < 0] <- 0

    zalpha <- rbind(qnorm(alpha, mean = ncp[1,], lower.tail = FALSE),
                    qnorm(alpha, mean = ncp[2,], lower.tail = TRUE))


  } else {

    stop("not a valid hypothesis type", call. = FALSE)

  }

  if(plot) {

    .plot.t.t1t2(ncp = ncp, df = Inf, alpha = alpha, alternative = alternative,
                 plot.main = plot.main, plot.sub = plot.sub)

  }

  if(verbose) {

    if(alternative == "equivalent") {
      if(any(c(is.matrix(df), is.matrix(alpha))))
        stop("'verbose' argument is not available for matrix type input", call. = FALSE)
      ncp.alt <- 0
      ncp.null <- round(ncp, 3)
      print(
        data.frame(power = power, ncp.alt = ncp.alt,
                   ncp.null.1 = ncp.null[1,], ncp.null.2 = ncp.null[2,],
                   alpha = alpha,
                   z.crit.1 = zalpha[1,], z.crit.2 = zalpha[2,]),
        row.names = FALSE)
    } else {
      if(any(c(is.matrix(ncp), is.matrix(df), is.matrix(alpha))))
        stop("'verbose' argument is not available for matrix type input", call. = FALSE)
      ncp.alt <- round(ncp, 3)
      ncp.null <- 0
      ifelse(alternative == "not equal",
             print(
               data.frame(power = power, ncp.alt = ncp.alt, ncp.null = ncp.null,
                          alpha = alpha,
                          z.crit.1 = zalpha[1,], z.crit.2 = zalpha[2,]),
               row.names = FALSE),
             print(
               data.frame(power = power, ncp.alt = ncp.alt, ncp.null = ncp.null,
                          alpha = alpha, z.crit = zalpha),
               row.names = FALSE))
    }

  } # end of verbose

  return(invisible(power))

} # end of plot.z.test()



############################
# generic F test functions #
############################

# type = 1 for light red, 2 for light blue
.plot.f.dist <- function(ncp = 0, df1, df2, xlim, type = 1) {


  plot.window.dim <- dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  ifelse(type == 1,
         color <- adjustcolor(2, alpha.f = 1),
         color <- adjustcolor(4, alpha.f = 1))

  # non-central f function
  funf <- function(x){
    df(x, df1 = df1, df2 = df2, ncp = ncp)
  }

  # mod of f dist
  ifelse(df1 > 2, mod.x <- (df2 / df1) * (df1 - 2) / (df2 + 2), mod.x <- 0)
  ifelse(df1 > 2, mod.y <- max(df(mod.x, df1 = df1, df2 = df2, ncp = 0),
                               df(mod.x, df1 = df1, df2 = df2, ncp = ncp)),
         mod.y <- 0.60)

  # plot central f distribution
  ylim <- c(0, mod.y + 0.05)
  plot(funf, xlim = xlim, ylim = ylim,
       xaxs = "i", yaxs = "i", bty = "l",
       col = color, lwd = 2, lty = type,
       xlab = "", ylab = "", cex.axis = cex.axis)

}  # end of .plot.f.dist()


# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.f.dist<- function(ncp = 0, df1, df2, xlim, type = 1) {

  color <- switch (type,
                   `1` = adjustcolor(2, alpha.f = 0.3),
                   `2` = adjustcolor(4, alpha.f = 0.3),
                   `3` = adjustcolor(1, alpha.f = 0.3))


  # non-central f function
  funf <- function(x){
    df(x, df1 = df1, df2 = df2, ncp = ncp)
  }

  x <- seq(min(xlim), max(xlim), by = .001)
  y <- funf(x)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if(type == 1 | type == 2) {
    polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    polygon(x = xs, y = ys, col = color, density = 20, angle = 45, border = NA)
  }

  prob <- pf(max(xlim), df1 = df1, df2 = df2, ncp = ncp, lower.tail = TRUE) -
    pf(min(xlim), df1 = df1, df2 = df2, ncp = ncp, lower.tail = TRUE)

  return(invisible(prob))

} # end of .paint.f.dist()


.plot.f.t1t2 <- function(ncp, df1, df2, alpha,
                         plot.main = NULL, plot.sub = NULL) {

  if(length(ncp) > 1) stop("not a valid plotting option", call. = FALSE)

  falpha <- qf(alpha, df1 = df1, df2 = df2, ncp = 0, lower.tail = FALSE)
  yfalpha <- df(falpha, df1 = df1, df2 = df2, ncp = 0)

  # x-axis limits
  ifelse(df1 < 2, prob.lower <- 0.30, prob.lower <- 0.001)
  lower <- qf(prob.lower, df1 = df1, df2 = df2,  ncp = 0)
  upper <- qf(0.999, df1 = df1, df2 = df2,  ncp = ncp)
  xlim <- c(lower, upper)

  plot.window.dim <- dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  .plot.f.dist(ncp = ncp, df1 = df1, df2 = df2, xlim = xlim, type = 2)
  par(new = TRUE)
  .plot.f.dist(ncp = 0, df1 = df1, df2 = df2, xlim = xlim, type = 1)

  # draw vertical lines for critical region
  segments(x0 = falpha, y0 = rep(0, length(falpha)),
           x1 = falpha, y1 = yfalpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  .paint.f.dist(ncp = 0, df1 = df1, df2 = df2, xlim = c(falpha, max(xlim)), type = 1)
  .paint.f.dist(ncp = ncp, df1 = df1, df2 = df2, xlim = c(lower, falpha), type = 2)
  power <- .paint.f.dist(ncp = ncp, df1 = df1, df2 = df2, xlim = c(falpha, max(xlim)), type = 3)
  # end of paint regions

  # axes labels and subtitle
  title(main = plot.main, line = 2, cex.main = cex.title)
  title(sub = plot.sub, line = 3, cex.sub = cex.title)
  title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))
  title(xlab = paste0("F Value (df1 = ", round(df1, digits = 2), ", df2 = ", round(df2, digits = 2), ")"),
        line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))

  alpha <- round(alpha, 2)
  beta <- round(1 - power, 2)
  power <- round(power, 2)

  legend("topright", cex = cex.legend,
         c(as.expression(bquote(Power == .(power))),
           as.expression(bquote(alpha == .(alpha))),
           as.expression(bquote(beta == .(beta)))),
         fill = c(adjustcolor(1, alpha.f = 0.3),
                  adjustcolor(2, alpha.f = 0.3),
                  adjustcolor(4, alpha.f = 0.3)),
         border = c(adjustcolor(1, alpha.f = 0.15),
                    adjustcolor(2, alpha.f = 0.15),
                    adjustcolor(4, alpha.f = 0.15)),
         bg = adjustcolor(1, alpha.f = 0.08),
         box.col = adjustcolor(1, alpha.f = 0),
         density = c(30, NA, NA),
         angle = c(45, NA, NA))

} # end of .plot.f.t1t2()


# power for the generic f test with (optional) type I and type II error plots
power.f.test <- function(ncp, df1, df2, alpha = 0.05,
                         plot = TRUE, plot.main = NULL, plot.sub = NULL,
                         verbose = TRUE) {

  if(sum(c(length(ncp), length(df1), length(df2), length(alpha)) > 1) > 1)
    stop("only one of the 'ncp', 'df1', 'df2', and 'alpha' arguments can take multiple values", call. = FALSE)

  if(any(c(length(ncp),  length(df1), length(df2), length(alpha)) > 1) & isTRUE(plot))
    stop("plotting is not available for multiple values", call. = FALSE)

  if(any(df1 < 0)) stop("'df1' sould be positive", call. = FALSE)
  if(any(df2 < 0)) stop("'df2' sould be positive", call. = FALSE)
  if(any(alpha < 0) | any(alpha > .999)) stop("'alpha' sould be between 0 and 0.999", call. = FALSE)
  if(any(ncp < 0)) stop("'ncp' sould be positive", call. = FALSE)

  falpha <- qf(alpha, df1 = df1, df2 = df2, ncp = 0, lower.tail = FALSE)

  power <- pf(qf(alpha, df1 = df1, df2 = df2, lower.tail = FALSE),
              df1 = df1, df2 = df2, ncp = ncp, lower.tail = FALSE)

  if(plot) {

    .plot.f.t1t2(ncp = ncp, df1 = df1, df2 = df2, alpha = alpha,
                 plot.main = plot.main, plot.sub = plot.sub)

  }

  if(verbose) {

    if(any(c(is.matrix(ncp), is.matrix(df1), is.matrix(df2), is.matrix(alpha))))
      stop("'verbose' argument is not available for matrix type input", call. = FALSE)

    ncp.alt <- round(ncp, 3)
    ncp.null <- 0
    print(data.frame(power = power, ncp.alt = ncp.alt, ncp.null = ncp.null,
                     alpha = alpha, df1 = df1, df2 = df2, f.crit = falpha), row.names = FALSE)

  } # end of verbose

  return(invisible(power))

} # end of plot.f.test()

#####################################
# generic Chi-square test functions #
#####################################

# type = 1 for light red, 2 for light blue
.plot.chisq.dist <- function(ncp = 0, df, xlim, type = 1) {


  plot.window.dim <- dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  ifelse(type == 1,
         color <- adjustcolor(2, alpha.f = 1),
         color <- adjustcolor(4, alpha.f = 1))

  # non-central f function
  funchisq <- function(x){
    dchisq(x, df = df, ncp = ncp)
  }

  # plot central t distribution
  ylim <- c(0, 0.25)
  plot(funchisq, xlim = xlim, ylim = ylim,
       xaxs = "i", yaxs = "i", bty = "l",
       col = color, lwd = 2, lty = type,
       xlab = "", ylab = "", cex.axis = cex.axis)

} # end of .plot.chisq.dist()


# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.chisq.dist<- function(ncp = 0, df, xlim, type = 1) {

  color <- switch (type,
                   `1` = adjustcolor(2, alpha.f = 0.3),
                   `2` = adjustcolor(4, alpha.f = 0.3),
                   `3` = adjustcolor(1, alpha.f = 0.3))


  # non-central f function
  funchisq <- function(x){
    dchisq(x, df = df, ncp = ncp)
  }

  x <- seq(min(xlim), max(xlim), by = .001)
  y <- funchisq(x)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if(type == 1 | type == 2) {
    polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    polygon(x = xs, y = ys, col = color, density = 20, angle = 45, border = NA)
  }

  prob <- pchisq(max(xlim), df = df, ncp = ncp, lower.tail = TRUE) -
    pchisq(min(xlim), df = df, ncp = ncp, lower.tail = TRUE)

  return(invisible(prob))

} # end of .paint.chisq.dist()



.plot.chisq.t1t2 <- function(ncp, df, alpha,
                             plot.main = NULL, plot.sub = NULL) {

  if(length(ncp) > 1) stop("not a valid plotting option", call. = FALSE)

  chisq.alpha <- qchisq(alpha, df = df, ncp = 0, lower.tail = FALSE)
  y.chisq.alpha <- dchisq(chisq.alpha, df = df, ncp = 0)

  # x-axis limits
  ifelse(df < 2, prob.lower <- 0.30, prob.lower <- 0.001)
  lower <- qchisq(prob.lower, df = df, ncp = 0)
  upper <- qchisq(0.999, df = df, ncp = ncp)
  xlim <- c(lower, upper)

  plot.window.dim <- dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  .plot.chisq.dist(ncp = ncp, df = df, xlim = xlim, type = 2)
  par(new = TRUE)
  .plot.chisq.dist(ncp = 0, df = df, xlim = xlim, type = 1)

  # draw vertical lines for critical region
  segments(x0 = chisq.alpha, y0 = rep(0, length(chisq.alpha)),
           x1 = chisq.alpha, y1 = y.chisq.alpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  .paint.chisq.dist(ncp = 0, df = df,  xlim = c(chisq.alpha, max(xlim)), type = 1)
  .paint.chisq.dist(ncp = ncp, df = df, xlim = c(lower, chisq.alpha), type = 2)
  power <- .paint.chisq.dist(ncp = ncp, df = df, xlim = c(chisq.alpha, max(xlim)), type = 3)
  # end of paint regions

  # axes labels and subtitle
  title(main = plot.main, line = 2, cex.main = cex.title)
  title(sub = plot.sub, line = 3, cex.sub = cex.title)
  title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))
  title(xlab = as.expression(bquote(chi[2] ~ "Value (" * df==.(df) * ")")),
        line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))

  alpha <- round(alpha, 2)
  beta <- round(1 - power, 2)
  power <- round(power, 2)

  legend("topright", cex = cex.legend,
         c(as.expression(bquote(Power == .(power))),
           as.expression(bquote(alpha == .(alpha))),
           as.expression(bquote(beta == .(beta)))),
         fill = c(adjustcolor(1, alpha.f = 0.3),
                  adjustcolor(2, alpha.f = 0.3),
                  adjustcolor(4, alpha.f = 0.3)),
         border = c(adjustcolor(1, alpha.f = 0.15),
                    adjustcolor(2, alpha.f = 0.15),
                    adjustcolor(4, alpha.f = 0.15)),
         bg = adjustcolor(1, alpha.f = 0.08),
         box.col = adjustcolor(1, alpha.f = 0),
         density = c(30, NA, NA),
         angle = c(45, NA, NA))

} # end of .plot.f.t1t2()


# power for the generic Chi-square test with (optional) type I and type II error plots
power.chisq.test <- function(ncp, df, alpha = 0.05,
                             plot = TRUE, plot.main = NULL, plot.sub = NULL,
                             verbose = TRUE) {

  if(sum(c(length(ncp), length(df), length(alpha)) > 1) > 1)
    stop("only one of the 'ncp', 'df', and 'alpha' arguments can take multiple values", call. = FALSE)

  if(any(c(length(ncp), length(df), length(alpha)) > 1) & isTRUE(plot))
    stop("plotting is not available for multiple values", call. = FALSE)

  if(any(df < 0)) stop("'df' sould be positive", call. = FALSE)
  if(any(alpha < 0) | any(alpha > .999)) stop("'alpha' sould be between 0 and 0.999", call. = FALSE)
  if(any(ncp < 0)) stop("'ncp' sould be positive", call. = FALSE)

  chisq.alpha <- qchisq(alpha, df = df, ncp = 0, lower.tail = FALSE)

  power <- pchisq(chisq.alpha, df = df, ncp = ncp, lower.tail = FALSE)

  if(plot) {

    .plot.chisq.t1t2(ncp = ncp, df = df, alpha = alpha,
                     plot.main = plot.main, plot.sub = plot.sub)

  }

  if(verbose) {

    if(any(c(is.matrix(ncp), is.matrix(df), is.matrix(alpha))))
      stop("'verbose' argument is not available for matrix type input", call. = FALSE)

    ncp.alt <- round(ncp, 3)
    ncp.null <- 0
    print(data.frame(power = power, ncp.alt = ncp.alt, ncp.null = ncp.null,
                     alpha = alpha, df = df, chisq.crit = chisq.alpha), row.names = FALSE)

  } # end of verbose

  return(invisible(power))

} # end of plot.chisq.test()

