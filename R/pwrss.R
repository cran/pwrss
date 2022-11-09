pwrss.z.prop <- function (p, p0 = 0, margin = 0, arcsin.trans = TRUE,
                          alpha = 0.05,
                          alternative = c("not equal", "greater", "less",
                                          "equivalent", "non-inferior", "superior"),
                          n = NULL, power = NULL, verbose = TRUE)
{

  alternative <- match.arg(alternative)

  if (is.null(n) & is.null(power)) stop("`n` and `power` cannot be `NULL` at the same time", call. = FALSE)
  if (!is.null(n) & !is.null(power)) stop("one of the `n` or `power` should be `NULL`", call. = FALSE)
  if ((alternative == "not equal" | alternative == "not equal" | alternative == "not equal") & margin != 0) warning("`margin` argument is ignored", call. = FALSE)
  if (alternative == "superior" & margin < 0) warning("expecting `margin > 0`", call. = FALSE)
  if (alternative == "non-inferior" & margin > 0) warning("expecting `margin < 0`", call. = FALSE)
  if (alternative == "greater" & (p < p0)) stop("alternative = 'greater' but p < p0", call. = FALSE)
  if (alternative == "less" & (p > p0)) stop("alternative = 'less' but p > p0", call. = FALSE)

  if(arcsin.trans) {
    var.num <- 1
    ifelse(margin < 0, sign <- -1, sign <- 1)
    margin <- abs(margin)
    h <- switch(alternative,
                `not equal` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)),
                `greater` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)),
                `less` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)),
                `non-inferior` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)) - sign*2*asin(sqrt(margin)),
                `superior` = 2*asin(sqrt(p)) - 2*asin(sqrt(p0)) - sign*2*asin(sqrt(margin)),
                `equivalent` = abs(2*asin(sqrt(p)) - 2*asin(sqrt(p0))) - sign*2*asin(sqrt(margin)))
    if(verbose) cat(" Approach: Arcsine transformation \n")
  } else {
    var.num <- p * (1 - p)
    h <- switch(alternative,
                `not equal` = p - p0,
                `greater` = p - p0,
                `less` = p - p0,
                `non-inferior` = p - p0 - margin,
                `superior` = p - p0 - margin,
                `equivalent` = abs(p - p0) - margin)
    if(verbose) cat(" Approach: Normal approximation \n")
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
      n <- M^2 * var.num / h^2
      lambda <- h / sqrt(var.num / n)

    }
    if (is.null(power)) {
      lambda <- h / sqrt(var.num / n)
      power <- 2*(1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }

  } else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }

  if(verbose) {
    cat(" One proportion compared to a constant (one sample z test) \n",
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

  invisible(structure(list(parms = list(p = p, p0 = p0, arcsin.trans = arcsin.trans,
                                        alpha = alpha, margin = margin, alternative = alternative, verbose = verbose),
                           test = "z",
                           ncp = lambda,
                           power = power,
                           n = n),
                      class = c("pwrss", "z")))
}

#################################################################################

pwrss.z.2props <- function (p1, p2, margin = 0, arcsin.trans = TRUE,
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
  if (alternative == "superior" & margin < 0) warning("expecting `margin > 0`", call. = FALSE)
  if (alternative == "non-inferior" & margin > 0) warning("expecting `margin < 0`", call. = FALSE)
  if (alternative == "greater" & (p1 < p2)) stop("alternative = 'greater' but p1 < p2", call. = FALSE)
  if (alternative == "less" & (p1 > p2)) stop("alternative = 'less' but p1 > p2", call. = FALSE)

  if(arcsin.trans) {
    ifelse(margin < 0, sign <- -1, sign <- 1)
    margin <- abs(margin)
    h <- switch(alternative,
                `not equal` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)),
                `greater` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)),
                `less` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)),
                `non-inferior` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)) - sign*2*asin(sqrt(margin)),
                `superior` = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)) - sign*2*asin(sqrt(margin)),
                `equivalent` = abs(2*asin(sqrt(p1)) - 2*asin(sqrt(p2))) - sign*2*asin(sqrt(margin)))
    if(verbose) cat(" Approach: Arcsine transformation \n")
  } else {
    h <- switch(alternative,
                `not equal` = p1 - p2,
                `greater` = p1 - p2,
                `less` = p1 - p2,
                `non-inferior` = p1 - p2 - margin,
                `superior` = p1 - p2 - margin,
                `equivalent` = abs(p1 - p2) - margin)
    if(verbose) cat(" Approach: Normal approximation \n")
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

      power <- 1 - pnorm(qnorm(alpha / 2, lower.tail = FALSE), lambda) +
        pnorm(-qnorm(alpha / 2, lower.tail = FALSE), lambda)
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

      if(alternative == "less") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
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

      if(alternative == "non-inferior") lambda <- abs(lambda)
      power <- 1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)
    }

  }
  else if (alternative == "equivalent") {

    HA_H0 <- abs(p1 - p2) - margin
    if (is.null(n2)) {
      beta <- 1 - power
      M <- qnorm(alpha, lower.tail = FALSE) + qnorm(beta / 2, lower.tail = FALSE)
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

      power <- 2*(1 - pnorm(qnorm(alpha, lower.tail = FALSE), lambda)) - 1
      if(power < 0) stop("design is not feasible", call. = FALSE)
    }

  }
  else {
    stop("`alternative` should be in c('not equal', 'greater', 'less', 'equivalent', 'non-inferior', 'superior')",
         call. = FALSE)
  }

  if(verbose) {
    cat(" Difference between two proportions (independent samples z test) \n",
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
                      class = c("pwrss", "z")))
}

#################################################################################

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
    if (alternative == "superior" & margin < 0)
      warning("expecting 'margin > 0' when mu - mu0 > 0", call. = FALSE)
    if (alternative == "non-inferior" & margin > 0)
      warning("expecting 'margin < 0' when mu - mu0 > 0", call. = FALSE)

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
    cat(" One mean compared to a constant (one sample z test) \n",
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
                      class = c("pwrss", "z")))
}

#################################################################################

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
    cat(" Difference between two means (independent samples z test) \n",
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
                      class = c("pwrss", "z")))
}

#################################################################################

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
      (sd1^2 / (n1^2 * (n1 - 1)) + sd2^2 / (n2^2 * (n2 - 1)))
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
      warning("`margin` argument is ignored")
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
    if (alternative == "superior" & margin < 0)
      warning("expecting 'margin > 0' when mu1 - mu2 > 0", call. = FALSE)
    if (alternative == "non-inferior" & margin > 0)
      warning("expecting 'margin < 0' when mu1 - mu2 > 0", call. = FALSE)

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
    ncp <- (HA_H0) / sqrt((sd1^2 + sd2^2 - 2 * sd1 * sd2 * paired.r) / n2)
  } else {
    ncp <- (HA_H0) / sqrt(sd1^2 / n1 + sd2^2 / n2)
  }
  ifelse(paired, n <- n2, n <- c(n1 = n1, n2 = n2))

  hypothesis <- alternative

  if(verbose) {
    cat(ifelse(paired,
               " Difference between two means (paired samples t test) \n",
               " Difference between two means (independent samples t test) \n"),
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
                      class = c("pwrss", "t")))
}

#################################################################################

pwrss.z.corr <- function (r = 0.50, r0 = 0, alpha = 0.05,
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
    cat(" One correlation compared to a constant (one sample z test) \n",
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
                      class = c("pwrss", "z")))

}

#################################################################################

pwrss.z.2corrs <- function (r1 = 0.50, r2 = 0.30,
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
    cat(" Difference between two correlations (independent samples z test) \n",
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
                      class = c("pwrss", "z")))

}

#################################################################################

# when k = 1 and predictor is binary
# d <- 0.20
# r2 <- d^2 / (d^2 + 4)
# f2 <- r2 /(1 - r2)
# results will be same as pwrss.t.2means(mu1 = d,...)
# specify k = m (the default) to get the test of r2 difference from zero
# specify k > m to get the test of r2 change from a change of zero
pwrss.f.reg <- function (r2 = 0.10, f2 = r2 /(1 - r2),
                         k = 4, m = k, alpha = 0.05,
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
               " R-squared compared to 0 in linear regression (F test) \n",
               " R-squared change in hierarchical linear regression (F test) \n"),
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

#################################################################################

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
    cat(" ", switch(n.way,
                    `1` = "One",
                    `2` = "Two",
                    `3` = "Three"),
        ifelse(n.covariates > 0,
               "-way Analysis of Covariance (ANCOVA) \n ",
               "-way Analysis of Variance (ANOVA) \n "),
        " H0: 'eta2' or 'f2' = 0 \n  HA: 'eta2' or 'f2' > 0 \n ------------------------------------\n",
        switch(n.way,
               `1` = c(" Factor A: ", n.levels, " levels \n"),
               `2` = c(" Factor A: ", n.levels[1], " levels \n", " Factor B: ", n.levels[2], " levels \n"),
               `3` = c(" Factor A: ", n.levels[1], " levels \n", " Factor B: ", n.levels[2], " levels \n", " Factor C: ", n.levels[3], " levels \n")),
        " ------------------------------------\n",
        sep = "")
  }

  if( n.way == 1) {

    df1 <- n.levels[1] - 1
    if(requested == "n") {
      n <- ss(df1 = df1, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha, power = power)
      df2 <- n - n.groups - n.covariates
      ncp <- n*f2

      if(verbose) {
        cat(" Given eta2 =", round(eta2, 3), "or f2 =", round(f2, 3), "\n",
            " Total n =", ceiling(n), "\n ------------------------------------\n",
            "Numerator degrees of freedom = ", round(df1, 2), "\n",
            "Denominator degrees of freedom =", round(df2, 2),"\n",
            "Non-centrality parameter =", round(ncp, 2),"\n")
      }

    } else if (requested == "power") {
      power <- pwr(df1 = df1, n = n, n.groups = n.groups, n.covariates = n.covariates, f2 = f2, alpha = alpha)
      df2 <- n - n.groups - n.covariates
      ncp <- n*f2

      if(verbose) {
        cat(" Given eta2 =", round(eta2, 3), "or f2 =", round(f2, 3), "\n",
            " Statistical power =", round(power, 3), "\n ------------------------------------\n",
            "Numerator degrees of freedom = ", round(df1, 2), "\n",
            "Denominator degrees of freedom =", round(df2, 2),"\n",
            "Non-centrality parameter =", round(ncp, 2),"\n")
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
        cat(" Given eta2 =", round(eta2, 3), "or f2 =", round(f2, 3), "\n",
            " Total n =", ceiling(ss.f1), "(for A) \n",
            " Total n =", ceiling(ss.f2), "(for B) \n",
            " Total n =", ceiling(ss.f1f2), "(for A x B) \n ------------------------------------\n",
            "Numerator degrees of freedom = ", round(df1, 2), "\n",
            "Denominator degrees of freedom =", round(df2, 2),"\n",
            "Non-centrality parameter =", round(ncp, 2),"\n")
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
        cat(" Given eta2 =", round(eta2, 3), "or f2 =", round(f2, 3), "\n",
            " Statistical power =", round(pwr.f1, 3), "(for A) \n",
            " Statistical power =", round(pwr.f2, 3), "(for B) \n",
            " Statistical power =", round(pwr.f1f2, 3), "(for A x B) \n ------------------------------------\n",
            "Numerator degrees of freedom = ", round(df1, 2), "\n",
            "Denominator degrees of freedom =", round(df2, 2),"\n",
            "Non-centrality parameter =", round(ncp, 2),"\n")
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
        cat(" Given eta2 =", round(eta2, 3), "or f2 =", round(f2, 3), "\n",
            " Total n =", ceiling(ss.f1), "(for A) \n",
            " Total n =", ceiling(ss.f2), "(for B) \n",
            " Total n =", ceiling(ss.f3), "(for C) \n",
            " Total n =", ceiling(ss.f1f2), "(for A x B) \n",
            " Total n =", ceiling(ss.f1f3), "(for A x C) \n",
            " Total n =", ceiling(ss.f2f3), "(for B x C) \n",
            " Total n =", ceiling(ss.f1f2f3), "(for A x B x C) \n ------------------------------------\n",
            "Numerator degrees of freedom = ", round(df1, 2), "\n",
            "Denominator degrees of freedom =", round(df2, 2),"\n",
            "Non-centrality parameter =", round(ncp, 2),"\n")
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
        cat(" Given eta2 =", round(eta2, 3), "or f2 =", round(f2, 3), "\n",
            " Statistical power =", round(pwr.f1, 3), "(for A) \n",
            " Statistical power =", round(pwr.f2, 3), "(for B) \n",
            " Statistical power =", round(pwr.f3, 3), "(for C) \n",
            " Statistical power =", round(pwr.f1f2, 3), "(for A x B) \n",
            " Statistical power =", round(pwr.f1f3, 3), "(for A x C) \n",
            " Statistical power =", round(pwr.f2f3, 3), "(for B x C) \n",
            " Statistical power =", round(pwr.f1f2f3, 3), "(for A x B x C) \n ------------------------------------\n",
            "Numerator degrees of freedom = ", round(df1, 2), "\n",
            "Denominator degrees of freedom =", round(df2, 2),"\n",
            "Non-centrality parameter =", round(ncp, 2),"\n")
      }

    } else {

      stop("Invalid 'n' or 'power'", call. = FALSE)

    }

  } else {

    stop("More than three-way ANOVA or ANCOVA is not allowed")

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
}

#################################################################################

# nonsphericity.cf is epsilon in SPSS output
# 1 / (n.measurements - 1) < nonsphericity.cf < 1 when spechiricty assumptions does not hold (or when compound symmetry is violated)
# Greenhouse and Geisser (1959), the other is by Huynh and Feldt (1976).
pwrss.f.rmanova <- function (eta2 = 0.10, f2 = eta2/(1 - eta2),
                             repmeasures.r = 0.50, n.levels = 2, n.measurements = 2,
                             epsilon = 1, alpha = 0.05,
                             type = c("between", "within", "interaction"),
                             n = NULL, power = NULL, verbose = TRUE)
{
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

  ss <- function(f2, n.levels, n.measurements, epsilon, alpha, power, type) {
    n.min <- n.levels + 1
    n <- try(silent = TRUE,
             uniroot(function(n) {
               if(type == 0) {
                 df1 <- n.levels - 1
                 df2 <- n - n.levels
               } else if(type == 1) {
                 df1 <- (n.measurements - 1) * epsilon
                 df2 <- (n - n.levels) * (n.measurements - 1) * epsilon
               } else if(type == 2) {
                 df1 <- (n.levels - 1) * (n.measurements - 1) * epsilon
                 df2 <- (n - n.levels) * (n.measurements - 1) * epsilon
               } else {
                 stop("Unknown type of effect", call. = FALSE)
               }
               u <- df1
               v <- df2
               if(u < 1 | v < 1) stop("design is not feasible", call. = FALSE)
               lambda <- f2 * n * epsilon
               power - pf(qf(alpha, df1 = u, df2 = v, lower.tail = FALSE),
                          df1 = u, df2 = v, ncp = lambda, lower.tail = FALSE)
             }, interval = c(n.min, 1e10))$root
    ) # try
    if(inherits(n, "try-error") | n == 1e+10) stop("design is not feasible", call. = FALSE)
    n
  }

  pwr <- function(f2, n, n.levels, n.measurements, epsilon, alpha, type) {
    if(type == 0) {
      df1 <- n.levels - 1
      df2 <- n - n.levels
    } else if(type == 1) {
      df1 <- (n.measurements - 1) * epsilon
      df2 <- (n - n.levels) * df1
    } else if(type == 2) {
      df1 <- (n.levels - 1) * (n.measurements - 1) * epsilon
      df2 <- (n - n.levels) * (n.measurements - 1) * epsilon
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
    f2 <- f2  * (n.measurements / (1 + (n.measurements - 1) * repmeasures.r))
  } else {
    f2 <- f2  * (n.measurements / (1 - repmeasures.r))
  }

  if (is.null(n)) {
    n <- ss(f2 = f2, n.levels = n.levels, n.measurements = n.measurements,
            epsilon = epsilon, alpha = alpha, power = power, type = type)
  }

  if (is.null(power)) {
    power <- pwr(f2 = f2, n = n, n.levels = n.levels, n.measurements = n.measurements,
                 epsilon = epsilon, alpha = alpha, type = type)
  }

  if(type == 0) {
    df1 <- n.levels - 1
    df2 <- n - n.levels
  } else if(type == 1) {
    df1 <- (n.measurements - 1) * epsilon
    df2 <- (n - n.levels) * df1
  } else if(type == 2) {
    df1 <- (n.levels - 1) * (n.measurements - 1) * epsilon
    df2 <- (n - n.levels) * (n.measurements - 1) * epsilon
  }

  ncp <- f2 * n * epsilon

  if(verbose) {
    cat(" One-way repeated measures analysis of variance (F test) \n",
        "H0: eta2 = 0 (or f2 = 0) \n HA: eta2 > 0 (or f2 > 0) \n",
        "------------------------------ \n",
        "Number of levels (groups) =", n.levels, "\n",
        "Number of measurement time points =",n.measurements, "\n",
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
                                        n.measurements = n.measurements, repmeasures.r = repmeasures.r,
                                        epsilon = epsilon, alpha = alpha, verbose = verbose),
                           test = "F",
                           df1 = df1,
                           df2 = df2,
                           ncp = f2 * n * epsilon,
                           power = power,
                           n = n),
                      class = c("pwrss", "f", "rmanova")))

}

#################################################################################

# power for the generic t test wit type I and type II erro plots
power.t.test <- function(ncp, df, alpha, alternative,
                         plot = TRUE, plot.main = NULL, plot.sub = NULL){

  alternative <- switch(alternative,
                        `not equal` = 0,
                        `greater` = 1,
                        `less` = 2,
                        `non-inferior` = 3,
                        `superior` = 4,
                        `equivalent` = 5)

  if(df < 0) stop("'df' sould be positive", call. = FALSE)
  if(alpha < 0 | alpha > .999) stop("'alpha' sould be between 0 and 0.999", call. = FALSE)
  if(ncp < 0) stop("'ncp' sould be positive - provide the absolute value", call. = FALSE)

  if(alternative == 0) {
    talpha <- qnorm(alpha / 2, lower.tail = FALSE)
  } else {
    talpha <- qnorm(alpha, lower.tail = FALSE)
  }

  if(alternative == 0) {
    power <- 1 - pt(qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = ncp) +
      pt(-qt(alpha / 2, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = ncp)
  } else if(alternative %in% 1:4) {
    power <- 1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = ncp)
  } else if(alternative == 5) {
    power <- 2 * (1 - pt(qt(alpha, df = df, ncp = 0, lower.tail = FALSE), df = df, ncp = ncp)) - 1
  } else {
    stop("not a valid test type", call. = FALSE)
  }


  if(plot) {

    alpha <- round(alpha, 3)
    talpha <- round(talpha, 2)
    ncp <- round(ncp, 2)
    beta <- round(1 - power, 3)


    # central t function
    funt0 <- function(x){
      dt(x, df = df, ncp = 0)
    }
    # non-central t function
    funt1 <- function(x){
      dt(x, df = df, ncp = ncp)
    }

    # plot central t distribution
    plot(funt0, xlim = c(-4, ncp + 4), ylim = c(0, 0.5),
         yaxs = "i", xaxs = "i", bty = "l",
         col = adjustcolor(2, alpha.f = 1), lwd = 2,
         # sub = paste("Type I Error Rate = ", round(alpha ,digits = 3), ",",
         #             "Type II Error Rate = ", round(beta, digits = 3), "\n",
         #             "Non-centrality Parameter (NCP) =  ", round(ncp, digits = 3)),
         xlab = "", ylab = "")

    par(new = TRUE)

    # plot non-central t distribution
    plot(funt1, xlim = c(-4, ncp + 4), ylim = c(0, 0.5),
         col = adjustcolor(4, alpha.f = 1), lwd = 2, lty = 2,
         yaxs = "i", xaxs = "i", yaxt = "n", xaxt = "n",
         bty = "l", xlab = "", ylab = "")
    legend("topright",
           c(as.expression(bquote(alpha == .(alpha))),
             as.expression(bquote(beta == .(beta))),
             ifelse(alternative == 0,
                    as.expression(bquote(t[alpha/2] == .(talpha))),
                    as.expression(bquote(t[alpha] == .(talpha)))
             ),
             as.expression(bquote(lambda == .(ncp)))),
           pch = c(19, 19, NA, NA), lty = c(NA, NA, 2, 2), cex = 1,
           col = c(adjustcolor(2, alpha.f = 0.3), adjustcolor(4, alpha.f = 0.3), 2, 4))

    # axes labels and subtitle
    title(main = plot.main, line = 2)
    title(sub = plot.sub, line = 3)
    title(ylab = "Probability Density", line = 2)
    title(xlab = paste0(expression(t), "(df = ", round(df, digits = 2), ")"), line = 2)

    # draw vertical lines
    # abline(v = 0, lty = 2, col = 4) # mean of central t in dashed blue line
    # abline(v = ncp, lty = 2, col = 4) # mean of non-central t in dashed blue line
    abline(v = talpha, lty = 2, col = 2) # t-critical in dashed red line

    # shaded area in light red for alpha
    xalpha0 <- seq(from = talpha, to = ncp + 4, by = .001)
    yalpha0 <- rep(NA, length(xalpha0))
    for(i in 1:length(xalpha0)){
      yalpha0[i] <- funt0(xalpha0[i])
    }

    xalpha <- c(xalpha0, rev(xalpha0))
    yalpha <- c(yalpha0, rep(0, length(yalpha0)))
    polygon(x = xalpha, y = yalpha, col = adjustcolor(2, alpha.f = 0.3), border = NA)

    # shaded area in light blue for beta
    xbeta0 <- seq(from = -4, to = talpha, by = .001)
    ybeta0 <- rep(NA, length(xbeta0))
    for(i in 1:length(xbeta0)){
      ybeta0[i] <- funt1(xbeta0[i])
    }
    xbeta <- c(xbeta0, rev(xbeta0))
    ybeta <- c(ybeta0, rep(0, length(ybeta0)))
    polygon(x = xbeta, y = ybeta, col = adjustcolor(4, alpha.f = 0.3), border = NA)

  }

  return(power)

}

#################################################################################

# power for the generic z test wit type I and type II erro plots
power.z.test <- function(ncp, alpha, alternative,
                         plot = TRUE, plot.main = NULL, plot.sub = NULL){

  alternative <- switch(alternative,
                        `not equal` = 0,
                        `greater` = 1,
                        `less` = 2,
                        `non-inferior` = 3,
                        `superior` = 4,
                        `equivalent` = 5)

  if(alpha < 0 | alpha > .999) stop("'alpha' sould be between 0 and 0.999", call. = FALSE)
  if(ncp < 0) stop("'ncp' sould be positive - provide the absolute value", call. = FALSE)

  if(alternative == 0) {
    zalpha <- qnorm(alpha / 2, lower.tail = FALSE)
  } else {
    zalpha <- qnorm(alpha, lower.tail = FALSE)
  }

  if(alternative == 0) {
    power <- 1 - pnorm(qnorm(alpha / 2, mean = 0, lower.tail = FALSE), mean = ncp) +
      pnorm(-qnorm(alpha / 2, mean = 0, lower.tail = FALSE), mean = ncp)
  } else if(alternative %in% 1:4) {
    power <- 1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = ncp)
  } else if(alternative == 5) {
    power <- 2 * (1 - pnorm(qnorm(alpha, mean = 0, lower.tail = FALSE), mean = ncp)) - 1
  } else {
    stop("not a valid test type", call. = FALSE)
  }

  if(plot) {

    alpha <- round(alpha, 3)
    zalpha <- round(zalpha, 2)
    ncp <- round(ncp, 2)
    beta <- round(1 - power, 3)

    # central t function
    funz0 <- function(x){
      dnorm(x, mean = 0)
    }
    # non-central t function
    funz1 <- function(x){
      dnorm(x, mean = ncp)
    }

    # plot central t distribution
    plot(funz0, xlim = c(-4, ncp + 4), ylim = c(0, 0.5),
         yaxs = "i", xaxs = "i", bty = "l",
         col = adjustcolor(2, alpha.f = 1), lwd = 2,
         # sub = paste("Type I Error Rate = ", round(alpha ,digits = 3), ",",
         #             "Type II Error Rate = ", round(beta, digits = 3), "\n",
         #             "Non-centrality Parameter (NCP) =  ", round(ncp, digits = 3)),
         xlab = "", ylab = "")

    par(new = TRUE)

    # plot non-central t distribution
    plot(funz1, xlim = c(-4, ncp + 4), ylim = c(0, 0.5),
         col = adjustcolor(4, alpha.f = 1), lwd = 2, lty = 2,
         yaxs = "i", xaxs = "i", yaxt = "n", xaxt = "n",
         bty = "l", xlab = "", ylab = "")
    legend("topright",
           c(as.expression(bquote(alpha == .(alpha))),
             as.expression(bquote(beta == .(beta))),
             ifelse(alternative == 0,
                    as.expression(bquote(z[alpha/2] == .(zalpha))),
                    as.expression(bquote(z[alpha] == .(zalpha)))
             ),
             as.expression(bquote(lambda == .(ncp)))),
           pch = c(19, 19, NA, NA), lty = c(NA, NA, 2, 2), cex = 1,
           col = c(adjustcolor(2, alpha.f = 0.3), adjustcolor(4, alpha.f = 0.3), 2, 4))

    # axes labels and subtitle
    title(main = plot.main, line = 2)
    title(sub = plot.sub, line = 3)
    title(ylab = "Probability Density", line = 2)
    title(xlab = paste0(expression(z)), line = 2)

    # draw vertical lines
    # abline(v = 0, lty = 2, col = 4) # mean of central t in dashed blue line
    # abline(v = ncp, lty = 2, col = 4) # mean of non-central t in dashed blue line
    abline(v = zalpha, lty = 2, col = 2) # t-critical in dashed red line

    # shaded area in light red for alpha
    xalpha0 <- seq(from = zalpha, to = ncp + 4, by = .001)
    yalpha0 <- rep(NA, length(xalpha0))
    for(i in 1:length(xalpha0)){
      yalpha0[i] <- funz0(xalpha0[i])
    }

    xalpha <- c(xalpha0, rev(xalpha0))
    yalpha <- c(yalpha0, rep(0, length(yalpha0)))
    polygon(x = xalpha, y = yalpha, col = adjustcolor(2, alpha.f = 0.3), border = NA)

    # shaded area in light blue for beta
    xbeta0 <- seq(from = -4, to = zalpha, by = .001)
    ybeta0 <- rep(NA, length(xbeta0))
    for(i in 1:length(xbeta0)){
      ybeta0[i] <- funz1(xbeta0[i])
    }
    xbeta <- c(xbeta0, rev(xbeta0))
    ybeta <- c(ybeta0, rep(0, length(ybeta0)))
    polygon(x = xbeta, y = ybeta, col = adjustcolor(4, alpha.f = 0.3), border = NA)

  }

  return(power)
}

#################################################################################

# power for the generic f test wit type I and type II erro plots
power.f.test <- function(ncp, df1, df2, alpha,
                         plot = TRUE, plot.main = NULL, plot.sub = NULL) {

  if(df1 < 0) stop("'df1' sould be positive", call. = FALSE)
  if(df2 < 0) stop("'df2' sould be positive", call. = FALSE)
  if(alpha < 0 | alpha > .999) stop("'alpha' sould be between 0 and 0.999", call. = FALSE)
  if(ncp < 0) stop("'ncp' sould be positive - provide the absolute value", call. = FALSE)

  falpha <- qf(alpha, df1 = df1, df2 = df2, ncp = 0, lower.tail = FALSE)

  power <- pf(qf(alpha, df1 = df1, df2 = df2, lower.tail = FALSE),
              df1 = df1, df2 = df2, ncp = ncp, lower.tail = FALSE)

  if(plot) {

    alpha <- round(alpha, 3)
    falpha <- round(falpha, 2)
    ncp <- round(ncp, 2)
    beta <- round(1 - power, 3)

    # central t function
    funf0 <- function(x){
      df(x, df1 = df1, df2 = df2, ncp = 0)
    }
    # non-central t function
    funf1 <- function(x){
      df(x, df1 = df1, df2 = df2, ncp = ncp)
    }

    # plot central t distribution
    plot(funf0, xlim = c(0, ncp + 4), ylim = c(0, 0.90),
         col = adjustcolor(2, alpha.f = 1), lwd = 2,
         yaxs = "i", xaxs = "i", bty = "l",
         # sub = paste("Type I Error Rate = ", round(alpha ,digits = 3), ",",
         #             "Type II Error Rate = ", round(beta, digits = 3), "\n",
         #             "Non-centrality Parameter (NCP) =  ", round(ncp, digits = 3)),
         xlab = "", ylab = "")

    par(new = TRUE)

    # plot non-central t distribution
    plot(funf1, xlim = c(0, ncp + 4), ylim = c(0, 0.90),
         col = adjustcolor(4, alpha.f = 1), lwd = 2, lty = 2,
         yaxs = "i", xaxs = "i", yaxt = "n", xaxt = "n",
         bty = "l", xlab = "", ylab = "")
    legend("topright",
           c(as.expression(bquote(alpha == .(alpha))),
             as.expression(bquote(beta == .(beta))),
             as.expression(bquote(F[alpha] == .(falpha))),
             as.expression(bquote(lambda == .(ncp)))),
           pch = c(19, 19, NA, NA), lty = c(NA, NA, 2, 2), cex = 1,
           col = c(adjustcolor(2, alpha.f = 0.3), adjustcolor(4, alpha.f = 0.3), 2, 4))

    # axes labels and subtitle
    title(main = plot.main, line = 2)
    title(sub = plot.sub, line = 3)
    title(ylab = "Probability Density", line = 2)
    title(xlab = paste0(expression(F),
                        "(df1 = ", round(df1, digits = 2),
                        ", df2 = ", round(df2, digits = 2), ")"), line = 2)

    # draw vertical lines
    # abline(v = 0, lty = 2, col = 4) # mean of central f in dashed blue line
    abline(v = falpha, lty = 2, col = 2) # f-critical in dashed red line

    # shaded area in light red for alpha
    xalpha0 <- seq(from = falpha, to = ncp + 4, by = .001)
    yalpha0 <- rep(NA, length(xalpha0))
    for(i in 1:length(xalpha0)){
      yalpha0[i] <- funf0(xalpha0[i])
    }

    xalpha <- c(xalpha0, rev(xalpha0))
    yalpha <- c(yalpha0, rep(0, length(yalpha0)))
    polygon(x = xalpha, y = yalpha, col = adjustcolor(2, alpha.f = 0.3), border = NA)

    # shaded area in light blue for beta
    xbeta0 <- seq(from = 0, to = falpha, by = .001)
    ybeta0 <- rep(NA, length(xbeta0))
    for(i in 1:length(xbeta0)){
      ybeta0[i] <- funf1(xbeta0[i])
    }
    xbeta <- c(xbeta0, rev(xbeta0))
    ybeta <- c(ybeta0, rep(0, length(ybeta0)))
    polygon(x = xbeta, y = ybeta, col = adjustcolor(4, alpha.f = 0.3), border = NA)

  }

  return(power)

}

#################################################################################

plot.pwrss <- function(x, ...) {

   if(all(c("pwrss","t") %in% class(x))) {
    power.t.test(ncp = abs(x$ncp),
                 df = x$df,
                 alpha = x$parms$alpha,
                 alternative = x$parms$alternative)
  } else if(all(c("pwrss","z") %in% class(x))) {
    power.z.test(ncp = abs(x$ncp),
                 alpha = x$parms$alpha,
                 alternative = x$parms$alternative)
  } else if(all(c("pwrss","f") %in% class(x))) {

    if("reg" %in% class(x)) {
      power.f.test(ncp = abs(x$ncp),
                   df1 = x$df1,
                   df2 = x$df2,
                   alpha = x$parms$alpha)
    } else if("ancova" %in% class(x)) {
      if(x$parms$n.way == 1) {
        power.f.test(ncp = abs(x$ncp),
                     df1 = x$df1,
                     df2 = x$df2,
                     alpha = x$parms$alpha)
      } else if(x$parms$n.way == 2) {

        layout(matrix(c(1,3,
                        2,3), nrow = 2, ncol = 2))
        default.pars <- par(no.readonly = TRUE)
        on.exit(par(default.pars))

        power.f.test(ncp = abs(x$ncp["A"]), plot.main = "A",
                     df1 = x$df1["A"],
                     df2 = x$df2["A"],
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["B"]), plot.main = "B",
                     df1 = x$df1["B"],
                     df2 = x$df2["B"],
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["AxB"]), plot.main = "A x B",
                     df1 = x$df1["AxB"],
                     df2 = x$df2["AxB"],
                     alpha = x$parms$alpha)

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
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["B"]), plot.main = "B",
                     df1 = x$df1["B"],
                     df2 = x$df2["B"],
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["C"]), plot.main = "C",
                     df1 = x$df1["C"],
                     df2 = x$df2["C"],
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["AxB"]), plot.main = "A x B",
                     df1 = x$df1["AxB"],
                     df2 = x$df2["AxB"],
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["AxC"]), plot.main = "A x C",
                     df1 = x$df1["AxC"],
                     df2 = x$df2["AxC"],
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["BxC"]), plot.main = "B x C",
                     df1 = x$df1["BxC"],
                     df2 = x$df2["BxC"],
                     alpha = x$parms$alpha)
        power.f.test(ncp = abs(x$ncp["AxBxC"]), plot.main = "A x B x C",
                     df1 = x$df1["AxBxC"],
                     df2 = x$df2["AxBxC"],
                     alpha = x$parms$alpha)

        par(mfrow = c(1,1))

      }

    } else if("rmanova" %in% class(x)) {
      power.f.test(ncp = abs(x$ncp),
                   df1 = x$df1,
                   df2 = x$df2,
                   alpha = x$parms$alpha)
    }

  } else {
    stop("not an object of the type 'pwrss'", call. = FALSE)
  }

}
