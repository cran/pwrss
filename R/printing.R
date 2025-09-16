.print.ascii.pwrss.logistic <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")
  cat("  Method          : ", x$method, "\n", sep = "")
  cat("  Predictor Dist. : ", x$dist, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    cat("  H0 (Null Claim) : Odds Ratio = 1\n")
    cat("  H1 (Alt. Claim) : Odds Ratio != 1\n\n")
  } else if (x$alt == "one.sided" && x$odds.ratio > 1) {
    cat("  H0 (Null Claim) : Odds Ratio <= 1\n")
    cat("  H1 (Alt. Claim) : Odds Ratio >  1\n\n")
  } else if (x$alt == "one.sided" && x$odds.ratio < 1) {
    cat("  H0 (Null Claim) : Odds Ratio >= 1\n")
    cat("  H1 (Alt. Claim) : Odds Ratio <  1\n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  Base Probability     = %.*f \n", digits, x$base.prob))
    cat(sprintf("  Odds Ratio           = %.*f \n", digits, x$odds.ratio))
    cat(sprintf("  Var. Corr. Factor    = %.*f \n", digits, x$vcf))
    cat(sprintf("  Mean of Alt.         = %.*f \n", digits, x$mean.alternative))
    cat(sprintf("  SD of Alt.           = %.*f \n", digits, x$sd.alternative))
    cat(sprintf("  Mean of Null         = %.*f \n", digits, x$mean.null))
    cat(sprintf("  SD of Null           = %.*f \n", digits, x$sd.null))
    cat(sprintf("  Critical Value       = %s \n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "  <<\n"))
  } else {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))

  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  Odds Ratio = [prob/(1-prob)] / [base.prob/(1-base.prob)] \n")
    cat("  prob       : Base probability when predictor = 0 \n")
    cat("  base.prob  : Probability when predictor = 1 \n")
    cat("  beta1      = log(Odds Ratio) \n")
    cat("  beta0      = log[base.prob/(1-base.prob)] \n\n")
  }

} # .print.ascii.pwrss.logistic()


.print.ascii.pwrss.poisson <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Method          : ", x$method, "\n", sep = "")
  cat("  Predictor Dist. : ", x$dist, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt ==  "two.sided") {
    cat("  H0 (Null Claim) : Rate Ratio = 1 \n")
    cat("  H1 (Alt. Claim) : Rate Ratio != 1 \n\n")
  } else if (x$alt == "one.sided" && x$rate.ratio > 1) {
    cat("  H0 (Null Claim) : Rate Ratio <= 1 \n")
    cat("  H1 (Alt. Claim) : Rate Ratio >  1 \n\n")
  } else if (x$alt == "one.sided" && x$rate.ratio < 1) {
    cat("  H0 (Null Claim) : Rate Ratio >= 1 \n")
    cat("  H1 (Alt. Claim) : Rate Ratio <  1 \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Base Rate            = %.*f \n", digits, x$base.rate))
    cat(sprintf("  Rate Ratio           = %.*f \n", digits, x$rate.ratio))
    cat(sprintf("  Var. Corr. Factor    = %.*f\n", digits, x$vcf))
    cat(sprintf("  Mean of Alt.         = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  SD of Alt.           = %.*f\n", digits, x$sd.alternative))
    cat(sprintf("  Mean of Null         = %.*f\n", digits, x$mean.null))
    cat(sprintf("  SD of Null           = %.*f\n", digits, x$sd.null))
    cat(sprintf("  Critical Value       = %s\n\n",
                paste(round(x$z.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "  <<\n"))
  } else {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "\n"))
  }
  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  Base Rate  = exp(beta0) \n")
    cat("  Rate Ratio = exp(beta1) \n\n")
  }

} # .print.ascii.pwrss.poisson()



.print.ascii.pwrss.regression <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt ==  "two.sided") {
    if (x$margin == 0) {
      cat("  H0 (Null Claim) : beta - null.beta = 0 \n")
      cat("  H1 (Alt. Claim) : beta - null.beta != 0 \n\n")
    } else {
      cat("  H0 (Null Claim) : beta - null.beta = margin \n")
      cat("  H1 (Alt. Claim) : beta - null.beta != margin \n\n")
    }
  } else if (x$alt == "one.sided" && x$ncp.alternative > x$ncp.null) {
    if (x$margin == 0) {
      cat("  H0 (Null Claim) : beta - null.beta <= 0 \n")
      cat("  H1 (Alt. Claim) : beta - null.beta >  0 \n\n")
    } else {
      cat("  H0 (Null Claim) : beta - null.beta <= margin \n")
      cat("  H1 (Alt. Claim) : beta - null.beta >  margin \n\n")
    }
  } else if (x$alt == "one.sided" && x$ncp.alternative < x$ncp.null) {
    if (x$margin == 0) {
      cat("  H0 (Null Claim) : beta - null.beta >= 0 \n")
      cat("  H1 (Alt. Claim) : beta - null.beta <  0 \n\n")
    } else {
      cat("  H0 (Null Claim) : beta - null.beta >= margin \n")
      cat("  H1 (Alt. Claim) : beta - null.beta <  margin \n\n")
    }
  } else if (x$alt == "two.one.sided" && (x$ncp.alternative > min(x$ncp.null) && x$ncp.alternative < max(x$ncp.null))) {
    cat("  H0 (Null Claim) : beta - null.beta <= min(margin) or \n                    beta - null.beta >= max(margin) \n")
    cat("  H1 (Alt. Claim) : beta - null.beta > min(margin) and \n                    beta - null.beta < max(margin)\n\n")
  } else if (x$alt == "two.one.sided" && (x$ncp.alternative < min(x$ncp.null) || x$ncp.alternative > max(x$ncp.null))) {
    cat("  H0 (Null Claim) : beta - null.beta >= min(margin) and \n                    beta - null.beta <= max(margin) \n")
    cat("  H1 (Alt. Claim) : beta - null.beta < min(margin) or \n                    beta - null.beta > max(margin) \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Std. Beta Under Alt.   = %.*f \n", digits, x$std.beta))
    cat(sprintf("  Std. Beta Under Null   = %.*f \n", digits, x$std.null.beta))
    cat(sprintf("  Std. Margin            = %s \n", paste(round(x$std.margin, digits), collapse = " and ")))
    cat(sprintf("  Degrees of Freedom     = %.*f \n", 0, x$df))
    cat(sprintf("  Non-centrality of Alt. = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  Non-centrality of Null = %s \n", paste(round(x$ncp.null, digits), collapse = " and ")))
    cat(sprintf("  Critical Value         = %s \n\n", paste(round(x$t.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0("  Sample Size            = ", round(x$n, digits = digits), "  <<\n"))
  } else {
    cat(paste0("  Sample Size            = ", round(x$n, digits = digits), "\n"))
  }
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error           = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  beta                 : Regression coefficient under alt. \n")
    cat("  null.beta            : Regression coefficient under null \n")
    cat("  margin               : Smallest beta - null.beta difference that matters \n\n")
    cat("  Std. Beta Under Alt. = beta * [SD(X) / SD(Y)]\n")
    cat("  Std. Beta Under Null = null.beta * [SD(X) / SD(Y)]\n")
    cat("  Std. Margin          = margin * [SD(X) / SD(Y)] \n\n")
    cat("  SD(X)                : Standard deviation of the predictor \n")
    cat("  SD(Y)                : Standard deviation of the outcome \n\n")
  }

} # .print.ascii.pwrss.regression()


.print.ascii.pwrss.f.regression <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$k.tested < x$k.total) {
    if (x$margin == 0) {
      cat("  H0 (Null Claim) : Change in R-squared = 0 \n")
      cat("  H1 (Alt. Claim) : Change in R-squared > 0 \n\n")
    } else {
      cat("  H0 (Null Claim) : 0 <= Change in R-squared <= margin \n")
      cat("  H1 (Alt. Claim) : Change in R-squared > margin \n\n")
    }
  } else {
    if (x$margin == 0) {
      cat("  H0 (Null Claim) : R-squared = 0 \n")
      cat("  H1 (Alt. Claim) : R-squared > 0 \n\n")
    } else {
      cat("  H0 (Null Claim) : 0 <= R-squared <= margin \n")
      cat("  H1 (Alt. Claim) : R-squared > margin \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    if (x$k.tested < x$k.total) {
      cat(sprintf("  Change in R-squared    = %.*f \n", digits, x$r.squared.change))
    } else {
      cat(sprintf("  R-squared              = %.*f \n", digits, x$r.squared.change))
    }
    cat(sprintf("  Margin                 = %.*f \n", digits, x$margin))
    cat(sprintf("  Num. Deg. of Freedom   = %.*f\n", 0, x$df1))
    cat(sprintf("  Denom. Deg. of Freedom = %.*f\n", 0, x$df2))
    cat(sprintf("  Non-centrality of Alt. = %.*f \n", digits, x$ncp.alternative))
    cat(sprintf("  Non-centrality of Null = %.*f \n", digits, x$ncp.null))
    cat(sprintf("  Critical Value         = %s \n\n", paste(round(x$f.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "  <<\n"))
  } else {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "\n"))
  }
  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    if (x$k.tested < x$k.total) {
      cat("  Margin : Smallest change in R-squared that matters \n\n")
    } else {
      cat("  Margin : Smallest R-squared that matters \n\n")
    }
  }

} # .print.ascii.pwrss.f.regression()



.print.ascii.pwrss.med <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `sobel` = "Sobel",
                   `aroian` = "Aroian",
                   `goodman` = "Goodman",
                   `joint` = "Joint",
                   `monte.carlo` = "Monte Carlo")

  cat(x$test, "\n\n", sep = "")
  cat("  Method            : ", method, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    cat("  H0 (Null Claim)   : beta[a*b] = 0 \n")
    cat("  H1 (Alt. Claim)   : beta[a*b] != 0 \n\n")
  } else if (x$alt == "one.sided") {
    if (x$std.beta.indirect < 0) {
      cat("  H0 (Null Claim)   : beta[a*b] >= 0 \n")
      cat("  H1 (Alt. Claim)   : beta[a*b] < 0 \n\n")
    } else {
      cat("  H0 (Null Claim)   : beta[a*b] <= 0 \n")
      cat("  H1 (Alt. Claim)   : beta[a*b] > 0 \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Std. beta[a]         = %.*f \n", digits, x$std.beta.a))
    cat(sprintf("  Std. beta[b]         = %.*f \n", digits, x$std.beta.b))
    if (x$method %in% c("sobel", "aorian", "goodman")) {
      cat(sprintf("  Std. beta[a*b]       = %.*f \n", digits, x$std.beta.indirect))
      cat(sprintf("  Mean of Alt.         = %.*f\n", digits, x$mean.alternative))
      cat(sprintf("  Mean of Null.        = %s \n", paste(round(x$mean.null, digits), collapse = " and ")))
      cat(sprintf("  Critical Value       = %s \n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
    } else {
      cat(sprintf("  Std. beta[a*b]       = %.*f \n\n", digits, x$std.beta.indirect))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "  <<\n"))
  } else {
    cat(paste0("  Sample Size          = ", round(x$n, digits = digits), "\n"))
  }
  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error         = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  beta[a]        : Regression coefficient for path a\n")
    cat("  beta[b]        : Regression coefficient for path b\n")
    cat("  beta[a*b]      : Coefficient for the indirect path a*b\n\n")
    cat("  Std. beta[a]   = beta[a] * [SD(predictor) / SD(mediator)]\n")
    cat("  Std. beta[b]   = beta[b] * [SD(mediator) / SD(outcome)]\n")
    cat("  Std. beta[a*b] = Std. beta[a] * Std. beta[b]\n\n")
    cat("  SD(predictor)  : Standard deviation of the predictor\n")
    cat("  SD(mediator)   : Standard deviation of the mediator\n")
    cat("  SD(outcome)    : Standard deviation of the outcome\n\n")
  }

} # .print.ascii.pwrss.med()



.print.ascii.pwrss.student <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    if (any(x$margin == 0)) {
      cat("  H0 (Null Claim) : d - null.d = 0 \n")
      cat("  H1 (Alt. Claim) : d - null.d != 0 \n\n")
    } else {
      cat("  H0 (Null Claim) : d - null.d = margin \n")
      cat("  H1 (Alt. Claim) : d - null.d != margin \n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$ncp.alternative < x$ncp.null) {
      if (x$margin == 0) {
        cat("  H0 (Null Claim) : d - null.d >= 0 \n")
        cat("  H1 (Alt. Claim) : d - null.d < 0 \n\n")
      } else {
        cat("  H0 (Null Claim) : d - null.d >= margin \n")
        cat("  H1 (Alt. Claim) : d - null.d < margin \n\n")
      }
    } else {
      if (any(x$margin == 0)) {
        cat("  H0 (Null Claim) : d - null.d <= 0 \n")
        cat("  H1 (Alt. Claim) : d - null.d > 0 \n\n")
      } else {
        cat("  H0 (Null Claim) : d - null.d <= margin \n")
        cat("  H1 (Alt. Claim) : d - null.d > margin \n\n")
      }
    }
  } else {
    if (x$ncp.alternative > min(x$ncp.null) && x$ncp.alternative < max(x$ncp.null)) {
      cat("  H0 (Null Claim) : d - null.d <= min(margin) or \n                    d - null.d >= max(margin) \n")
      cat("  H1 (Alt. Claim) : d - null.d > min(margin) and \n                    d - null.d < max(margin) \n\n")
    } else {
      cat("  H0 (Null Claim  : d - null.d >= min(margin) and \n                    d - null.d <= max(margin) \n")
      cat("  H1 (Alt. Claim) : d - null.d < min(margin) or \n                    d - null.d > max(margin) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Cohen's d              = %.*f\n", digits, x$d))
    cat(sprintf("  Cohen's d Under Null   = %.*f\n", digits, x$null.d))
    cat(sprintf("  Margin                 = %s\n", paste(round(x$margin, digits), collapse = " and ")))
    cat(sprintf("  Degrees of Freedom     = %.*f\n", 0, x$df))
    cat(sprintf("  Non-centrality of Alt. = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  Non-centrality of Null = %s\n", paste(round(x$ncp.null, digits), collapse = " and ")))
    cat(sprintf("  Critical Value         = %s \n\n", paste(round(x$t.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("  Sample Size            = ", n.text, "  <<\n"))
  } else {
    cat(paste0("  Sample Size            = ", n.text, "\n"))
  }
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  Margin : Smallest d - null.d difference that matters \n\n")
  }

} # .print.ascii.pwrss.student()



.print.ascii.pwrss.wilcoxon <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  dist <- switch(x$dist,
                 `normal` = "Normal",
                 `uniform` = "Uniform",
                 `double.exponential` = "Double Exponential",
                 `laplace` = "Laplace",
                 `logistic` = "Logistic")
  method <- switch(x$method,
                   `guenther` = "Guenther",
                   `noether` = "Noether")

  cat(x$test, "\n\n", sep = "")
  cat("  Method       : ", method, "\n", sep = "")
  cat("  Distribution : ", dist, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    if (x$margin == 0) {
      cat("  H0 (Null Claim) : d - null.d = 0 \n")
      cat("  H1 (Alt. Claim) : d - null.d != 0 \n\n")
    } else {
      cat("  H0 (Null Claim) : d - null.d = margin \n")
      cat("  H1 (Alt. Claim) : d - null.d != margin \n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$method == "guenther") {
      is.less <- x$ncp < x$null.ncp
    } else {
      is.less <- x$mean < x$null.mean
    }
    if (is.less) {
      if (x$margin == 0) {
        cat("  H0 (Null Claim) : d - null.d >= 0 \n")
        cat("  H1 (Alt. Claim) : d - null.d < 0 \n\n")
      } else {
        cat("  H0 (Null Claim) : d - null.d >= margin \n")
        cat("  H1 (Alt. Claim) : d - null.d < margin \n\n")
      }
    } else {
      if (x$margin == 0) {
        cat("  H0 (Null Claim) : d - null.d <= 0 \n")
        cat("  H1 (Alt. Claim) : d - null.d > 0 \n\n")
      } else {
        cat("  H0 (Null Claim) : d - null.d <= margin \n")
        cat("  H1 (Alt. Claim) : d - null.d > margin \n\n")
      }
    }
  } else {
    if (x$method == "guenther") {
      is.equivalent <- x$ncp > min(x$null.ncp) && x$ncp < max(x$null.ncp)
    } else {
      is.equivalent <- x$mean > min(x$null.mean) && x$mean < max(x$null.mean)
    }

    if (is.equivalent) {
      cat("  H0 (Null Claim) : d - null.d <= min(margin) or \n                    d - null.d >= max(margin) \n")
      cat("  H1 (Alt. Claim) : d - null.d > min(margin) and \n                    d - null.d < max(margin) \n\n")
    } else {
      cat("  H0 (Null Claim  : d - null.d >= min(margin) and \n                    d - null.d <= max(margin) \n")
      cat("  H1 (Alt. Claim) : d - null.d < min(margin) or \n                    d - null.d > max(margin) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Cohen's d              = %.*f\n", digits, x$d))
    cat(sprintf("  Cohen's d Under Null   = %.*f\n", digits, x$null.d))
    cat(sprintf("  Margin                 = %s\n", paste(round(x$margin, digits), collapse = " and ")))
    if (x$method == "guenther") {
      cat(sprintf("  Degrees of Freedom     = %.*f\n", 0, x$df))
      cat(sprintf("  Non-centrality of Alt. = %.*f\n", digits, x$ncp))
      cat(sprintf("  Non-centrality of Null = %s\n", paste(round(x$null.ncp, digits), collapse = " and ")))
      cat(sprintf("  Critical Value         = %s\n\n", paste(round(x$t.alpha, digits), collapse = " and ")))
    } else {
      cat(sprintf("  Mean of Alt.           = %.*f\n", digits, x$mean))
      cat(sprintf("  Mean of Null           = %s\n", paste(round(x$null.mean, digits), collapse = " and ")))
      cat(sprintf("  Critical Value         = %s\n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("  Sample Size            = ", n.text, "  <<\n"))
  } else {
    cat(paste0("  Sample Size            = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  Margin : Smallest d - null.d difference that matters \n\n")
  }

} # .print.ascii.pwrss.wilcoxon()



.print.ascii.pwrss.gof <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  cat("  H0 (Null Claim)   : P[i,j] = P0[i,j] for all (i,j) \n")
  cat("  H1 (Alt. Claim)   : P[i,j] != P0[i,j] for some (i,j)\n\n")

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Degrees of Freedom     = %.*f \n", 0, x$df))
    cat(sprintf("  Non-centrality of Alt. = %.*f \n", digits, x$ncp.alternative))
    cat(sprintf("  Non-centrality of Null = %.*f \n", digits, x$ncp.null))
    cat(sprintf("  Critical Value         = %s\n\n",
                paste(round(x$chisq.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0("  Total Sample Size      = ", round(x$n, digits = digits), "  << \n"))
  } else {
    cat(paste0("  Total Sample Size      = ", round(x$n, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  For goodness-of-fit, comparisons are P[i] vs P0[i] \n")
    cat("  For independence, comparisons are P[i,j] vs P0[i,j] \n")
    cat("  Independence implies (default) P0[i,j] = P[i,.] * P[.,j] \n\n")
    cat("  P[i,j] : Joint probability for cell (i,j) \n")
    cat("  P[i,.] : Marginal probability for row i (sum over j) \n")
    cat("  P[.,j] : Marginal probability for column j (sum over i) \n\n")
  }

} # .print.ascii.pwrss.gof()



.print.ascii.pwrss.chisq <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$ncp.null > 0) {
    cat("  H0 (Null Claim)   : 0 <= ncp <= null.ncp \n")
    cat("  H1 (Alt. Claim)   : ncp > null.ncp \n\n")
  } else {
    cat("  H0 (Null Claim)   : ncp = null.ncp \n")
    cat("  H1 (Alt. Claim)   : ncp > null.ncp \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Degrees of Freedom     = %.*f\n", 0, x$df))
    cat(sprintf("  Non-centrality of Alt. = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  Non-centrality of Null = %.*f\n", digits, x$ncp.null))
    cat(sprintf("  Critical Value         = %s\n\n",
                paste(round(x$chisq.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  ncp      : Non-centrality parameter of alt. \n")
    cat("  null.ncp : Non-centrality parameter of null \n\n")
  }

} # .print.ascii.pwrss.chisq()



.print.ascii.pwrss.t <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    cat("  H0 (Null Claim) : ncp = null.ncp \n")
    cat("  H1 (Alt. Claim) : ncp != null.ncp \n\n")
  } else if (x$alt == "one.sided") {
    if (x$ncp.alternative < x$ncp.null) {
      cat("  H0 (Null Claim) : ncp >= null.ncp \n")
      cat("  H1 (Alt. Claim) : ncp < null.ncp \n\n")
    } else {
      cat("  H0 (Null Claim) : ncp <= null.ncp \n")
      cat("  H1 (Alt. Claim) : ncp > null.ncp \n\n")
    }
  } else {
    if (x$ncp.alternative > min(x$ncp.null) && x$ncp.alternative < max(x$ncp.null)) {
      cat("  H0 (Null Claim) : ncp <= min(null.ncp) or \n                    ncp >= max(null.ncp) \n")
      cat("  H1 (Alt. Claim) : ncp > min(null.ncp) and \n                    ncp < max(null.ncp) \n\n")
    } else {
      cat("  H0 (Null Claim) : ncp >= min(null.ncp) and \n                    ncp <= max(null.ncp) \n")
      cat("  H1 (Alt. Claim) : ncp < min(null.ncp) or \n                    ncp > max(null.ncp) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    if (isFALSE(is.infinite(x$df))) {
      cat(sprintf("  Degrees of Freedom     = %.*f\n", 0, x$df))
    }
    cat(sprintf("  Non-centrality of Alt. = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  Non-centrality of Null = %s\n", paste(round(x$ncp.null, digits), collapse = " and ")))
    cat(sprintf("  Critical Value         = %s\n\n", paste(round(x$t.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "  << \n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  ncp      : Non-centrality parameter of Alt. \n")
    cat("  null.ncp : Non-centrality parameter of Null \n\n")
  }

} # .print.ascii.pwrss.t()



.print.ascii.pwrss.z <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)


  if (x$alt == "two.sided") {
    cat("  H0 (Null Claim) : mean = null.mean \n")
    cat("  H1 (Alt. Claim) : mean != null.mean \n\n")
  } else if (x$alt == "one.sided") {
    if (x$mean.alternative < x$mean.null) {
      cat("  H0 (Null Claim) : mean >= null.mean \n")
      cat("  H1 (Alt. Claim) : mean < null.mean \n\n")
    } else {
      cat("  H0 (Null Claim) : mean <= null.mean \n")
      cat("  H1 (Alt. Claim) : mean > null.mean \n\n")
    }
  } else {
    if (x$mean.alternative > min(x$mean.null) && x$mean.alternative < max(x$mean.null)) {
      cat("  H0 (Null Claim) : mean <= min(null.mean) or \n                    mean >= max(null.mean) \n")
      cat("  H1 (Alt. Claim) : mean > min(null.mean) and \n                    mean < max(null.mean) \n\n")
    } else {
      cat("  H0 (Null Claim) : mean >= min(null.mean) and \n                    mean <= max(null.mean) \n")
      cat("  H1 (Alt. Claim) : mean < min(null.mean) or \n                    mean > max(null.mean) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Mean of Alt.         = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  Mean of Null         = %s\n", paste(round(x$mean.null, digits), collapse = " and ")))
    cat(sprintf("  Critical Value       = %s\n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  mean      : Mean of alt. \n")
    cat("  null.mean : Mean of null \n\n")
  }

} # .print.ascii.pwrss.z()



.print.ascii.pwrss.binom <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    cat("  H0 (Null Claim) : prob = null.prob \n")
    cat("  H1 (Alt. Claim) : prob != null.prob \n\n")
  } else if (x$alt == "one.sided") {
    if (x$prob.alt < x$prob.null) {
      cat("  H0 (Null Claim) : prob >= null.prob \n")
      cat("  H1 (Alt. Claim) : prob < null.prob \n\n")
    } else {
      cat("  H0 (Null Claim) : prob <= null.prob \n")
      cat("  H1 (Alt. Claim) : prob > null.prob \n\n")
    }
  } else {
    if (x$prob.alt  > min(x$prob.null) && x$prob.alt < max(x$prob.null)) {
      cat("  H0 (Null Claim) : prob <= min(null.prob) or \n                    prob >= max(null.prob) \n")
      cat("  H1 (Alt. Claim) : prob > min(null.prob) and \n                    prob < max(null.prob) \n\n")
    } else {
      cat("  H0 (Null Claim) : prob >= min(null.prob) and \n                    prob <= max(null.prob) \n")
      cat("  H1 (Alt. Claim) : prob < min(null.prob) or \n                    prob > max(null.prob) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Size                   = %.*f\n", 0, x$size))
    cat(sprintf("  Probability Under Alt. = %.*f\n", digits, x$prob.alt))
    cat(sprintf("  Probability Under Null = %s\n", paste(round(x$prob.null, digits), collapse = " and ")))
    cat(sprintf("  Critical Value         = %s\n\n", paste(round(x$binom.alpha, 0), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  prob      : Probability under alt. \n")
    cat("  null.prob : Probability under null \n\n")
  }

} # .print.ascii.pwrss.binom()



.print.ascii.pwrss.f <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$ncp.null > 0) {
    cat("  H0 (Null Claim) : 0 <= ncp <= null.ncp \n")
    cat("  H1 (Alt. Claim) : ncp > null.ncp \n\n")
  } else {
    cat("  H0 (Null Claim) : ncp = null.ncp \n")
    cat("  H1 (Alt. Claim) : ncp > null.ncp \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Num. Deg. of Freedom   = %.*f \n", 0, x$df1))
    cat(sprintf("  Denom. Deg. of Freedom = %.*f \n", 0, x$df2))
    cat(sprintf("  Non-centrality of Alt. = %.*f \n", digits, x$ncp.alternative))
    cat(sprintf("  Non-centrality of Null = %.*f \n", digits, x$ncp.null))
    cat(sprintf("  Critical Value         = %s \n\n", paste(round(x$f.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  ncp      : Non-centrality parameter of alt. \n")
    cat("  null.ncp : Non-centrality parameter of null \n\n")
  }

} # .print.ascii.pwrss.f()




.print.ascii.pwrss.ancova <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n.total",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$null.ncp == 0) {
    cat("  H0 (Null Claim) : eta.squared = 0 \n")
    cat("  H1 (Alt. Claim) : eta.squared > 0 \n\n")
  } else {
    cat("  H0 (Null Claim) : 0 <= eta.squared <= null.eta.squared \n")
    cat("  H1 (Alt. Claim) : eta.squared > null.eta.squared\n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Design                 = %s \n", x$effect))
    cat(sprintf("  Num. Deg. of Freedom   = %.*f \n", 0, x$df1))
    cat(sprintf("  Denom. Deg. of Freedom = %.*f \n", 0, x$df2))
    cat(sprintf("  Non-centrality of Alt. = %.*f \n", digits, x$ncp))
    cat(sprintf("  Non-centrality of Null = %.*f \n", digits, x$null.ncp))
    cat(sprintf("  Critical Value         = %s \n\n", paste(round(x$f.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n.total") {
    cat(paste0("  Total Sample Size      = ", round(x$n.total, digits), "  <<\n"))
  } else {
    cat(paste0("  Total Sample Size      = ", round(x$n.total, digits), "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)    = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power      = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    if (x$null.ncp != 0) cat("  null.eta.squared : (Partial) Eta-squared under null \n")
    cat("  eta.squared      : (Partial) Eta-squared under alt. \n\n")
  }

} # .print.ascii.pwrss.ancova()




.print.ascii.pwrss.contrast <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n.total",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  cat("  H0 (Null Claim) : psi = 0 \n")
  cat("  H1 (Alt. Claim) : psi != 0 \n\n")

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Contrast Est. (psi)    = %.*f \n", digits, x$psi))
    cat(sprintf("  Standardized psi (d)   = %.*f \n", digits, x$d))
    cat(sprintf("  Degrees of Freedom     = %.*f \n", 0, x$df))
    cat(sprintf("  Non-centrality of Alt. = %.*f \n", digits, x$ncp))
    cat(sprintf("  Non-centrality of Null = %.*f \n", digits, x$null.ncp))
    cat(sprintf("  Critical Value         = %s \n\n", paste(round(x$t.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n.total") {
    cat(paste0("  Total Sample Size     = ", round(x$n.total, digits), "  <<\n"))
  } else {
    cat(paste0("  Total Sample Size     = ", round(x$n.total, digits), "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)   = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power     = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power     = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  psi : Contrast estimate, sum(contrast[i] * mu[i]) \n")
    cat("  d   : Standardized contrast estimate \n\n")
  }

} # .print.ascii.pwrss.contrast()



.print.ascii.pwrss.contrasts <- function(x, data = NULL, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n.total",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  cat("  H0 (Null Claim) : psi = 0 \n")
  cat("  H1 (Alt. Claim) : psi != 0 \n\n")

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    adjust.alpha <- switch(x$adjust.alpha,
                           `none` = "None",
                           `fdr` = "False Discovery Rate",
                           `hochberg` = "Hochberg (1988)",
                           `BH` = "Benjamini & Hochberg (1995)",
                           `BY` = "Benjamini & Yakutieli (2001)",
                           `holm` = "Holm (1979)",
                           `hommel` = "Hommel (1988)",
                           `bonferroni` = "Bonferroni",
                           `tukey` = "Tukey")
    cat(sprintf("  Alpha Adjustment       = %s \n", adjust.alpha))
    cat(sprintf("  Adjusted Alpha         = %.*f \n", digits, x$alpha))
    cat(sprintf("  Non-centrality of Null = %.*f \n\n", digits, x$null.ncp))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (!is.null(x$data)) {
    print(x$data, row.names = FALSE)
    cat("\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  psi : Contrast estimate, sum(contrast[i] * mu[i])\n")
    cat("  d   : Standardized contrast estimate \n")
    cat("  ncp : Non-centrality parameter under alt. \n\n")
  }

} # .print.ascii.pwrss.contrasts()




.print.ascii.pwrss.fisher <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `z` = "Normal Approximation",
                   `exact` = "Fisher's Exact")

  cat(x$test, "\n\n", sep = "")
  cat("  Method          : ", method, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    if (any(x$margin == 0)) {
      cat("  H0 (Null Claim) : prob1 - prob2 = 0 \n")
      cat("  H1 (Alt. Claim) : prob1 - prob2 != 0 \n\n")
    } else {
      cat("  H0 (Null Claim) : prob1 - prob2 = margin \n")
      cat("  H1 (Alt. Claim) : prob1 - prob2 != margin \n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      if (any(x$margin == 0)) {
        cat("  H0 (Null Claim) : prob1 - prob2 >= 0 \n")
        cat("  H1 (Alt. Claim) : prob1 - prob2 < 0 \n\n")
      } else {
        cat("  H0 (Null Claim) : prob1 - prob2 >= margin \n")
        cat("  H1 (Alt. Claim) : prob1 - prob2 < margin \n\n")
      }
    } else {
      if (any(x$margin == 0)) {
        cat("  H0 (Null Claim) : prob1 - prob2 <= 0 \n")
        cat("  H1 (Alt. Claim) : prob1 - prob2 > 0 \n\n")
      } else {
        cat("  H0 (Null Claim) : prob1 - prob2 <= margin \n")
        cat("  H1 (Alt. Claim) : prob1 - prob2 > margin \n\n")
      }
    }
  } else {
    if (x$delta > min(x$margin) && x$delta < max(x$margin)) {
      cat("  H0 (Null Claim) : prob1 - prob2 <= min(margin) or \n                    prob1 - prob2 >= max(margin) \n")
      cat("  H1 (Alt. Claim) : prob1 - prob2 > min(margin) and \n                    prob1 - prob2 < max(margin) \n\n")
    } else {
      cat("  H0 (Null Claim) : prob1 - prob2 >= min(margin) and \n                    prob1 - prob2 <= max(margin) \n")
      cat("  H1 (Alt. Claim) : prob1 - prob2 < min(margin) or \n                    prob1 - prob2 > max(margin) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  prob1 - prob2        = %.*f\n", digits, x$delta))
    cat(sprintf("  Odds Ratio           = %.*f\n", digits, x$odds.ratio))
    if (x$method == "z") {
      cat(sprintf("  Mean of Alt.         = %.*f\n", digits, x$mean.alternative))
      cat(sprintf("  Mean of Null         = %s\n", paste(round(x$mean.null, digits), collapse = " and ")))
      cat(sprintf("  Critical Value       = %s\n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
    } else {
      cat("\n")
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)

  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("  Sample Size          = ", n.text, "  <<\n"))
  } else {
    cat(paste0("  Sample Size          = ", n.text, "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))

  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  prob1       : Probability of success in the first group \n")
    cat("  prob2       : Probability of success in the second group \n")
    if (any(x$margin != 0))
      cat("  margin      : Smallest prob1 - prob2 that matters \n")
    cat("  Odds Ratio  : Odds(prob1) / Odds(prob2) \n")
    cat("  Odds(prob1) : prob1 / (1 - prob1) \n")
    cat("  Odds(prob2) : prob2 / (1 - prob2) \n\n")
  }

} # .print.ascii.pwrss.fisher()




.print.ascii.pwrss.mcnemar <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `z` = "Normal Approximation",
                   `exact` = "McNemar's Exact")

  cat(x$test, "\n\n", sep = "")
  cat("  Method          : ", method, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    cat("  H0 (Null Claim) : prob10 - prob01 = 0\n")
    cat("  H1 (Alt. Claim) : prob10 - prob01 != 0\n\n")
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      cat("  H0 (Null Claim) : prob10 - prob01 >= 0\n")
      cat("  H1 (Alt. Claim) : prob10 - prob01 < 0\n\n")
    } else {
      cat("  H0 (Null Claim) : prob10 - prob01 <= 0\n")
      cat("  H1 (Alt. Claim) : prob10 - prob01 > 0\n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    cat(sprintf("  Odds Ratio           = %.*f\n", digits, x$odds.ratio))
    cat(sprintf("  prob10 - prob01      = %.*f\n", digits, x$delta))
    if (x$method == "exact") {
      cat(sprintf("  Size of Disc. Pairs  = %.*f \n", 0, x$size))
      cat(sprintf("  prob10 || DP for Alt. = %.*f \n", digits, x$prob.alternative))
      cat(sprintf("  prob10 || DP for Null = %s \n", paste(round(x$prob.null, digits), collapse = " and ")))
      cat(sprintf("  Critical Value       = %s \n\n", paste(round(x$binom.alpha, digits), collapse = " and ")))
    } else {
      cat(sprintf("  Mean of Alt.        = %.*f\n", digits, x$mean.alternative))
      cat(sprintf("  Mean of Null        = %s\n", paste(round(x$mean.null, digits), collapse = " and ")))
      cat(sprintf("  Critical Value      = %s\n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)

  n.text <- paste(round(x$n.paired, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("  Paired Sample Size   = ", n.text, "  <<\n"))
  } else {
    cat(paste0("  Paired Sample Size   = ", n.text, "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))

  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  prob10      : Joint prob. of observing {1,0} \n")
    cat("  prob01      : Joint prob. of observing {0,1} \n")
    cat("  Odds Ratio  : prob10 / prob01 \n")
    cat("  prob10 | DP : Conditional prob. of observing {1,0} \n                among DP, prob10 / (prob10 + prob01) \n")
    cat("  DP          : Discordant pairs \n\n")
  }

} # .print.ascii.pwrss.mcnemar()




.print.ascii.pwrss.oneprop <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `z` = "Normal Approximation",
                   `exact` = "Exact")

  cat(x$test, "\n\n", sep = "")

  if (x$method == "exact") {
    cat("  Method                 : ", method, "\n\n", sep = "")
  } else {
    stderr <- switch(x$std.err,
                     `alternative` = "Alternative",
                     `null` = "Null")
    cat("  Method                 : ", method, "\n", sep = "")
    cat("  Continuity Correction  : ", x$correct, "\n", sep = "")
    cat("  Arcsine Transformation : ", x$arcsine, "\n", sep = "")
    cat("  Standard Error         : Calculated From ", stderr, "\n\n", sep = "")
  }

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    cat("  H0 (Null Claim)        : prob - null.prob = 0\n")
    cat("  H1 (Alt. Claim)        : prob - null.prob != 0\n\n")
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      cat("  H0 (Null Claim)        : prob - null.prob >= 0\n")
      cat("  H1 (Alt. Claim)        : prob - null.prob < 0\n\n")
    } else {
      cat("  H0 (Null Claim)        : prob - null.prob <= 0\n")
      cat("  H1 (Alt. Claim)        : prob - null.prob > 0\n\n")
    }
  } else if (x$alt == "two.one.sided") {
    if (x$delta[1] > 0 && x$delta[2] < 0) {
      cat("  H0 (Null Claim)        : prob - min(null.prob) <= 0 or \n                           prob - max(null.prob) >= 0\n")
      cat("  H1 (Alt. Claim)        : prob - min(null.prob) > 0 and \n                           prob - max(null.prob) < 0\n\n")
    } else {
      cat("  H0 (Null Claim)        : prob - min(null.prob) >= 0 and \n                           prob - max(null.prob) <= 0\n")
      cat("  H1 (Alt. Claim)        : prob - min(null.prob) < 0 or \n                           prob - max(null.prob) > 0\n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    delta.text <- paste(round(x$delta, digits), collapse = " and ")
    cat(sprintf("  prob - null.prob      = %s\n", delta.text))
    or.text <- paste(round(x$odds.ratio, digits), collapse = " and ")
    cat(sprintf("  Odds Ratio            = %s\n", or.text))
    if (x$method == "exact") {
      cat(sprintf("  Size                  = %.*f\n", 0, x$size))
      cat(sprintf("  Prob. Under Alt       = %.*f\n", digits, x$prob.alternative))
      cat(sprintf("  Prob. Under Null      = %s\n", paste(round(x$prob.null, digits), collapse = " and ")))
      cat(sprintf("  Critical Value        = %s\n\n", paste(round(x$binom.alpha, digits), collapse = " and ")))
    } else {
      cat(sprintf("  Mean of Alt.          = %.*f\n", digits, x$mean.alternative))
      cat(sprintf("  Mean of Null          = %s\n", paste(round(x$mean.null, digits), collapse = " and ")))
      cat(sprintf("  Critical Value        = %s\n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)

  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("  Sample Size           = ", n.text, "  <<\n"))
  } else {
    cat(paste0("  Sample Size           = ", n.text, "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)   = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power     = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power     = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    cat("  Odds Ratio      : Odds(prob) / Odds(null.prob) \n")
    cat("  Odds(prob)      : prob / (1 - prob) \n")
    cat("  Odds(null.prob) : null.prob / (1 - null.prob) \n\n")
  }

} # .print.ascii.pwrss.oneprop()




.print.ascii.pwrss.steiger <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n", sep = "")
  cat("  Common Index    : ", x$common, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    if (x$common) {
      cat("  H0 (Null Claim) : rho12 - rho13 = 0\n")
      cat("  H1 (Alt. Claim) : rho12 - rho13 != 0\n\n")
    } else {
      cat("  H0 (Null Claim) : rho12 - rho34 = 0\n")
      cat("  H1 (Alt. Claim) : rho12 - rho34 != 0\n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      if (x$common) {
        cat("  H0 (Null Claim) : rho12 - rho13 >= 0\n")
        cat("  H1 (Alt. Claim) : rho12 - rho13 < 0\n\n")
      } else {
        cat("  H0 (Null Claim) : rho12 - rho34 >= 0\n")
        cat("  H1 (Alt. Claim) : rho12 - rho34 < 0\n\n")
      }
    } else {
      if (x$common) {
        cat("  H0 (Null Claim) : rho12 - rho13 <= 0\n")
        cat("  H1 (Alt. Claim) : rho12 - rho13 > 0\n\n")
      } else {
        cat("  H0 (Null Claim) : rho12 - rho34 <= 0\n")
        cat("  H1 (Alt. Claim) : rho12 - rho34 > 0\n\n")
      }
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    if (x$common) {
      cat(sprintf("  rho12 - rho13        = %.*f\n", digits, x$delta))
    } else {
      cat(sprintf("  rho12 - rho34        = %.*f\n", digits, x$delta))
    }
    cat(sprintf("  Cohen's q            = %.*f\n", digits, x$q))
    cat(sprintf("  Mean of Alt.         = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  SD of Alt.           = %.*f\n", digits, x$sd.alternative))
    cat(sprintf("  Mean of Null         = %.*f\n", digits, x$mean.null))
    cat(sprintf("  SD of Null           = %.*f\n", digits, x$sd.null))
    cat(sprintf("  Critical Value       = %s\n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)

  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("  Sample Size          = ", n.text, "  <<\n"))
  } else {
    cat(paste0("  Sample Size          = ", n.text, "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    if (x$common) {
      cat("  rho12 : Correlation between variable V1 and V2 \n")
      cat("  rho13 : Correlation between variable V1 and V3 \n\n")
    } else {
      cat("  rho12 : Correlation between variable V1 and V2 \n")
      cat("  rho34 : Correlation between variable V3 and V4 \n\n")
    }
  }

} # .print.ascii.pwrss.steiger()




.print.ascii.pwrss.twocors <- function(x, digits = 3, verbose = 1, ...) {

  UL <- UR <- LL <- LR <- "+"
  HL <- "-"
  VL <- "|"
  HR <- strrep(HL, 50)
  line <- paste0("-", strrep("-", 50), "\n")

  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "             SAMPLE SIZE CALCULATION              ",
                 "                POWER CALCULATION                 "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")

  cat(x$test, "\n\n")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    if (x$design %in% c("independent", "paired")) {
      cat("  H0 (Null Claim) : rho1 - rho2 = 0\n")
      cat("  H1 (Alt. Claim) : rho1 - rho2 != 0\n\n")
    } else {
      cat("  H0 (Null Claim) : rho -  null.rho = 0\n")
      cat("  H1 (Alt. Claim) : rho -  null.rho != 0\n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      if (x$design %in% c("independent", "paired")) {
        cat("  H0 (Null Claim) : rho1 - rho2 >= 0\n")
        cat("  H1 (Alt. Claim) : rho1 - rho2 < 0\n\n")
      } else {
        cat("  H0 (Null Claim) : rho -  null.rho >= 0\n")
        cat("  H1 (Alt. Claim) : rho -  null.rho < 0\n\n")
      }
    } else {
      if (x$design %in% c("independent", "paired")) {
        cat("  H0 (Null Claim) : rho1 - rho2 <= 0\n")
        cat("  H1 (Alt. Claim) : rho1 - rho2 > 0\n\n")
      } else {
        cat("  H0 (Null Claim) : rho -  null.rho <= 0\n")
        cat("  H1 (Alt. Claim) : rho -  null.rho > 0\n\n")
      }
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters\n")
    cat(line)
    if (x$design %in% c("independent", "paired")) {
      cat(sprintf("  rho1 - rho2          = %.*f\n", digits, x$delta))
    } else {
      cat(sprintf("  rho - null.rho       = %.*f\n", digits, x$delta))
    }
    cat(sprintf("  Cohen's q            = %.*f\n", digits, x$q))
    cat(sprintf("  Mean of Alt.         = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  SD of Alt.           = %.*f\n", digits, x$sd.alternative))
    cat(sprintf("  Mean of Null         = %.*f\n", digits, x$mean.null))
    cat(sprintf("  SD of Null           = %.*f\n", digits, x$sd.null))
    cat(sprintf("  Critical Value       = %s\n\n", paste(round(x$z.alpha, digits), collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)

  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("  Sample Size          = ", n.text, "  <<\n"))
  } else {
    cat(paste0("  Sample Size          = ", n.text, "\n"))
  }

  cat(sprintf("  Type 1 Error (alpha) = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (beta)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "  <<\n\n"))
  } else {
    cat(paste0("  Statistical Power    = ", round(x$power, digits), "\n\n"))
  }

  if (verbose == 2) {
    cat(line)
    cat("Definitions\n")
    cat(line)
    if (x$design %in% c("independent", "paired")) {
      cat("  rho1 : Correlation in group 1 \n")
      cat("  rho2 : Correlation in group 2 \n\n")
    } else {
      cat("  rho      : Correlation under alt. \n")
      cat("  null.rho : Correlation under null \n\n")
    }
  }

} # .print.ascii.pwrss.twocors()
