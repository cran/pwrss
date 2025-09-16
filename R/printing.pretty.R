
.print.pwrss.logistic <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Method           : ", x$method, "\n", sep = "")
  cat("  Predictor Dist.  : ", x$dist, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    cat("  H\u2080 (Null)        : Odds Ratio (OR) = 1\n")
    cat("  H\u2081 (Alternative) : Odds Ratio (OR) \u2260 1\n\n")
  } else if (x$alt == "one.sided" && x$odds.ratio > 1) {
    cat("  H\u2080 (Null)        : Odds Ratio (OR) \u2264 1\n")
    cat("  H\u2081 (Alternative) : Odds Ratio (OR) > 1\n\n")
  } else if (x$alt == "one.sided" && x$odds.ratio < 1) {
    cat("  H\u2080 (Null)        : Odds Ratio (OR) \u2265 1\n")
    cat("  H\u2081 (Alternative) : Odds Ratio (OR) < 1\n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  Base Probability  = %.*f \n", digits, x$base.prob))
    cat(sprintf("  Odds Ratio (OR)   = %.*f \n", digits, x$odds.ratio))
    cat(sprintf("  Var. Corr. Factor = %.*f\n", digits, x$vcf))
    cat(sprintf("  \u03BC (Alternative)   = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  \u03C3 (Alternative)   = %.*f\n", digits, x$sd.alternative))
    cat(sprintf("  \u03BC\u2080 (Null)         \u2009= %.*f \n", digits, x$mean.null))
    cat(sprintf("  \u03C3\u2080 (Null)         \u2009= %.*f \n", digits, x$sd.null))
    cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)     = %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0(" \033[34m Sample Size       = ",  round(x$n, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  round(x$n, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", " OR = [p / (1 \u2212 p)] / [p\u2080 / (1 \u2212 p\u2080)]\n", "\033[0m")
    cat("\033[36m", "\u03B2\u2081 \u2009= log(OR)\n", "\033[0m")
    cat("\033[36m", "\u03B2\u2080 \u2009= log[p\u2080 / (1 \u2212 p\u2080)]\n", "\033[0m")
    cat("\033[36m", "p\u2080 \u2009: Base probability when predictor = 0 \n", "\033[0m")
    cat("\033[36m", "p  : Probability when predictor = 1 \n", "\033[0m")
    cat("\033[36m", " \u03BC  : Mean \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3  : Standard deviation \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.logistic()




.print.pwrss.poisson <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Method           : ", x$method, "\n", sep = "")
  cat("  Predictor Dist.  : ", x$dist, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt ==  "two.sided") {
    cat("  H\u2080 (Null)        : Rate Ratio (RR) = 1\n")
    cat("  H\u2081 (Alternative) : Rate Ratio (RR) \u2260 1\n\n")
  } else if (x$alt == "one.sided" && x$rate.ratio > 1) {
    cat("  H\u2080 (Null)        : Rate Ratio (RR) \u2264 1\n")
    cat("  H\u2081 (Alternative) : Rate Ratio (RR) > 1\n\n")
  } else if (x$alt == "one.sided" && x$rate.ratio < 1) {
    cat("  H\u2080 (Null)        : Rate Ratio (RR) \u2265 1\n")
    cat("  H\u2081 (Alternative) : Rate Ratio (RR) < 1\n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  Base Rate         = %.*f \n", digits, x$base.rate))
    cat(sprintf("  Rate Ratio (RR)   = %.*f \n", digits, x$rate.ratio))
    cat(sprintf("  Var. Corr. Factor = %.*f\n", digits, x$vcf))
    cat(sprintf("  \u03BC (Alternative)   = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  \u03C3 (Alternative)   = %.*f\n", digits, x$sd.alternative))
    cat(sprintf("  \u03BC\u2080 (Null)         \u2009= %.*f \n", digits, x$mean.null))
    cat(sprintf("  \u03C3\u2080 (Null)         \u2009= %.*f \n", digits, x$sd.null))
    cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)     = %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0(" \033[34m Sample Size       = ",  round(x$n, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  round(x$n, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  Base Rate       = exp(\u03B2\u2080)\n", "\033[0m", sep = "")
    cat("\033[36m", "  Rate Ratio (RR) = exp(\u03B2\u2081)\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BC               : Mean \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3               : Standard deviation \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.poisson()


.print.pwrss.regression <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt ==  "two.sided") {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)        : \u03B2 - \u03B2\u2080 = 0 \n")
      cat("  H\u2081 (Alternative) : \u03B2 - \u03B2\u2080 \u2260 0 \n\n")
    } else {
      cat("  H\u2080 (Null)        : \u03B2 - \u03B2\u2080 = \u03B4 \n")
      cat("  H\u2081 (Alternative) : \u03B2 - \u03B2\u2080 \u2260 \u03B4 \n\n")
    }
  } else if (x$alt == "one.sided" && x$ncp.alternative > x$ncp.null) {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)        : \u03B2 - \u03B2\u2080 \u2264 0 \n")
      cat("  H\u2081 (Alternative) : \u03B2 - \u03B2\u2080 > 0 \n\n")
    } else {
      cat("  H\u2080 (Null)        : \u03B2 - \u03B2\u2080 \u2264 \u03B4 \n")
      cat("  H\u2081 (Alternative) : \u03B2 - \u03B2\u2080 > \u03B4 \n\n")
    }
  } else if (x$alt == "one.sided" && x$ncp.alternative < x$ncp.null) {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)        : \u03B2 - \u03B2\u2080 \u2265 0 \n")
      cat("  H\u2081 (Alternative) : \u03B2 - \u03B2\u2080 < 0 \n\n")
    } else {
      cat("  H\u2080 (Null)        : \u03B2 - \u03B2\u2080 \u2265 \u03B4 \n")
      cat("  H\u2081 (Alternative) : \u03B2 - \u03B2\u2080 < \u03B4 \n\n")
    }
  } else if (x$alt == "two.one.sided" && (x$ncp.alternative > min(x$ncp.null) && x$ncp.alternative < max(x$ncp.null))) {
    cat("  H\u2080 (Null)        : \u03B2 - \u03B2\u2080 \u2264 min(\u03B4) \u222A \u03B2 - \u03B2\u2080 \u2265 max(\u03B4) \n")
    cat("  H\u2081 (Alternative) : min(\u03B4) < \u03B2 - \u03B2\u2080 < max(\u03B4) \n\n")
  } else if (x$alt == "two.one.sided" && (x$ncp.alternative < min(x$ncp.null) || x$ncp.alternative > max(x$ncp.null))) {
    cat("  H\u2080 (Null)        : min(\u03B4) \u2264 \u03B2 - \u03B2\u2080 \u2264 max(\u03B4) \n")
    cat("  H\u2081 (Alternative) : \u03B2 - \u03B2\u2080 < min(\u03B4) \u222A \u03B2 - \u03B2\u2080 > max(\u03B4) \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  Std. \u03B2            = %.*f \n", digits, x$std.beta))
    cat(sprintf("  Std. \u03B2\u2080           \u2009= %.*f \n", digits, x$std.null.beta))
    cat(sprintf("  Std. \u03B4            = %s \n", paste(round(x$std.margin, digits),  collapse = " and ")))
    cat(sprintf("  df                = %.*f \n", 0, x$df))
    cat(sprintf("  \u03BB (Alternative)   = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  \u03BB\u2080 (Null)         \u2009= %s \n", paste(round(x$ncp.null, digits),  collapse = " and ")))
    cat(sprintf("  T\u207B\u00B9(\u03B1, \u03BB\u2080)        \u2009\u2009= %s \n\n", paste(round(x$t.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0(" \033[34m Sample Size       = ",  round(x$n, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  round(x$n, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error      = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03B2      : Regression coefficient under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03B2\u2080     \u2009: Regression coefficient under null \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03B4      : Margin(s) - ignorable \u03B2 - \u03B2\u2080 difference \n\n", "\033[0m", sep = "")
    cat("\033[36m", "  Std. \u03B2 = \u03B2 * [\u03C3(X) / \u03C3(Y)] \n", "\033[0m", sep = "")
    cat("\033[36m", "  Std. \u03B2\u2080\u2009= \u03B2\u2080 * [\u03C3(X) / \u03C3(Y)] \n", "\033[0m", sep = "")
    cat("\033[36m", "  Std. \u03B4 = \u03B4 * [\u03C3(X) / \u03C3(Y)] \n\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3(X)   : Standard devition of the predictor \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3(Y)   : Standard devition of the outcome \n\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB      : Non-centrality parameter under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080     \u2009: Non-centrality parameter under null \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.regression()


.print.pwrss.f.regression <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$k.tested < x$k.total) {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)        : \u0394R\u00B2 = 0 \n")
      cat("  H\u2081 (Alternative) : \u0394R\u00B2 > 0 \n\n")
    } else {
      cat("  H\u2080 (Null)        : 0 \u2264 \u0394R\u00B2 \u2264 \u03B4 \n")
      cat("  H\u2081 (Alternative) : \u0394R\u00B2 > \u03B4 \n\n")
    }
  } else {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)        : R\u00B2 = 0 \n")
      cat("  H\u2081 (Alternative) : R\u00B2 > 0 \n\n")
    } else {
      cat("  H\u2080 (Null)        : 0 \u2264 R\u00B2 \u2264 \u03B4 \n")
      cat("  H\u2081 (Alternative) : R\u00B2 > \u03B4 \n\n")
    }

  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    if (x$k.tested < x$k.total) {
      cat(sprintf("  \u0394R\u00B2               = %.*f\n", digits, x$r.squared.change))
    } else {
      cat(sprintf("  R\u00B2                = %.*f\n", digits, x$r.squared.change))
    }
    cat(sprintf("  \u03B4                 = %.*f\n", digits, x$margin))
    cat(sprintf("  df1               = %.*f\n", 0, x$df1))
    cat(sprintf("  df2               = %.*f\n", 0, x$df2))
    cat(sprintf("  \u03BB                 = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  \u03BB\u2080                \u2009= %.*f\n", digits, x$ncp.null))
    cat(sprintf("  F\u207B\u00B9(\u03B1, \u03BB\u2080)        \u2009\u2009= %s \n\n", paste(round(x$f.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0(" \033[34m Sample Size       = ",  round(x$n, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  round(x$n, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03B4 \u2009\u2009: Margin - ignorable R\u00B2 or \u0394R\u00B2 \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB \u2009\u2009: Non-centrality parameter under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080 : Non-centrality parameter under null \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.f.regression()


.print.pwrss.student <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)         : d - d\u2080 = 0 \n")
      cat("  H\u2081 (Alternative)  : d - d\u2080 \u2260 0 \n\n")
    } else {
      cat("  H\u2080 (Null)         : d - d\u2080 = \u03B4 \n")
      cat("  H\u2081 (Alternative)  : d - d\u2080 \u2260 \u03B4 \n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$ncp.alternative < x$ncp.null) {
      if (x$margin == 0) {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2265 0 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 < 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2265 \u03B4 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 < \u03B4 \n\n")
      }
    } else {
      if (x$margin == 0) {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2264 0 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 > 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2264 \u03B4 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 > \u03B4 \n\n")
      }
    }
  } else {
    if (x$ncp.alternative > min(x$ncp.null) && x$ncp.alternative < max(x$ncp.null)) {
      cat("  H\u2080 (Null)         : d - d\u2080 \u2264 min(\u03B4) \u222A d - d\u2080 \u2265 max(\u03B4) \n")
      cat("  H\u2081 (Alternative)  : min(\u03B4) < d - d\u2080 < max(\u03B4) \n\n")
    } else {
      cat("  H\u2080 (Null)         : min(\u03B4) \u2264 d - d\u2080 \u2264 max(\u03B4) \n")
      cat("  H\u2081 (Alternative)  : d - d\u2080 < min(\u03B4) \u222A d - d\u2080 > max(\u03B4) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  d                 = %.*f\n", digits, x$d))
    cat(sprintf("  d\u2080                \u2009= %.*f\n", digits, x$null.d))
    cat(sprintf("  \u03B4                 = %s \n", paste(round(x$margin, digits),  collapse = " and ")))
    cat(sprintf("  df                = %.*f\n", 0, x$df))
    cat(sprintf("  \u03BB                 = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  \u03BB\u2080                \u2009= %s \n", paste(round(x$ncp.null, digits),  collapse = " and ")))
    cat(sprintf("  T\u207B\u00B9(\u03B1, \u03BB\u2080)        \u2009\u2009= %s \n\n", paste(round(x$t.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("\033[34m  Sample Size       = ", n.text, "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  d \u2009\u2009: Cohen's d under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  d\u2080 : Cohen's d under null \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03B4 \u2009\u2009: Margin - ignorable d - d\u2080 difference \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB \u2009\u2009: Non-centrality parameter under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080 : Non-centrality parameter under null \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.student()


.print.pwrss.wilcoxon <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
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
  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Method           \u2009\u2009: ", method, "\n", sep = "")
  cat("  Distribution     \u2009\u2009: ", dist, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)         : d - d\u2080 = 0 \n")
      cat("  H\u2081 (Alternative)  : d - d\u2080 \u2260 0 \n\n")
    } else {
      cat("  H\u2080 (Null)         : d - d\u2080 = \u03B4 \n")
      cat("  H\u2081 (Alternative)  : d - d\u2080 \u2260 \u03B4 \n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$method == "guenther") {
      is.less <- x$ncp < x$null.ncp
    } else {
      is.less <- x$mean < x$null.mean
    }
    if (is.less) {
      if (x$margin == 0) {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2265 0 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 < 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2265 \u03B4 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 < \u03B4 \n\n")
      }
    } else {
      if (x$margin == 0) {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2264 0 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 > 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : d - d\u2080 \u2264 \u03B4 \n")
        cat("  H\u2081 (Alternative)  : d - d\u2080 > \u03B4 \n\n")
      }
    }
  } else {
    if (x$method == "guenther") {
      is.equivalent <- x$ncp > min(x$null.ncp) && x$ncp < max(x$null.ncp)
    } else {
      is.equivalent <- x$mean > min(x$null.mean) && x$mean < max(x$null.mean)
    }
    if (is.equivalent) {
      cat("  H\u2080 (Null)         : d - d\u2080 \u2264 min(\u03B4) \u222A d - d\u2080 \u2265 max(\u03B4) \n")
      cat("  H\u2081 (Alternative)  : min(\u03B4) < d - d\u2080 < max(\u03B4) \n\n")
    } else {
      cat("  H\u2080 (Null)         : min(\u03B4) \u2264 d - d\u2080 \u2264 max(\u03B4) \n")
      cat("  H\u2081 (Alternative)  : d - d\u2080 < min(\u03B4) \u222A d - d\u2080 > max(\u03B4) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  d                 = %.*f\n", digits, x$d))
    cat(sprintf("  d\u2080                \u2009= %.*f\n", digits, x$null.d))
    cat(sprintf("  \u03B4                 = %s \n", paste(round(x$margin, digits),  collapse = " and ")))
    if (x$method == "guenther") {
      cat(sprintf("  \u03BB                 = %.*f\n", digits, x$ncp))
      cat(sprintf("  df                = %.*f\n", 0, x$df))
      cat(sprintf("  \u03BB\u2080                \u2009= %s \n", paste(round(x$null.ncp, digits),  collapse = " and ")))
      cat(sprintf("  T\u207B\u00B9(\u03B1, \u03BB\u2080)        \u2009\u2009= %s \n\n", paste(round(x$t.alpha, digits),  collapse = " and ")))
    } else {
      cat(sprintf("  \u03BC                 = %.*f\n", digits, x$mean))
      cat(sprintf("  \u03BC\u2080                \u2009= %s \n", paste(round(x$null.mean, digits),  collapse = " and ")))
      cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080)        \u2009\u2009= %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("\033[34m  Sample Size       = ", n.text, "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  d \u2009\u2009: Cohen's d under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  d\u2080 : Cohen's d under null \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03B4 \u2009\u2009: Margin - ignorable d - d\u2080 difference \n", "\033[0m", sep = "")
    if (x$method == "guenther") {
      cat("\033[36m", "  \u03BB \u2009\u2009: Non-centrality parameter under alternative \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03BB\u2080 : Non-centrality parameter under null \n\n", "\033[0m", sep = "")
    } else {
      cat("\033[36m", "  \u03BC \u2009\u2009: Mean of the alternative distribution \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03BC\u2080 : Mean of the null distribution. \n\n", "\033[0m", sep = "")
    }
  }

} # .print.pwrss.wilcoxon()


.print.pwrss.med <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `sobel` = "Sobel",
                   `aroian` = "Aroian",
                   `goodman` = "Goodman",
                   `joint` = "Joint",
                   `monte.carlo` = "Monte Carlo")
  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Method            : ", method, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    cat("  H\u2080 (Null)         : \u03B2[a*b] = 0 \n")
    cat("  H\u2081 (Alternative)  : \u03B2[a*b] \u2260 0 \n\n")
  } else if (x$alt == "one.sided") {
    if (x$std.beta.indirect < 0) {
      cat("  H\u2080 (Null)         : \u03B2[a*b] \u2265 0 \n")
      cat("  H\u2081 (Alternative)  : \u03B2[a*b] < 0 \n\n")
    } else {
      cat("  H\u2080 (Null)         : \u03B2[a*b] \u2264 0 \n")
      cat("  H\u2081 (Alternative)  : \u03B2[a*b] > 0 \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  Std. \u03B2[a]         = %.*f \n", digits, x$std.beta.a))
    cat(sprintf("  Std. \u03B2[b]         = %.*f \n", digits, x$std.beta.b))
    if (x$method %in% c("sobel", "aorian", "goodman")) {
      cat(sprintf("  Std. \u03B2[a*b]       = %.*f \n", digits, x$std.beta.indirect))
      cat(sprintf("  \u03BC (Alternative)   = %.*f\n", digits, x$mean.alternative))
      cat(sprintf("  \u03BC\u2080 (Null)         \u2009= %s \n", paste(round(x$mean.null, digits),  collapse = " and ")))
      cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080)        \u2009\u2009= %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
    } else {
      cat(sprintf("  Std. \u03B2[a*b]       = %.*f \n\n", digits, x$std.beta.indirect))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0(" \033[34m Sample Size       = ",  round(x$n, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  round(x$n, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error      = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03B2[a]         : Regression coefficient for path `a` \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03B2[b]         : Regression coefficient for path `b` \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03B2[a*b]       : Coefficient for the indirect path `a*b` \n\n", "\033[0m", sep = "")
    cat("\033[36m", "  Std. \u03B2[a]    = \u03B2[a] * [\u03C3(predictor) / \u03C3(mediator)] \n", "\033[0m", sep = "")
    cat("\033[36m", "  Std. \u03B2[b]    = \u03B2[b] * [\u03C3(mediator) / \u03C3(outcome)] \n", "\033[0m", sep = "")
    cat("\033[36m", "  Std. \u03B2[a*b]  = Std. \u03B2[a] * Std. \u03B2[b] \n\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3(predictor) : Standard devition of the predictor \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3(mediator)  : Standard devition of the mediator \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3(outcome)   : Standard devition of the outcome \n\n", "\033[0m", sep = "")
    if (x$method %in% c("sobel", "aorian", "goodman")) {
      cat("\033[36m", "  \u03BC            : Mean of the alternative distribution \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03BC\u2080           \u2009: Mean of the null distribution \n\n", "\033[0m", sep = "")
    }
  }

} # .print.pwrss.med()


.print.pwrss.gof <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  cat("  H\u2080 (Null)         : P[i,j] = P\u2080[i,j] for \u2200(i,j) \n")
  cat("  H\u2081 (Alternative)  : P[i,j] \u2260 P\u2080[i,j] for \u2203(i,j)\n\n")

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  df                = %.*f\n", 0, x$df))
    cat(sprintf("  \u03BB                 = %.*f \n", digits, x$ncp.alternative))
    cat(sprintf("  \u03BB\u2080                \u2009= %.*f \n\n", digits, x$ncp.null))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n") {
    cat(paste0(" \033[34m Total Sample Size = ",  round(x$n, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Total Sample Size = ",  round(x$n, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  Independence implies P\u2080[i,j] = P[i,.] * P[.,j] \n", "\033[0m", sep = "")
    cat("\033[36m", "  P[i,j] : Joint prob. for cell (i,j) \n", "\033[0m", sep = "")
    cat("\033[36m", "  P[i,.] : Marginal prob. for row i (sum over j) \n", "\033[0m", sep = "")
    cat("\033[36m", "  P[.,j] : Marginal prob. for column j (sum over i) \n", "\033[0m", sep = "")
    cat("\033[36m", "  For goodness-of-fit test, it is P[i] vs P\u2080[i] \n\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB      : Non-centrality parameter under alternative\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080     \u2009: Non-centrality parameter under null\n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.gof()


.print.pwrss.chisq <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$ncp.null > 0) {
    cat("  H\u2080 (Null)         : 0 \u2264 \u03BB \u2264 \u03BB\u2080 \n")
    cat("  H\u2081 (Alternative)  : \u03BB > \u03BB\u2080 \n\n")
  } else {
    cat("  H\u2080 (Null)         : \u03BB = \u03BB\u2080 \n")
    cat("  H\u2081 (Alternative)  : \u03BB > \u03BB\u2080 \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  df                = %.*f\n", 0, x$df))
    cat(sprintf("  \u03BB                 = %.*f \n", digits, x$ncp.alternative))
    cat(sprintf("  \u03BB\u2080                \u2009= %.*f \n", digits, x$ncp.null))
    cat(sprintf("  Inv-\u03C7\u00B2(\u03B1, \u03BB\u2080)     \u2009= %s \n\n", paste(round(x$chisq.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB \u2009\u2009: Non-centrality parameter under alternative\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080 : Non-centrality parameter under null\n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.chisq()




.print.pwrss.t <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    cat("  H\u2080 (Null)         : \u03BB = \u03BB\u2080 \n")
    cat("  H\u2081 (Alternative)  : \u03BB \u2260 \u03BB\u2080 \n\n")
  } else if (x$alt == "one.sided") {
    if (x$ncp.alternative < x$ncp.null) {
      cat("  H\u2080 (Null)         : \u03BB \u2265 \u03BB\u2080 \n")
      cat("  H\u2081 (Alternative)  : \u03BB < \u03BB\u2080 \n\n")
    } else {
      cat("  H\u2080 (Null)         : \u03BB \u2264 \u03BB\u2080 \n")
      cat("  H\u2081 (Alternative)  : \u03BB > \u03BB\u2080 \n\n")
    }
  } else {
    if (x$ncp.alternative > min(x$ncp.null) && x$ncp.alternative < max(x$ncp.null)) {
      cat("  H\u2080 (Null)         : \u03BB \u2264 min(\u03BB\u2080) \u222A \u03BB \u2265 max(\u03BB\u2080) \n")
      cat("  H\u2081 (Alternative)  : min(\u03BB\u2080) < \u03BB < max(\u03BB\u2080) \n\n")
    } else {
      cat("  H\u2080 (Null)         : min(\u03BB\u2080) \u2264 \u03BB \u2264 max(\u03BB\u2080) \n")
      cat("  H\u2081 (Alternative)  : \u03BB < min(\u03BB\u2080) \u222A \u03BB > max(\u03BB\u2080) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    if (isFALSE(is.infinite(x$df))) {
      cat(sprintf("  df                = %.*f \n", 0, x$df))
    }
    cat(sprintf("  \u03BB                 = %.*f \n", digits, x$ncp.alternative))
    cat(sprintf("  \u03BB\u2080                \u2009= %s \n", paste(round(x$ncp.null, digits),  collapse = " and ")))
    cat(sprintf("  T\u207B\u00B9(\u03B1, \u03BB\u2080)        \u2009\u2009= %s \n\n", paste(round(x$t.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results \n")
  cat(line)
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB \u2009\u2009: Non-centrality parameter under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080 : Non-centrality parameter under null \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.t()


.print.pwrss.z <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    cat("  H\u2080 (Null)         : \u03BC = \u03BC\u2080 \n")
    cat("  H\u2081 (Alternative)  : \u03BC \u2260 \u03BC\u2080 \n\n")
  } else if (x$alt == "one.sided") {
    if (x$mean.alternative < x$mean.null) {
      cat("  H\u2080 (Null)         : \u03BC \u2265 \u03BC\u2080 \n")
      cat("  H\u2081 (Alternative)  : \u03BC < \u03BC\u2080 \n\n")
    } else {
      cat("  H\u2080 (Null)         : \u03BC \u2264 \u03BC\u2080 \n")
      cat("  H\u2081 (Alternative)  : \u03BC > \u03BC\u2080 \n\n")
    }
  } else {
    if (x$mean.alternative > min(x$mean.null) && x$mean.alternative < max(x$mean.null)) {
      cat("  H\u2080 (Null)         : \u03BC \u2264 min(\u03BC\u2080) \u222A \u03BC \u2265 max(\u03BC\u2080) \n")
      cat("  H\u2081 (Alternative)  : min(\u03BC\u2080) < \u03BC < max(\u03BC\u2080) \n\n")
    } else {
      cat("  H\u2080 (Null)         : min(\u03BC\u2080) \u2264 \u03BC \u2264 max(\u03BC\u2080) \n")
      cat("  H\u2081 (Alternative)  : \u03BC < min(\u03BC\u2080) \u222A \u03BC > max(\u03BC\u2080) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  \u03BC                 = %.*f \n", digits, x$mean.alternative))
    cat(sprintf("  \u03BC\u2080                \u2009= %s \n", paste(round(x$mean.null, digits),  collapse = " and ")))
    cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080)        \u2009\u2009= %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03BC \u2009\u2009: Mean under alternative\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BC\u2080 : Mean under null\n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.z()


.print.pwrss.binom <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    cat("  H\u2080 (Null)         : P = P\u2080 \n")
    cat("  H\u2081 (Alternative)  : P \u2260 P\u2080 \n\n")
  } else if (x$alt == "one.sided") {
    if (x$prob.alternative < x$prob.null) {
      cat("  H\u2080 (Null)         : P \u2265 P\u2080 \n")
      cat("  H\u2081 (Alternative)  : P < P\u2080 \n\n")
    } else {
      cat("  H\u2080 (Null)         : P \u2264 P\u2080 \n")
      cat("  H\u2081 (Alternative)  : P > P\u2080 \n\n")
    }
  } else {
    if (x$prob.alternative > min(x$prob.null) && x$prob.alternative < max(x$prob.null)) {
      cat("  H\u2080 (Null)         : P \u2264 min(P\u2080) \u222A P \u2265 max(P\u2080) \n")
      cat("  H\u2081 (Alternative)  : min(P\u2080) < P < max(P\u2080) \n\n")
    } else {
      cat("  H\u2080 (Null)         : min(P\u2080) \u2264 P \u2264 max(P\u2080) \n")
      cat("  H\u2081 (Alternative)  : P < min(P\u2080) \u222A P > max(P\u2080) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  Size              = %.*f \n", 0, x$size))
    cat(sprintf("  P                 = %.*f \n", digits, x$prob.alternative))
    cat(sprintf("  P\u2080                \u2009= %s \n", paste(round(x$prob.null, digits),  collapse = " and ")))
    cat(sprintf("  Bin\u207B\u00B9(\u03B1, P\u2080)      \u2009\u2009= %s \n\n", paste(round(x$binom.alpha, 0),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  P \u2009\u2009: Probability under alternative\n", "\033[0m", sep = "")
    cat("\033[36m", "  P\u2080 : Probability under null\n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.binom()


.print.pwrss.f <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$ncp.null > 0) {
    cat("  H\u2080 (Null)         : 0 \u2264 \u03BB \u2264 \u03BB\u2080 \n")
    cat("  H\u2081 (Alternative)  : \u03BB > \u03BB\u2080 \n\n")
  } else {
    cat("  H\u2080 (Null)         : \u03BB = \u03BB\u2080 \n")
    cat("  H\u2081 (Alternative)  : \u03BB > \u03BB\u2080 \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  df1               = %.*f\n", 0, x$df1))
    cat(sprintf("  df2               = %.*f\n", 0, x$df2))
    cat(sprintf("  \u03BB                 = %.*f\n", digits, x$ncp.alternative))
    cat(sprintf("  \u03BB\u2080                \u2009= %.*f\n", digits, x$ncp.null))
    cat(sprintf("  F\u207B\u00B9(\u03B1, \u03BB\u2080)        \u2009\u2009= %s \n\n", paste(round(x$f.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results \n")
  cat(line)
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB \u2009\u2009: Non-centrality parameter under alternative\n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080 : Non-centrality parameter under null\n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.f()

.print.pwrss.ancova <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n.total",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")
  cat(LL, HR, LR, "\n\n", sep = "")
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$null.ncp == 0) {
    cat("  H\u2080 (Null)             : \u03B7\u00B2 = 0 \n")
    cat("  H\u2081 (Alternative)      : \u03B7\u00B2 > 0 \n\n")
  } else {
    cat("  H\u2080 (Null)             : 0 \u2264 \u03B7\u00B2 \u2264 \u03B7\u2080\u00B2 \n")
    cat("  H\u2081 (Alternative)      : \u03B7\u00B2 > \u03B7\u2080\u00B2 \n\n")
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  Design                = %s \n", x$effect))
    cat(sprintf("  df1                   = %.*f\n", 0, x$df1))
    cat(sprintf("  df2                   = %.*f\n", 0, x$df2))
    cat(sprintf("  \u03BB                     = %.*f\n", digits, x$ncp))
    cat(sprintf("  \u03BB\u2080                    \u2009= %.*f\n", digits, x$null.ncp))
    cat(sprintf("  F\u207B\u00B9(\u03B1, \u03BB\u2080)            \u2009\u2009= %s \n\n", paste(round(x$f.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n.total") {
    cat(paste0("\033[34m  Total Sample Size     = ",  round(x$n.total, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Total Sample Size     = ",  round(x$n.total, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)      = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)      = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power     = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power     = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03B7\u00B2 \u2009\u2009: (Partial) Eta-squared under alternative \n", "\033[0m", sep = "")
    if (x$null.ncp != 0) cat("\033[36m", "  \u03B7\u2080\u00B2 : (Partial) Eta-squared under null \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB  \u2009\u2009: Non-centrality parameter under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080  : Non-centrality parameter under null \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.ancova()


.print.pwrss.contrast <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n.total",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  cat("  H\u2080 (Null)             : \u03C8 = 0 \n")
  cat("  H\u2081 (Alternative)      : \u03C8 \u2260 0 \n\n")

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  \u03C8                     = %.*f\n", digits, x$psi))
    cat(sprintf("  d                     = %.*f\n", digits, x$d))
    cat(sprintf("  df                    = %.*f\n", 0, x$df))
    cat(sprintf("  \u03BB                     = %.*f\n", digits, x$ncp))
    cat(sprintf("  \u03BB\u2080                    \u2009= %.*f\n", digits, x$null.ncp))
    cat(sprintf("  T\u207B\u00B9(\u03B1, \u03BB\u2080)            \u2009\u2009= %s \n\n", paste(round(x$t.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (x$requested == "n.total") {
    cat(paste0("\033[34m  Total Sample Size     = ",  round(x$n.total, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Total Sample Size     = ",  round(x$n.total, digits = digits)), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)      = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)      = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power     = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power     = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03C8  \u2009\u2009: Contrast est. defined as \u2211(contrast\u1D62 * \u03BC\u1D62) \n", "\033[0m", sep = "")
    cat("\033[36m", "  d  \u2009\u2009: Standardized contrast estimate \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB  \u2009\u2009: Non-centrality parameter under alternative \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03BB\u2080  : Non-centrality parameter under null \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.contrast()


.print.pwrss.contrasts <- function(x, data = NULL, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n.total",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  cat("  H\u2080 (Null)                 : \u03C8 = 0 \n")
  cat("  H\u2081 (Alternative)          : \u03C8 \u2260 0 \n\n")

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
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
    cat(sprintf("  \u03BB\u2080                    \u2009= %.*f\n", digits, x$null.ncp))
    cat(sprintf("  Adjusted Type 1 Error (\u03B1) = %.*f \n", digits, x$alpha))
    cat(sprintf("  \u03B1 Adjustment              = %s \n\n", adjust.alpha))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  if (!is.null(x$data)) {
    print(x$data, row.names = FALSE)
    cat("\n")
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  \u03C8 (psi) \u2009\u2009: Contrast est. defined as \u2211(contrast\u1D62 * \u03BC\u1D62) \n", "\033[0m", sep = "")
    cat("\033[36m", "  d       \u2009\u2009: Standardized contrast estimate \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.contrasts()


.print.pwrss.fisher <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `z` = "Normal Approximation",
                   `exact` = "Fisher's Exact")
  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Method           \u2009: ", method, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    if (x$margin == 0) {
      cat("  H\u2080 (Null)         : P\u2081 - P\u2082 = 0 \n")
      cat("  H\u2081 (Alternative)  : P\u2081 - P\u2082 \u2260 0 \n\n")
    } else {
      cat("  H\u2080 (Null)         : P\u2081 - P\u2082 = \u03B4 \n")
      cat("  H\u2081 (Alternative)  : P\u2081 - P\u2082 \u2260 \u03B4 \n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      if (x$margin == 0) {
        cat("  H\u2080 (Null)         : P\u2081 - P\u2082 \u2265 0 \n")
        cat("  H\u2081 (Alternative)  : P\u2081 - P\u2082 < 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : P\u2081 - P\u2082 \u2265 \u03B4 \n")
        cat("  H\u2081 (Alternative)  : P\u2081 - P\u2082 < \u03B4 \n\n")
      }
    } else {
      if (x$margin == 0) {
        cat("  H\u2080 (Null)         : P\u2081 - P\u2082 \u2264 0 \n")
        cat("  H\u2081 (Alternative)  : P\u2081 - P\u2082 > 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : P\u2081 - P\u2082 \u2264 \u03B4 \n")
        cat("  H\u2081 (Alternative)  : P\u2081 - P\u2082 > \u03B4 \n\n")
      }
    }
  } else {
    if (x$delta > min(x$margin) && x$delta < max(x$margin)) {
      cat("  H\u2080 (Null)         : P\u2081 - P\u2082 \u2264 min(\u03B4) \u222A P\u2081 - P\u2082 \u2265 max(\u03B4) \n")
      cat("  H\u2081 (Alternative)  : min(\u03B4) < P\u2081 - P\u2082 < max(\u03B4) \n\n")
    } else {
      cat("  H\u2080 (Null)         : min(\u03B4) \u2264 P\u2081 - P\u2082 \u2264 max(\u03B4) \n")
      cat("  H\u2081 (Alternative)  : P\u2081 - P\u2082 < min(\u03B4) \u222A P\u2081 - P\u2082 > max(\u03B4) \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    if (x$method == "z") {
      cat(sprintf("  P\u2081 - P\u2082           \u2009\u2009= %.*f \n", digits, x$delta))
      cat(sprintf("  Odds Ratio (OR)   = %.*f \n", digits, x$odds.ratio))
      cat(sprintf("  \u03BC                 = %.*f \n", digits, x$mean.alternative))
      cat(sprintf("  \u03BC\u2080                \u2009= %s \n", paste(round(x$mean.null, digits),  collapse = " and ")))
      cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080)        \u2009\u2009= %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
    } else {
      cat(sprintf("  P\u2081 - P\u2082           \u2009\u2009= %.*f \n", digits, x$delta))
      cat(sprintf("  Odds Ratio (OR)   = %.*f \n\n", digits, x$odds.ratio))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("\033[34m  Sample Size       = ", n.text, "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  P\u2081       : Probability of success in the first group \n", "\033[0m", sep = "")
    cat("\033[36m", "  P\u2082       : Probability of success in the second group \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03B4       \u2009\u2009: Margin - ignorable P\u2081 - P\u2082 difference \n", "\033[0m", sep = "")
    cat("\033[36m", "  OR      \u2009\u2009: Odds(P\u2081) / Odds(P\u2082) \n", "\033[0m", sep = "")
    cat("\033[36m", "  Odds(P\u2081) : P\u2081 / (1 - P\u2081) \n", "\033[0m", sep = "")
    if (x$method == "z") {
      cat("\033[36m", "  Odds(P\u2082) : P\u2082 / (1 - P\u2082) \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03BC       \u2009\u2009: Mean of the alternative distribution \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03BC\u2080       : Mean of the null distribution \n\n", "\033[0m", sep = "")
    } else {
      cat("\033[36m", "  Odds(P\u2082) : P\u2082 / (1 - P\u2082) \n\n", "\033[0m", sep = "")
    }
  }

} # .print.pwrss.fisher()


.print.pwrss.mcnemar <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `z` = "Normal Approximation",
                   `exact` = "McNemar's Exact")
  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Method            \u2009\u2009: ", method, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)

  if (x$alt == "two.sided") {
    cat("  H\u2080 (Null)          : P\u2081\u2080 - P\u2080\u2081 = 0 \n")
    cat("  H\u2081 (Alternative)   : P\u2081\u2080 - P\u2080\u2081 \u2260 0 \n\n")
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      cat("  H\u2080 (Null)          : P\u2081\u2080 - P\u2080\u2081 \u2265 0 \n")
      cat("  H\u2081 (Alternative)   : P\u2081\u2080 - P\u2080\u2081 < 0 \n\n")
    } else {
      cat("  H\u2080 (Null)          : P\u2081\u2080 - P\u2080\u2081 \u2264 0 \n")
      cat("  H\u2081 (Alternative)   : P\u2081\u2080 - P\u2080\u2081 > 0 \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    cat(sprintf("  P\u2081\u2080 - P\u2080\u2081           \u200A= %.*f \n", digits, x$delta))
    cat(sprintf("  Odds Ratio (OR)    = %.*f \n", digits, x$odds.ratio))
    if (x$method == "exact") {
      cat(sprintf("  Size (Discordant)  = %.*f \n", 0, x$size))
      cat(sprintf("  P\u2081                 \u2009= %.*f \n", digits, x$prob.alternative))
      cat(sprintf("  P\u2080                 \u2009= %s \n", paste(round(x$prob.null, digits),  collapse = " and ")))
      cat(sprintf("  Bin\u207B\u00B9(\u03B1, P\u2080)       \u2009\u2009= %s \n\n", paste(round(x$binom.alpha, digits),  collapse = " and ")))
    } else {
      cat(sprintf("  \u03BC\u2081                 \u2009= %.*f \n", digits, x$mean.alternative))
      cat(sprintf("  \u03BC\u2080                 \u2009= %s \n", paste(round(x$mean.null, digits),  collapse = " and ")))
      cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080)         \u2009\u2009= %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n.paired, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("\033[34m  Paired Sample Size = ", n.text, "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Paired Sample Size = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)   = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power  = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power  = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  OR : P\u2081\u2080 / P\u2080\u2081 \n", "\033[0m", sep = "")
    if (x$method == "exact") {
      cat("\033[36m", "  P\u2081 \u2009: Prob. of {1,0} among discordant pairs under alt. \n", "\033[0m", sep = "")
      cat("\033[36m", "  P\u2080 \u2009: Prob. of {1,0} among discordant pairs under null \n\n", "\033[0m", sep = "")
    } else {
      cat("\033[36m", "  \u03BC\u2081 \u2009: Mean of the alternative distribution \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03BC\u2080 \u2009: Mean of the null distribution \n\n", "\033[0m", sep = "")
    }
  }

} # .print.pwrss.mcnemar()


.print.pwrss.oneprop <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  method <- switch(x$method,
                   `z` = "Normal Approximation",
                   `exact` = "Exact")

  # Header
  cat(x$test, "\n\n", sep = "")

  if (x$method == "exact") {
    cat("  Method            \u2009\u2009: ", method, "\n\n", sep = "")
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
    cat("  H\u2080 (Null)          : P - P\u2080 = 0 \n")
    cat("  H\u2081 (Alternative)   : P - P\u2080 \u2260 0 \n\n")
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      cat("  H\u2080 (Null)          : P - P\u2080 \u2265 0 \n")
      cat("  H\u2081 (Alternative)   : P - P\u2080 < 0 \n\n")
    } else {
      cat("  H\u2080 (Null)          : P - P\u2080 \u2264 0 \n")
      cat("  H\u2081 (Alternative)   : P - P\u2080 > 0 \n\n")
    }
  } else if (x$alt == "two.one.sided") {
    if (x$delta[1] > 0 && x$delta[2] < 0) {
      cat("  H\u2080 (Null)          : P - min(P\u2080) \u2264 0 \u222A P - max(P\u2080) \u2265 0 \n")
      cat("  H\u2081 (Alternative)   : P - min(P\u2080) > 0 \u2229 P - max(P\u2080) < 0 \n\n")
    } else {
      cat("  H\u2080 (Null)          : P - min(P\u2080) \u2265 0 \u2229 P - max(P\u2080) \u2264 0 \n")
      cat("  H\u2081 (Alternative)   : P - min(P\u2080) < 0 \u222A P - max(P\u2080) > 0 \n\n")
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    delta.text <- paste(round(x$delta, digits), collapse = " and ")
    cat(sprintf("  P - P\u2080             \u2009= %s \n", delta.text))
    or.text <- paste(round(x$odds.ratio, digits), collapse = " and ")
    cat(sprintf("  Odds Ratio (OR)    = %s \n", or.text))
    if (x$method == "exact") {
      cat(sprintf("  Size               = %.*f \n", 0, x$size))
      cat(sprintf("  P                  = %.*f \n", digits, x$prob.alternative))
      cat(sprintf("  P\u2080                 \u2009= %s \n", paste(round(x$prob.null, digits),  collapse = " and ")))
      cat(sprintf("  Bin\u207B\u00B9(\u03B1, P\u2080)       \u2009\u2009= %s \n\n", paste(round(x$binom.alpha, digits),  collapse = " and ")))
    } else {
      cat(sprintf("  \u03BC\u2081                 \u2009= %.*f \n", digits, x$mean.alternative))
      cat(sprintf("  \u03BC\u2080                 \u2009= %s \n", paste(round(x$mean.null, digits),  collapse = " and ")))
      cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080)         \u2009\u2009= %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
    }
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("\033[34m  Paired Sample Size = ", n.text, "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size        = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)   = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)   = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power  = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power  = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "  OR       : Odds(P) / Odds(P\u2080) \n", "\033[0m", sep = "")
    cat("\033[36m", "  Odds(P)  : P / (1 - P) \n", "\033[0m", sep = "")
    cat("\033[36m", "  Odds(P\u2080) \u2009: P\u2080 / (1 - P\u2080) \n", "\033[0m", sep = "")
    if (x$method == "exact") {
      cat("\033[36m", "  P        : Probability of success under alternative \n", "\033[0m", sep = "")
      cat("\033[36m", "  P\u2080       \u2009: Probability of success under null \n\n", "\033[0m", sep = "")
    } else {
      cat("\033[36m", "  \u03BC\u2081       \u2009: Mean of the alternative distribution \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03BC\u2080       \u2009: Mean of the null distribution \n\n", "\033[0m", sep = "")
    }
  }

} # .print.pwrss.oneprop()


.print.pwrss.steiger <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")
  cat("  Common Index       \u2009: ", x$common, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    if (x$common) {
      cat("  H\u2080 (Null)         : \u03C1\u2081\u2082 - \u03C1\u2081\u2083 = 0 \n")
      cat("  H\u2081 (Alternative)  : \u03C1\u2081\u2082 - \u03C1\u2081\u2083 \u2260 0 \n\n")
    } else {
      cat("  H\u2080 (Null)         : \u03C1\u2081\u2082 - \u03C1\u2083\u2084 = 0 \n")
      cat("  H\u2081 (Alternative)  : \u03C1\u2081\u2082 - \u03C1\u2083\u2084 \u2260 0 \n\n")
    }

  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      if (x$common) {
        cat("  H\u2080 (Null)         : \u03C1\u2081\u2082 - \u03C1\u2081\u2083 \u2265 0 \n")
        cat("  H\u2081 (Alternative)  :\u03C1\u2081\u2082 - \u03C1\u2081\u2083 < 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : \u03C1\u2081\u2082 - \u03C1\u2083\u2084 \u2265 0 \n")
        cat("  H\u2081 (Alternative)  : \u03C1\u2081\u2082 - \u03C1\u2083\u2084 < 0 \n\n")
      }

    } else {
      if (x$common) {
        cat("  H\u2080 (Null)         : \u03C1\u2081\u2082 - \u03C1\u2081\u2083 \u2264 0 \n")
        cat("  H\u2081 (Alternative)  : \u03C1\u2081\u2082 - \u03C1\u2081\u2083 > 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : \u03C1\u2081\u2082 - \u03C1\u2083\u2084 \u2264 0 \n")
        cat("  H\u2081 (Alternative)  : \u03C1\u2081\u2082 - \u03C1\u2083\u2084 > 0 \n\n")
      }

    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    if (x$common) {
      cat(sprintf("  \u03C1\u2081\u2082 - \u03C1\u2081\u2083          \u2009= %.*f \n", digits, x$delta))
    } else {
      cat(sprintf("  \u03C1\u2081\u2082 - \u03C1\u2083\u2084          \u2009= %.*f \n", digits, x$delta))
    }
    cat(sprintf("  Cohen's q         = %.*f \n", digits, x$q))
    cat(sprintf("  \u03BC (Alternative)   = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  \u03C3 (Alternative)   = %.*f\n", digits, x$sd.alternative))
    cat(sprintf("  \u03BC\u2080 (Null)         \u2009= %.*f \n", digits, x$mean.null))
    cat(sprintf("  \u03C3\u2080 (Null)         \u2009= %.*f \n", digits, x$sd.null))
    cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)     = %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("\033[34m  Sample Size       = ", n.text, "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    if (x$common) {
      cat("\033[36m", "  \u03C1\u2081\u2082 : Correlation between variable V1 and V2 \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03C1\u2081\u2083 : Correlation between variable V1 and V3 \n", "\033[0m", sep = "")
    } else {
      cat("\033[36m", "  \u03C1\u2081\u2082 : Correlation between variable V1 and V2 \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03C1\u2083\u2084 : Correlation between variable V3 and V4 \n", "\033[0m", sep = "")
    }
    cat("\033[36m", "  \u03BC  \u2009: Mean \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3  \u2009: Standard deviation \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.stegier()


.print.pwrss.twocors <- function(x, digits = 3, verbose = 1, ...) {

  UL <- "\u2554"
  UR <- "\u2557"
  LL <- "\u255A"
  LR <- "\u255D"
  HL <- "\u2550"
  VL <- "\u2551"
  HR <- strrep(HL, 50)
  line <- paste0("\u2500", strrep("\u2500", 50), "\n")
  cat(UL, HR, UR, "\n", sep = "")
  cat(VL, ifelse(x$requested == "n",
                 "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 "               \033[34m POWER CALCULATION \033[0m                "),
      VL, "\n", sep = "")

  cat(LL, HR, LR, "\n\n", sep = "")

  # Header
  cat(x$test, "\n\n", sep = "")

  cat(line)
  cat("Hypotheses\n")
  cat(line)
  if (x$alt == "two.sided") {
    if (x$design %in% c("independent", "paired")) {
      cat("  H\u2080 (Null)         : \u03C1\u2081 - \u03C1\u2082 = 0 \n")
      cat("  H\u2081 (Alternative)  : \u03C1\u2081 - \u03C1\u2082 \u2260 0 \n\n")
    } else {
      cat("  H\u2080 (Null)         : \u03C1\u2081 - \u03C1\u2080 = 0 \n")
      cat("  H\u2081 (Alternative)  : \u03C1\u2081 - \u03C1\u2080 \u2260 0 \n\n")
    }
  } else if (x$alt == "one.sided") {
    if (x$delta < 0) {
      if (x$design %in% c("independent", "paired")) {
        cat("  H\u2080 (Null)         : \u03C1\u2081 - \u03C1\u2082 \u2265 0 \n")
        cat("  H\u2081 (Alternative)  : \u03C1\u2081 - \u03C1\u2082 < 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : \u03C1\u2081 - \u03C1\u2080 \u2265 0 \n")
        cat("  H\u2081 (Alternative)  : \u03C1\u2081 - \u03C1\u2080 < 0 \n\n")
      }
    } else {
      if (x$design %in% c("independent", "paired")) {
        cat("  H\u2080 (Null)         : \u03C1\u2081 - \u03C1\u2082 \u2264 0 \n")
        cat("  H\u2081 (Alternative)  : \u03C1\u2081 - \u03C1\u2082 > 0 \n\n")
      } else {
        cat("  H\u2080 (Null)         : \u03C1\u2081 - \u03C1\u2080 \u2264 0 \n")
        cat("  H\u2081 (Alternative)  : \u03C1\u2081 - \u03C1\u2080 > 0 \n\n")
      }
    }
  }

  if (verbose == 2) {
    cat(line)
    cat("Key Parameters \n")
    cat(line)
    if (x$design %in% c("independent", "paired")) {
      cat(sprintf("  \u03C1\u2081 - \u03C1\u2082           \u2009\u2009= %.*f \n", digits, x$delta))
    } else {
      cat(sprintf("  \u03C1\u2081 - \u03C1\u2080           \u2009\u2009= %.*f \n", digits, x$delta))
    }
    cat(sprintf("  Cohen's q         = %.*f \n", digits, x$q))
    cat(sprintf("  \u03BC (Alternative)   = %.*f\n", digits, x$mean.alternative))
    cat(sprintf("  \u03C3 (Alternative)   = %.*f\n", digits, x$sd.alternative))
    cat(sprintf("  \u03BC\u2080 (Null)         \u2009= %.*f \n", digits, x$mean.null))
    cat(sprintf("  \u03C3\u2080 (Null)         \u2009= %.*f \n", digits, x$sd.null))
    cat(sprintf("  Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)     = %s \n\n", paste(round(x$z.alpha, digits),  collapse = " and ")))
  }

  cat(line)
  cat("Results\n")
  cat(line)
  n.text <- paste(round(x$n, digits), collapse = " and ")
  if (x$requested == "n") {
    cat(paste0("\033[34m  Sample Size       = ", n.text, "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n", "\033[0m"))
  } else {
    cat(paste0("  Sample Size       = ",  n.text), "\n")
  }
  cat(sprintf("  Type 1 Error (\u03B1)  = %.*f\n", digits, x$alpha))
  cat(sprintf("  Type 2 Error (\u03B2)  = %.*f\n", digits, 1 - x$power))
  if (x$requested == "power") {
    cat(paste0(" \033[34m Statistical Power = ", round(x$power, digits = digits), "\033[0m", "\033[1;35m", "  \U00025C4\U00025C4 \n\n", "\033[0m"))
  } else {
    cat(paste0("  Statistical Power = ", round(x$power, digits = digits), "\n\n"))
  }

  if (verbose == 2) {
    cat("\033[36m", line, "\033[0m", sep = "")
    cat("\033[36m", "Definitions\n", "\033[0m", sep = "")
    cat("\033[36m", line, "\033[0m", sep = "")
    if (x$design %in% c("independent", "paired")) {
      cat("\033[36m", "  \u03C1\u2081 : Correlation (for some V1 ~ V2) in the first group \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03C1\u2082 : Correlation (for some V1 ~ V2) in the second group \n", "\033[0m", sep = "")
    } else {
      cat("\033[36m", "  \u03C1\u2081 : Correlation (for some V1 ~ V2) under alternative \n", "\033[0m", sep = "")
      cat("\033[36m", "  \u03C1\u2080 : Correlation (for some V1 ~ V2) under null \n", "\033[0m", sep = "")
    }
    cat("\033[36m", "  \u03BC \u2009\u2009: Mean \n", "\033[0m", sep = "")
    cat("\033[36m", "  \u03C3 \u2009\u2009: Standard deviation \n\n", "\033[0m", sep = "")
  }

} # .print.pwrss.twocors()
