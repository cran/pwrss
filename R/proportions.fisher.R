# type 1 and type 2 error plots are not available for this function
# n.ratio = n1 / n2
power.exact.fisher <- function(prob1, prob2,
                               n2 = NULL, n.ratio = 1,
                               alpha = 0.05, power = NULL,
                               alternative = c("two.sided", "one.sided"),
                               method = c("exact", "approximate"),
                               ceiling = TRUE, verbose = TRUE, pretty = FALSE) {

  check.positive(n.ratio)
  check.proportion(prob1, prob2, alpha)
  check.logical(ceiling)
  if (!is.null(power)) check.proportion(power)
  if (!is.null(n2)) check.sample.size(n2)

  alternative <- tolower(match.arg(alternative))
  method <- tolower(match.arg(method))

  ifelse(is.null(power),
         requested <- "power",
         requested <- "n")

  pwr.approx <- function(prob1, prob2, n2, n.ratio,
                         alpha, alternative,
                         # correct.continuity = FALSE,
                         pooled.stderr = FALSE) {

    # Gpower
    # sigma0 <- sqrt(((n1 * (1 - prob1) + n2 * (1 - prob2)) / (n1 * n2)) * ((n1 * prob1 + n2 *prob2) / (n1 + n2)))
    # stderr <- (1 / sigma0) * sqrt((prob1 * (1 - prob1)) / n1 + (prob2 * (1 - prob2)) / n2)
    # delta <- (1 / sigma0) * (prob1 - prob2 - (k / 2) * (1 / n1 + 1 / n2)) # w/ continuity correction

    n1 <- n.ratio * n2

    delta <- prob1 - prob2

    if (pooled.stderr) {

      p.bar <- (n1 * prob1 + n2 * prob2) / (n1 + n2)
      stderr <- sqrt(p.bar * (1 - p.bar) * (1 / n1 + 1 / n2))

    } else {

      stderr <- sqrt((prob1 * (1 - prob1)) / n1 + (prob2 * (1 - prob2)) / n2)

    }

    # if (correct.continuity) {
    #
    #   ifelse(prob1 < prob2, k <- -1, k <- 1)
    #   if (alternative %in% c("not equal", "two.sided")) k <- c(-1, 1)
    #   delta <- (prob1 - prob2 - (k / 2) * (1 / n1 + 1 / n2))
    #
    # } else {
    #
    #   delta <- prob1 - prob2
    #
    # }

    pwr.obj <- power.z.test(mean = delta / stderr, sd = 1, null.mean = 0, null.sd = 1,
                            alpha = alpha, alternative = alternative,
                            plot = FALSE, verbose = FALSE)

    power <- pwr.obj$power
    mean.alternative <- pwr.obj$mean.alternative
    sd.alternative <- pwr.obj$sd.alternative
    mean.null <- pwr.obj$mean.null
    sd.null <- pwr.obj$sd.null
    z.alpha <- pwr.obj$z.alpha

    return(list(power = power, mean.alternative =  mean.alternative,
                sd.alternative = sd.alternative, mean.null = mean.null,
                sd.null = sd.null, z.alpha = z.alpha))

  } # pwr.approx


  ss.approx <- function(prob1, prob2, power, n.ratio,
                        alpha, alternative,
                        pooled.stderr = FALSE) {

    n2 <- uniroot(function(n2) {
      power - pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                         alpha = alpha, alternative = alternative,
                         pooled.stderr = pooled.stderr)$power
    }, interval = c(2, 1e10))$root

    return(n2)

  } # ss.approx()

  pwr.exact <- function(prob1, prob2, n2, n.ratio, alpha, alternative) {

    if (n2 > 2000) {

      power <- pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                      alpha = alpha, alternative = alternative)$power

      n.total <- 2000 + ceiling(n.ratio * 2000)

      stop(paste(" Consider method = 'approximate' for total sample size >", n.total), call. = FALSE)

    } else {

      eps <- 1e-10

      n1 <- round(n.ratio * n2)
      x1.seq <- rep(0:n1, times = n2 + 1)
      x2.seq <- rep(0:n2, each = n1 + 1)
      m <- x1.seq + x2.seq
      n.total <- n1 + n2
      k <- n1

      if (alternative == "one.sided") {

        ifelse(prob1 < prob2,
               one.sided.less <- TRUE,
               one.sided.less <- FALSE)

        joint.probs <- dbinom(x1.seq, n1, prob1) * dbinom(x2.seq, n2, prob2)
        p.values <- phyper(x1.seq - ifelse(one.sided.less, 0, 1), m, n.total - m, k, lower.tail = one.sided.less)
        reject <- !any(is.na(p.values)) && any(p.values <= alpha)
        power <- sum(joint.probs[reject])

      } else if (alternative == "two.sided") {

        joint.probs <- dbinom(x1.seq, n1, prob1) * dbinom(x2.seq, n2, prob2)
        valid <- joint.probs > 0
        x1.seq <- x1.seq[valid]
        x2.seq <- x2.seq[valid]
        joint.probs <- joint.probs[valid]
        m.seq <- x1.seq + x2.seq

        m.unique <- unique(m.seq)
        pmf.lookup <- lapply(m.unique, function(m) {
          support <- max(0, m - n2):min(n1, m)
          probs <- dhyper(support, m, n.total - m, n1)
          list(support = support, probs = probs)
        })
        names(pmf.lookup) <- as.character(m.unique)

        passed <- vapply(seq_along(x1.seq), function(i) {
          x1 <- x1.seq[i]
          m  <- m.seq[i]
          key <- as.character(m)
          lookup <- pmf.lookup[[key]]
          p.obs <- dhyper(x1, m, n.total - m, n1)
          p.value <- sum(lookup$probs[lookup$probs <= (p.obs + eps)])
          !any(is.na(p.value)) && any(p.value <= alpha)
        }, logical(1))

        power <- sum(joint.probs[passed])

      } else {

        stop("Hypothesis type not supported.", call. = FALSE)

      }

    } # if n2 > 2000

    power

  } # pwr()


  ss.exact <- function(prob1, prob2, power, n.ratio, alpha, alternative) {

    n2 <- ss.approx(prob1 = prob1, prob2 = prob2, power = power, n.ratio = n.ratio,
                    alpha = alpha, alternative = alternative)
    n2 <- ceiling(n2)

    if (n2 > 500) {

      n.total <- 500 + ceiling(n.ratio * 500)

      stop(paste(" Consider method = 'approximate' for total sample size >", n.total), call. = FALSE)

    } else {


      if (n2 > 200) {

        steprob20 <- 20
        achieved.power <- 0

        while (achieved.power < power) {

          achieved.power <- pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                                      alpha = alpha, alternative = alternative)

          if (achieved.power < power) n2 <- n2 + steprob20

        } # while

        n2 <- n2 - steprob20

      } # n2 > 200

      if (n2 > 50) {

        step5 <- 5
        achieved.power <- 0

        while (achieved.power < power) {

          achieved.power <- pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                                      alpha = alpha, alternative = alternative)

          if (achieved.power < power) n2 <- n2 + step5

        } # while

        n2 <- n2 - step5

      } # n2 > 50

      steprob1 <- 1
      achieved.power <- 0

      while (achieved.power < power) {

        achieved.power <- pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                                    alpha = alpha, alternative = alternative)

        if (achieved.power < power) n2 <- n2 + steprob1

      } # n2 < 50

    } # if (n2 > 1000)

    n2

  } #  ss.exact()



  # method
  if (method == "exact") {

    if (is.null(power)) {

      power <- pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                         alpha = alpha, alternative = alternative)

      mean.alternative <- NA
      sd.alternative <- NA
      mean.null <- NA
      sd.null <- NA
      z.alpha <- NA

      n1 <- n.ratio * n2
      n.total <- n1 + n2

    }

    if (is.null(n2)) {

      n2 <- ss.exact(prob1 = prob1, prob2 = prob2, power = power, n.ratio = n.ratio,
                     alpha = alpha, alternative = alternative)

      n1 <- n.ratio * n2
      n.total <- n1 + n2

      if (ceiling) {

        n1 <- ceiling(n1)
        n2 <- ceiling(n2)
        n.ratio <- n1 / n2
        n.total <- n1 + n2

      }

      power <- pwr.exact(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                         alpha = alpha, alternative = alternative)

      mean.alternative <- NA
      sd.alternative <- NA
      mean.null <- NA
      sd.null <- NA
      z.alpha <- NA

    }

  }  else if (method == "approximate") {

    if (is.null(power)) {

      pwr.obj <- pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                            alpha = alpha, alternative = alternative)

      power <- pwr.obj$power
      mean.alternative <- pwr.obj$mean.alternative
      sd.alternative <- pwr.obj$sd.alternative
      mean.null <- pwr.obj$mean.null
      sd.null <- pwr.obj$sd.null
      z.alpha <- pwr.obj$z.alpha

      n1 <- n.ratio * n2
      n.total <- n1 + n2

    }

    if (is.null(n2)) {

      n2 <- ss.approx(prob1 = prob1, prob2 = prob2, power = power, n.ratio = n.ratio,
                     alpha = alpha, alternative = alternative)

      n1 <- n.ratio * n2
      n.total <- n1 + n2

      if (ceiling) {

        n1 <- ceiling(n1)
        n2 <- ceiling(n2)
        n.ratio <- n1 / n2
        n.total <- n1 + n2

      }

      pwr.obj <- pwr.approx(prob1 = prob1, prob2 = prob2, n2 = n2, n.ratio = n.ratio,
                            alpha = alpha, alternative = alternative)

      power <- pwr.obj$power
      mean.alternative <- pwr.obj$mean.alternative
      sd.alternative <- pwr.obj$sd.alternative
      mean.null <- pwr.obj$mean.null
      sd.null <- pwr.obj$sd.null
      z.alpha <- pwr.obj$z.alpha

    }

  }  # method

  ifelse(method == "exact",
         class <- c("pwrss", "exact", "fisher"),
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
                      test = "Independent Proportions",
                      alpha = alpha,
                      alt = alternative,
                      method = ifelse(method == "exact", "exact", "z"),
                      delta = prob1 - prob2,
                      margin = 0,
                      odds.ratio = (prob1 / (1 - prob1)) /  (prob2 / (1 - prob2)),
                      mean.alternative = mean.alternative,
                      sd.alternative = sd.alternative,
                      mean.null = mean.null,
                      sd.null = sd.null,
                      z.alpha = z.alpha,
                      power = power,
                      n = c(n1 = n1, n2 = n2),
                      n.total = n.total)

    if (pretty) {
      .print.pwrss.fisher(print.obj, verbose = verbose)
    } else {
      .print.ascii.pwrss.fisher(print.obj, verbose = verbose)
    }

  } # verbose


  invisible(structure(list(parms = list(prob1 = prob1, prob2 = prob2,
                                        alpha = alpha, method = method,
                                        alternative = alternative, ceiling = ceiling,
                                        verbose = verbose),
                           test = ifelse(method == "exact", "exact", "z"),
                           delta = prob1 - prob2,
                           odds.ratio = (prob1 / (1 - prob1)) /  (prob2 / (1 - prob2)),
                           mean = mean.alternative,
                           sd = sd.alternative,
                           null.mean = mean.null,
                           null.sd = sd.null,
                           alternative = alternative,
                           z.alpha = z.alpha,
                           power = power,
                           n = c(n1 = n1, n2 = n2),
                           n.total = n.total),
                      class = class))

} # power.exact.fisher()

power.exact.twoprops.fisher <- power.exact.fisher
