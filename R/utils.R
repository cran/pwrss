inflate.sample <- function(n, rate = 0.05,
                           ceiling = TRUE, verbose = TRUE) {
  check.sample.size(n)
  n.adj <- n / (1 - rate)
  if (ceiling) n.adj <- ceiling(n.adj)
  if (verbose) cat(n.adj)
  return(invisible(n.adj))
} # inflate.sample

etasq.to.f <- function(eta.squared, verbose = TRUE) {
  check.nonnegative(eta.squared)
  f.squared <- eta.squared / (1 - eta.squared)
  if (verbose) print(c(f.squared = f.squared, f = sqrt(f.squared), eta.squared = eta.squared))
  invisible(list(f.squared = f.squared, f = sqrt(f.squared), eta.squared = eta.squared))
} # etasq.to.f


f.to.etasq <- function(f, verbose = TRUE) {
  check.nonnegative(f)
  f.squared <- f ^ 2
  eta.squared <- f.squared / (1 + f.squared)
  if (verbose) print(c(eta.squared = eta.squared, f.squared = f.squared, f = sqrt(f.squared)))
  invisible(list(eta.squared = eta.squared, f.squared = f.squared, f = sqrt(f.squared)))
} # f.to.etasq

# Fisher's z transformation
cor.to.z <- function(rho, verbose = TRUE) {

  z <- numeric(length = length(rho))
  for (i in seq_along(rho)) {
    check.correlation(rho[i])
    z[i] <- 0.50 * log((1 + rho[i]) / (1 - rho[i]))
  }

  if (verbose) print(c(z = z, rho = rho))
  invisible(list(z = z, rho = rho))

} # cor.to.z()

z.to.cor <- function(z, verbose = TRUE) {

  rho <- numeric(length = length(z))
  for (i in seq_along(z)) {
    if (!is.numeric(z) && !is.finite(z))
      stop("Incorrect value for z", call. = FALSE)
    rho[i] <- (exp(2 * z[i]) - 1) / (exp(2 * z[i]) + 1)
  }

  if (verbose) print(c(rho = rho, z = z))
  invisible(list(rho = rho, z = z))

} # z.to.cor()


# Cohen (1988, S. 109)
# q <.1: no effect;
# .1 <= q < .3: small effect;
# .3 <= q < .5: intermediate effect;
# q =>.5: large effect

cors.to.q <- function(rho1, rho2, verbose = TRUE) {

  check.correlation(rho1, rho2)

  q <- cor.to.z(rho1, FALSE)$z - cor.to.z(rho2, FALSE)$z

  if (verbose) print(c(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))
  invisible(list(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))

} # cors.to.q()


q.to.cors <- function(q, rho1 = NULL, rho2 = NULL, verbose = TRUE) {

  if (!is.null(rho1)) check.correlation(rho1)
  if (!is.null(rho2)) check.correlation(rho2)
  check.logical(verbose)
  check.numeric(q)

  if (is.null(rho1) && is.null(rho2))
    stop("Both `rho1` and `rho2` cannot be NULL.", call. = FALSE)
  if (!is.null(rho1) && !is.null(rho2))
    stop("Exactly one of the `rho1` or `rho2` should be NULL.", call. = FALSE)

  if (is.null(rho1) && !is.null(rho2)) {
    z1 <- cor.to.z(rho2, FALSE)$z - q
    rho1 <- (exp(2 * z1) - 1) / (exp(2 * z1) + 1)
  }

  if (!is.null(rho1) && is.null(rho2)) {
    z2 <- cor.to.z(rho1, FALSE)$z + q
    rho2 <- (exp(2 * z2) - 1) / (exp(2 * z2) + 1)
  }

  if (verbose) print(c(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))
  invisible(list(q = q, delta = rho1 - rho2, rho1 = rho1, rho2 = rho2))

} # q.to.cors()

d.to.cles <- function(d, design = c("independent", "paired", "one.sample"), verbose = TRUE) {

  check.numeric(d)
  check.logical(verbose)

  design <- tolower(match.arg(design))

  if (design == "independent") {
    prob <- pnorm(d / sqrt(2))
  } else {
    prob <- pnorm(d)
  }

  if (verbose) print(c(cles = prob, d = d))
  invisible(list(cles = prob, d = d))

} # d.to.cles


cles.to.d <- function(cles, design = c("independent", "paired", "one.sample"), verbose = TRUE) {

  check.proportion(cles)

  design <- tolower(match.arg(design))

  if (design == "independent") {
    d <- sqrt(2) * qnorm(cles)
  } else {
    d <- qnorm(cles)
  }

  if (verbose) print(c(d = d, cles = cles))
  invisible(list(d = d, cles = cles))

} # cles.to.d


means.to.d <- function(mu1, mu2 = 0,
                       sd1 = 1, sd2 = 1,
                       n2, n.ratio = 1,
                       paired = FALSE,
                       rho.paired = 0.50,
                       verbose = TRUE) {

  check.logical(paired)
  check.numeric(mu1, mu2)
  check.correlation(rho.paired)
  check.positive(sd1, sd2, n.ratio)
  check.sample.size(n2)

  n1 <- n.ratio * n2

  alpha <- 0.05
  var.ratio <- sd1 ^ 2 / sd2 ^ 2
  f.crit.lower <- qf(alpha, df1 = n1 - 1, df2 = n2 - 1, lower.tail = TRUE)
  f.crit.upper <- qf(1 - alpha, df1 = n1 - 1, df2 = n2 - 1, lower.tail = TRUE)
  if (var.ratio <= f.crit.lower || var.ratio >= f.crit.upper)
    warning("Interpretation of Cohen's d may no longer be valid when variances differ beyond sampling error.\n", call. = FALSE)

  if (paired) {
    pooled.sd <- sqrt(sd1 ^ 2 + sd2 ^ 2 - 2 * sd1 * sd2 * rho.paired)
  } else {
    pooled.sd <- sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
  }

  d <- (mu1 - mu2) / pooled.sd

  if (verbose) print(c(d = d))
  invisible(list(d = d, mu1 = mu1, mu2 = mu2,
                 sd1 = sd1, sd2 = sd2,
                 pooled.sd = pooled.sd,
                 var.ratio = var.ratio,
                 n1 = n1, n2 = n2,
                 n.ratio = n.ratio,
                 paired = paired, rho.paired = rho.paired,
                 verbose = verbose))

} # means.to.d

probs.to.h <- function(prob1, prob2 = 0.50, verbose = TRUE) {

  check.proportion(prob1, prob2)
  check.logical(verbose)

  h <- 2 * asin(sqrt(prob1)) - 2 * asin(sqrt(prob2))

  if (verbose) print(c(h = h, prob1 = prob1, prob2 = prob2))

  invisible(list(h = h, prob1 = prob1, prob2 = prob2))

} # probs.to.h

# Zhang, Cao, and Ahn (2017)
joint.probs.2x2 <- function(prob1, prob2, rho = 0.50, verbose = TRUE) {

  check.proportion(prob1, prob2)
  check.logical(verbose)

  if (!is.numeric(rho) && rho < -1 && rho > 1) stop("Incorrect value for `rho`", call. = FALSE)

  rho.min <- max(
    -sqrt(prob1 * prob2 / ((1 - prob1) * (1 - prob2))),
    -sqrt((1 - prob1) * (1 - prob2) / (prob2 * prob2))
  )

  rho.max <- min(
    sqrt(prob1 * (1 - prob2) / (prob2 * (1 - prob1))),
    sqrt(prob2 * (1 - prob1) / (prob1 * (1 - prob2)))
  )

  if (rho < rho.min || rho > rho.max) {

    stop(paste("Combination of 'prob1', 'prob2' and 'rho' is not feasible. \n 'rho' should be between", round(rho.min, 3), "and", round(rho.max, 3)), call. = FALSE)

  }

  prob11 <- rho * sqrt(prob1 * (1 - prob1) * prob2 * (1 - prob2)) + prob1 * prob2
  prob10 <- prob1 - prob11
  prob01 <- prob2 - prob11
  prob00 <- 1 - (prob11 + prob10 + prob01)

  if (verbose) {

    print(c(prob11 = prob11, prob10 = prob10, prob01 = prob01, prob00 = prob00))

  }

  return(invisible(list(parms = list(prob1 = prob1, prob2 = prob2, rho = rho,
                                     rho.min = rho.min, rho.max = rho.max),
                        prob11 = prob11, prob10 = prob10, prob01 = prob01, prob00 = prob00)))

} # joint.probs.2x2


marginal.probs.2x2 <- function(prob11, prob10, prob01, prob00, verbose = TRUE) {

  check.proportion(prob11, prob10, prob01, prob00)
  check.logical(verbose)

  total <- prob11 + prob10 + prob01 + prob00

  if (abs(total - 1) > 1e-6) {
    stop("Joint probabilities must sum to 1.", call. = FALSE)
  }

  prob1 <- prob11 + prob10
  prob2 <- prob11 + prob01

  if (prob1 == 1 || prob1 == 0 || prob2 == 1 || prob2 == 0) {

    rho <- NA
    warning("Undefined correlation: division by zero in denominator.")

  } else {

    rho <- (prob11 - prob1 * prob2) / sqrt(prob1 * (1 - prob1) * prob2 * (1 - prob2))

  }

  if (verbose) {

    print(c(prob1 = prob1, prob2 = prob2, rho = rho))

  }

  return(invisible(list(parms = list(prob11 = prob11, prob10 = prob10, prob01 = prob01, prob00 = prob00, rho = rho),
                        prob1 = prob1, prob2 = prob2, rho = rho)))

} # marginal.probs.2x2


# internal function to get some chisq stat
probs.to.w <- function(prob.matrix, null.prob.matrix = NULL, verbose = TRUE) {

  if (any(prob.matrix < 0) || any(prob.matrix > 1))
    stop("Matrix elements outside of [0, 1] range.", call. = FALSE)

  if (is.null(null.prob.matrix)) {
    ifelse(is.vector(prob.matrix),
           null.prob.matrix <- rep(1 / length(prob.matrix), length(prob.matrix)),
           null.prob.matrix <- outer(rowSums(prob.matrix), colSums(prob.matrix)) / sum(prob.matrix))
  } else {
    if (any(null.prob.matrix < 0) || any(null.prob.matrix > 1))
      stop("Matrix elements outside of [0, 1] range.", call. = FALSE)
    if (isFALSE(all(dim(prob.matrix) == dim(null.prob.matrix))))
      stop("Dimensions should match `prob.matrix`.", call. = FALSE)
  }

  if (is.vector(prob.matrix)) {
    if (length(prob.matrix) != length(null.prob.matrix))
      stop("Length of 'prob.matrix' and 'null.prob.matrix' should match.", call. = FALSE)
    if (sum(prob.matrix) != 1 || sum(null.prob.matrix) != 1)
      stop("Cell probabilities should sum to 1.", call. = FALSE)
  } else if (is.matrix(prob.matrix)) {
    if (any(dim(prob.matrix) != dim(null.prob.matrix)))
      stop("Dimensions for 'prob.matrix' and 'null.prob.matrix' do not match.", call. = FALSE)
  } else {
    stop("Not a valid 'prob.matrix'.", call. = FALSE)
  }

  ifelse(is.vector(prob.matrix),
         df <- length(prob.matrix) - 1,
         df <- (nrow(prob.matrix) - 1) * (ncol(prob.matrix) - 1))

  chisq <- sum((prob.matrix - null.prob.matrix) ^ 2 / null.prob.matrix)

  ifelse(is.vector(prob.matrix),
         mdf <- 1,
         mdf <- min(nrow(prob.matrix) - 1, ncol(prob.matrix) - 1))

  w <- sqrt(chisq / (sum(prob.matrix) * mdf))

  if (w > 1) warning("w > 1 is unrealistic, please check your input.", call. = FALSE)
  if (verbose) print(c(w = w, df = df))
  invisible(list(w = w, df = df, prob.matrix = prob.matrix, null.prob.matrix = null.prob.matrix))

} # probs.to.w
