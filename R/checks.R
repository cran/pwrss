# internal check functions
check.proportion <- function(...) {

  dots <- substitute(list(...))[-1]
  args <- list(...)

  arg.names <- vapply(dots, function(expr) paste0("`", deparse(expr, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) {
    is.numeric(x) && length(x) == 1 && x >= 0 && x < 1
  }, logical(1))

  bad.names <- arg.names[!check]

  if (length(bad.names) > 0) {
    stop(
      sprintf(
        "Argument%s %s %s not have valid value%s (must be numeric, length 1, > 0 and < 1)",
        if (length(bad.names) > 1) "s" else "",
        paste(bad.names, collapse = if (length(bad.names) > 2) ", " else " and "),
        if (length(bad.names) > 1) "do" else "does",
        if (length(bad.names) > 1) "s" else ""
      ),
      call. = FALSE
    )
  } # bad.names

} # check.probs()

# internal check functions
check.correlation <- function(...) {

  dots <- substitute(list(...))[-1]
  args <- list(...)

  arg.names <- vapply(dots, function(expr) paste0("`", deparse(expr, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) {
    is.numeric(x) && length(x) == 1 && x > -1 && x < 1
  }, logical(1))

  bad.names <- arg.names[!check]

  if (length(bad.names) > 0) {
    stop(
      sprintf(
        "Argument%s %s %s not have valid value%s (must be numeric, length 1, > -1 and < 1)",
        if (length(bad.names) > 1) "s" else "",
        paste(bad.names, collapse = if (length(bad.names) > 2) ", " else " and "),
        if (length(bad.names) > 1) "do" else "does",
        if (length(bad.names) > 1) "s" else ""
      ),
      call. = FALSE
    )
  } # bad.names

} # check.correlations()

check.logical <- function(...) {

  dots <- substitute(list(...))[-1]
  args <- list(...)

  arg.names <- vapply(dots, function(expr) paste0("`", deparse(expr, nlines = 1), "`"), character(1))

  # Check: is logical and of length 1
  check <- vapply(args, function(x) {
    is.logical(x) && length(x) == 1 && !is.na(x)
  }, logical(1))

  bad.names <- arg.names[!check]

  if (length(bad.names) > 0) {
    stop(
      sprintf(
        "Argument%s %s %s not have valid logical value%s (must be TRUE or FALSE only)",
        if (length(bad.names) > 1) "s" else "",
        paste(bad.names, collapse = if (length(bad.names) > 2) ", " else " and "),
        if (length(bad.names) > 1) "do" else "does",
        if (length(bad.names) > 1) "s" else ""
      ),
      call. = FALSE
    )
  } # bad.names

} # check.logical()


check.sample.size <- function(...) {

  dots <- substitute(list(...))[-1]
  args <- list(...)

  arg.names <- vapply(dots, function(expr) paste0("`", deparse(expr, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) {
    is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 0 && abs(x - round(x)) < .Machine$double.eps^0.5
  }, logical(1))

  bad.names <- arg.names[!check]

  if (length(bad.names) > 0) {
    stop(
      sprintf(
        "Argument%s %s %s not have valid sample size value%s (must be integer-like, > 1, and finite)",
        if (length(bad.names) > 1) "s" else "",
        paste(bad.names, collapse = if (length(bad.names) > 2) ", " else " and "),
        if (length(bad.names) > 1) "do" else "does",
        if (length(bad.names) > 1) "s" else ""
      ),
      call. = FALSE
    )
  } # bad.names

} # check.sample.size

check.nonnegative <- function(...) {

  dots <- substitute(list(...))[-1]
  args <- list(...)

  arg.names <- vapply(dots, function(expr) paste0("`", deparse(expr, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) {
    is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 0
  }, logical(1))

  bad.names <- arg.names[!check]

  if (length(bad.names) > 0) {
    stop(sprintf(
      "Argument%s %s %s not have valid non-negative value%s (must be numeric, >= 0, finite)",
      if (length(bad.names) > 1) "s" else "",
      paste(bad.names, collapse = if (length(bad.names) > 2) ", " else " and "),
      if (length(bad.names) > 1) "do" else "does",
      if (length(bad.names) > 1) "s" else ""
    ), call. = FALSE)
  } # bad.names

} # check.nonnegative

check.positive <- function(...) {

  dots <- substitute(list(...))[-1]
  args <- list(...)

  arg.names <- vapply(dots, function(expr) paste0("`", deparse(expr, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) {
    is.numeric(x) && length(x) == 1 && is.finite(x) && x > 0
  }, logical(1))

  bad.names <- arg.names[!check]

  if (length(bad.names) > 0) {
    stop(sprintf(
      "Argument%s %s %s not have valid non-negative value%s (must be numeric, > 0, finite)",
      if (length(bad.names) > 1) "s" else "",
      paste(bad.names, collapse = if (length(bad.names) > 2) ", " else " and "),
      if (length(bad.names) > 1) "do" else "does",
      if (length(bad.names) > 1) "s" else ""
    ), call. = FALSE)
  } # bad.names

} # check.positive

check.numeric <- function(...) {

  dots <- substitute(list(...))[-1]
  args <- list(...)

  arg.names <- vapply(dots, function(expr) paste0("`", deparse(expr, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) {
    is.numeric(x) && length(x) == 1 && is.finite(x)
  }, logical(1))

  bad.names <- arg.names[!check]

  if (length(bad.names) > 0) {
    stop(sprintf(
      "Argument%s %s %s not have valid numeric value%s (must be numeric, scalar, and finite)",
      if (length(bad.names) > 1) "s" else "",
      paste(bad.names, collapse = if (length(bad.names) > 2) ", " else " and "),
      if (length(bad.names) > 1) "do" else "does",
      if (length(bad.names) > 1) "s" else ""
    ), call. = FALSE)
  } # bad.names

} # check.positive

check.correlation.matrix <- function(x) {

  is.square <- nrow(x) == ncol(x)
  if (!is.square) stop("Correlation matrix is not square", call. = FALSE)

  is.symmetric <- isSymmetric.matrix(x)
  if (!is.symmetric) stop("Correlation matrix is not symmetric", call. = FALSE)

  is.positive.definite <- all(eigen(x, symmetric = TRUE)$values > 0)
  if (!is.positive.definite) stop("Correlation matrix is not positive definite", call. = FALSE)

  is.invertible <- det(x) > 0
  if (!is.invertible) stop("Correlation matrix is not invertible", call. = FALSE)

  is.well.conditioned <- kappa(x) < 1000
  if (!is.well.conditioned) stop("Correlation matrix is not well-conditioned", call. = FALSE)

} # check.correlation.matrix()

# check.logical(correct, paired)
# check.proportions(p1, p2, alpha, r2, eta2)
# check.sample.size(n2, n2, k.covariates)
# check.nonnegative(sd1, sd2, f2)
# check.numeric(mu1, mu2, beta0, beta1)
