# internal check functions
check.proportion <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))
  fnc.name <- as.character(match.call()[[1]])
  val.rng <- c(ifelse(fnc.name == "check.power", 0.01, 0), ifelse(fnc.name == "check.power", 0.99, 1))

  check <- vapply(args, function(x) length(x) == 1 && is.numeric(x) && x >= val.rng[1] && x <= val.rng[2], logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], fnc.name), call. = FALSE)

} # check.proportion()

check.power <- check.proportion

check.correlation <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) length(x) == 1 && is.numeric(x) && x >= -1 && x <= 1, logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], as.character(match.call()[[1]])), call. = FALSE)

} # check.correlation()

check.logical <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) length(x) == 1 && is.logical(x) && !is.na(x), logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], as.character(match.call()[[1]])), call. = FALSE)

} # check.logical()


check.sample.size <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) length(x) == 1 && isInt(x) && is.finite(x) && x >= 2, logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], as.character(match.call()[[1]])), call. = FALSE)

} # check.sample.size

check.factor.level <- check.sample.size # check.factor.level

check.size <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) length(x) == 1 && isInt(x) && is.finite(x) && x >= 0, logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], as.character(match.call()[[1]])), call. = FALSE)

} # check.size

check.nonnegative <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) length(x) == 1 && is.numeric(x) && is.finite(x) && x >= 0, logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], as.character(match.call()[[1]])), call. = FALSE)

} # check.nonnegative

check.positive <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) length(x) == 1 && is.numeric(x) && is.finite(x) && x > 0, logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], as.character(match.call()[[1]])), call. = FALSE)

} # check.positive

check.numeric <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))

  check <- vapply(args, function(x) length(x) == 1 && is.numeric(x) && is.finite(x), logical(1))

  if (!all(check))
    stop(format_errmsg(arg.names[!check], as.character(match.call()[[1]])), call. = FALSE)

} # check.positive

check.vector <- function(x, fnc, min.length = 2, max.length = Inf) {

  var.name <- deparse(substitute(x), nlines = 1)
  if (!is.vector(x) || length(x) < min.length || length(x) > max.length) {
    length.msg <- paste0(ifelse(min.length != max.length, "at least ", ""), sprintf("%d", min.length),
                         ifelse(max.length < Inf && min.length != max.length, sprintf(" and maximally %d", max.length), ""))
    stop(sprintf("`%s` neeeds to be a vector with a length of %s.", var.name, length.msg), call. = FALSE)
  }

  fnc.name <- deparse(substitute(fnc), nlines = 1)
  err.msg <- sprintf("All elements of `%s` need to be valid %s values (%s)", var.name, fnc2type(fnc.name), valid.cond(fnc.name))
  tryCatch(invisible(unlist(lapply(x, fnc))), # returns NULL if successful
           error = function(e) stop(err.msg, call. = FALSE))

} # check.vector

check.same.lengths <- function(...) {

  args <- list(...)
  arg.names <- vapply(substitute(list(...))[-1], function(fn) paste0("`", deparse(fn, nlines = 1), "`"), character(1))
  if (is.null(args[[1]]))
    stop(sprintf("To use `check.same.lengths`, %s can not be NULL.", arg.names[[1]]), call. = FALSE)

  check <- rep(TRUE, length(args))
  check[c(FALSE, !unlist(lapply(args, is.null))[-1])] <-
    unlist(lapply(args, length)[c(FALSE, !unlist(lapply(args, is.null))[-1])]) == length(args[[1]])

  if (!all(check))
    stop(
      sprintf("%s should have the same length as %s.",
              paste(arg.names[!check], collapse = ifelse(sum(!check) == 2, " and ", ", ")),
              arg.names[1]),
      call. = FALSE
    )

} # check.same.lengths

check.margins <- function(x, fnc, alternative = "two.sided") {

  var.name <- deparse(substitute(x), nlines = 1)

  fnc.name <- deparse(substitute(fnc), nlines = 1)
  err.msg <- sprintf("All elements of `%s` need to be valid %s values (%s)", var.name, fnc2type(fnc.name), valid.cond(fnc.name))
  tryCatch(invisible(unlist(lapply(x, fnc))),
           error = function(e) stop(err.msg, call. = FALSE))

  # if a margin needs order (typically the case for proportions), it can't be duplicated and the order needs checking
  needs_order <- fnc.name == "check.proportion"

  # ensure correct length and correct values of the variable
  if (alternative %in% c("two.sided", "one.sided") && length(x) != 1) {
    stop("If `alternative` is \"two.sided\" or \"one.sided\", `", var.name, "` must be of length one.", call. = FALSE)
  } else if (alternative == "two.one.sided" && !needs_order && length(x) == 1 && x != 0) {
    x <- abs(x) * c(-1, 1) # assumes equidistance from zero; if x is 0, a scalar is generated (thus: `x != 0` above)
  } else if (alternative == "two.one.sided" && (length(x) != 2 || (needs_order && x[1] > x[2]))) {
    stop("If `alternative` is \"two.one.sided\", `", var.name, "` must be of ",
         ifelse(needs_order, "", "length one (absolute value, that must be different from 0), or "), "length two ",
         "(lower and upper bounds",
         ifelse(needs_order, ", with the upper bound being larger than the lower bound", ""), ").",
         call. = FALSE)
  }

  x
}

check.correlation.matrix <- function(x) {

  var.name <- deparse(substitute(x), nlines = 1)

  is.square <- nrow(x) == ncol(x)
  if (!is.square)
    stop(sprintf("Correlation matrix `%s` is not square", var.name), call. = FALSE)

  is.symmetric <- isSymmetric.matrix(x)
  if (!is.symmetric)
    stop(sprintf("Correlation matrix `%s` is not symmetric", var.name), call. = FALSE)

  # check that all values (only the upper triangular, as the matrix is symmetric) are valid correlations,
  # either all values are valid or an error is thrown
  errmsg <- sprintf("The values in the correlation matrix (`%s`) must be numeric, >= -1 and <= 1", var.name)
  tryCatch(invisible(unlist(lapply(x[upper.tri(x)], check.correlation))),
           error = function(e) stop(errmsg, call. = FALSE))

  correct.diagonal <- all(diag(x) == 1)
  if (!correct.diagonal)
    stop(sprintf("All values in the main diagonal of the correlation matrix (`%s`) must be 1", var.name), call. = FALSE)

  is.positive.definite <- all(eigen(x, symmetric = TRUE)$values > 0)
  if (!is.positive.definite)
    stop(sprintf("Correlation matrix `%s` is not positive definite", var.name), call. = FALSE)

  is.invertible <- det(x) > 0
  if (!is.invertible)
    stop(sprintf("Correlation matrix `%s` is not invertible", var.name), call. = FALSE)

  is.well.conditioned <- kappa(x) < 1000
  if (!is.well.conditioned)
    stop(sprintf("Correlation matrix `%s` is not well-conditioned", var.name), call. = FALSE)

} # check.correlation.matrix()

check.na <- function(...) {

  args <- list(...)

  vapply(args, function(a) !is.null(a) && all(is.na(a)), logical(1))

} # check.na

check.null <- function(...) {

  args <- list(...)

  vapply(args, is.null, logical(1))

} # check.null

check.not_null <- function(...) !check.null(...) # check.not_null

check.pos_sign <- function(req.sign) {
  if (req.sign %in% c("+", 1, "1", "+1", "positive", "pozitive")) {
    TRUE
  } else if (req.sign %in% c("-", -1, "-1", "negative")) {
    FALSE
  } else {
    stop("`req.sign` can only be `+` and `-` for this function.", call. = FALSE)
  }
} # check.pos_sign

check.null_sign <- function(req.sign, alternative) {
  if (req.sign %in% c(" ", 0, "0", "", "null")) {
    if (alternative != "two.one.sided")
      stop("`req.sign` can not be `0` for 'one.sided' or 'two.sided' hypothesis tests.", call. = FALSE)

    TRUE
  } else if (req.sign %in% c("+", 1, "1", "+1", "positive", "pozitive", "-", -1, "-1", "negative")) {
    FALSE
  } else {
    stop("`req.sign` can only be `+`, `-` and ` ` for this function.", call. = FALSE)
  }
} # check.null_sign

# helper function(s) ---------------------------------------------------------------------------------------------------
format_errmsg <- function(names = NULL, fnc.name = NULL) {
  sprintf("Argument%s %s %s not have %svalid %s value%s (must be length 1, %s).",
          ifelse(length(names) > 1, "s", ""),
          paste(names, collapse = ifelse(length(names) > 2, ", ", " and ")),
          ifelse(length(names) > 1, "do", "does"),
          ifelse(length(names) > 1, "", "a "),
          fnc2type(fnc.name),
          ifelse(length(names) > 1, "s", ""),
          valid.cond(fnc.name))
}

fnc2type <- function(fnc.name = "") {
  gsub("factor.level", "factor level", fixed = TRUE,
       gsub("nonnegative", "non-negative", fixed = TRUE,
       gsub("sample.size", "sample size", fixed = TRUE,
       gsub("check.", "", fnc.name, fixed = TRUE))))
}

valid.cond <- function(fnc.name = "") {
  if (fnc.name == "check.proportion") {
    "numeric, >= 0, and <= 1"
  } else if (fnc.name == "check.power") {
    "numeric, >= 0.01, and <= 0.99"
  } else if (fnc.name == "check.correlation") {
    "numeric, >= -1, and <= 1"
  } else if (fnc.name == "check.logical") {
    "TRUE or FALSE"
  } else if (fnc.name %in% c("check.sample.size", "check.factor.level")) {
    "integer-like, > 1, and finite"
  } else if (fnc.name == "check.size") {
    "integer-like, >= 0, and finite"
  } else if (fnc.name == "check.nonnegative") {
    "numeric, >= 0, and finite"
  } else if (fnc.name == "check.positive") {
    "numeric, > 0, and finite"
  } else if (fnc.name == "check.numeric") {
    "numeric, and finite"
  } else {
    stop(sprintf("%s is not a valid check-function", fnc.name), call. = FALSE)
  }
}

# check.logical(correct, paired)
# check.proportions(p1, p2, alpha, r2, eta2)
# check.sample.size(n2, n2)
# check.nonnegative(sd1, sd2, f2)
# check.positive(k.covariates)
# check.numeric(mu1, mu2, beta0, beta1)
