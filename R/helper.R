# ensure valid values for verbose: logical values are converted to 0/1, numerical values are kept as long as they
# have the length 1, are integer and have the value 0, 1 or 2, in any other case, 1 (the default) is returned
ensure.verbose <- function(verbose = NULL) {
  if (is.logical(verbose)) {
    as.integer(verbose)
  } else if (is.numeric(verbose)) {
    ifelse(length(verbose) == 1 && verbose %% 1 == 0 && verbose %in% c(-1, 0, 1, 2), verbose, 1)
  } else {
    1
  }
}

# check the input parameters n, power and es, and return which calculation is requested
get.requested <- function(es = NULL, n = NULL, power = NULL) {

  es_vars <- gsub("list", "", deparse(substitute(es), nlines = 1))
  n_vars <- deparse(substitute(n), nlines = 1)

  if (is.list(es)) {
    if        (sum(unlist(lapply(es, is.null))) == 2) {
      stop(sprintf("Exactly one element / entry of `%s` can be NULL, not both.",
                   paste(strsplit(gsub("[() ]", "", es_vars), ",")[[1]], collapse = "` or `")), call. = FALSE)
    } else if (sum(unlist(lapply(es, is.null))) == 1) {
      es <- NULL
    }
  }
  es <- unlist(es)[1] # works for vectors too (then unlist has no effect)

  is_na <- check.na(es, n, power)
  if (sum(is_na) > 1)
    stop(sprintf("Maximally one input parameter (`%s`, `%s`, or `power`) can be excluded from checking.", es_vars, n_vars),
         call. = FALSE)

  if (sum(check.not_null(n, power, es)) != 2) {
    n.parms <- ifelse(!any(is_na), "two", "one")
    s.parms <- do.call(sprintf,
                     list(ifelse(!any(is_na), "%s`, `%s`, or `%s", "%s` or `%s"), es_vars, n_vars, "power")[c(TRUE, !is_na)])
    stop(sprintf("Exactly %s of the parameters `%s` must be given, one has to be NULL.", n.parms, s.parms), call. = FALSE)
  }

  invisible(c("es", "n", "power")[check.null(es, n, power)]) # return what is requested / to be calculated

} # get.requested

get.interval <- function(null.ncp, distribution = c("z", "t", "lp", "binom"), alpha = 0.05, alternative = "two.sided",
                         req.sign = "+", sd = NULL, df = NULL) {

    if (distribution == "z") {

      min.alt <- stats::qnorm(1e-10,     mean = stats::qnorm(alpha,     mean = min(null.ncp), sd = sd), sd = sd)
      max.alt <- stats::qnorm(1 - 1e-10, mean = stats::qnorm(1 - alpha, mean = max(null.ncp), sd = sd), sd = sd)

    } else if (distribution == "t") {

      suppressWarnings({
        min.alt <- stats::qt(1e-10,     ncp = stats::qt(alpha,     ncp = min(null.ncp), df = df), df = df)
        max.alt <- stats::qt(1 - 1e-10, ncp = stats::qt(1 - alpha, ncp = max(null.ncp), df = df), df = df)
      })

    } else if (distribution == "lp") {

      suppressMessages({
        min.alt <- sadists::qlambdap(1e-10,     t = sadists::qlambdap(alpha,     t = min(null.ncp), df = df), df = df)
        max.alt <- sadists::qlambdap(1 - 1e-10, t = sadists::qlambdap(1 - alpha, t = max(null.ncp), df = df), df = df)
      })

    } else if (distribution == "binom") {

        min.alt <- 1e-4
        max.alt <- 1 - 1e-4

    }

    if (check.null_sign(req.sign, alternative)) {
      sort(null.ncp)
    } else if (check.pos_sign(req.sign)) {
      c(max(null.ncp), max.alt)
    } else { # req.sign neither null nor positive
      c(min.alt, min(null.ncp))
    }

} # get.interval

isInt <- function(x) is.numeric(x) && !any(abs(x - round(x)) > .Machine$double.eps ^ 2 / 3)

# lenInt <- function(n) ifelse(n <= 1, 1, ceiling(log10(abs(n))) + as.integer(n %% 10 == 0))

check.snap4plot <- function(snpFle = "", pltFnc = NULL, pltPrm = list(), pltWdt = 800, pltHgh = 800) {
  if (nchar(Sys.getenv("GITHUB_ACTIONS")) == 0) { # ensures that the code only runs on a local machine, not as GitHub action
    tmpFle <- tempfile(fileext = ".png")
    addPrm <- list(alpha = 0.05, verbose = 0)[c("alpha", "verbose") %in% names(formals(pltFnc))]
    testthat::announce_snapshot_file(name = snpFle)
    grDevices::png(tmpFle, width = pltWdt, height = pltHgh)
    if (any(c("plot", "plot.main") %in% names(formals(pltFnc)))) {
      do.call(pltFnc, c(pltPrm, addPrm))
    } else {
      plot(do.call(pltFnc, c(pltPrm, addPrm)))
    }
    grDevices::dev.off()
    testthat::expect_snapshot_file(path = tmpFle, name = snpFle, variant = Sys.info()[["sysname"]])
    unlink(tmpFle)
  }
}

uniroot_break <- function(ur.obj) {
  !inherits(ur.obj, "try-error") ||
    (inherits(ur.obj, "try-error") && attr(ur.obj, "condition")[["message"]] != "f() values at end points not of opposite sign")
}
