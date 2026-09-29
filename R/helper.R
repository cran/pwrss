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

  es_vars <- gsub("list", "", deparse(substitute(es), nlines = 1), fixed = TRUE)
  n_vars <- deparse(substitute(n), nlines = 1)

  if (is.list(es)) {
    if        (sum(unlist(lapply(es, is.null))) == 2) {
      stop(sprintf("Exactly one element / entry of `%s` can be NULL, not both.",
                   paste(strsplit(gsub("[() ]", "", es_vars, fixed = TRUE), ",", fixed = TRUE)[[1]], collapse = "` or `")),
           call. = FALSE)
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
    n.parms <- ifelse(any(is_na), "one", "two")
    s.parms <- do.call(sprintf,
                     list(ifelse(any(is_na), "%s` or `%s", "%s`, `%s`, or `%s"), es_vars, n_vars, "power")[c(TRUE, !is_na)])
    stop(sprintf("Exactly %s of the parameters `%s` must be given, one has to be NULL.", n.parms, s.parms), call. = FALSE)
  }

  invisible(c("es", "n", "power")[check.null(es, n, power)]) # return what is requested / to be calculated

} # get.requested

get.interval <- function(null.ncp, req.sign, distribution = c("z", "t", "lp", "binom"), alpha = 0.05,
                         alternative = c("two.sided", "one.sided", "two.one.sided"),
                         sd = NULL, df = NULL) {

  distribution <- match.arg(distribution)
  alternative  <- match.arg(alternative)
  if (distribution %in% c("lp", "t")) {
    check.positive(df)
  } else if (distribution == "z") {
    check.positive(sd)
  }

  if (check.null_sign(req.sign, alternative)) { # req.sign is "0" (or equivalent)

    sort(null.ncp)

  } else if (check.pos_sign(req.sign)) { # req.sign is "+" (or equivalent)

    alpha.upr <- ifelse(alternative == "two.sided", 1 - alpha / 2, 1 - alpha)
    alt.upr <- switch(distribution,
                      "z" = stats::qnorm(1 - 1e-10, mean = stats::qnorm(alpha.upr, mean = max(null.ncp), sd = sd), sd = sd),
                      "t" = suppressWarnings(stats::qt(1 - 1e-10, ncp = stats::qt(alpha.upr, ncp = max(null.ncp), df = df), df = df)),
                      "lp" = qlambdap(1 - 1e-10, ncp = qlambdap(alpha.upr, ncp = max(null.ncp), df = df), df = df),
                      "binom" = 1 - 1e-4)
    c(max(null.ncp), alt.upr)

  } else { # req.sign is "-" (or equivalent)

    alpha.lwr <- ifelse(alternative == "two.sided", alpha / 2, alpha)
    alt.lwr <- switch(distribution,
                      "z" = stats::qnorm(1e-10, mean = stats::qnorm(alpha.lwr, mean = min(null.ncp), sd = sd), sd = sd),
                      "t" = suppressWarnings(stats::qt(1e-10, ncp = stats::qt(alpha.lwr, ncp = min(null.ncp), df = df), df = df)),
                      "lp" = qlambdap(1e-10, ncp = qlambdap(alpha.lwr, ncp = min(null.ncp), df = df), df = df),
                      "binom" = 1e-4)
    c(alt.lwr, min(null.ncp))

  }

} # get.interval

isInt <- function(x) is.numeric(x) && !any(abs(x - round(x)) > .Machine$double.eps ^ (2 / 3))

# lenInt <- function(n) ifelse(n <= 1, 1, ceiling(log10(abs(n))) + as.integer(n %% 10 == 0))

check.snap4plot <- function(snpFle = "", pltFnc = NULL, pltPrm = list(), pltWdt = 800, pltHgh = 800) {
  # ensures that the code only runs on a local machine, not as GitHub action or in a Docker
  if (nchar(Sys.getenv("GITHUB_ACTIONS")) > 0 || any(file.exists(c("/.dockerenv", "/run/.containerenv")))) {
    testthat::announce_snapshot_file(name = snpFle)
  } else {
    tmpFle <- tempfile(fileext = ".png")
    addPrm <- list(alpha = 0.05, verbose = 0)[c("alpha", "verbose") %in% names(formals(pltFnc))]
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
