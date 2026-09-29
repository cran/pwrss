# helper and formatting functions --------------------------------------------------------------------------------------
.header <- function(requested, ncp = FALSE, utf = FALSE) {

  # eff.dsc <- ifelse(ncp, "NCP", "EFFECT")

  if (utf) {

    RC <- switch(requested,
                 `n` =        "           \033[34m SAMPLE SIZE CALCULATION \033[0m              ",
                 `power` =    "               \033[34m POWER CALCULATION \033[0m                ",
                 `es` = ifelse(ncp,
                               paste("      \033[34m MINIMUM DETECTABLE NCP CALCULATION \033[0m        "),
                               paste("     \033[34m MINIMUM DETECTABLE EFFECT CALCULATION \033[0m      ")))

    paste0(paste0("\u2554", strrep("\u2550", 50), "\u2557", "\n"),
           paste0("\u2551", RC,                   "\u2551", "\n"),
           paste0("\u255A", strrep("\u2550", 50), "\u255D", "\n\n"))

  } else {

    RC <- switch(requested,
                 `n` =        "             SAMPLE SIZE CALCULATION              ",
                 `power` =    "                POWER CALCULATION                 ",
                 `es` = ifelse(ncp,
                               paste("       MINIMUM DETECTABLE NCP CALCULATION         "),
                               paste("      MINIMUM DETECTABLE EFFECT CALCULATION       ")))

    paste0(paste0("+", strrep("-", 50), "+", "\n"),
           paste0("|", RC,              "|", "\n"),
           paste0("+", strrep("-", 50), "+", "\n\n"))

  }

}

.topic <- function(topic, utf = FALSE, emphasize = FALSE) {

  tline    <- strrep(ifelse(utf, "\u2500", "-"), 52)
  emph_beg <- ifelse(utf && emphasize, "\033[36m", "")
  emph_end <- ifelse(utf && emphasize, "\033[0m",  "")

  paste0(paste0(emph_beg, tline, emph_end, "\n"),
         paste0(emph_beg, topic, emph_end, "\n"),
         paste0(emph_beg, tline, emph_end, "\n"))

}

.fmt_val <- function(val = NA, digits = 3, spf_neg = "") {
  if (is.numeric(val)) {
    fmt <- paste(rep(paste0("%", spf_neg, paste0(".", ifelse(isInt(val), 0, digits), "f")), length(val)),
                 collapse = ifelse(length(val) < 3, " and ", ", "))
  } else if (is.character(val) || is.na(val)) {
    fmt <- paste0(strrep(" ", nchar(spf_neg)), "%s")
  }

  if (spf_neg != "")
    gsub("^\\+| \\+", " ", do.call(sprintf, c(fmt = fmt, as.list(val))))
  else
    do.call(sprintf, c(fmt = fmt, as.list(val)))
}

.fmt_utf <- function(n) {
  rpl_list <- list(`0` = "\u2080", `1` = "\u2081", `2` = "\u2082", `3` = "\u2083", `4` = "\u2084", prob = "P",
                   beta = "\u03B2", rho = "\u03C1", eta = "\u03B7")

  for (r in names(rpl_list)) n <- gsub(r, rpl_list[[r]], n)

  n
}

.fmt_adsc <- function(f, utf = FALSE) {
  f <- gsub("^std.", "Std. ", gsub("odds.ratio", "Odds Ratio", gsub("rate.ratio", "Rate Ratio", f)))
  if (utf) {
    .fmt_utf(gsub("mean", "\u03BC", gsub("lambda", "\u03BB", gsub(".squared|.squared.change", "\u00b2",
               gsub("^r\\.", "R.", f)), fixed = TRUE), fixed = TRUE))
  } else {
    gsub(".squared|.squared.change", "-squared", gsub("^r\\.", "R.", f))
  }
}

.nspacer <- function(x) floor(length(gregexpr("\u2009", x, fixed = TRUE)[[1]]) / 2)

.pad <- function(dsc, maxlen) strrep(" ", maxlen + .nspacer(dsc) - nchar(dsc))

.a_pad <- function(n, utf) {
  if (any(grepl("^tgt.effect$", n))) {
    if (any(grepl("^rho\\d{2}",   n))) return(as.integer(utf) + 1)
    if (any(grepl("^prob\\d{1}$", n))) return(as.integer(utf) + 1)
    if (any(grepl("^prob\\d{2}$", n))) return(2)
    if (any(grepl("^std.beta$",   n))) return(5)
    as.integer(utf)
  } else if (any(grepl("^df$", n)) && !any(grepl("^n.total$", n))) {
    as.integer(utf)
  } else if (any(grepl("^n.pres$", n))) {
    as.integer(utf)
  } else if (any(grepl("^eta.squared", n))) {
    as.integer(!utf) * 5
  } else if (any(grepl("^rate.ratio", n))) {
    as.integer(utf) * 3 + 4
  } else if (any(grepl("^odds.ratio", n)) && !any(grepl("^tgt.effect",  n))) {
    as.integer(utf) * 3 + 4
  } else if (any(grepl("^r.squared", n))) {
    as.integer(!utf) * 3
  } else {
    0
  }
}

.keyparms <- function(x, parms_mtx, utf = FALSE, digits = 3) {
  # ensures that the variables in the parms_mtx exist in the printing object and remove those that are NA
  parms_mtx <- parms_mtx[vapply(parms_mtx[, "var.name"], function(n) utils::hasName(x, n), logical(1)), ]
  parms_mtx <- parms_mtx[vapply(parms_mtx[, "var.name"], function(n) !all(is.na(x[[n]])),  logical(1)), ]
  spf_neg   <- ifelse(any(stats::na.omit(suppressWarnings(as.numeric(unlist(x[parms_mtx[, 1]])))) < 0), "+", "")
  parms_col <- ifelse(utf, 3, 2)
  parms_max <- max(vapply(parms_mtx[, parms_col], function(x) nchar(x) - .nspacer(x), numeric(1)))
  parms_out <- .topic("Key Parameters", utf)
  for (r in seq_len(nrow(parms_mtx))) {
    dsc <- parms_mtx[[r, parms_col]]
    parms_out <- c(parms_out,
                   sprintf("  %s%s = %s\n", dsc, .pad(dsc, parms_max), .fmt_val(x[[parms_mtx[r, 1]]], digits, spf_neg)))
  }

  paste0(paste(parms_out, collapse = ""), "\n")
}

.eq <- function(utf = FALSE) ifelse(utf, "=",      " =")
.le <- function(utf = FALSE) ifelse(utf, "\u2264", "<=")
.ls <- function(utf = FALSE) ifelse(utf, "<",      " <")
.ge <- function(utf = FALSE) ifelse(utf, "\u2265", ">=")
.gt <- function(utf = FALSE) ifelse(utf, ">",      " >")
.ne <- function(utf = FALSE) ifelse(utf, "\u2260", "!=")

.h0_sign <- function(alt = c("two.sided", "one.sided"), less = FALSE, utf = FALSE) {
  alt <- match.arg(alt)

  if (alt == "two.sided") .eq(utf) else if (alt == "one.sided" && less) .ge(utf) else if (alt == "one.sided" && !less) .le(utf)
}

.h1_sign <- function(alt = c("two.sided", "one.sided"), less = FALSE, utf = FALSE) {
  alt <- match.arg(alt)

  if (alt == "two.sided") .ne(utf) else if (alt == "one.sided" && less) .ls(utf) else if (alt == "one.sided" && !less) .gt(utf)
}

.h0_twoone <- function(tgt, mrg, utf, alt = "two.sided", val.alt = NA, val.null = NA, val.mrg = NULL) {
  if (alt %in% c("one.sided", "two.sided")) {
    paste(tgt, .h0_sign(alt, val.alt < val.null, utf), ifelse(!is.null(val.mrg) && val.mrg == 0, "0", mrg))
  } else if (alt == "two.one.sided" && (val.alt > min(val.null) && val.alt < max(val.null))) {
    if (utf)
      c(paste0(tgt, " ", .le(utf), " min(", mrg, ") \u222A ", tgt, " ", .ge(utf), " max(", mrg, ")"))
    else
      c(paste0(tgt, " ", .le(utf), " min(", mrg, ") or"),  paste0(tgt, " ", .ge(utf), " max(", mrg, ")"))
  } else if (alt == "two.one.sided" && (val.alt < min(val.null) || val.alt > max(val.null))) {
    if (utf)
      c(paste0("min(", mrg, ") ", .le(utf), " ", tgt, " ", .le(utf), " max(", mrg, ")"))
    else
      c(paste0(tgt, " ", .ge(utf), " min(", mrg, ") and"), paste0(tgt, " ", .le(utf), " max(", mrg, ")"))
  }
}

.h1_twoone <- function(tgt, mrg, utf, alt = "two.sided", val.alt = NA, val.null = NA, val.mrg = NULL) {
  if (alt %in% c("one.sided", "two.sided")) {
    paste(tgt, .h1_sign(alt, val.alt < val.null, utf), ifelse(!is.null(val.mrg) && val.mrg == 0, "0", mrg))
  } else if (alt == "two.one.sided" && (val.alt > min(val.null) && val.alt < max(val.null))) {
    if (utf)
      c(paste0("min(", mrg, ") < ", tgt, " < max(", mrg, ")"))
    else
      c(paste0(tgt, " ", .gt(utf), " min(", mrg, ") and"), paste0(tgt, " ", .ls(utf), " max(", mrg, ")"))
  } else if (alt == "two.one.sided" && (val.alt < min(val.null) || val.alt > max(val.null))) {
    if (utf)
      c(paste0(tgt, " ", .ls(utf), " min(", mrg, ") \u222A ", tgt, " ", .gt(utf), " max(", mrg, ")"))
    else
      c(paste0(tgt, " ", .ls(utf), " min(", mrg, ") or"),  paste0(tgt, " ", .gt(utf), " max(", mrg, ")"))
  }
}

.hypotheses <- function(h0_text, h1_text, utf = FALSE) {
  if (length(h0_text) > 1) h0_text <- paste(h0_text, collapse = paste0("\n", .pad("", 21)))
  if (length(h1_text) > 1) h1_text <- paste(h1_text, collapse = paste0("\n", .pad("", 21)))

  paste0(.topic("Hypotheses", utf),
         sprintf("  %s %s : %s\n",   ifelse(utf, "H\u2080", "H0"), "(Null)       ", h0_text),
         sprintf("  %s %s : %s\n\n", ifelse(utf, "H\u2081", "H1"), "(Alternative)", h1_text))

}

.c_a <- function(mark = FALSE, utf = FALSE) {
  # if the parameter is requested (mark = TRUE), c_a[1 / 2] and c_a[3 / 5] are used for coloring (in utf)
  # and empty for ascii; c_a[4] is used for <<
  if (mark && utf)
    c("\033[34m", "\033[0m", "  \033[1;35m", "\u25C4\u25C4", "\033[0m")
  else if (mark && !utf)
    c("",         "",        "  ",           "<<",           "")
  else
    rep("", 5)
}

# assembles / formats the "Effect Size / Target Effect / Target NCP" line
.esline <- function(x, utf = FALSE, digits = 3, a_pad = 0) {
  if (utils::hasName(x, "tgt.effect")) {
    es_sel <- names(x) %in% x[["tgt.effect"]]
    es_nme <- "Target Effect"
  } else if (utils::hasName(x, "tgt.ncp")) {
    es_sel <- names(x) %in% x[["tgt.ncp"]]
    es_nme <- "Target NCP"
  } else {
    es_sel <- names(x) %in% c("d", "eta.squared", "odds.ratio", "r.squared", "r.squared.change", "rate.ratio", "w")
    es_nme <- "Effect Size"
  }

  if (sum(es_sel) == 1) {
    # colors and arrows when requested, otherwise ""
    c_a <- .c_a(x$requested == "es", utf)

    # determine effect size and null
    es_fld  <- names(x)[es_sel]
    es_adsc <- .fmt_adsc(es_fld, utf)
    es_pad  <- .pad(es_adsc, ifelse(utf, 14, 17) + a_pad - nchar(es_nme))
    es_val  <- sprintf("%s", .fmt_val(x[[es_fld]], digits))
    if (paste0("null.", es_fld) %in% names(x)) {
      es_ndsc <- ifelse(utf, paste0(es_adsc, "\u2080"), paste0("null.", es_adsc))
      es_val  <- sprintf("%s (vs. %s = %s)", es_val, es_ndsc, .fmt_val(x[[paste0("null.", es_fld)]], digits))
    } else if (es_fld %in% c("prob1", "prob2", "prob10", "prob01", "rho1", "rho2", "rho12", "rho13", "rho34")) {
      if (es_fld %in% c("prob1", "prob2")) {
        es_nfld <- setdiff(c("prob1", "prob2"), es_fld)
      } else if (es_fld %in% c("prob10", "prob01")) {
        es_nfld <- setdiff(c("prob10", "prob01"), es_fld)
      } else if (es_fld %in% c("rho1", "rho2")) {
        es_nfld <- setdiff(c("rho1", "rho2"), es_fld)
      } else if (es_fld %in% c("rho12", "rho13", "rho34")) {
        es_nfld <- if (es_fld != "rho12") "rho12" else if (x$common) "rho13" else "rho34"
      }
      es_val  <- sprintf("%s (vs. %s = %s)", es_val, ifelse(utf, .fmt_utf(es_nfld), es_nfld), .fmt_val(x[[es_nfld]], digits))
    }
    sprintf("  %s%s (%s)%s = %s%s%s%s%s\n", c_a[1], es_nme, es_adsc, es_pad, es_val, c_a[2], c_a[3], c_a[4], c_a[5])
  } else if (sum(es_sel) == 0) {
    if (x$test != "Indirect Effect in a Mediation Model")
      stop("printing: no effect size")
    # return an empty string (for the generic functions, that only calculate power)
    ""
  } else {
    print(names(x))
    stop("`print.obj` contains more than one effect size.", call. = FALSE)
  }
}

# assembles "Sample Size" line
.nline <- function(x, utf = FALSE, digits = 3, a_pad = 0) {
  if (any(c("n", "n.total", "n.paired", "n.pres", "df", "size") %in% names(x))) {
    # colors and arrows when requested, otherwise ""
    c_a <- .c_a(x$requested == "n", utf)
    # assembles / formats the sample size line
    if (utils::hasName(x, "n.paired")) {
      n_name <- "Paired Sample Size"
      n <- x[["n.paired"]]
    } else if (utils::hasName(x, "n")) {
      n_name <- "Sample Size"
      n <- x[["n"]]
    } else if (utils::hasName(x, "n.total")) {
      n_name <- "Total Sample Size"
      n <- x[["n.total"]]
    } else if (utils::hasName(x, "n.pres")) {
      n_name <- "Presumed Sample S."
      n <- x[["n.pres"]]
    } else if (utils::hasName(x, "df")) {
      if (!is.finite(x[["df"]])) return("")
      n_name <- "Degrees of Freedom"
      n <- x[["df"]]
    } else if (utils::hasName(x, "size")) {
      n_name <- "Number of Trials"
      n <- x[["size"]]
    }
    n_pad    <- .pad(n_name, ifelse(utf, 17, 20) + a_pad)
    n_text   <- paste(round(n, digits), collapse = ifelse(length(n) == 2, " and ", ", "))
    sprintf("  %s%s%s = %s%s%s%s%s\n", c_a[1], n_name, n_pad, n_text, c_a[2], c_a[3], c_a[4], c_a[5])
  } else {
    # return an empty string (for the generic functions, that only calculate power)
    ""
  }
}

# assembles / formats the "Statistical Power" line
.pline <- function(x, utf = FALSE, digits = 3, a_pad = 0) {
  # colors and arrows when requested, otherwise ""
  c_a <- .c_a(x$requested == "power", utf)
  sprintf("  %sStatistical Power%s = %.*f%s%s%s%s\n\n", c_a[1], .pad("", ifelse(utf, 0, 3) + a_pad), digits, x$power, c_a[2], c_a[3], c_a[4], c_a[5])
}

.results <- function(x, utf = FALSE, digits = 3) {
  a_pad <- .a_pad(names(x), utf)
  paste0(.topic("Results", utf),
         .esline(x, utf, digits, a_pad),
         .nline(x,  utf, digits, a_pad),
         sprintf("  Type 1 Error %s%s = %.*f\n", ifelse(utf, "(\u03B1) ", "(alpha)"), .pad("", a_pad), digits, x$alpha),
         sprintf("  Type 2 Error %s%s = %.*f\n", ifelse(utf, "(\u03B2) ", "(beta) "), .pad("", a_pad), digits, 1 - x$power),
         .pline(x,  utf, digits, a_pad))
}

.defs <- function(defs_mtx, utf = FALSE) {
  defs_col <- ifelse(utf, 2, 1)
  defs_out <- .topic("Definitions", utf, emphasize = utf)
  emph_beg <- ifelse(utf, "\033[36m", "")
  emph_end <- ifelse(utf, "\033[0m",  "")

  for (r in seq_len(nrow(defs_mtx))) {
    dsc <- defs_mtx[[r, defs_col]]
    if (dsc == "") next
    defs_out <- c(defs_out, sprintf("  %s%s%s%s", emph_beg, gsub("\n$", "", dsc), emph_end, ifelse(grepl("\n$", dsc), "\n\n", "\n")))
  }

  paste0(paste(defs_out, collapse = ""), "\n")
}


# print functions ------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.ancova <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  if (utf) {
    h0_text <- ifelse(x$null.ncp == 0, "\u03B7\u00B2 = 0", "0 \u2264 \u03B7\u00B2 \u2264 \u03B7\u00B2\u2080")
    h1_text <- ifelse(x$null.ncp == 0, "\u03B7\u00B2 > 0", "\u03B7\u00B2 > \u03B7\u00B2\u2080")
  } else {
    h0_text <- ifelse(x$null.ncp == 0, "eta.squared = 0", "0 <= eta.squared <= null.eta.squared")
    h1_text <- ifelse(x$null.ncp == 0, "eta.squared > 0", "eta.squared > null.eta.squared")
  }
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.n. |  ascii           |  utf
    parms_mtx <- t(matrix(c("effect",   "Design",                 "Design",
                            "df1",      "Num. Deg. of Freedom",   "df1",
                            "df2",      "Denom. Deg. of Freedom", "df2",
                            "ncp",      "Non-centrality of Alt.", "\u03BB",
                            "null.ncp", "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "f.alpha",  "Critical Value",         "F\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii (utf at next line indented)
    defs_mtx <- t(matrix(c("null.eta.squared : (Partial) eta-squared under null",
                             "\u03B7\u00B2\u2080 : (Partial) eta-squared under null",
                           "eta.squared      : (Partial) eta-squared under alt.",
                             "\u03B7\u00B2\u2009 : (Partial) eta-squared under alternative",
                           "",
                             "\u03BB\u2009  : Non-centrality parameter under alternative",
                           "",
                             "\u03BB\u2080  : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))
    if (x$null.ncp == 0) defs_mtx <- defs_mtx[-1, ]

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.ancova() ----------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.contrast <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  h0_text <- ifelse(utf, "\u03C8 = 0",      "psi  = 0")
  h1_text <- ifelse(utf, "\u03C8 \u2260 0", "psi != 0")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name        |  ascii           |  utf
    parms_mtx <- t(matrix(c("psi",      "Contrast Est. (psi)",    "\u03C8",
                            "d",        "Standardized psi (d)",   "d",
                            "df",       "Degrees of Freedom",     "df",
                            "ncp",      "Non-centrality of Alt.", "\u03BB",
                            "null.ncp", "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "t.alpha",  "Critical Value",         "T\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii (utf at next line indented)
    defs_mtx <- t(matrix(c("psi : Contrast estimate, sum(contrast[i] * mu[i])",
                             "\u03C8\u2009 : Contrast est. defined as \u2211(contrast\u1D62 * \u03BC\u1D62)",
                           "d   : Standardized contrast estimate",
                             "d\u2009 : Standardized contrast estimate",
                           "",
                             "\u03BB\u2009 : Non-centrality parameter under alternative",
                           "",
                             "\u03BB\u2080 : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.contrast() --------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.contrasts <- function(x, data = NULL, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  h0_text <- ifelse(utf, "\u03C8 = 0",      "psi  = 0")
  h1_text <- ifelse(utf, "\u03C8 \u2260 0", "psi != 0")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    x$adjust.alpha <- switch(x$adjust.alpha,
                             `none` = "None", `fdr` = "False Discovery Rate", `hochberg` = "Hochberg (1988)",
                             `BH` = "Benjamini & Hochberg (1995)", `BY` = "Benjamini & Yakutieli (2001)",
                             `holm` = "Holm (1979)", `hommel` = "Hommel (1988)", `bonferroni` = "Bonferroni",
                             `tukey` = "Tukey")

    #                      | var.name        |  ascii                  |  utf
    parms_mtx <- t(matrix(c("null.ncp",     "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "alpha",        "Adjusted Alpha",         "Adjusted Type 1 Error (\u03B1)",
                            "adjust.alpha", "Alpha Adjustment",       "\u03B1 Adjustment"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.topic("Results", utf))
  if (!is.null(x$data)) {
    print(x$data, row.names = FALSE)
    cat("\n")
  }

  if (verbose == 2) {
    #                     | ascii (utf at next line indented)
    defs_mtx <- t(matrix(c("psi : Contrast estimate, sum(contrast[i] * mu[i])",
                             "\u03C8 (psi) : Contrast est. defined as \u2211(contrast\u1D62 * \u03BC\u1D62)",
                           "d   : Standardized contrast estimate",
                             "d       : Standardized contrast estimate",
                           "ncp : Non-centrality parameter under alt.",
                             "\u03BB\u2080      : Non-centrality parameter under alt."),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.contrasts() -------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.steiger <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")
  cat("  Common Index : ", x$common, "\n\n", sep = "")

  if (utf)
    tgt <- paste("\u03C1\u2081\u2082 - ", ifelse(x$common, "\u03C1\u2081\u2083", "\u03C1\u2083\u2084"))
  else
    tgt <- paste("rho12 -",               ifelse(x$common, "rho13",              "rho34"))
  h0_text <- paste(tgt, .h0_sign(x$alternative, x$delta < 0, utf), "0")
  h1_text <- paste(tgt, .h1_sign(x$alternative, x$delta < 0, utf), "0")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name  |  ascii          |  utf
    parms_mtx <- t(matrix(c("delta",     "rho12 - rho13",  "\u03C1\u2081\u2082 - \u03C1\u2081\u2083\u2009",
                            "delta",     "rho12 - rho34",  "\u03C1\u2081\u2082 - \u03C1\u2083\u2084\u2009",
                            "q",         "Cohen's q",      "Cohen's q",
                            "mean",      "Mean of Alt.",   "\u03BC (Alternative)",
                            "sd",        "SD of Alt.",     "\u03C3 (Alternative)",
                            "null.mean", "Mean of Null",   "\u03BC\u2080 (Null)\u2009",
                            "null.sd",   "SD of Null",     "\u03C3\u2080 (Null)\u2009",
                            "z.alpha",   "Critical Value", "Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx[-1 - as.integer(x$common), ], utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                          |  utf
    defs_mtx <- t(matrix(c("rho12 : Correlation between variable V1 and V2", "\u03C1\u2081\u2082 : Correlation between variable V1 and V2",
                           "rho13 : Correlation between variable V1 and V3", "\u03C1\u2081\u2083 : Correlation between variable V1 and V3",
                           "rho34 : Correlation between variable V3 and V4", "\u03C1\u2083\u2084 : Correlation between variable V3 and V4",
                           "",                                               "\u03BC\u2009  : Mean",
                           "",                                               "\u03C3\u2009  : Standard deviation"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx[-2 - as.integer(x$common), ], utf))
  }

} # .print.pwrss.steiger() ---------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.twocors <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {
  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  if (utf)
    tgt <- ifelse(x$design == "one.sample", "\u03C1 - \u03C1\u2080", "\u03C1\u2081 - \u03C1\u2082")
  else
    tgt <- ifelse(x$design == "one.sample", "rho - null.rho",        "rho1 - rho2")
  h0_text <- paste(tgt, .h0_sign(x$alternative, x$delta < 0, utf), "0")
  h1_text <- paste(tgt, .h1_sign(x$alternative, x$delta < 0, utf), "0")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name  |  ascii          |  utf
    parms_mtx <- t(matrix(c("delta",     "rho1 - rho2",    "\u03C1\u2081 - \u03C1\u2082\u2009",
                            "delta",     "rho - null.rho", "\u03C1\u2081 - \u03C1\u2080\u2009",
                            "q",         "Cohen's q",      "Cohen's q",
                            "mean",      "Mean of Alt.",   "\u03BC (Alternative)",
                            "sd",        "SD of Alt.",     "\u03C3 (Alternative)",
                            "null.mean", "Mean of Null",   "\u03BC\u2080 (Null)\u2009",
                            "null.sd",   "SD of Null",     "\u03C3\u2080 (Null)\u2009",
                            "z.alpha",   "Critical Value", "Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx[ifelse(x$design == "one.sample", -1, -2), ], utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                                |  utf
    defs_mtx <- t(matrix(c("rho1 : Correlation (for some V1 ~ V2) in group 1",     "\u03C1\u2081 : Correlation (for some V1 ~ V2) in group 1",
                           "rho2 : Correlation (for some V1 ~ V2) in group 2",     "\u03C1\u2082 : Correlation (for some V1 ~ V2) in group 2",
                           "rho      : Correlation (for some V1 ~ V2) under alt.", "\u03C1\u2081 : Correlation (for some V1 ~ V2) under alternative",
                           "null.rho : Correlation (for some V1 ~ V2) under null", "\u03C1\u2080 : Correlation (for some V1 ~ V2) under null",
                           "",                                                     "\u03BC\u2009 : Mean",
                           "",                                                     "\u03C3\u2009 : Standard deviation"),

                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx[ifelse(x$design == "one.sample", -1, -2) * 2 + c(1, 0), ], utf))
  }

} # .print.pwrss.twocors() ---------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.binom <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, TRUE, utf))
  cat(x$test, "\n\n", sep = "")

  tgt <- ifelse(utf, "P",       "prob")
  mrg <- ifelse(utf, "P\u2080", "null.prob")
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$prob, val.null = x$null.prob)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$prob, val.null = x$null.prob)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name    |  ascii                  |  utf
    parms_mtx <- t(matrix(c("size",        "Size",                   "Size",
                            "prob",        "Probability Under Alt.", "P",
                            "null.prob",   "Probability Under Null", "P\u2080\u2009",
                            "binom.alpha", "Critical Value",         "Bin\u207B\u00B9(\u03B1, P\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                        |  utf
    defs_mtx <- t(matrix(c("prob      : Probability under alt.", "P\u2009 : Probability under alternative",
                           "null.prob : Probability under null", "P\u2080 : Probability under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.binom() -----------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.chisq <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, TRUE, utf))
  cat(x$test, "\n\n", sep = "")

  if (utf) {
    h0_text <- ifelse(x$null.lambda == 0, "\u03BB = \u03BB\u2080", "0 \u2264 \u03BB \u2264 \u03BB\u2080")
    h1_text <- ifelse(x$null.lambda == 0, "\u03BB > \u03BB\u2080", "\u03BB > \u03BB\u2080")
  } else {
    h0_text <- ifelse(x$null.lambda == 0, "lambda = null.lambda",        "0 <= lambda <= null.lambda")
    h1_text <- ifelse(x$null.lambda == 0, "lambda > null.lambda",        "lambda > null.lambda")
  }
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name    |  ascii                  |  utf
    parms_mtx <- t(matrix(c("df",          "Degrees of Freedom",     "df",
                            "lambda",      "Non-centrality of Alt.", "\u03BB",
                            "null.lambda", "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "chisq.alpha", "Critical Value",         "Inv-\u03C7\u00B2(\u03B1, \u03BB\u2080)\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                          |  utf
    defs_mtx <- t(matrix(c("lambda      : Non-centrality parameter of alt.", "\u03BB\u2009 : Non-centrality parameter under alternative",
                           "null.lambda : Non-centrality parameter of null", "\u03BB\u2080 : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.chisq() -----------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.f <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, TRUE, utf))
  cat(x$test, "\n\n", sep = "")

  if (utf) {
    h0_text <- ifelse(x$null.lambda == 0, "\u03BB = \u03BB\u2080", "0 \u2264 \u03BB \u2264 \u03BB\u2080")
    h1_text <- ifelse(x$null.lambda == 0, "\u03BB > \u03BB\u2080", "\u03BB > \u03BB\u2080")
  } else {
    h0_text <- ifelse(x$null.lambda == 0, "lambda = 0", "0 <= lambda <= null.lambda")
    h1_text <- ifelse(x$null.lambda == 0, "lambda > 0", "lambda > null.lambda")
  }
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name    |  ascii                  |  utf
    parms_mtx <- t(matrix(c("df1",         "Num. Deg. of Freedom",   "df1",
                            "df2",         "Denom. Deg. of Freedom", "df2",
                            "lambda",      "Non-centrality of Alt.", "\u03BB",
                            "null.lambda", "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "f.alpha",     "Critical Value",         "F\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                       |  utf
    defs_mtx <- t(matrix(c("lambda      : Non-centrality parameter of alt.", "\u03BB\u2009 : Non-centrality parameter under alternative",
                           "null.lambda : Non-centrality parameter of null", "\u03BB\u2080 : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.f() ---------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.t <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, TRUE, utf))
  cat(x$test, "\n\n", sep = "")

  tgt <- ifelse(utf, "\u03BB",       "lambda")
  mrg <- ifelse(utf, "\u03BB\u2080", "null.lambda")
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$lambda, val.null = x$null.lambda)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$lambda, val.null = x$null.lambda)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name    |  ascii                  |  utf
    parms_mtx <- t(matrix(c("df",          "Degrees of Freedom",     "df",
                            "lambda",      "Non-centrality of Alt.", "\u03BB",
                            "null.lambda", "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "t.alpha",     "Critical Value",         "T\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))
    if (is.infinite(x$df)) parms_mtx <- parms_mtx[-1, ]

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                       |  utf
    defs_mtx <- t(matrix(c("lambda      : Non-centrality parameter of alt.", "\u03BB\u2009 : Non-centrality parameter under alternative",
                           "null.lambda : Non-centrality parameter of null", "\u03BB\u2080 : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.t() ---------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.z <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, TRUE, utf))
  cat(x$test, "\n\n", sep = "")

  tgt <- ifelse(utf, "\u03BC",       "mean")
  mrg <- ifelse(utf, "\u03BC\u2080", "null.mean")
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$mean, val.null = x$null.mean)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$mean, val.null = x$null.mean)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name  |  ascii                  |  utf
    parms_mtx <- t(matrix(c("mean",      "Mean of Alt.",   "\u03BC",
                            "null.mean", "Mean of Null",   "\u03BC\u2080\u2009",
                            "z.alpha",   "Critical Value", "Z\u207B\u00B9(\u03B1, \u03BC\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                        |  utf
    defs_mtx <- t(matrix(c("mean      : Mean of alt.", "\u03BC\u2009 : Mean under alternative",
                           "null.mean : Mean of null", "\u03BC\u2080 : Mean under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.z() ---------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.ttest <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {


  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  tgt <- ifelse(utf, "d - d\u2080", "d - null.d")
  mrg <- ifelse(utf, "\u03B4",      "margin")
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$ncp, val.null = x$null.ncp, val.mrg = x$margin)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$ncp, val.null = x$null.ncp, val.mrg = x$margin)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name |  ascii                  |  utf
    parms_mtx <- t(matrix(c("margin",   "Margin",                 "\u03B4",
                            "df",       "Degrees of Freedom",     "df",
                            "ncp",      "Non-centrality of Alt.", "\u03BB",
                            "null.ncp", "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "t.alpha",  "Critical Value",         "T\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                                |  utf
    defs_mtx <- t(matrix(c("",                                                     "d\u2009 : Cohen's d under alternative",
                           "",                                                     "d\u2080 : Cohen's d under null",
                           "Margin : Smallest d - null.d difference that matters", "\u03B4\u2009 : Margin - ignorable d - d\u2080 difference",
                           "",                                                     "\u03BB\u2009 : Non-centrality parameter under alternative",
                           "",                                                     "\u03BB\u2080 : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.ttest() -----------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.wilcoxon <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  dist <- switch(x$dist, `normal` = "Normal", `uniform` = "Uniform", `double.exponential` = "Double Exponential",
                         `laplace` = "Laplace", `logistic` = "Logistic")
  method <- switch(x$method, `guenther` = "Guenther", `noether` = "Noether")
  cat(x$test, "\n\n", sep = "")
  cat("  Method       : ", method, "\n", sep = "")
  cat("  Distribution : ", dist, "\n\n", sep = "")

  tgt <- ifelse(utf, "d - d\u2080", "d - null.d")
  mrg <- ifelse(utf, "\u03B4",      "margin")
  val.alt  <- ifelse(x$method == "guenther", x$ncp,      x$mean)
  val.null <- ifelse(x$method == "guenther", x$null.ncp, x$null.mean)
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = val.alt, val.null = val.null, val.mrg = x$margin)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = val.alt, val.null = val.null, val.mrg = x$margin)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #               | var.name  |  ascii                  |  utf
    parms_vec   <- c("margin",    "Margin",                 "\u03B4")
    if (x$method == "guenther") {
      #             | var.name  |  ascii                  |  utf
      parms_vec <- c(parms_vec,
                     "df",        "Degrees of Freedom",     "df",
                     "ncp",       "Non-centrality of Alt.", "\u03BB",
                     "null.ncp",  "Non-centrality of Null", "\u03BB\u2080\u2009",
                     "t.alpha",   "Critical Value",         "T\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009")
    } else if (x$method == "noether") {
      #             | var.name  |  ascii          |  utf
      parms_vec <- c(parms_vec,
                     "mean",      "Mean of Alt.",   "\u03BC",
                     "null.mean", "Mean of Null",   "\u03BC\u2080\u2009",
                     "z.alpha",   "Critical Value", "Z\u207B\u00B9(\u03B1, \u03BC\u2080)\u2009\u2009")
    }

    parms_mtx <- t(matrix(parms_vec, nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))
    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                         |  utf
    defs_vec   <- c("",                                                     "d\u2009 : Cohen's d under alternative",
                    "",                                                     "d\u2080 : Cohen's d under null",
                    "Margin : Smallest d - null.d difference that matters", "\u03B4\u2009 : Margin - ignorable d - d\u2080 difference")
    if (x$method == "guenther") {
      #            | utf (ascii is "" to keep the correct shape)
      defs_vec <- c(defs_vec,
                    "", "\u03BB\u2009 : Non-centrality parameter under alternative",
                    "", "\u03BB\u2080 : Non-centrality parameter under null")
    } else if (x$method == "noether") {
      #            | utf (ascii is "" to keep the correct shape)
      defs_vec <- c(defs_vec,
                    "", "\u03BC\u2009 : Mean of the alternative distribution",
                    "", "\u03BC\u2080 : Mean of the null distribution")
    }

    defs_mtx <- t(matrix(defs_vec, nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))
    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.wilcoxon() --------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.fisher <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")
  cat("  Method : ", switch(x$method, `z` = "Normal Approximation", `exact` = "Fisher's Exact"), "\n\n", sep = "")

  tgt <- ifelse(utf, "P\u2081 - P\u2082", "prob1 - prob2")
  mrg <- ifelse(utf, "\u03B4",            "margin")
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$delta, val.null = x$margin, val.mrg = x$margin)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$delta, val.null = x$margin, val.mrg = x$margin)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #               | var.name   |  ascii            |  utf
    parms_vec <-   c("delta",       "prob1 - prob2",   "P\u2081 - P\u2082",
                     "odds.ratio",  "Odds Ratio (OR)", "Odds Ratio (OR)")
    if (x$method == "z") {
      #             | var.name  |   ascii            |  utf
      parms_vec <- c(parms_vec,
                     "mean",       "Mean of Alt.",     "\u03BC",
                     "null.mean",  "Mean of Null",     "\u03BC\u2080\u2009",
                     "z.alpha",    "Critical Value",   "Z\u207B\u00B9(\u03B1, \u03BC\u2080)\u2009")
    }
    parms_mtx <- t(matrix(parms_vec, nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #              | ascii (utf at next line indented)
    defs_vec <-   c("prob1       : Probability of success in the first group",
                      "P\u2081       : Probability of success in the first group",
                    "prob2       : Probability of success in the second group",
                      "P\u2082       : Probability of success in the second group",
                    "margin      : Smallest prob1 - prob2 that matters",
                      "\u03B4\u2009\u2009       : Margin - ignorable P\u2081 - P\u2082 difference",
                    "Odds Ratio  : Odds(prob1) / Odds(prob2)",
                      "OR       : Odds(P\u2081) / Odds(P\u2082)",
                    "Odds(prob1) : prob1 / (1 - prob1)",
                      "Odds(P\u2081) : P\u2081 / (1 - P\u2081)",
                    "Odds(prob2) : prob2 / (1 - prob2)",
                      "Odds(P\u2082) : P\u2082 / (1 - P\u2082)")
    if (x$method == "z" && utf) {
      #            |     utf (ascii is not required and thus empty - "")
      defs_vec <- c(defs_vec,
                    "", "\u03BC\u2009       : Mean of the alternative distribution",
                    "", "\u03BC\u2080       : Mean of the null distribution")
    }
    defs_mtx <- t(matrix(defs_vec, nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))
    if (all(x$margin == 0)) defs_mtx <- defs_mtx[-3, ]

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.fisher() ----------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.gof <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  h0_text <- ifelse(utf,  "P[i,j] = P\u2080[i,j] for \u2200(i,j)",      "P[i,j]  = P0[i,j] for all (i,j)")
  h1_text <- ifelse(utf,  "P[i,j] \u2260 P\u2080[i,j] for \u2203(i,j)", "P[i,j] != P0[i,j] for some (i,j)")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name    |  ascii                  |  utf
    parms_mtx <- t(matrix(c("df",          "Degrees of Freedom",     "df",
                            "ncp",         "Non-centrality of Alt.", "\u03BB",
                            "null.ncp",    "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "chisq.alpha", "Critical Value",         "Inv-\u03C7\u00B2(\u03B1, \u03BB\u2080)\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                                      |  utf
    defs_mtx <- t(matrix(c("For goodness-of-fit, comparisons are P[i] vs P0[i]",         "For goodness-of-fit test, it is P[i] vs P\u2080[i]",
                           "For independence, comparisons are P[i,j] vs P0[i,j]",        "",
                           "Independence implies (default) P0[i,j] = P[i,.] * P[.,j]\n", "Independence implies P\u2080[i,j] = P[i,.] * P[.,j]",
                           "P[i,j] : Joint probability for cell (i,j)",                  "P[i,j] : Joint prob. for cell (i,j)",
                           "P[i,.] : Marginal probability for row i (sum over j)",       "P[i,.] : Marginal prob. for row i (sum over j)",
                           "P[.,j] : Marginal probability for column j (sum over i)",    "P[.,j] : Marginal prob. for column j (sum over i)\n",
                           "",                                                           "\u03BB      : Non-centrality parameter under alternative",
                           "",                                                           "\u03BB\u2080\u2009    : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.gof() -------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.mcnemar <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")
  cat("  Method : ", switch(x$method, `z` = "Normal Approximation", `exact` = "McNemar's Exact"), "\n\n", sep = "")

  tgt <- ifelse(utf, "P\u2081\u2080 - P\u2080\u2081", "prob10 - prob01")
  h0_text <- paste(tgt, .h0_sign(x$alternative, x$delta < 0, utf), "0")
  h1_text <- paste(tgt, .h1_sign(x$alternative, x$delta < 0, utf), "0")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #               | var.name   |   ascii                |  utf
    parms_vec <-   c("odds.ratio",  "Odds Ratio",           "Odds Ratio (OR)",
                     "delta",       "prob10 - prob01",      "P\u2081\u2080 - P\u2080\u2081")
    if (x$method == "exact") {
      #             | var.name    |  ascii                |  utf
      parms_vec <- c(parms_vec,
                     "size",        "Size of Disc. Pairs",  "Size (Discordant)",
                     "prob",        "prob10 | DP for Alt.", "P\u2081\u2009",
                     "null.prob",   "prob10 | DP for Null", "P\u2080\u2009",
                     "binom.alpha", "Critical Value",       "Bin\u207B\u00B9(\u03B1, P\u2080)\u2009\u2009")
    } else {
      #             | var.name  |    ascii                |  utf
      parms_vec <- c(parms_vec,
                     "mean",        "Mean of Alt.",         "\u03BC\u2081\u2009",
                     "null.mean",   "Mean of Null",         "\u03BC\u2080\u2009",
                     "z.alpha",     "Critical Value",       "Z\u207B\u00B9(\u03B1, \u03BC\u2080)\u2009\u2009")
    }
    parms_mtx <- t(matrix(parms_vec, nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #              | ascii                                              | utf
    defs_vec <-   c("Odds Ratio  : prob10 / prob01",                      "OR : P\u2081\u2080 / P\u2080\u2081",
                    "prob10      : Joint prob. of observing {1,0}",       "",
                    "prob01      : Joint prob. of observing {0,1}",       "",
                    "prob10 | DP : Conditional prob. of observing {1,0}", "",
                    "              among DP, prob10 / (prob10 + prob01)", "",
                    "DP          : Discordant pairs", "")
    if (x$method == "exact") {
      #            | ascii | utf
      defs_vec <- c(defs_vec,
                    "",     "P\u2081\u2009 : Prob. of {1,0} among discordant pairs under alt.",
                    "",     "P\u2080\u2009 : Prob. of {1,0} among discordant pairs under null")
    } else {
      #            | ascii | utf
      defs_vec <- c(defs_vec,
                    "",     "\u03BC\u2081\u2009 : Mean of the alternative distribution",
                    "",     "\u03BC\u2080\u2009 : Mean of the null distribution")
    }
    defs_mtx <- t(matrix(defs_vec, nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.mcnemar() ---------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.oneprop <- function(x, digits = 3, verbose = 1, utf = TRUE, ...) {

  cat(.header(x$requested, FALSE, utf))
  method <- switch(x$method, `z` = "Normal Approximation", `exact` = "Exact")
  cat(x$test, "\n\n", sep = "")
  if (x$method == "exact") {
    cat("  Method                 : ", method, "\n\n", sep = "")
  } else {
    stderr <- switch(x$std.err, `alternative` = "Alternative", `null` = "Null")
    cat("  Method                 : ", method, "\n", sep = "")
    cat("  Continuity Correction  : ", x$correct, "\n", sep = "")
    cat("  Arcsine Transformation : ", x$arcsine, "\n", sep = "")
    cat("  Standard Error         : Calculated From ", stderr, "\n\n", sep = "")
  }

  tgt <- ifelse(utf, "P",       "prob")
  mrg <- ifelse(utf, "P\u2080", "null.prob")
  val.alt  <- ifelse(x$method == "exact", x$prob,      x$mean)
  val.null <- ifelse(x$method == "exact", x$null.prob, x$null.mean)
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = val.alt, val.null = val.null)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = val.alt, val.null = val.null)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #               | var.name   |   ascii            |  utf
    parms_vec <-   c("delta",       "prob - null.prob", "P - P\u2080\u2009",
                     "odds.ratio",  "Odds Ratio (OR)",  "Odds Ratio (OR)")
    if (x$method == "exact") {
      #             | var.name    |  ascii            |  utf
      parms_vec <- c(parms_vec,
                     "size",             "Size",             "Size",
                     "prob",        "Prob. Under Alt",  "P",
                     "null.prob",   "Prob. Under Null", "P\u2080\u2009",
                     "binom.alpha", "Critical Value",   "Bin\u207B\u00B9(\u03B1, P\u2080)")
    } else {
      #             | var.name  |    ascii          |    utf
      parms_vec <- c(parms_vec,
                     "mean",        "Mean of Alt.",     "\u03BC\u2081\u2009",
                     "null.mean",   "Mean of Null",     "\u03BC\u2080\u2009",
                     "z.alpha",     "Critical Value",   "Z\u207B\u00B9(\u03B1, \u03BC\u2080)")
    }
    parms_mtx <- t(matrix(parms_vec, nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #              | ascii                                          |  utf
    defs_vec <-   c("Odds Ratio      : Odds(prob) / Odds(null.prob)", "OR       : Odds(P) / Odds(P\u2080)",
                    "Odds(prob)      : prob / (1 - prob)",            "Odds(P)  : P / (1 - P)",
                    "Odds(null.prob) : null.prob / (1 - null.prob)",  "Odds(P\u2080)\u2009: P\u2080 / (1 - P\u2080)")
    if (x$method == "exact") {
      #            | ascii | utf
      defs_vec <- c(defs_vec,
                    "",     "P        : Probability of success under alternative",
                    "",     "P\u2080\u2009      : Probability of success under null")
    } else {
      #            | ascii | utf
      defs_vec <- c(defs_vec,
                    "",     "\u03BC\u2081\u2009      : Mean of the alternative distribution",
                    "",     "\u03BC\u2080\u2009      : Mean of the null distribution")
    }
    defs_mtx <- t(matrix(defs_vec, nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.oneprop() ---------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.f.regression <- function(x, digits = 3, verbose = 1, utf = FALSE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  rsq <- if (x$k.tested < x$k.total) ifelse(utf, "\u0394R\u00B2", "Change in R-squared") else ifelse(utf, "R\u00B2", "R-squared")
  h0_text <- ifelse(x$margin == 0, sprintf("%s = 0", rsq), sprintf("0 %s %s %s margin", .le(utf), rsq, .le(utf)))
  h1_text <- ifelse(x$margin == 0, sprintf("%s > 0", rsq), sprintf("%s > margin",                    rsq))
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name          | ascii                  | utf
    parms_mtx <- t(matrix(c("r.squared.change", rsq,                      rsq,
                            "margin",           "Margin",                 "\u03B4",
                            "df1",              "Num. Deg. of Freedom",   "df1",
                            "df2",              "Denom. Deg. of Freedom", "df2",
                            "ncp",              "Non-centrality of Alt.", "\u03BB",
                            "null.ncp",         "Non-centrality of Null", "\u03BB\u2080\u2009",
                            "f.alpha",          "Critical Value",         "F\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    mrg_def_ascii <- sprintf("Margin : Smallest %s that matters", gsub("Change", "change", rsq, fixed = TRUE))
    #                     | ascii       |  utf
    defs_mtx <- t(matrix(c(mrg_def_ascii, "\u03B4\u2009 : Margin - ignorable R\u00B2 or \u0394R\u00B2",
                           "",            "\u03BB\u2009 : Non-centrality parameter under alternative",
                           "",            "\u03BB\u2080 : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.f.regression() ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.t.regression <- function(x, digits = 3, verbose = 1, utf = FALSE, ...) {

  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")

  tgt <- ifelse(utf, "\u03B2 - \u03B2\u2080", "beta - null.beta")
  mrg <- ifelse(utf, "\u03B4",                "margin")
  h0_text <- .h0_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$ncp, val.null = x$null.ncp, val.mrg = x$margin)
  h1_text <- .h1_twoone(tgt, mrg, utf, alt = x$alternative, val.alt = x$ncp, val.null = x$null.ncp, val.mrg = x$margin)
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name      | ascii                 | utf
    parms_mtx <- t(matrix(c("std.beta",      "Std. Beta Under Alt.",   "Std. \u03B2",
                            "std.null.beta", "Std. Beta Under Null",   "Std. \u03B2\u2080\u2009",
                            "std.margin",    "Std. Margin",            "Std. \u03B4",
                            "df",            "Degrees of Freedom",     "df",
                            "ncp",           "Non-centrality of Alt.", "\u03BB (Alternative)",
                            "null.ncp",      "Non-centrality of Null", "\u03BB\u2080 (Null)\u2009",
                            "t.alpha",       "Critical Value",         "T\u207B\u00B9(\u03B1, \u03BB\u2080)\u2009\u2009"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii (utf at next line indented)
    defs_mtx <- t(matrix(c("beta                 : Regression coefficient under alt.",
                             "\u03B2       : Regression coefficient under alternative",
                           "null.beta            : Regression coefficient under null",
                             "\u03B2\u2080\u2009     : Regression coefficient under null",
                           "margin               : Smallest beta - null.beta difference that matters\n",
                             "\u03B4       : Margin(s) - ignorable \u03B2 - \u03B2\u2080 difference\n",
                           "Std. Beta Under Alt. = beta * [SD(X) / SD(Y)]",
                             "Std. \u03B2  = \u03B2  * [\u03C3(X) / \u03C3(Y)]",
                           "Std. Beta Under Null = null.beta * [SD(X) / SD(Y)]",
                             "Std. \u03B2\u2080\u2009= \u03B2\u2080 * [\u03C3(X) / \u03C3(Y)]",
                           "Std. Margin          = margin * [SD(X) / SD(Y)]\n",
                             "Std. \u03B4  = \u03B4  * [\u03C3(X) / \u03C3(Y)]\n",
                           "SD(X)                : Standard deviation of the predictor",
                             "\u03C3(X)    : Standard devition of the predictor",
                           "SD(Y)                : Standard deviation of the outcome",
                             "\u03C3(Y)    : Standard devition of the outcome\n",
                           "",
                             "\u03BB       : Non-centrality parameter under alternative",
                           "",
                             "\u03BB\u2080\u2009     : Non-centrality parameter under null"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.t.regression() ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.logistic <- function(x, digits = 3, verbose = 1, utf = FALSE, ...) {
  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")
  cat("  Method          : ", x$method, "\n", sep = "")
  cat("  Predictor Dist. : ", x$dist, "\n\n", sep = "")

  h0_text <- paste("Odds Ratio (OR)", .h0_sign(x$alternative, x$odds.ratio < 1, utf), "1")
  h1_text <- paste("Odds Ratio (OR)", .h1_sign(x$alternative, x$odds.ratio < 1, utf), "1")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name   |  ascii              |  utf
    parms_mtx <- t(matrix(c("base.prob",  "Base Probability",   "Base Probability",
                            "odds.ratio", "Odds Ratio (OR)",    "Odds Ratio (OR)",
                            "vcf",        "Var. Corr. Factor",  "Var. Corr. Factor",
                            "mean",       "Mean (Alternative)", "\u03BC (Alternative)",
                            "sd",         "SD (Alternative)",   "\u03C3 (Alternative)",
                            "null.mean",  "Mean (Null)",        "\u03BC\u2080 (Null)\u2009",
                            "null.sd",    "SD (Null)",          "\u03C3\u2080 (Null)\u2009",
                            "z.alpha",    "Critical Value",     "Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                                                    |  utf
    defs_mtx <- t(matrix(c("Odds Ratio = [prob/(1-prob)] / [base.prob/(1-base.prob)]", "OR = [p / (1 \u2212 p)] / [p\u2080 / (1 \u2212 p\u2080)]",
                           "prob       : Base probability when predictor = 0",         "p\u2080\u2009 : Base probability when predictor = 0",
                           "base.prob  : Probability when predictor = 1",              "p  : Probability when predictor = 1",
                           "beta1      = log(Odds Ratio)",                             "\u03B2\u2081\u2009 = log(OR)",
                           "beta0      = log[base.prob/(1-base.prob)]",                "\u03B2\u2080\u2009 = log[p\u2080 / (1 \u2212 p\u2080)]",
                           "",                                                         "\u03BC  : Mean",
                           "",                                                         "\u03C3  : Standard deviation"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.ascii.pwrss.logistic() --------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.med <- function(x, digits = 3, verbose = 1, utf = FALSE, ...) {

  cat(.header(x$requested, FALSE, utf))
  method <- switch(x$method, `sobel` = "Sobel", `aroian` = "Aroian", `goodman` = "Goodman", `joint` = "Joint", `monte.carlo` = "Monte Carlo")
  cat(x$test, "\n\n", sep = "")
  cat("  Method : ", method, "\n\n", sep = "")

  tgt <- ifelse(utf, "\u03B2[a*b]", "beta[a*b]")
  h0_text <- paste(tgt, .h0_sign(x$alternative, x$std.beta.indirect < 0, utf), "0")
  h1_text <- paste(tgt, .h1_sign(x$alternative, x$std.beta.indirect < 0, utf), "0")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #               | var.name           | ascii           | utf
    parms_vec   <- c("std.beta.a",        "Std. beta[a]",   "Std. \u03B2[a]",
                     "std.beta.b",        "Std. beta[b]",   "Std. \u03B2[b]",
                     "std.beta.indirect", "Std. beta[a*b]", "Std. \u03B2[a*b]")
    if (x$method %in% c("sobel", "aorian", "goodman")) {
      #             | var.name           | ascii           | utf
      parms_vec <- c(parms_vec,
                     "mean",              "Mean of Alt.",   "\u03BC (Alternative)",
                     "null.mean",         "Mean of Null.",  "\u03BC\u2080 (Null)\u2009",
                     "z.alpha",           "Critical Value", "Z\u207B\u00B9(\u03B1, \u03BC\u2080)\u2009\u2009")
    }

    parms_mtx <- t(matrix(parms_vec, nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))
    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #              | ascii (utf at next line indented)
    defs_vec   <- c("beta[a]        : Regression coefficient for path `a`",
                      "\u03B2[a]         : Regression coefficient for path `a`",
                    "beta[b]        : Regression coefficient for path `b`",
                      "\u03B2[b]         : Regression coefficient for path `b`",
                    "beta[a*b]      : Coefficient for the indirect path `a*b`\n",
                      "\u03B2[a*b]       : Coefficient for the indirect path `a*b`\n",
                    "Std. beta[a]   = beta[a] * [SD(predictor) / SD(mediator)]",
                      "Std. \u03B2[a]    = \u03B2[a] * [\u03C3(predictor) / \u03C3(mediator)]",
                    "Std. beta[b]   = beta[b] * [SD(mediator) / SD(outcome)]",
                      "Std. \u03B2[b]    = \u03B2[b] * [\u03C3(mediator) / \u03C3(outcome)]",
                    "Std. beta[a*b] = Std. beta[a] * Std. beta[b]\n",
                      "Std. \u03B2[a*b]  = Std. \u03B2[a] * Std. \u03B2[b]\n",
                    "SD(predictor)  : Standard deviation of the predictor",
                      "\u03C3(predictor) : Standard devition of the predictor",
                    "SD(mediator)   : Standard deviation of the mediator",
                      "\u03C3(mediator)  : Standard devition of the mediator",
                    "SD(outcome)    : Standard deviation of the outcome",
                      "\u03C3(outcome)   : Standard devition of the outcome")

    if (x$method %in% c("sobel", "aorian", "goodman") && utf) {
      #            | utf (ascii is "" to keep the correct shape)
      defs_vec <- c(defs_vec[seq(length(defs_vec) - 1)], paste0(defs_vec[length(defs_vec)], "\n"),
                    "", "\u03BC            : Mean of the alternative distribution",
                    "", "\u03BC\u2080\u2009          : Mean of the null distribution")
    }

    defs_mtx <- t(matrix(defs_vec, nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))
    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.med() -------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
.print.pwrss.poisson <- function(x, digits = 3, verbose = 1, utf = FALSE, ...) {
  cat(.header(x$requested, FALSE, utf))
  cat(x$test, "\n\n", sep = "")
  cat("  Method          : ", x$method, "\n", sep = "")
  cat("  Predictor Dist. : ", x$dist, "\n\n", sep = "")

  h0_text <- paste("Rate Ratio (RR)", .h0_sign(x$alternative, x$rate.ratio < 1, utf), "1")
  h1_text <- paste("Rate Ratio (RR)", .h1_sign(x$alternative, x$rate.ratio < 1, utf), "1")
  cat(.hypotheses(h0_text, h1_text, utf))

  if (verbose == 2) {
    #                      | var.name   |  ascii               | utf
    parms_mtx <- t(matrix(c("base.rate",  "Base Rate",          "Base Rate",
                            "rate.ratio", "Rate Ratio (RR)",    "Rate Ratio (RR)",
                            "vcf",        "Var. Corr. Factor",  "Var. Corr. Factor",
                            "mean",       "Mean (Alternative)", "\u03BC (Alternative)",
                            "sd",         "SD (Alternative)",   "\u03C3 (Alternative)",
                            "null.mean",  "Mean (Null)",        "\u03BC\u2080 (Null)\u2009",
                            "null.sd",    "SD (Null)",          "\u03C3\u2080 (Null)\u2009",
                            "z.alpha",    "Critical Value",     "Z\u207B\u00B9(\u03B1, \u03BC\u2080, \u03C3\u2080)"),
                          nrow = 3, dimnames = list(c("var.name", "ascii", "utf"), NULL)))

    cat(.keyparms(x, parms_mtx, utf, digits))
  }

  cat(.results(x, utf, digits))

  if (verbose == 2) {
    #                     | ascii                        |  utf
    defs_mtx <- t(matrix(c("Base Rate       = exp(beta0)", "Base Rate       = exp(\u03B2\u2080)",
                           "Rate Ratio (RR) = exp(beta1)", "Rate Ratio (RR) = exp(\u03B2\u2081)",
                           "",                             "\u03BC               : Mean",
                           "",                             "\u03C3               : Standard deviation"),
                  nrow = 2, dimnames = list(c("ascii", "utf"), NULL)))

    cat(.defs(defs_mtx, utf))
  }

} # .print.pwrss.poisson() ---------------------------------------------------------------------------------------------
