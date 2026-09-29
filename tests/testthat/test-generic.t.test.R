test_that("generic.t.test.R works", {
    # power.t.test -----------------------------------------------------------------------------------------------------
    expect_equal(power.t.test(ncp = 1.96, df = 100, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.49254726, ncp = 1.96, null.ncp = 0, df = 100, alpha = 0.05,
                                alternative = "two.sided", t.alpha = 1.98397152 * c(-1, 1), beta = 0.50745274,
                                type.s = 0.0000970993, type.m = 1.427979969),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(ncp = 1.96, df = 120, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.619466, ncp = 1.96, null.ncp = 0, df = 120, alpha = 0.05, alternative = "one.sided",
                                t.alpha = 1.6576509, beta = 0.38053400, type.s = 0, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(ncp = 0, null.ncp = c(-2, 2), df = 100, alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.277047649, ncp = 0, null.ncp = c(-2, 2), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 0.35551537 * c(-1, 1), beta = 0.722952351,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(ncp = 0, null.ncp = 2, df = 100, alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.277047649, ncp = 0, null.ncp = c(-2, 2), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 0.35551537 * c(-1, 1), beta = 0.722952351,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(ncp = 2, null.ncp = c(-1, 1), df = 400, alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.167168042, ncp = 2, null.ncp = c(-1, 1), df = 400, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 2.97258983 * c(-1, 1), beta = 0.832831958,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(ncp = 2, null.ncp = -1, df = 400, alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.167168042, ncp = 2, null.ncp = c(-1, 1), df = 400, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 2.97258983 * c(-1, 1), beta = 0.832831958,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    # example 25.3 from the GPower manual
    expect_equal(power.t.test(ncp = -1.25, null.ncp = 0, df = 24, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.22452404, ncp = -1.25, null.ncp = 0, df = 24, alpha = 0.05, alternative = "two.sided",
                                t.alpha = 2.063898562 * c(-1, 1), beta = 0.775475958, type.s = 0.0034900921, type.m = 2.1835685),
                           class = c("pwrss", "generic", "t")))
    # result is identical: power ~ 0.224524/5

    expect_equal(power.t.test(power = 0.80, df = 100, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 2.8288472, null.ncp = 0, df = 100, alpha = 0.05,
                                alternative = "two.sided", t.alpha = 1.98397152 * c(-1, 1), beta = 0.200000001,
                                type.s = 0.00000119776692, type.m = 1.1334521),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(power = 0.80, df = 100, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 2.50347, null.ncp = 0, df = 100, alpha = 0.05, alternative = "one.sided",
                                t.alpha = 1.66023433, beta = 0.2, type.s = 0, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(power = 0.80, req.sign = "-", df = 100, alpha = 0.05, alternative = "one.sided",
                              plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = -2.50347, null.ncp = 0, df = 100, alpha = 0.05, alternative = "one.sided",
                                t.alpha = -1.66023433, beta = 0.2, type.s = 0, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(suppressWarnings(power.t.test(power = 0.80, req.sign = "0", null.ncp = c(-3, 3), df = 100, alpha = 0.05,
                                               alternative = "two.one.sided", plot = FALSE, verbose = 0)),
                 structure(list(power = 0.82026598, ncp = 0, null.ncp = c(-3, 3), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 1.35103226 * c(-1, 1), beta = 0.17973402,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(power = 0.80, req.sign = "+", null.ncp = c(-1, 1), df = 100, alpha = 0.05,
                              alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 3.8638821, null.ncp = c(-1, 1), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 3.0109934 * c(-1, 1), beta = 0.2,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(power = 0.80, ncp = 3, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 3, null.ncp = 0, df = 14.8520054, alpha = 0.05, alternative = "two.sided",
                                t.alpha = 2.1333009 * c(-1, 1), beta = 0.199999989, type.s = 0.00000106692622,
                                type.m = 1.186699976),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(power = 0.80, ncp = 3, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.800000662, ncp = 3, null.ncp = 0, df = 4.2402269, alpha = 0.05,
                                alternative = "one.sided", t.alpha = 2.09761361, beta = 0.199999338, type.s = 0, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(power = 0.80, ncp = 0, null.ncp = c(-2, 2), alpha = 0.05,
                              alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.277520063, ncp = 0, null.ncp = c(-2, 2), df = 9955311578.2, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 0.355146373 * c(-1, 1), beta = 0.722479937,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))
    expect_equal(power.t.test(power = 0.80, ncp = 4, null.ncp = c(-1, 1), alpha = 0.05,
                              alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 4, null.ncp = c(-1, 1), df = 32.409657, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 3.12206879 * c(-1, 1), beta = 0.199999993,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "t")))

    expect_error(power.t.test(power = 0.80, alpha = 0.05, alternative = "two.sided"),
                 "Exactly two of the parameters `ncp`, `df`, or `power` must be given, one has to be NULL.")
    expect_error(power.t.test(power = 0.8, df = 2, alpha = 0.05, alternative = "two.sided"),
                 "Degrees of freedom can not be smaller than 3.")
    expect_error(power.t.test(ncp = 2, null.ncp = 0, df = 0, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 "Argument `df` does not have a valid positive value \\(must be length 1, numeric, > 0, and finite\\).")
    expect_error(power.t.test(ncp = 2, null.ncp = c(-1, 1), df = 100, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 "If `alternative` is \"two.sided\" or \"one.sided\", `null.ncp` must be of length one.")
    expect_error(power.t.test(ncp = 2, null.ncp = c(-1, 1, 1), df = 100, alpha = 0.05, alternative = "two.one.sided",
                              plot = FALSE, verbose = 0),
                 paste("If `alternative` is \"two.one.sided\", `null.ncp` must be of length one \\(absolute value,",
                       "that must be different from 0\\), or length two \\(lower and upper bounds\\)."))
    expect_error(power.t.test(ncp = 2, alpha = 0.05, df = 100, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 paste("If `alternative` is \"two.one.sided\", `null.ncp` must be of length one \\(absolute value,",
                       "that must be different from 0\\), or length two \\(lower and upper bounds\\)."))
    expect_warning(power.t.test(power = 0.80, req.sign = "0", null.ncp = c(-2, 2), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 "The target power rate cannot be achieved within the null bounds.")
})
