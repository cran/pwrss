test_that("generic.lp.test.R works", {
    # power.lp.test ----------------------------------------------------------------------------------------------------
    expect_equal(power.lp.test(ncp = 1.96, df = 100, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.49812069, ncp = 1.96, null.ncp = 0, df = 100, alpha = 0.05,
                                alternative = "two.sided", t.alpha = 1.959964 * c(-1, 1), beta = 0.50187931,
                                type.s = 0.000104670166, type.m = 1.410026951),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(ncp = 1.96, df = 100, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.62069048, ncp = 1.96, null.ncp = 0, df = 100, alpha = 0.05,
                                alternative = "one.sided", t.alpha = 1.64485363, beta = 0.37930952, type.s = 0,
                                type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(ncp = 0, null.ncp = c(-2, 2), df = 100, alpha = 0.05, alternative = "two.one.sided",
                               plot = FALSE, verbose = 0),
                 structure(list(power = 0.261531695, ncp = 0, null.ncp = c(-2, 2), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 0.33388239 * c(-1, 1), beta = 0.73846831,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(ncp = 2, null.ncp = c(-1, 1), df = 100, alpha = 0.05, alternative = "two.one.sided",
                               plot = FALSE, verbose = 0),
                 structure(list(power = 0.16906976, ncp = 2, null.ncp = c(-1, 1), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 2.96236048 * c(-1, 1), beta = 0.83093024,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(power = 0.80, df = 100, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 2.82524975, null.ncp = 0, df = 100, alpha = 0.05,
                                alternative = "two.sided", t.alpha = 1.959964 * c(-1, 1), beta = 0.2,
                                type.s = 0.00000135394609, type.m = 1.12380131),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(power = 0.80, df = 100, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 2.505827, null.ncp = 0, df = 100, alpha = 0.05,
                                alternative = "one.sided", t.alpha = 1.64485363, beta = 0.2, type.s = 0,
                                type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(suppressWarnings(power.lp.test(power = 0.80, req.sign = "0", null.ncp = c(-3, 3), df = 100, alpha = 0.05,
                                                alternative = "two.one.sided", plot = FALSE, verbose = 0)),
                 structure(list(power = 0.810253357, ncp = 0, null.ncp = c(-3, 3), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 1.311328974 * c(-1, 1), beta = 0.189746643,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(power = 0.80, req.sign = "+", null.ncp = c(-1, 1), df = 100, alpha = 0.05,
                               alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 3.8441069, null.ncp = c(-1, 1), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 2.96236048 * c(-1, 1), beta = 0.2,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(power = 0.80, req.sign = "-", null.ncp = c(-1, 1), df = 100, alpha = 0.05,
                               alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = -3.8441069, null.ncp = c(-1, 1), df = 100, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 2.96236048 * c(-1, 1), beta = 0.2,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(ncp = 2.82524975, power = 0.80, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 2.82524975, null.ncp = 0, df = 100, alpha = 0.05,
                                alternative = "two.sided", t.alpha = 1.959964 * c(-1, 1), beta = 0.2,
                                type.s = 0.00000135394605, type.m = 1.12380131),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(ncp = 2.505827, power = 0.80, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 2.505827, null.ncp = 0, df = 100.000057, alpha = 0.05,
                                alternative = "one.sided", t.alpha = 1.64485363, beta = 0.2, type.s = 0,
                                type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(suppressWarnings(power.lp.test(ncp = 0, power = 0.810253357, req.sign = "0", null.ncp = c(-3, 3),
                                                alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0)),
                 structure(list(power = 0.810253357, ncp = 0, null.ncp = c(-3, 3), df = 99.9999944, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 1.311328974 * c(-1, 1), beta = 0.189746643,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(ncp = 3.8441069, power = 0.80, req.sign = "+", null.ncp = c(-1, 1), alpha = 0.05,
                               alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 3.8441069, null.ncp = c(-1, 1), df = 99.999972, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 2.96236048 * c(-1, 1), beta = 0.2,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))
    expect_equal(power.lp.test(ncp = -3.8441069, power = 0.80, req.sign = "-", null.ncp = c(-1, 1), alpha = 0.05,
                               alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = -3.8441069, null.ncp = c(-1, 1), df = 99.999972, alpha = 0.05,
                                alternative = "two.one.sided", t.alpha = 2.96236048 * c(-1, 1), beta = 0.2,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "lp")))

    expect_error(power.lp.test(power = 0.80),
                 "Exactly two of the parameters `ncp`, `df`, or `power` must be given, one has to be NULL.")
    expect_error(power.lp.test(power = 0.80, df = 2), "`df` can not be smaller than 3.")
    expect_warning(power.lp.test(ncp = 36, df = 100, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                   "Consider using a z-test. Lambda-prime distribution with a large non-centrality parameter can be unreliable.")
    expect_warning(power.lp.test(power = 0.80, req.sign = "0", null.ncp = c(-2, 2), df = 100, alpha = 0.05,
                                 alternative = "two.one.sided", plot = FALSE, verbose = 0),
                   "The target power rate cannot be achieved within the null bounds.")
})
