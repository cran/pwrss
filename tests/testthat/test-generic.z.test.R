test_that("generic.z.test.R works", {
    # power.z.test -----------------------------------------------------------------------------------------------------
    expect_equal(power.z.test(mean = 1.96, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.500058649, mean = 1.96, sd = 1, null.mean = 0, null.sd = 1, alpha = 0.05,
                                alternative = "two.sided", z.alpha = 1.959964 * c(-1, 1), beta = 0.49994135,
                                type.s = 0.000088551816, type.m = 1.40704658),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(mean = 1.96, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.62367474, mean = 1.96, sd = 1, null.mean = 0, null.sd = 1, alpha = 0.05,
                                alternative = "one.sided", z.alpha = 1.64485363, beta = 0.37632526, type.s = 0, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(mean = 0, null.mean = c(-2, 2), alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.277520063, mean = 0, sd = 1, null.mean = c(-2, 2), null.sd = 1, alpha = 0.05,
                                alternative = "two.one.sided", z.alpha = 0.355146373 * c(-1, 1), beta = 0.72247994,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(mean = 0, null.mean = 2, alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.277520063, mean = 0, sd = 1, null.mean = c(-2, 2), null.sd = 1, alpha = 0.05,
                                alternative = "two.one.sided", z.alpha = 0.355146373 * c(-1, 1), beta = 0.72247994,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(mean = 2, null.mean = c(-1, 1), alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.168537023, mean = 2, sd = 1, null.mean = c(-1, 1), null.sd = 1, alpha = 0.05,
                                alternative = "two.one.sided", z.alpha = 2.959964 * c(-1, 1), beta = 0.83146298,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(mean = 2, null.mean = -1, alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.168537023, mean = 2, sd = 1, null.mean = c(-1, 1), null.sd = 1, alpha = 0.05,
                                alternative = "two.one.sided", z.alpha = 2.959964 * c(-1, 1), beta = 0.83146298,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(power = 0.80, alpha = 0.05, alternative = "two.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, mean = 2.80158179, sd = 1, null.mean = 0, null.sd = 1, alpha = 0.05,
                                alternative = "two.sided", z.alpha = 1.959964 * c(-1, 1), beta = 0.2,
                                type.s = 0.00000120072319, type.m = 1.124912517),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(power = 0.80, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, mean = 2.486474861, sd = 1, null.mean = 0, null.sd = 1, alpha = 0.05,
                                alternative = "one.sided", z.alpha = 1.644853627, beta = 0.2, type.s = 0, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(suppressWarnings(power.z.test(power = 0.80, req.sign = "0", null.mean = c(-3, 3), alpha = 0.05,
                                               alternative = "two.one.sided", plot = FALSE, verbose = 0)),
                 structure(list(power = 0.824629074, mean = 0, sd = 1, null.mean = c(-3, 3), null.sd = 1, alpha = 0.05,
                                alternative = "two.one.sided", z.alpha = 1.355146373 * c(-1, 1), beta = 0.175370926,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(power = 0.80, req.sign = "+", null.mean = c(-1, 1), alpha = 0.05,
                              alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, mean = 3.80158522, sd = 1, null.mean = c(-1, 1), null.sd = 1, alpha = 0.05,
                                alternative = "two.one.sided", z.alpha = 2.959964 * c(-1, 1), beta = 0.2,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "z")))
    expect_equal(power.z.test(power = 0.80, req.sign = "-", null.mean = c(-1, 1), alpha = 0.05,
                              alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, mean = -3.80158522, sd = 1, null.mean = c(-1, 1), null.sd = 1, alpha = 0.05,
                                alternative = "two.one.sided", z.alpha = 2.959964 * c(-1, 1), beta = 0.2,
                                type.s = NA, type.m = NA),
                           class = c("pwrss", "generic", "z")))

    expect_error(power.z.test(), "Exactly one of the parameters `mean` or `power` must be given, one has to be NULL.")
    expect_error(power.z.test(mean = 2, null.mean = c(-1, 1), alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 "If `alternative` is \"two.sided\" or \"one.sided\", `null.mean` must be of length one.")
    expect_error(power.z.test(mean = 2, null.mean = c(-1, 1, 1), alpha = 0.05, alternative = "two.one.sided",
                              plot = FALSE, verbose = 0),
                 paste("If `alternative` is \"two.one.sided\", `null.mean` must be of length one \\(absolute value,",
                       "that must be different from 0\\), or length two \\(lower and upper bounds)."))
    expect_error(power.z.test(mean = 2, alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 paste("If `alternative` is \"two.one.sided\", `null.mean` must be of length one \\(absolute value,",
                       "that must be different from 0\\), or length two \\(lower and upper bounds\\)."))
    expect_warning(power.z.test(power = 0.80, req.sign = "0", null.mean = c(-2, 2), alpha = 0.05,
                                alternative = "two.one.sided", plot = FALSE, verbose = 0),
                   "The target power rate cannot be achieved within the null bounds.")
})
