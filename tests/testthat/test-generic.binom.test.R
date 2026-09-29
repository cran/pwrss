test_that("generic.binom.test.R works", {
    # power.binom.test -------------------------------------------------------------------------------------------------
    expect_equal(power.binom.test(size = 200, prob = 0.6, null.prob = 0.5, alpha = 0.05, alternative = "one.sided",
                                  plot = FALSE, verbose = 0),
                 structure(list(power = 0.860335667, size = 200, prob = 0.6, null.prob = 0.5, alpha = 0.038418816,
                                alternative = "one.sided", binom.alpha = 112),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(size = 200, prob = 0.4, null.prob = 0.5, alpha = 0.05, alternative = "two.sided",
                                  plot = FALSE, verbose = 0),
                  structure(list(power = 0.786848265, size = 200, prob = 0.4, null.prob = 0.5, alpha = 0.040037192,
                                 alternative = "two.sided", binom.alpha = c(85, 114)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(size = 200, prob = 0.5, null.prob = c(0.4, 0.6), alpha = 0.05, alternative = "two.one.sided",
                                  plot = FALSE, verbose = 0),
                 structure(list(power = 0.7707534, size = 200, prob = 0.5, null.prob = c(0.4, 0.6), alpha = 0.049184711,
                                alternative = "two.one.sided", binom.alpha = c(91, 108)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(size = 200, prob = 0.7, null.prob = c(0.4, 0.6), alpha = 0.05, alternative = "two.one.sided",
                                  plot = FALSE, verbose = 0),
                 structure(list(power = 0.84208844, size = 200, prob = 0.7, null.prob = c(0.4, 0.6), alpha = 0.049432239,
                                alternative = "two.one.sided", binom.alpha = c(66, 133)),
                           class = c("pwrss", "generic", "binom")))

    expect_equal(power.binom.test(size = 200, power = 0.80, null.prob = 0.5, alpha = 0.05, alternative = "two.sided",
                                  plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, size = 200, prob = 0.60158815, null.prob = 0.5, alpha = 0.040037192,
                                alternative = "two.sided", binom.alpha = c(85, 114)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(size = 200, power = 0.80, null.prob = 0.5, alpha = 0.05, alternative = "one.sided",
                                  plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, size = 200, prob = 0.591709046, null.prob = 0.5, alpha = 0.038418816,
                                alternative = "one.sided", binom.alpha = 112),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(suppressWarnings(power.binom.test(size = 200, req.sign = "0", power = 0.80, null.prob = c(0.4, 0.6),
                                                   alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0)),
                 structure(list(power = 0.7707534, size = 200, prob = 0.5, null.prob = c(0.4, 0.6), alpha = 0.049184711,
                                alternative = "two.one.sided", binom.alpha = c(91, 108)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(size = 200, req.sign = "+", power = 0.80, null.prob = c(0.4, 0.6), alpha = 0.05,
                                  alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.79999998, size = 200, prob = 0.694808861, null.prob = c(0.4, 0.6),
                                alpha = 0.049432239, alternative = "two.one.sided", binom.alpha = c(66, 133)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(size = 200, req.sign = "-", power = 0.80, null.prob = c(0.4, 0.6), alpha = 0.05,
                                 alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, size = 200, prob = 0.305191136, null.prob = c(0.4, 0.6), alpha = 0.049432239,
                                alternative = "two.one.sided", binom.alpha = c(66, 133)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(power = 0.80, prob = 0.6, null.prob = 0.5, alpha = 0.05, alternative = "two.sided",
                                  plot = FALSE, verbose = 0),
                 structure(list(power = 0.8103826, size = 201, prob = 0.6, null.prob = 0.5, alpha = 0.048001466,
                                alternative = "two.sided", binom.alpha = c(86, 114)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(power = 0.80, prob = 0.6, null.prob = 0.5, alpha = 0.05, alternative = "one.sided",
                                  plot = FALSE, verbose = 0),
                 structure(list(power = 0.81290017, size = 160, prob = 0.6, null.prob = 0.5, alpha = 0.048284677,
                                alternative = "one.sided", binom.alpha = 90),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(power = 0.80, prob = 0.5, null.prob = c(0.4, 0.6), alpha = 0.05,
                                  alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.80196258, size = 218, prob = 0.5, null.prob = c(0.4, 0.6), alpha = 0.045213449,
                                alternative = "two.one.sided", binom.alpha = c(99, 118)),
                           class = c("pwrss", "generic", "binom")))
    expect_equal(power.binom.test(power = 0.80, prob = 0.7, null.prob = c(0.4, 0.6), alpha = 0.05,
                                  alternative = "two.one.sided", plot = FALSE, verbose = 0),
                 structure(list(power = 0.79516075, size = 191, prob = 0.7, null.prob = c(0.4, 0.6), alpha = 0.038224051,
                                alternative = "two.one.sided", binom.alpha = c(62, 128)),
                           class = c("pwrss", "generic", "binom")))

    expect_error(power.binom.test(size = 200, prob = 0.6, null.prob = c(0.5, 0.7), alpha = 0.05, alternative = "one.sided",
                                  plot = FALSE, verbose = 0),
                 "If `alternative` is \"two.sided\" or \"one.sided\", `null.prob` must be of length one.")
    expect_error(power.binom.test(size = -1, prob = 0.6, null.prob = 0.5, alpha = 0.05, alternative = "one.sided", plot = FALSE, verbose = 0),
                 "All elements of `size` need to be valid size values \\(integer-like, >= 0, and finite\\)")
    expect_error(power.binom.test(size = 2, prob = 0.6, null.prob = 0.5, alpha = 0.05, alternative = "two.one.sided", plot = FALSE),
                 paste("If `alternative` is \"two.one.sided\", `null.prob` must be of length two \\(lower and upper bounds,",
                       "with the upper bound being larger than the lower bound\\)."))
    expect_error(power.binom.test(size = 2, prob = 0.6, null.prob = c(0.5, 0.3), alpha = 0.05, alternative = "two.one.sided", plot = FALSE),
                 paste("If `alternative` is \"two.one.sided\", `null.prob` must be of length two \\(lower and upper bounds,",
                       "with the upper bound being larger than the lower bound\\)."))
    expect_warning(power.binom.test(size = 200, req.sign = "0", power = 0.80, null.prob = c(0.4, 0.6),
                                    alpha = 0.05, alternative = "two.one.sided", plot = FALSE, verbose = 0),
                   "The target power rate cannot be achieved within the null bounds.")
})
