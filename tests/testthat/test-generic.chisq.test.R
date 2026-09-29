test_that("generic.chisq.test.R works", {
    # power.chisq.test -------------------------------------------------------------------------------------------------
    expect_equal(power.chisq.test(ncp = 20, df = 100, alpha = 0.05, plot = FALSE, verbose = 0),
                 structure(list(power = 0.381376391, ncp = 20, null.ncp = 0, df = 100, alpha = 0.05, chisq.alpha = 124.342113),
                           class = c("pwrss", "generic", "chisq")))
    expect_equal(power.chisq.test(ncp = 20, null.ncp = 10, df = 10, alpha = 0.05, plot = FALSE, verbose = 0),
                 structure(list(power = 0.31365895, ncp = 20, null.ncp = 10, df = 10, alpha = 0.05, chisq.alpha = 34.0886349),
                           class = c("pwrss", "generic", "chisq")))
    expect_equal(power.chisq.test(ncp = 20, df = 17, alpha = 0.05, plot = FALSE, verbose = 0),
                 structure(list(power = 0.80744863, ncp = 20, null.ncp = 0, df = 17, alpha = 0.05, chisq.alpha = 27.5871116),
                           class = c("pwrss", "generic", "chisq")))
    expect_equal(power.chisq.test(power = 0.8, ncp = 40, alpha = 0.05, plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 40, null.ncp = 0, df = 96.872062, alpha = 0.05, chisq.alpha = 120.846553),
                           class = c("pwrss", "generic", "chisq")))
    expect_equal(power.chisq.test(power = 0.8, df = 100, alpha = 0.05, plot = FALSE, verbose = 0),
                 structure(list(power = 0.8, ncp = 40.5563986, null.ncp = 0, df = 100, alpha = 0.05, chisq.alpha = 124.342113),
                           class = c("pwrss", "generic", "chisq")))

    expect_error(power.chisq.test(ncp = NULL, df = NULL, alpha = 0.05),
                 "Exactly two of the parameters `ncp`, `df`, or `power` must be given, one has to be NULL.")
    expect_error(power.chisq.test(ncp = 1, null.ncp = 11, df = 17, alpha = 0.05, plot = FALSE, verbose = 0),
                 "`ncp` should be greater than or equal to `null.ncp`.")
    expect_error(power.chisq.test(power = 0.8, ncp = NULL, df = 0.1, alpha = 0.05),
                 "Degrees of freedom can not be smaller than 1.")
})
