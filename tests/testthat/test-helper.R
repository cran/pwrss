# ensure.verbose -------------------------------------------------------------------------------------------------------
test_that("ensure.verbose works", {
    expect_equal(ensure.verbose(), 1)
    expect_equal(ensure.verbose("A"), 1)
    expect_equal(ensure.verbose(complex(1)), 1)
    expect_equal(ensure.verbose(-2), 1)
    expect_equal(ensure.verbose(-1), -1)
    expect_equal(ensure.verbose(0), 0)
    expect_equal(ensure.verbose(1), 1)
    expect_equal(ensure.verbose(2), 2)
    expect_equal(ensure.verbose(3), 1)
    expect_equal(ensure.verbose(-2L), 1)
    expect_equal(ensure.verbose(-1L), -1)
    expect_equal(ensure.verbose(0L), 0)
    expect_equal(ensure.verbose(1L), 1)
    expect_equal(ensure.verbose(2L), 2)
    expect_equal(ensure.verbose(3L), 1)
    expect_equal(ensure.verbose(FALSE), 0)
    expect_equal(ensure.verbose(TRUE), 1)
})

# get.requested --------------------------------------------------------------------------------------------------------
test_that("get.requested works", {
    probs <- list(NULL, NULL)
    expect_error(get.requested(es = probs, n = NULL, power = NULL),
                 "Exactly one element / entry of `probs` can be NULL, not both.")

    d <- NA
    n <- NA
    power <- NULL
    expect_error(get.requested(es = d, n = n, power = power),
                 "Maximally one input parameter \\(`d`, `n`, or `power`\\) can be excluded from checking.")

    d <- NULL
    n <- NULL
    power <- NULL
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly two of the parameters `d`, `n`, or `power` must be given, one has to be NULL.")
    d <- 0.20
    n <- NULL
    power <- NULL
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly two of the parameters `d`, `n`, or `power` must be given, one has to be NULL.")
    d <- NULL
    n <- 100
    power <- NULL
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly two of the parameters `d`, `n`, or `power` must be given, one has to be NULL.")
    d <- NULL
    n <- NULL
    power <- 0.80
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly two of the parameters `d`, `n`, or `power` must be given, one has to be NULL.")
    d <- 0.20
    n <- 200
    power <- 0.8
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly two of the parameters `d`, `n`, or `power` must be given, one has to be NULL.")
    d <- NA
    n <- 200
    power <- 0.8
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly one of the parameters `n` or `power` must be given, one has to be NULL.")
    d <- NA
    n <- NULL
    power <- NULL
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly one of the parameters `n` or `power` must be given, one has to be NULL.")

    d <- 0.2
    n <- NA
    power <- 0.8
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly one of the parameters `d` or `power` must be given, one has to be NULL.")
    d <- NULL
    n <- NA
    power <- NULL
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly one of the parameters `d` or `power` must be given, one has to be NULL.")

    d <- 0.2
    n <- 200
    power <- NA
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly one of the parameters `d` or `n` must be given, one has to be NULL.")
    d <- NULL
    n <- NULL
    power <- NA
    expect_error(get.requested(es = d, n = n, power = power),
                 "Exactly one of the parameters `d` or `n` must be given, one has to be NULL.")

    d <- NULL
    n <- 1000
    power <- 0.80
    expect_equal(get.requested(es = d, n = n, power = power), "es")

    d <- 0.20
    n <- NULL
    power <- 0.80
    expect_equal(get.requested(es = d, n = n, power = power), "n")

    d <- 0.20
    n <- 1000
    power <- NULL
    expect_equal(get.requested(es = d, n = n, power = power), "power")

    d <- NA
    n <- NULL
    power <- 0.80
    expect_equal(get.requested(es = d, n = n, power = power), "n")

    d <- NA
    n <- 1000
    power <- NULL
    expect_equal(get.requested(es = d, n = n, power = power), "power")

    d <- NULL
    n <- NA
    power <- 0.8
    expect_equal(get.requested(es = d, n = n, power = power), "es")

    d <- 0.2
    n <- NA
    power <- NULL
    expect_equal(get.requested(es = d, n = n, power = power), "power")

    d <- NULL
    n <- 1000
    power <- NA
    expect_equal(get.requested(es = d, n = n, power = power), "es")

    d <- 0.2
    n <- NULL
    power <- NA
    expect_equal(get.requested(es = d, n = n, power = power), "n")
})

# get.interval ---------------------------------------------------------------------------------------------------------
test_that("get.interval works", {
    # when using the distribution values with rather large values (e.g., 8331.498), MacOS and Windows need lower
    # tolerances when comparing to the values determined on Linux; on MacOS or with ARM processors, this already occurs
    # for the t-distribution functions that are part of base R
    isARM <- (Sys.info()["sysname"] == "Darwin" || Sys.info()[["machine"]] %in% c("arm64", "aarch64"))
    isWin <- (Sys.info()["sysname"] == "Windows")
    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "z", alternative = "two.sided", sd = 1), c(0, +8.32130487))
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "z", alternative = "two.sided", sd = 1), c(-8.32130487, 0))
    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "t", alternative = "two.sided", df = 3),
                 c(0, +8331.4980), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "t", alternative = "two.sided", df = 3),
                 c(-8331.4993, 0), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "t", alternative = "two.sided", df = 3e4),
                 c(0, +8.32545167), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "t", alternative = "two.sided", df = 3e4),
                 c(-8.32545167, 0), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "lp", alternative = "two.sided", df = 3e4),
                 c(0, +8.32167150), tolerance = if (isARM) 1e-2 else if (isWin) 1e-4 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "lp", alternative = "two.sided", df = 3e4),
                 c(-8.32167055, 0), tolerance = if (isARM) 1e-2 else if (isWin) 1e-4 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0.5, req.sign = "+", distribution = "binom", alternative = "two.sided"), c(0.5, 0.9999))
    expect_equal(get.interval(null.ncp = 0.5, req.sign = "-", distribution = "binom", alternative = "two.sided"), c(0.0001, 0.5))

    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "z", alternative = "one.sided", sd = 1), c(0, +8.00619452))
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "z", alternative = "one.sided", sd = 1), c(-8.00619452, 0))
    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "t", alternative = "one.sided", df = 3),
                 c(0, +6529.29741), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "t", alternative = "one.sided", df = 3),
                 c(-6529.29700, 0), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "t", alternative = "one.sided", df = 3e4),
                 c(0, +8.00984221), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "t", alternative = "one.sided", df = 3e4),
                 c(-8.00984221, 0), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "+", distribution = "lp", alternative = "one.sided", df = 3e4),
                 c(0, +8.0088033), tolerance = if (isARM) 1e-2 else if (isWin) 1e-3 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0, req.sign = "-", distribution = "lp", alternative = "one.sided", df = 3e4),
                 c(-8.0088029, 0), tolerance = if (isARM) 1e-2 else if (isWin) 1e-3 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = 0.5, req.sign = "+", distribution = "binom", alternative = "one.sided"), c(0.5, 0.9999))
    expect_equal(get.interval(null.ncp = 0.5, req.sign = "-", distribution = "binom", alternative = "one.sided"), c(0.0001, 0.5))

    expect_equal(get.interval(null.ncp = c(1, -1), req.sign = "0", distribution = "binom", alternative = "two.one.sided"), c(-1, 1))
    expect_equal(get.interval(null.ncp = c(2, -2), req.sign = "0", distribution = "z", alternative = "two.one.sided", sd = 1), c(-2, 2))
    expect_equal(get.interval(null.ncp = c(2, -2), req.sign = "0", distribution = "t", alternative = "two.one.sided", df = 2), c(-2, 2))
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "+", distribution = "z", alternative = "two.one.sided", sd = 1), c(2, 10.0061945))
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "-", distribution = "z", alternative = "two.one.sided", sd = 1), c(-10.0061945, -2))
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "+", distribution = "t", alternative = "two.one.sided", df = 3),
                 c(2, 16810.1936),  tolerance = if (isARM) 2e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "-", distribution = "t", alternative = "two.one.sided", df = 3),
                 c(-16810.206, -2), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "+", distribution = "t", alternative = "two.one.sided", df = 3e4),
                 c(2, 10.0120826), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "-", distribution = "t", alternative = "two.one.sided", df = 3e4),
                 c(-10.01208226, -2), tolerance = if (isARM) 1e-6 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "+", distribution = "lp", alternative = "two.one.sided", df = 3e4),
                 c(2, 10.00847078), tolerance = if (isARM) 1e-2 else if (isWin) 1e-4 else testthat_tolerance())
    expect_equal(get.interval(null.ncp = c(-2, 2), req.sign = "-", distribution = "lp", alternative = "two.one.sided", df = 3e4),
                 c(-10.0084683, -2), tolerance = if (isARM) 1e-2 else if (isWin) 1e-4 else testthat_tolerance())
})
