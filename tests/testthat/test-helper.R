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
