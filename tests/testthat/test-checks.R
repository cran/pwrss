test_that("checks.R works", {
    # check.proportion -------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(seq(0, 10) / 10, check.proportion)))
    expect_error(check.proportion(-1, 2, 1.1, 0, 1),
                 "Arguments `-1`, `2`, `1.1` do not have valid proportion values \\(must be length 1, numeric, >= 0, and <= 1\\)")

    # check.power ------------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(seq(1, 99, 14) / 100, check.power)))
    expect_error(check.power(-1, 2, 0, 1, 0.01, 0.99),
                 "Arguments `-1`, `2`, `0`, `1` do not have valid power values \\(must be length 1, numeric, >= 0.01, and <= 0.99\\)")

    # check.correlation ------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(seq(-10, 10) / 10, check.correlation)))
    expect_error(check.correlation(-1.1, 2, 1.1, 0, 1, -1),
                 "Arguments `-1.1`, `2`, `1.1` do not have valid correlation values \\(must be length 1, numeric, >= -1, and <= 1\\)")

    # check.logical ----------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(sample(c(TRUE, FALSE), 20, replace = TRUE), check.logical)))
    expect_error(check.logical(1, 0, NA, TRUE, FALSE),
                 "Arguments `1`, `0`, `NA` do not have valid logical values \\(must be length 1, TRUE or FALSE\\)")

    # check.sample.size ------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(c(2, sample(seq(2, 10000), 20)), check.sample.size)))
    expect_error(check.sample.size(1, 0, 1.1, NA, Inf, 2, 100),
                 "Arguments `1`, `0`, `1.1`, `NA`, `Inf` do not have valid sample size values \\(must be length 1, integer-like, > 1, and finite\\)")

    # check.factor.level -----------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(c(2, sample(seq(2, 10000), 20)), check.factor.level)))
    expect_error(check.factor.level(1, 0, 1.1, NA, Inf, 2, 100),
                 "Arguments `1`, `0`, `1.1`, `NA`, `Inf` do not have valid factor level values \\(must be length 1, integer-like, > 1, and finite\\)")

    # check.nonnegative ------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(seq(0, 100) / 10, check.nonnegative)))
    expect_error(check.nonnegative(-1, NA, Inf, 0, 1, 1.1, 10),
                 "Arguments `-1`, `NA`, `Inf` do not have valid non-negative values \\(must be length 1, numeric, >= 0, and finite\\)")

    # check.positive ---------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(seq(1, 100) / 10, check.positive)))
    expect_error(check.positive(-1, 0, NA, Inf, 1, 1.1, 10),
                 "Arguments `-1`, `0`, `NA`, `Inf` do not have valid positive values \\(must be length 1, numeric, > 0, and finite\\)")

    # check.numeric ----------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(rnorm(100), check.numeric)))
    expect_null(unlist(lapply(sample(1e5, 100) - 5e4, check.numeric)))
    expect_error(check.numeric(NA, Inf, seq(3), as.complex(1), 1, 1.1, 10),
                 "Arguments `NA`, `Inf`, `seq\\(3\\)`, `as.complex\\(1\\)` do not have valid numeric values \\(must be length 1, numeric, and finite\\)")

    # check.vector -----------------------------------------------------------------------------------------------------
    expect_null(unlist(lapply(seq(100), function(n) check.vector(rnorm(sample(10:20)), check.numeric))))
    expect_null(unlist(lapply(seq(100), function(n) check.vector(runif(sample(10:20), 0, 1), check.proportion))))
    expect_null(unlist(lapply(seq(100), function(n) check.vector(runif(sample(10:20), -1, 1), check.correlation))))
    expect_error(check.vector(rnorm(1), check.numeric), "`rnorm\\(1\\)` neeeds to be a vector with a length of at least 2.")
    expect_error(check.vector(1, check.numeric),
                 "`1` neeeds to be a vector with a length of at least 2.")
    expect_error(check.vector(sample(100, 5), check.numeric, 6),
                 "`sample\\(100, 5\\)` neeeds to be a vector with a length of at least 6.")
    expect_error(check.vector(sample(100, 5), check.numeric, 6, 10),
                 "`sample\\(100, 5\\)` neeeds to be a vector with a length of at least 6 and maximally 10.")
    expect_error(check.vector(sample(100, 5), check.numeric, 6, 6),
                 "`sample\\(100, 5\\)` neeeds to be a vector with a length of 6.")
    expect_error(check.vector(matrix(rnorm(100), ncol = 10), check.numeric),
                 "`matrix\\(rnorm\\(100\\), ncol = 10\\)` neeeds to be a vector with a length of at least 2.")
    expect_error(check.vector(rnorm(100), check.logical),
                 "All elements of `rnorm\\(100\\)` need to be valid logical values \\(TRUE or FALSE\\)")
    expect_error(check.vector(runif(100, -2, 2), check.correlation),
                 "All elements of `runif\\(100, -2, 2\\)` need to be valid correlation values \\(numeric, >= -1, and <= 1\\)")

    # check.same.lengths -----------------------------------------------------------------------------------------------
    expect_null(check.same.lengths(rnorm(100), runif(100), sample(100, 100), seq(100)))
    expect_null(check.same.lengths(rnorm(100), runif(100), NULL, seq(100), NULL))
    expect_error(check.same.lengths(NULL, runif(100), sample(100, 100), seq(100)),
                 "To use `check.same.lengths`, `NULL` can not be NULL.")
    expect_error(check.same.lengths(rnorm(100), runif(98), sample(100, 98), seq(100)),
                 "`runif\\(98\\)` and `sample\\(100, 98\\)` should have the same length as `rnorm\\(100\\)`.")
    expect_error(check.same.lengths(rnorm(100), NULL, runif(98), NULL, sample(100, 98), NULL, seq(100), NULL),
                 "`runif\\(98\\)` and `sample\\(100, 98\\)` should have the same length as `rnorm\\(100\\)`.")
    expect_error(check.same.lengths(rnorm(100), NULL, runif(98), NULL, sample(100, 100), NULL, NULL, seq(98)),
                 "`runif\\(98\\)` and `seq\\(98\\)` should have the same length as `rnorm\\(100\\)`.")

    # check.margins ----------------------------------------------------------------------------------------------------
    expect_equal(check.margins(0, check.nonnegative, "two.sided"), 0)
    expect_equal(check.margins(1, check.positive, "two.one.sided"), c(-1, 1))
    expect_equal(check.margins(0.2, check.positive, "two.one.sided"), c(-0.2, 0.2))
    expect_equal(check.margins(c(0.2, 0.4), check.proportion, "two.one.sided"), c(0.2, 0.4))
    expect_error(check.margins(0, check.positive, "two.one.sided"),
                 "All elements of `0` need to be valid positive values \\(numeric, > 0, and finite\\)")
    expect_error(check.margins(c(1, -1), check.proportion, "two.one.sided"),
                 "All elements of `c\\(1, -1\\)` need to be valid proportion values \\(numeric, >= 0, and <= 1\\)")
    expect_error(check.margins(c(0.2, 0.4), check.proportion, "one.sided"),
                 "If `alternative` is \"two.sided\" or \"one.sided\", `c\\(0.2, 0.4\\)` must be of length one.")
    expect_error(check.margins(c(0.4, 0.2), check.proportion, "two.one.sided"),
                 paste("If `alternative` is \"two.one.sided\", `c\\(0.4, 0.2\\)` must be of length two \\(lower and",
                       "upper bounds, with the upper bound being larger than the lower bound\\)."))
    expect_error(check.margins(c(2, 3, 4), check.positive, "two.one.sided"),
                 paste("If `alternative` is \"two.one.sided\", `c\\(2, 3, 4\\)` must be of length one \\(absolute value,",
                       "that must be different from 0\\), or length two \\(lower and upper bounds\\)."))
    expect_error(check.margins(0, check.numeric, "two.one.sided"),
                 paste("If `alternative` is \"two.one.sided\", `0` must be of length one \\(absolute value,",
                       "that must be different from 0\\), or length two \\(lower and upper bounds\\)."))

    # check.correlation.matrix  ----------------------------------------------------------------------------------------
    expect_null(check.correlation.matrix(cor(matrix(rnorm(1000), ncol = 10))))
    tstMtx <- diag(5)[-1, ]
    expect_error(check.correlation.matrix(tstMtx), "Correlation matrix `tstMtx` is not square")
    tstMtx <- rbind(diag(5)[-1, ], diag(5)[1, ])
    expect_error(check.correlation.matrix(tstMtx), "Correlation matrix `tstMtx` is not symmetric")
    tstMtx <- cor(matrix(rnorm(400), ncol = 4)) + matrix(rep(2, 16), nrow = 4) - (diag(4) * 2)
    expect_error(check.correlation.matrix(tstMtx), "The values in the correlation matrix \\(`tstMtx`\\) must be numeric, >= -1 and <= 1")
    tstMtx <- diag(5) + matrix(c(1, rep(0, 24)), nrow = 5)
    expect_error(check.correlation.matrix(tstMtx), "All values in the main diagonal of the correlation matrix \\(`tstMtx`\\) must be 1")
    tstMtx <- matrix(c(1, 0.7, 0.3, 0.7, 1, -0.5, 0.3, -0.5, 1), nrow = 3)
    expect_error(check.correlation.matrix(tstMtx), "Correlation matrix `tstMtx` is not positive definite")
    # neither I nor ChatGPT could not generate test cases for (1) a matrix that is not invertible (determinant of 0 or lower)
    # and (2) a matrix that is not well conditioned (kappa > 1000)

    # check.null  ------------------------------------------------------------------------------------------------------
    expect_equal(check.null(NULL, NULL, NULL, NULL), rep(TRUE, 4))
    expect_equal(check.null(NA, NULL, Inf, NULL), c(FALSE, TRUE, FALSE, TRUE))

    # check.not_null  --------------------------------------------------------------------------------------------------
    expect_equal(check.not_null(NULL, NULL, NULL, NULL), rep(FALSE, 4))
    expect_equal(check.not_null(NA, NULL, Inf, NULL), c(TRUE, FALSE, TRUE, FALSE))

    # check.pos_sign  --------------------------------------------------------------------------------------------------
    expect_true(check.pos_sign("+"))
    expect_true(check.pos_sign(1))
    expect_true(check.pos_sign("1"))
    expect_true(check.pos_sign("+1"))
    expect_true(check.pos_sign("positive"))
    expect_false(check.pos_sign("-"))
    expect_false(check.pos_sign(-1))
    expect_false(check.pos_sign("-1"))
    expect_false(check.pos_sign("negative"))
    expect_error(check.pos_sign(" "),       "`req.sign` can only be `\\+` and `-` for this function.")
    expect_error(check.pos_sign(""),        "`req.sign` can only be `\\+` and `-` for this function.")
    expect_error(check.pos_sign("A"),       "`req.sign` can only be `\\+` and `-` for this function.")
    expect_error(check.pos_sign(2),         "`req.sign` can only be `\\+` and `-` for this function.")
    expect_error(check.pos_sign(-2),        "`req.sign` can only be `\\+` and `-` for this function.")

    # check.null_sign  -------------------------------------------------------------------------------------------------
    expect_true(check.null_sign(" ",    "two.one.sided"))
    expect_true(check.null_sign(0,      "two.one.sided"))
    expect_true(check.null_sign("0",    "two.one.sided"))
    expect_true(check.null_sign("",     "two.one.sided"))
    expect_true(check.null_sign("null", "two.one.sided"))
    expect_false(check.null_sign(1,  "two.sided"))
    expect_false(check.null_sign(-1, "two.sided"))
    expect_false(check.null_sign(1,  "one.sided"))
    expect_false(check.null_sign(-1, "one.sided"))
    expect_error(check.null_sign(0, "two.sided"), "`req.sign` can not be `0` for 'one.sided' or 'two.sided' hypothesis tests.")
    expect_error(check.null_sign(0, "one.sided"), "`req.sign` can not be `0` for 'one.sided' or 'two.sided' hypothesis tests.")
    expect_error(check.null_sign(2,   "two.sided"),     "`req.sign` can only be `\\+`, `-` and ` ` for this function.")
    expect_error(check.null_sign(2,   "two.one.sided"), "`req.sign` can only be `\\+`, `-` and ` ` for this function.")
    expect_error(check.null_sign(-2,  "two.sided"),     "`req.sign` can only be `\\+`, `-` and ` ` for this function.")
    expect_error(check.null_sign(-2,  "two.one.sided"), "`req.sign` can only be `\\+`, `-` and ` ` for this function.")
    expect_error(check.null_sign("A", "two.sided"),     "`req.sign` can only be `\\+`, `-` and ` ` for this function.")
    expect_error(check.null_sign("A", "two.one.sided"), "`req.sign` can only be `\\+`, `-` and ` ` for this function.")

})
