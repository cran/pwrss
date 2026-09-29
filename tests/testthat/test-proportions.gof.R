test_that("proportions.gof.R works", {
    # power.chisq.gof (= pwrss.chisq.gofit) ----------------------------------------------------------------------------
    mtxW <- probs.to.w(c(0.28, 0.72), verbose = 0)
    crrRes <- power.chisq.gof(w = mtxW$w, df = mtxW$df, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = mtxW$w, null.w = 0, df = mtxW$df, n = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 7.9376, null.ncp = 0, chisq.alpha = 3.84145882, w = mtxW$w,
                      power = 0.8043919, n = 41))
    expect_equal(crrRes,     pwrss.chisq.gofit(w = mtxW$w, df = mtxW$df, power = 0.80, verbose = FALSE))
    expect_equal(crrRes,     pwrss.chisq.gofit(p1 = mtxW$prob.matrix, power = 0.80, verbose = FALSE))
    expect_equal(crrRes,     pwrss.chisq.gofit(p1 = mtxW$prob.matrix, p0 = mtxW$null.prob.matrix, power = 0.80, verbose = FALSE))
    expect_equal(crrRes[-1], pwrss.chisq.gofit(w = mtxW$w, df = mtxW$df, n = 41, verbose = FALSE)[-1])

    crrRes <- power.chisq.gof(w = mtxW$w, df = mtxW$df, n = 41, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = mtxW$w, null.w = 0, df = mtxW$df, n = 41, power = NULL, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 7.9376, null.ncp = 0,
                      chisq.alpha = 3.84145882, w = mtxW$w, power = 0.8043919, n = 41))

    crrRes <- power.chisq.gof(df = mtxW$df, n = 41, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = NULL, null.w = 0, df = mtxW$df, n = 41, power = 0.80, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 7.8486511, null.ncp = 0, chisq.alpha = 3.84145882,
                      w = 0.43752773, power = 0.79998954, n = 41))

    mtxW <- probs.to.w(rbind(c(0.056, 0.132), c(0.944, 0.868)), verbose = 0)
    crrRes <- power.chisq.gof(w = mtxW$w, df = mtxW$df, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = mtxW$w, null.w = 0, df = mtxW$df, n = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 7.8504063, null.ncp = 0, chisq.alpha = 3.84145882,
                      w = mtxW$w, power = 0.80007722, n = 463))

    crrRes <- power.chisq.gof(w = mtxW$w, df = mtxW$df, n = 463, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = mtxW$w, null.w = 0, df = mtxW$df, n = 463, power = NULL, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 7.8504063, null.ncp = 0, chisq.alpha = 3.84145882,
                      w = mtxW$w, power = 0.80007722, n = 463))

    crrRes <- power.chisq.gof(df = mtxW$df, n = 463, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = NULL, null.w = 0, df = mtxW$df, n = 463, power = 0.80, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 7.8461815, null.ncp = 0, chisq.alpha = 3.84145882,
                      w = 0.130178325, power = 0.7998661, n = 463))

    mtxW <- probs.to.w(cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078), c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101)), verbose = 0)
    crrRes <- power.chisq.gof(w = mtxW$w, df = mtxW$df, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = mtxW$w, null.w = 0, df = mtxW$df, n = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 11.9353027, null.ncp = 0, chisq.alpha = 9.48772904,
                      w = mtxW$w, power = 0.80000063, n = 13069))

    crrRes <- power.chisq.gof(w = mtxW$w, df = mtxW$df, n = 13069, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = mtxW$w, null.w = 0, df = mtxW$df, n = 13069, power = NULL, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 11.9353027, null.ncp = 0, chisq.alpha = 9.48772904,
                      w = mtxW$w, power = 0.80000063, n = 13069))

    crrRes <- power.chisq.gof(df = mtxW$df, n = 13069, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "chisq", "gof"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(w = NULL, null.w = 0, df = mtxW$df, n = 13069, power = 0.80, alpha = 0.05, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "chisq.alpha", "w", "power", "n")],
                 list(test = "chisq", df = mtxW$df, ncp = 11.9548133, null.ncp = 0, chisq.alpha = 9.48772904,
                      w = 0.030244765, power = 0.800734563, n = 13069))

    expect_error(power.chisq.gof(w = 0.01, null.w = 0.1, df = 1, power = 0.80, alpha = 0.05, verbose = 0),
                 "`w` should be greater than or equal to `null.w`.")
    expect_error(power.chisq.gof(w = 1e-6, df = 1, power = 0.99, alpha = 1e-6),
                 "Design is not feasible.")
    expect_error(power.chisq.gof(df = 1, n = 2, power = 0.80, alpha = 0.05, verbose = 0),
                 "Design is not feasible.")
    expect_error(pwrss.chisq.gofit(p1 = data.frame, power = 0.80, alpha = 0.05, verbose = FALSE),
                 "`p1` needs to be either a vector or a matrix.")
    expect_error(pwrss.chisq.gofit(p1 = array(1, c(2, 2, 2)), power = 0.80, alpha = 0.05, verbose = FALSE),
                 "`p1` needs to be either a vector or a matrix.")
    expect_error(pwrss.chisq.gofit(p1 = c(0.280, 0.721), power = 0.80, alpha = 0.05, verbose = FALSE),
                 "Cell probabilities in `p1` \\(and `p0` if given\\) should sum to 1.")
    expect_error(pwrss.chisq.gofit(p1 = matrix(c(0.2, 0.3, 0.4, 0.3, 0.5, 0.2), ncol = 2), power = 0.80),
                 "Cell probabilities \\(per column\\) in `p1` \\(and `p0` if given\\) should sum to 1.")
    expect_error(pwrss.chisq.gofit(p1 = matrix(c(0.2, 0.5, 0.3), ncol = 1), p0 = matrix(c(0.2, 0.5, 0.3), nrow = 1), power = 0.80),
                 "Dimensions of `p1` and `p0` do not match up.")
    expect_error(pwrss.chisq.gofit(w = 0.44, power = 0.80, alpha = 0.05, verbose = 0),
                 "You need to specify either `w` and `df` or `p1`.")
    expect_warning(pwrss.chisq.gofit(p1 = c(0.72, 0.28), w = 0.44, df = 1, power = 0.80, alpha = 0.05, verbose = 0),
                   "Ignoring any specifications to `p1`, or `p0`.")
})
