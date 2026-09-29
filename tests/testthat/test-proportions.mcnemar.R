test_that("proportions.mcnemar.R works", {
    # NB: the code is partially tested through power.exact.twoprops which serves as a wrapper to power.exact.fisher
    # (for independent tests) and power.exact.mcnemar (for paired tests)

    # power.exact.mcnemar ----------------------------------------------------------------------------------------------
    crrRes <- power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10, n.paired = 100, alpha = 0.05, alternative = "two.sided",
                                  method = "exact", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.20, prob01 = 0.10, req.sign = "+", n.paired = 100, power = NULL, alpha = 0.05,
                          alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.20, prob01 = 0.10, delta = 0.10, odds.ratio = 2, size = 30,
                      prob = 2 / 3, null.prob = 0.5, binom.alpha = c(9, 20), mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.042773945, power = 0.373077084, n.paired = 100))

    crrRes <- power.exact.mcnemar(prob10 = 0.20, req.sign = "-", n.paired = 100, power = 0.3731, alpha = 0.05,
                                  alternative = "two.sided", method = "exact", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.20, prob01 = NULL, req.sign = "-", n.paired = 100, power = 0.3731, alpha = 0.05,
                          alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.20, prob01 = 0.100000178, delta = 0.099999822, odds.ratio = 1.999996436,
                      size = 31, prob = 0.66666627, null.prob = 0.5, binom.alpha = c(9, 21), mean = NA, sd = NA,
                      null.mean = NA, null.sd = NA, z.alpha = NA, alpha = 0.0294493735, power = 0.373075611,
                      n.paired = 100))

    crrRes <- power.exact.mcnemar(prob01 = 0.10, req.sign = "+", n.paired = 100, power = 0.3731, alpha = 0.05,
                                  alternative = "two.sided", method = "exact", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = NULL, prob01 = 0.10, req.sign = "+", n.paired = 100, power = 0.3731, alpha = 0.05,
                          alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.199989317, prob01 = 0.10, delta = 0.099989317, odds.ratio = 1.999893165,
                      size = 30, prob = 0.6666548, null.prob = 0.5, binom.alpha = c(9, 20), mean = NA, sd = NA,
                      null.mean = NA, null.sd = NA, z.alpha = NA, alpha = 0.042773945, power = 0.37301145, n.paired = 100))

    crrRes <- power.exact.mcnemar(prob10 = 0.10, prob01 = 0.20, n.paired = 100, alpha = 0.05, alternative = "one.sided",
                                  method = "exact", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.10, prob01 = 0.20, req.sign = "+", n.paired = 100, power = NULL, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.10, prob01 = 0.20, delta = -0.10, odds.ratio = 0.5, size = 30,
                      prob = 1 / 3, null.prob = 0.5, binom.alpha = 10, mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.049368573, power = 0.514606143, n.paired = 100))

    crrRes <- power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10, n.paired = 100, alpha = 0.05, alternative = "one.sided",
                                  method = "exact", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.20, prob01 = 0.10, req.sign = "+", n.paired = 100, power = NULL, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.20, prob01 = 0.10, delta = 0.10, odds.ratio = 2, size = 30,
                      prob = 2 / 3, null.prob = 0.5, binom.alpha = 19, mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.049368573, power = 0.514606143, n.paired = 100))

    crrRes <- power.exact.mcnemar(prob10 = 0.10, prob01 = 0.20, n.paired = 500, alpha = 0.05, alternative = "two.sided",
                                  method = "approximate", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.10, prob01 = 0.20, req.sign = "+", n.paired = 500, power = NULL, alpha = 0.05,
                          alternative = "two.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "z", prob10 = 0.10, prob01 = 0.20, delta = -0.10, odds.ratio = 0.5, size = 150,
                      prob = 1 / 3, null.prob = 0.5, binom.alpha = c(62, 87), mean = -4.118767908, sd = 1,
                      null.mean = 0, null.sd = 1, z.alpha = 1.95996398 * c(-1, 1), alpha = 0.05, power = 0.98456731,
                      n.paired = 500))

    crrRes <- power.exact.mcnemar(prob10 = 0.10, req.sign = "+", n.paired = 500, power = 0.9846, alpha = 0.05,
                                  alternative = "two.sided", method = "approximate", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.10, prob01 = NULL, req.sign = "+", n.paired = 500, power = 0.9846, alpha = 0.05,
                          alternative = "two.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "z", prob10 = 0.10, prob01 = 0.200023938, delta = -0.100023938, odds.ratio = 0.49994016,
                      size = 151, prob = 0.33330674, null.prob = 0.5, binom.alpha = c(62, 88), mean = -4.119611, sd = 1,
                      null.mean = 0, null.sd = 1, z.alpha = 1.95996398 * c(-1, 1), alpha = 0.05, power = 0.9846,
                      n.paired = 500))

    crrRes <- power.exact.mcnemar(prob01 = 0.20, req.sign = "-", n.paired = 500, power = 0.9846, alpha = 0.05,
                                  alternative = "two.sided", method = "approximate", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = NULL, prob01 = 0.20, req.sign = "-", n.paired = 500, power = 0.9846, alpha = 0.05,
                          alternative = "two.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "z", prob10 = 0.099982902, prob01 = 0.20, delta = -0.100017098, odds.ratio = 0.49991451,
                      size = 150, prob = 0.33329534, null.prob = 0.5, binom.alpha = c(62, 87), mean = -4.11961108,
                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = 1.95996398 * c(-1, 1), alpha = 0.05, power = 0.9846,
                      n.paired = 500))

    crrRes <- power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10, n.paired = 500, alpha = 0.05, alternative = "one.sided",
                                  method = "approximate", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.20, prob01 = 0.10, req.sign = "+", n.paired = 500, power = NULL, alpha = 0.05,
                          alternative = "one.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "z", prob10 = 0.20, prob01 = 0.10, delta = 0.10, odds.ratio = 2, size = 150,
                      prob = 2 / 3, null.prob = 0.5, binom.alpha = 85, mean = 4.1241548, sd = 1, null.mean = 0, null.sd = 1,
                      z.alpha = 1.64485363, alpha = 0.05, power = 0.993418, n.paired = 500))

    crrRes <- power.exact.mcnemar(prob10 = 0.10, prob01 = 0.20, n.paired = 500, alpha = 0.05, alternative = "one.sided",
                                  method = "approximate", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.10, prob01 = 0.20, req.sign = "+", n.paired = 500, power = NULL, alpha = 0.05,
                          alternative = "one.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "z", prob10 = 0.10, prob01 = 0.20, delta = -0.10, odds.ratio = 0.5, size = 150,
                      prob = 1 / 3, null.prob = 0.5, binom.alpha = 64, mean = -4.1241548, sd = 1, null.mean = 0, null.sd = 1,
                      z.alpha = -1.64485363, alpha = 0.05, power = 0.993418, n.paired = 500))

    crrRes <- power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10, power = 0.80, alpha = 0.05, alternative = "two.sided",
                                  method = "exact", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.20, prob01 = 0.10, req.sign = "+", n.paired = NULL, power = 0.80, alpha = 0.05,
                          alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.20, prob01 = 0.10, delta = 0.10, odds.ratio = 2, size = 75,
                      prob = 2 / 3, null.prob = 0.5, binom.alpha = c(28, 46), mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.036954904, power = 0.801446921, n.paired = 249))

    crrRes <- power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10, power = 0.80, alpha = 0.05, alternative = "one.sided",
                                  method = "exact", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.20, prob01 = 0.10, req.sign = "+", n.paired = NULL, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.20, prob01 = 0.10, delta = 0.10, odds.ratio = 2, size = 60,
                      prob = 2 / 3, null.prob = 0.5, binom.alpha = 36, mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.04623049, power = 0.80298386, n.paired = 200))

    crrRes <- power.exact.mcnemar(prob10 = 0.08, prob01 = 0.32, n = 50, alpha = 0.05, alternative = "one.sided",
                                  method = "exact", verbose = 0) # example 5.3 from GPower
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes, c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob",
                           "binom.alpha", "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_identical(crrRes[["parms"]],
                     list(prob10 = 0.08, prob01 = 0.32, req.sign = "+", n.paired = 50, power = NULL, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.08, prob01 = 0.32, delta = -0.24, odds.ratio = 0.25, size = 20,
                      prob = 0.2, null.prob = 0.5, binom.alpha = 5, mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.0206947327, power = 0.839356985, n.paired = 50))
    # -> reasonable close for power = 0.839343 (GPower) / 0.839357 (pwrss); quite off for alpha = 0.032578 (GPower) / 0.020695 (pwrss)
    # likely due to missing alpha-balancing

    expect_error(power.exact.mcnemar(prob10 = 0.002, prob01 = 0.001, power = 0.99, alpha = 1e-6, alternative = "two.sided"),
                 "Sample size exceeds 100,000. Please check the assumptions.")
})
