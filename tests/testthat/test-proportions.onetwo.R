# power.exact.oneprop --------------------------------------------------------------------------------------------------
test_that("power.exact.oneprop works", {
    crrRes <- power.exact.oneprop(prob = 0.45, null.prob = 0.50, n = 500, alpha = 0.05, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = 500, power = NULL, alpha = 0.05,
                          alternative = "one.sided", verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n")],
                 list(test = "exact", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, size = 500,
                      binom.alpha = 231, alpha = 0.048945014, power = 0.7208035, n = 500))

    crrRes <- power.exact.oneprop(power = 0.80, prob = 0.45, null.prob = 0.50, alpha = 0.05, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = NULL, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n")],
                 list(test = "exact", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, size = 631,
                      binom.alpha = 294, alpha = 0.047223688, power = 0.800822354, n = 631))

    crrRes <- power.exact.oneprop(req.sign = "-", null.prob = 0.50, power = 0.80, n = 500, alpha = 0.05, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, req.sign = "-", null.prob = 0.50, n = 500, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n")],
                 list(test = "exact", prob = 0.44430789, null.prob = 0.50, delta = -0.055692113, odds.ratio = 0.79955766,
                      size = 500, binom.alpha = 231, alpha = 0.048945014, power = 0.8, n = 500))

    crrRes <- power.exact.oneprop(prob = 0.80, null.prob = 0.65, n = 20, alternative = "one.sided", verbose = 0) # example 4.3 from GPower
    expect_s3_class(crrRes, c("pwrss", "exact", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.80, req.sign = "+", null.prob = 0.65, n = 20, power = NULL, alpha = 0.05,
                          alternative = "one.sided", verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "size", "binom.alpha", "alpha", "power", "n")],
                 list(test = "exact", prob = 0.80, null.prob = 0.65, delta = 0.15, odds.ratio = 2.153846154, size = 20,
                      binom.alpha = 16, alpha = 0.044375603, power = 0.41144886, n = 20))
    # -> results are identical: power ~ 0.411449, actual alpha ~ 0.044376

    expect_error(power.exact.oneprop(prob = 0.45, null.prob = c(0.40, 0.50), alpha = 0.05, n = 500, alternative = "one.sided"),
                 "If `alternative` is \"two.sided\" or \"one.sided\", `null.prob` must be of length one.")
    expect_error(power.exact.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, n = 500, alternative = "two.one.sided"),
                 paste("If `alternative` is \"two.one.sided\", `null.prob` must be of length two \\(lower and upper bounds,",
                       "with the upper bound being larger than the lower bound\\)."))
    expect_error(power.exact.oneprop(prob = 0.45, null.prob = 1.2, alpha = 0.05, n = 500, alternative = "one.sided"),
                 "All elements of `null.prob` need to be valid proportion values \\(numeric, >= 0, and <= 1\\)")
})

# power.z.oneprop (= pwrss.z.prop) -------------------------------------------------------------------------------------
test_that("power.z.oneprop works", {
    crrRes <- power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, n = 500, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = 500, power = NULL, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = FALSE, correct = FALSE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, mean = -2.2473329,
                      sd = 1, null.mean = 0, null.sd = 1.00503782, z.alpha = -1.653140096, power = 0.7238084, n = 500))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, n = 500, alternative = "one.sided", arcsine = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = 500, power = NULL, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = TRUE, correct = FALSE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, mean = -2.23981163,
                      sd = 1, null.mean = 0, null.sd = 1.00503782, z.alpha = -1.653140096, power = 0.72128783, n = 500))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, n = 500, alternative = "one.sided", correct = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = 500, power = NULL, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = FALSE, correct = TRUE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, mean = -2.2023862,
                      sd = 1, null.mean = 0, null.sd = 1.00503782, z.alpha = -1.653140096, power = 0.708581722, n = 500))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = NULL, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = FALSE, correct = FALSE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, mean = -2.49646214,
                      sd = 1, null.mean = 0, null.sd = 1.00503782, z.alpha = -1.653140096, power = 0.800475822, n = 617))
    expect_equal(crrRes, pwrss.z.prop(p = 0.45, p0 = 0.50, alpha = 0.05, power = 0.80, alternative = "less", verbose = FALSE))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, power = 0.80, alternative = "one.sided", arcsine = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = NULL, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = TRUE, correct = FALSE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, mean = -2.49615927,
                      sd = 1, null.mean = 0, null.sd = 1.00503782, z.alpha = -1.653140096, power = 0.80039114, n = 621))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, power = 0.80, alternative = "one.sided", correct = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = NULL, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = FALSE, correct = TRUE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, mean = -2.49677972,
                      sd = 1, null.mean = 0, null.sd = 1.00503782, z.alpha = -1.653140096, power = 0.8005646, n = 637))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, power = 0.80, alternative = "two.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = 0.50, n = NULL, power = 0.80, alpha = 0.05,
                          alternative = "two.sided", std.error = "null", arcsine = FALSE, correct = FALSE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = 0.50, delta = -0.05, odds.ratio = 0.818181818, mean = -2.8123106,
                      sd = 1, null.mean = 0, null.sd = 1.00503782, z.alpha = 1.96983792100888 * c(-1, 1), power = 0.80023915,
                      n = 783))
    expect_equal(crrRes, pwrss.z.prop(p = 0.45, p0 = 0.50, alpha = 0.05, power = 0.80, alternative = "not equal", verbose = FALSE))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = c(0.40, 0.50), alpha = 0.05, power = 0.80, alternative = "two.one.sided",
                              std.error = "alternative", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = c(0.40, 0.50), n = NULL, power = 0.80, alpha = 0.05,
                          alternative = "two.one.sided", std.error = "alternative", arcsine = FALSE, correct = FALSE,
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = c(0.40, 0.50), delta = 0.05 * c(1, -1), odds.ratio = c(1.2272727, 0.818181818),
                      mean = 0, sd = 1, null.mean = 2.9267143 * c(-1, 1), null.sd = 1, z.alpha = c(-1.28186067, 1.28186067),
                      power = 0.800108473, n = 848))

    crrRes <- power.z.oneprop(prob = 0.45, null.prob = c(0.35, 0.40), alpha = 0.05, power = 0.80, alternative = "two.one.sided",
                              std.error = "alternative", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.45, req.sign = "+", null.prob = c(0.35, 0.40), n = NULL, power = 0.80, alpha = 0.05,
                          alternative = "two.one.sided", std.error = "alternative", arcsine = FALSE, correct = FALSE,
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                          "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.45, null.prob = c(0.35, 0.40), delta = 0.05 * c(2, 1),
                      odds.ratio = c(1.519480519, 1.227272727), mean = 0, sd = 1, null.mean = c(-5.606633877, -2.803316939),
                      null.sd = 1, z.alpha = c(-7.5665979, -0.8433530), power = 0.800484462, n = 778))
    expect_equal(crrRes, pwrss.z.prop(p = 0.45, p0 = c(0.35, 0.40), alpha = 0.05, power = 0.80, alternative = "equivalent", verbose = FALSE))

    crrRes <- power.z.oneprop(null.prob = 0.50, req.sign = "-", alpha = 0.05, power = 0.80, n = 617, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, req.sign = "-", null.prob = 0.50, n = 617, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = FALSE, correct = FALSE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.450033949, null.prob = 0.50, delta = -0.049966051, odds.ratio = 0.818294054,
                      mean = -2.49475, sd = 1, null.mean = 0, null.sd = 1.00503092, z.alpha = -1.653128762,
                      power = 0.8, n = 617))

    crrRes <- power.z.oneprop(null.prob = 0.50, req.sign = "-", alpha = 0.05, n = 621, power = 0.80,
                              alternative = "one.sided", arcsine = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, req.sign = "-", null.prob = 0.50, n = 621, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = TRUE, correct = FALSE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.4500281, null.prob = 0.50, delta = -0.049971904, odds.ratio = 0.8182747,
                      mean = -2.494751942, sd = 1, null.mean = 0, null.sd = 1.00503211, z.alpha = -1.6531307,
                      power = 0.8, n = 621))

    crrRes <- power.z.oneprop(null.prob = 0.50, req.sign = "-", alpha = 0.05, n = 637, power = 0.80,
                              alternative = "one.sided", correct = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, req.sign = "-", null.prob = 0.50, n = 637, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", std.error = "null", arcsine = FALSE, correct = TRUE, ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.450039651, null.prob = 0.50, delta = -0.049960349, odds.ratio = 0.818312904,
                      mean = -2.4947482, sd = 1, null.mean = 0, null.sd = 1.00502977, z.alpha = -1.653126859,
                      power = 0.800000026, n = 637))

    crrRes <- suppressWarnings(power.z.oneprop(null.prob = c(0.60, 0.70), req.sign = "0", alpha = 0.05, n = 800,
                                               power = 0.80, alternative = "two.one.sided", correct = TRUE,
                                               std.error = "alternative", verbose = 0))
    expect_s3_class(crrRes, c("pwrss", "z", "oneprop"))
    expect_named(crrRes,
                 c("parms", "test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, req.sign = "0", null.prob = c(0.60, 0.70), n = 800, power = 0.80, alpha = 0.05,
                          alternative = "two.one.sided", std.error = "alternative", arcsine = FALSE, correct = TRUE,
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob", "null.prob", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", prob = 0.65042911, null.prob = c(0.60, 0.70), delta = c(0.05042911, -0.04957089),
                      odds.ratio = c(1.240433395, 0.797421468), mean = 0, sd = 1, null.mean = c(-2.95421812, 2.90331135),
                      null.sd = 1, z.alpha = c(-1.30936449, 1.25845773), power = 0.8006814, n = 800))


    expect_warning(power.z.oneprop(prob = 0.45, null.prob = 0.50, n = 500, alternative = "one.sided", arcsine = TRUE,
                                   correct = TRUE, verbose = 0),
                   "Continuity correction does not apply to arcsine transformation approach.")
    expect_warning(power.z.oneprop(prob = 0.45, null.prob = c(0.40, 0.50), n = 500, alternative = "two.one.sided", verbose = 0),
                   "`std.error` = \"null\" is ignored. Using \"alternative\" for equivalence or minimal effect testing.")
    expect_warning(power.z.oneprop(null.prob = c(0.40, 0.60), req.sign = "0", alpha = 0.05, n = 100, power = 0.80,
                                   alternative = "two.one.sided", correct = TRUE, std.error = "alternative", verbose = 0),
                   "The target power rate cannot be achieved within the null bounds.")
    expect_error(power.z.oneprop(prob = 0.45, null.prob = c(0.40, 0.50), alpha = 0.05, n = 500, alternative = "one.sided"),
                 "If `alternative` is \"two.sided\" or \"one.sided\", `null.prob` must be of length one.")
    expect_error(power.z.oneprop(prob = 0.45, null.prob = 0.50, alpha = 0.05, n = 500, alternative = "two.one.sided"),
                 paste("If `alternative` is \"two.one.sided\", `null.prob` must be of length two \\(lower and upper bounds,",
                       "with the upper bound being larger than the lower bound\\)."))
    expect_error(power.z.oneprop(prob = 0.45, null.prob = 1.2, alpha = 0.05, n = 500, alternative = "one.sided"),
                 "All elements of `null.prob` need to be valid proportion values \\(numeric, >= 0, and <= 1\\)")
    expect_error(pwrss.z.prop(p = 0.45, p0 = 0.50, margin = 1, alpha = 0.05, power = 0.80, alternative = "less", verbose = FALSE),
                 "Provide a reasonable margin consistent with `p` - `p0`.")
})

# power.exact.twoprops -------------------------------------------------------------------------------------------------
test_that("power.exact.twoprops works", {
    crrRes <- power.exact.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, n2 = 500, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "fisher"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "alternative", "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.60, req.sign = "+", n.ratio = 1, n2 = 500, power = NULL, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "alternative",
                          "z.alpha", "power", "n", "n.total")],
                 list(test = "exact", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = NA, sd = NA,
                      null.mean = NA, null.sd = NA, alternative = "one.sided", z.alpha = NA, power = 0.466794272,
                      n = c(n1 = 500, n2 = 500), n.total = 1000))

    crrRes <- power.exact.twoprops(prob1 = 0.70, prob2 = 0.60, alpha = 0.05, n2 = 500, paired = TRUE, rho.paired = 0.7, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.12285039, prob01 = 0.02285039, req.sign = "+", n.paired = 500, power = NULL, alpha = 0.05,
                      alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", prob10 = 0.12285039, prob01 = 0.02285039, delta = 0.10, odds.ratio = 5.376293, size = 73,
                      prob = 0.843169063, null.prob = 0.5, binom.alpha = c(27, 45), mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.034415914, power = 0.999990832, n.paired = 500))

    crrRes <- power.exact.twoprops(prob1 = 0.70, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "fisher"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "alternative", "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.70, prob2 = 0.60, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "alternative", "z.alpha", "power", "n", "n.total")],
                 list(test = "exact", delta = 0.10, odds.ratio = 1.55555556, mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      alternative = "one.sided", z.alpha = NA, power = 0.801117239, n = c(n1 = 302, n2 = 302), n.total = 604))

    crrRes <- power.exact.twoprops(prob1 = 0.70, prob2 = 0.60, alpha = 0.05, power = 0.80, paired = TRUE, rho.paired = 0.7, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.12285039, prob01 = 0.02285039, req.sign = "+", n.paired = NULL, power = 0.80, alpha = 0.05,
                      alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha", "mean", "sd",
                          "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", delta = 0.10, odds.ratio = 5.376293, size = 18, prob = 0.843169063, null.prob = 0.5,
                      binom.alpha = c(4, 13), mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.030883789, power = 0.807733751, n.paired = 120))

    crrRes <- power.exact.twoprops(prob1 = 0.70, req.sign = "-", alpha = 0.05, n2 = 302, power = 0.80,
                                   alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "fisher"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "alternative", "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.70, prob2 = NULL, req.sign = "-", n.ratio = 1, n2 = 302, power = 0.80, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "alternative", "z.alpha", "power", "n", "n.total")],
                 list(test = "exact", delta = 0.099829897, odds.ratio = 1.55445335, mean = NA, sd = NA, null.mean = NA,
                      null.sd = NA, alternative = "one.sided", z.alpha = NA, power = 0.79992023, n = c(n1 = 302, n2 = 302),
                      n.total = 604))

    crrRes <- power.exact.twoprops(prob2 = 0.60, alpha = 0.05, n2 = 302, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "fisher"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "alternative", "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = NULL, prob2 = 0.60, req.sign = "+", n.ratio = 1, n2 = 302, power = 0.8, alpha = 0.05,
                          alternative = "one.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "alternative", "z.alpha", "power", "n", "n.total")],
                 list(test = "exact", delta = 0.09986048, odds.ratio = 1.55452255, mean = NA, sd = NA, null.mean = NA,
                      null.sd = NA, alternative = "one.sided", z.alpha = NA, power = 0.8000729, n = c(n1 = 302, n2 = 302),
                      n.total = 604))

    crrRes <- power.exact.twoprops(prob1 = 0.70, req.sign = "-", alpha = 0.05, n2 = 120, power = 0.80, paired = TRUE,
                                   rho.paired = 0.7, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.122314849, prob01 = 0.02315924041, req.sign = "+", n.paired = 120, power = NULL,
                      alpha = 0.05, alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha", "mean", "sd",
                          "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", delta = 0.099155608, odds.ratio = 5.28147066, size = 18, prob = 0.840801613,
                      null.prob = 0.5, binom.alpha = c(4, 13), mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.030883789, power = 0.79997889, n.paired = 120))

    crrRes <- power.exact.twoprops(prob2 = 0.60, req.sign = "+", alpha = 0.05, n2 = 120, power = 0.80, paired = TRUE,
                                   rho.paired = 0.7, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "exact", "mcnemar"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.12241839, prob01 = 0.023204822, req.sign = "+", n.paired = 120, power = NULL, alpha = 0.05,
                      alternative = "two.sided", method = "exact", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha", "mean", "sd",
                          "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired")],
                 list(test = "exact", delta = 0.099213568, odds.ratio = 5.27555827, size = 18, prob = 0.84065163,
                      null.prob = 0.5, binom.alpha = c(4, 13), mean = NA, sd = NA, null.mean = NA, null.sd = NA,
                      z.alpha = NA, alpha = 0.030883789, power = 0.800035, n.paired = 120))
})

# power.z.twoprops (= pwrss.z.2props) ------------------------------------------------------------------------------
test_that("power.z.twoprops works", {
    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, n2 = 500, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = 500, power = NULL,
                          alpha = 0.05, alternative = "one.sided", arcsine = FALSE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 1.63517485,
                      sd = 1, null.mean = 0, null.sd = 1.001336, z.alpha = 1.64705116, power = 0.495262149,
                      n = c(n1 = 500, n2 = 500), n.total = 1000))
    expect_equal(crrRes, pwrss.z.2props(p1 = 0.65, p2 = 0.60, alpha = 0.05, n2 = 500, alternative = "less", verbose = FALSE))
    expect_equal(crrRes, pwrss.z.2prop(p1  = 0.65, p2 = 0.60, alpha = 0.05, n2 = 500, alternative = "less", verbose = FALSE))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", arcsine = FALSE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.48955035,
                      sd = 1, null.mean = 0, null.sd = 1.001336, z.alpha = 1.64705116, power = 0.8002457,
                      n = c(n1 = 1159, n2 = 1159), n.total = 2318))
    expect_equal(crrRes, pwrss.z.2props(p1 = 0.65, p2 = 0.60, alpha = 0.05, power = 0.80, alternative = "less", verbose = FALSE))
    expect_equal(crrRes, pwrss.z.2prop(p1  = 0.65, p2 = 0.60, alpha = 0.05, power = 0.80, alternative = "less", verbose = FALSE))

    crrRes <- power.z.twoprops(prob1 = 0.65, req.sign = "-", n2 = 1159, alpha = 0.05, power = 0.80,
                               alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = NULL, req.sign = "-", margin = 0, n.ratio = 1, n2 = 1159, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", arcsine = FALSE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.600022015, delta = 0.049977985, odds.ratio = 1.23798167,
                      mean = 2.488465902, sd = 1, null.mean = 0, null.sd = 1.001334843, z.alpha = 1.647049248,
                      power = 0.79994272, n = c(n1 = 1159, n2 = 1159), n.total = 2318))

    crrRes <- power.z.twoprops(prob2 = 0.60, req.sign = "+", n2 = 1159, alpha = 0.05, power = 0.80,
                               alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = NULL, prob2 = 0.60, req.sign = "+", margin = 0, n.ratio = 1, n2 = 1159, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", arcsine = FALSE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.64998433, prob2 = 0.6, delta = 0.049984331, odds.ratio = 1.238009969,
                      mean = 2.4887577, sd = 1, null.mean = 0, null.sd = 1.00133516, z.alpha = 1.64704976,
                      power = 0.80002426, n = c(n1 = 1159, n2 = 1159), n.total = 2318))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "one.sided",
                               arcsine = TRUE, std.error = "unpooled", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", arcsine = TRUE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "unpooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.48648363,
                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = 1.64485363, power = 0.800002455,
                      n = c(n1 = 1158, n2 = 1158), n.total = 2316))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "two.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", arcsine = FALSE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.8046943,
                      sd = 1, null.mean = 0, null.sd = 1.001336, z.alpha = 1.96258250806517 * c(-1, 1),
                      power = 0.800138245, n = c(n1 = 1471, n2 = 1471), n.total = 2942))
    expect_equal(crrRes, pwrss.z.2props(p1 = 0.65, p2 = 0.60, alpha = 0.05, power = 0.80, alternative = "not equal", verbose = 0))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "two.sided",
                               arcsine = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", arcsine = TRUE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.80244863,
                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = 1.959964 * c(-1, 1), power = 0.8002426,
                      n = c(n1 = 1471, n2 = 1471), n.total = 2942))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "two.sided",
                               correct = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", arcsine = FALSE, correct = TRUE, paired = FALSE,
                          rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.8049467,
                      sd = 1, null.mean = 0, null.sd = 1.001336, z.alpha = 1.96258250806517 * c(-1, 1),
                      power = 0.80020886, n = c(n1 = 1511, n2 = 1511), n.total = 3022))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "two.sided",
                               std.error = "unpooled", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", arcsine = FALSE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "unpooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.80183286,
                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = 1.959964 * c(-1, 1), power = 0.80007028,
                      n = c(n1 = 1468, n2 = 1468), n.total = 2936))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, alpha = 0.05, power = 0.80, alternative = "two.sided",
                               arcsine = TRUE, std.error = "unpooled", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = 0, n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", arcsine = TRUE, correct = FALSE, paired = FALSE,
                          rho.paired = 0.5, std.error = "unpooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.8024486,
                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = 1.959964 * c(-1, 1), power = 0.8002426,
                      n = c(n1 = 1471, n2 = 1471), n.total = 2942))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, margin = c(0, 0.10), alpha = 0.05, power = 0.80,
                               alternative = "two.one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = c(0, 0.10), n.ratio = 1, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", arcsine = FALSE, correct = FALSE,
                          paired = FALSE, rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.92874378,
                      sd = 1, null.mean = c(0, 5.85748755251), null.sd = 1.001336, z.alpha = c(1.64705116, 4.21043639),
                      power = 0.80004950, n = c(n1 = 1604, n2 = 1604), n.total = 3208))

    crrRes <- power.z.twoprops(prob1 = 0.65, prob2 = 0.60, margin = c(0.10, 0.20), alpha = 0.05, power = 0.80,
                               alternative = "two.one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = c(0.10, 0.20), n.ratio = 1, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", arcsine = FALSE, correct = FALSE,
                          paired = FALSE, rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.8046943, sd = 1,
                      null.mean = c(5.60938862, 11.21877723), null.sd = 1.001336, z.alpha = c(3.64680611, 13.18135974),
                      power = 0.80013731, n = c(n1 = 1471, n2 = 1471), n.total = 2942))

    crrRes <- pwrss.z.2props(p1 = 0.65, p2 = 0.60, margin = 0.10, alpha = 0.05, power = 0.80, alternative = "equivalent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd",
                   "z.alpha", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(prob1 = 0.65, prob2 = 0.6, req.sign = "+", margin = c(-0.10, 0.10), n.ratio = 1, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", arcsine = FALSE, correct = FALSE,
                          paired = FALSE, rho.paired = 0.5, std.error = "pooled", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob1", "prob2", "delta", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "z.alpha",
                          "power", "n", "n.total")],
                 list(test = "z", prob1 = 0.65, prob2 = 0.60, delta = 0.05, odds.ratio = 1.23809524, mean = 2.48955035,
                      sd = 1, null.mean = 4.97910071 * c(-1, 1), null.sd = 1.001336, z.alpha = 3.33204955 * c(-1, 1),
                      power = 0.8002457, n = c(n1 = 1159, n2 = 1159), n.total = 2318))

    crrRes <- power.z.twoprops(prob1 = 0.55, prob2 = 0.45, rho.paired = 0.41414141, power = 0.8, paired = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.20, prob01 = 0.10, req.sign = "+", n.paired = NULL, power = 0.80, alpha = 0.05,
                      alternative = "two.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2")],
                 list(test = "z", prob10 = 0.2, prob01 = 0.1, delta = 0.1, odds.ratio = 2, size = 71, prob = 2 / 3,
                      null.prob = 0.5, binom.alpha = c(26, 44), mean = 2.81314954, sd = 1, null.mean = 0, null.sd = 1,
                      z.alpha = 1.959964 * c(-1, 1), alpha = 0.05, power = 0.8032218, n.paired = 235, prob1 = 0.55,
                      prob2 = 0.45))

    crrRes <- power.z.twoprops(prob1 = 0.55, prob2 = 0.45, rho.paired = 0.41414141, n2 = 235, paired = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.20, prob01 = 0.10, req.sign = "+", n.paired = 235, power = NULL, alpha = 0.05,
                      alternative = "two.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2")],
                 list(test = "z", prob10 = 0.2, prob01 = 0.1, delta = 0.1, odds.ratio = 2, size = 71, prob = 2 / 3,
                      null.prob = 0.5, binom.alpha = c(26, 44), mean = 2.81314954, sd = 1, null.mean = 0, null.sd = 1,
                      z.alpha = 1.959964 * c(-1, 1), alpha = 0.05, power = 0.8032218, n.paired = 235, prob1 = 0.55,
                      prob2 = 0.45))

    crrRes <- power.z.twoprops(prob1 = 0.55, req.sign = "-", rho.paired = 0.41414141, n2 = 235, power = 0.8, paired = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.199770384, prob01 = 0.100172781, req.sign = "+", n.paired = 235, power = NULL, alpha = 0.05,
                      alternative = "two.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2")],
                 list(test = "z", prob10 = 0.199770384, prob01 = 0.100172781, delta = 0.099597604, odds.ratio = 1.994258149,
                      size = 71, prob = 0.66602746, null.prob = 0.5, binom.alpha = c(26, 44), mean = 2.80184946, sd = 1,
                      null.mean = 0, null.sd = 1, z.alpha = 1.959964 * c(-1, 1), alpha = 0.05, power = 0.80007397,
                      n.paired = 235, prob1 = 0.55, prob2 = 0.4504024))

    crrRes <- power.z.twoprops(prob2 = 0.45, req.sign = "+", rho.paired = 0.41414141, n2 = 235, power = 0.8, paired = TRUE, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "twoprops"))
    expect_named(crrRes,
                 c("parms", "test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                   "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2"))
    expect_equal(crrRes[["parms"]],
                 list(prob10 = 0.199770384, prob01 = 0.100172781, req.sign = "+", n.paired = 235, power = NULL, alpha = 0.05,
                      alternative = "two.sided", method = "approximate", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "prob10", "prob01", "delta", "odds.ratio", "size", "prob", "null.prob", "binom.alpha",
                          "mean", "sd", "null.mean", "null.sd", "z.alpha", "alpha", "power", "n.paired", "prob1", "prob2")],
                 list(test = "z", prob10 = 0.199770384, prob01 = 0.100172781, delta = 0.099597604, odds.ratio = 1.994258149,
                      size = 71, prob = 0.66602746, null.prob = 0.5, binom.alpha = c(26, 44), mean = 2.80184946, sd = 1,
                      null.mean = 0, null.sd = 1, z.alpha = 1.959964 * c(-1, 1), alpha = 0.05, power = 0.80007397,
                      n.paired = 235, prob1 = 0.5495976, prob2 = 0.45))

    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, margin = 1, n2 = 500, alternative = "one.sided"),
                 "Provide a reasonable `margin` consistent with `prob1` - `prob2`.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, margin = 0.05, n2 = 500, alternative = "one.sided"),
                 "The value of margin should be different from the prob1 - prob2 difference.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, n2 = 500, alternative = "one.sided", correct = TRUE, paired = TRUE),
                 "Continuity correction is currently not available for paired proportions.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, n2 = 500, alternative = "one.sided", arcsine = TRUE, correct = TRUE),
                 "Continuity correction does not apply to arcsine transformation approach.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, n2 = 500, alternative = "one.sided", arcsine = TRUE, paired = TRUE),
                 "Arcsine transformation is currently not available for paired proportions.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, margin = 0.1, n2 = 500, alternative = "one.sided", arcsine = TRUE),
                 "Arcsine transformation is currently not available for non-zero null.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, n2 = 500, alternative = "two.one.sided", paired = TRUE),
                 "Two one-sided tests are currently not available for paired proportions.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, n2 = 500, alternative = "two.one.sided", arcsine = TRUE),
                 "Arcsine transformation is currently not available for two one-sided tests.")
    expect_error(power.z.twoprops(prob1 = 0.65, prob2 = 0.60, margin = 0.10, n2 = 500, alternative = "two.one.sided"),
                 "Provide margins in the form of margin = c\\(lower, upper\\).")
    expect_error(power.z.twoprops(prob1 = 0.55, prob2 = 0.45, margin = 0.10, power = 0.8, verbose = 0),
                 "The value of margin should be different from the prob1 - prob2 difference.")
    expect_warning(power.z.twoprops(prob1 = 0.65, margin = 0.10, alpha = 0.05, n2 = 200, power = 0.80, verbose = 0),
                   "`margin` argument is ignored.")
})
