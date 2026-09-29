test_that("regression.mediation.R works", {
    set.seed(1)

    # power.z.mediation (= pwrss.z.mediation)
    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 200, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 200,
                          power = NULL, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.1,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.0725, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = 2.54702338, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.721421379, n = 200))

    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.1,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.0725, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = 2.80172572, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80004029, n = 242))

#    crrRes <- power.z.mediation(beta.a = 0.25, beta.cp = 0.10, power = 0.80, n = 242, verbose = 0)
#    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
#    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
#                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
#                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
#    expect_identical(crrRes[["parms"]],
#                     list(beta.a = 0.25, beta.b = NULL, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
#                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 242,
#                          power = 0.80, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
#                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
#    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
#                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
#                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
#                 list(test = "z", beta.a = 0.25, beta.b = 0.24997071, beta.indirect = 0.062492677, beta.cp = 0.1,
#                      r.squared.mediator = 0.0625, r.squared.outcome = 0.072485355, std.beta.a = 0.25,
#                      std.beta.b = 0.24997071, std.beta.cp = 0.1, std.beta.indirect = 0.062492677, mean = 2.8015458,
#                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.799989922,
#                      n = 242))

#    crrRes <- power.z.mediation(beta.b = 0.25, beta.cp = 0.10, power = 0.80, n = 242, verbose = 0)
#    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
#    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
#                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
#                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
#    expect_identical(crrRes[["parms"]],
#                     list(beta.a = NULL, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
#                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 242,
#                          power = 0.80, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
#                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
#    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
#                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
#                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
#                 list(test = "z", beta.a = 0.24997071, beta.b = 0.25, beta.indirect = 0.062492677, beta.cp = 0.1,
#                      r.squared.mediator = 0.0625, r.squared.outcome = 0.072485355, std.beta.a = 0.25,
#                      std.beta.b = 0.24997071, std.beta.cp = 0.1, std.beta.indirect = 0.062492677, mean = 2.8015458,
#                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.799989922,
#                      n = 242))

    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 500, method = "sobel", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 500,
                          power = NULL, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.1,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.0725, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = 4.0271976, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.980643926, n = 500))

    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 500, method = "aroian", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 500,
                          power = NULL, alpha = 0.05, alternative = "two.sided", method = "aroian", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.1,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.0725, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = 3.9965349, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.97915347, n = 500))

    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 500, method = "goodman", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 500,
                          power = NULL, alpha = 0.05, alternative = "two.sided", method = "goodman", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.1,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.0725, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = 4.058577, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.98207449, n = 500))

    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 500, method = "joint", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 500,
                          power = NULL, alpha = 0.05, alternative = "two.sided", method = "joint", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "joint", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.10,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.0725, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = NA, sd = 1, null.mean = NA, null.sd = 1,
                      z.alpha = NA, power = 0.999805534, n = 500))

    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 500, method = "monte.carlo", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL,
                          n = 500, power = NULL, alpha = 0.05, alternative = "two.sided", method = "monte.carlo",
                          n.simulation = 1000, n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "monte.carlo", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.10,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.0725, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = NA, sd = 1, null.mean = NA, null.sd = 1,
                      z.alpha = NA, power = 0.999, n = 500))

    crrRes <- power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, r.squared.outcome = 0.50, power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.25, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 1,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = 0.50, n = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.25, beta.b = 0.25, beta.indirect = 0.0625, beta.cp = 0.10,
                      r.squared.mediator = 0.0625, r.squared.outcome = 0.5, std.beta.a = 0.25, std.beta.b = 0.25,
                      std.beta.cp = 0.1, std.beta.indirect = 0.0625, mean = 2.80378134, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.800615212, n = 185))

    crrRes <- power.z.mediation(beta.a = 0.40, beta.b = 0.25, beta.cp = 0.10, sd.predictor = 0.5, n = 200, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.40, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 0.50,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = 200,
                          power = NULL, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.40, beta.b = 0.25, beta.indirect = 0.10, beta.cp = 0.10,
                      r.squared.mediator = 0.04, r.squared.outcome = 0.065, std.beta.a = 0.20, std.beta.b = 0.25,
                      std.beta.cp = 0.05, std.beta.indirect = 0.05, mean = 2.247805948, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.613279038, n = 200))

    crrRes <- power.z.mediation(beta.a = 0.40, beta.b = 0.25, beta.cp = 0.10, sd.predictor = 0.5, power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.40, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 0.50,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = NULL, n = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.40, beta.b = 0.25, beta.indirect = 0.10, beta.cp = 0.10,
                      r.squared.mediator = 0.04, r.squared.outcome = 0.065, std.beta.a = 0.20, std.beta.b = 0.25,
                      std.beta.cp = 0.05, std.beta.indirect = 0.05, mean = 2.8030059, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.800398454, n = 311))

    crrRes <- power.z.mediation(beta.a = 0.40, beta.b = 0.25, beta.cp = 0.10, r.squared.outcome = 0.50, sd.predictor = 0.5,
                                power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.40, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1, sd.predictor = 0.50,
                          sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL, r.squared.outcome = 0.50, n = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", method = "sobel", n.simulation = 1000,
                          n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.40, beta.b = 0.25, beta.indirect = 0.10, beta.cp = 0.10,
                      r.squared.mediator = 0.04, r.squared.outcome = 0.5, std.beta.a = 0.20, std.beta.b = 0.25,
                      std.beta.cp = 0.05, std.beta.indirect = 0.05, mean = 2.8027968, sd = 1, null.mean = 0,
                      null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.800340, n = 254))

    crrRes <- power.z.mediation(beta.a = 0.40, beta.b = 0.25, beta.cp = 0.10, r.squared.outcome = 0.50,
                                sd.predictor = sqrt(2 / 3 * 1 / 3), power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "mediation"))
    expect_named(crrRes, c("parms", "test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                                  "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect",
                                  "mean", "sd", "null.mean", "null.sd", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(beta.a = 0.40, beta.b = 0.25, ab.ratio = 1, req.sign = "+", beta.cp = 0.1,
                          sd.predictor = sqrt(2 / 3 * 1 / 3), sd.mediator = 1, sd.outcome = 1, r.squared.mediator = NULL,
                          r.squared.outcome = 0.50, n = NULL, power = 0.80, alpha = 0.05, alternative = "two.sided",
                          method = "sobel", n.simulation = 1000, n.draws = 1000, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "beta.a", "beta.b", "beta.indirect", "beta.cp", "r.squared.mediator",
                          "r.squared.outcome", "std.beta.a", "std.beta.b", "std.beta.cp", "std.beta.indirect", "mean",
                          "sd", "null.mean", "null.sd", "z.alpha", "power", "n")],
                 list(test = "z", beta.a = 0.40, beta.b = 0.25, beta.indirect = 0.10, beta.cp = 0.10,
                      r.squared.mediator = 0.035555556, r.squared.outcome = 0.50, std.beta.a = 0.188561808,
                      std.beta.b = 0.25, std.beta.cp = 0.047140452, std.beta.indirect = 0.047140452, mean = 2.80658531,
                      sd = 1, null.mean = 0, null.sd = 1, z.alpha = c(-1.959964, 1.959964), power = 0.801397824,
                      n = 279))

    expect_error(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, verbose = 0),
                 "Exactly one of the parameters `n` or `power` must be given, one has to be NULL.")
    expect_error(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 200, power = 0.8, verbose = 0),
                 "Exactly one of the parameters `n` or `power` must be given, one has to be NULL.")
#    expect_error(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, verbose = 0),
#                 "Exactly two of the parameters `\\(beta.a, beta.b\\)`, `n`, or `power` must be given, one has to be NULL.")
#    expect_error(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, n = 200, power = 0.8, verbose = 0),
#                 "Exactly two of the parameters `\\(beta.a, beta.b\\)`, `n`, or `power` must be given, one has to be NULL.")
    expect_error(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, power = 0.8, method = "joint", verbose = 0),
                 "Sample size calculation not supported by this method")
    crrWrn <- capture_warnings(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, r.squared.outcome = 0, n = 200, verbose = 0, utf = FALSE))
    expect_match(crrWrn[1], "Ignoring any specification to `beta.cp`.")
    expect_match(crrWrn[2], "Specified `r.squared.outcome` is smaller than the base `r.squared.outcome`.")
    expect_warning(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, r.squared.mediator = 0.05, n = 200, verbose = 0),
                   "Specified `r.squared.mediator` is smaller than the base `r.squared.mediator`.")
    expect_warning(power.z.mediation(beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10, r.squared.outcome = 0.07, n = 200, verbose = 0),
                   "Specified `r.squared.outcome` is smaller than the base `r.squared.outcome`.")

    mthNme <- c("sobel", "aroian", "goodman", "joint", "mc")
    crrRes <- pwrss.z.mediation(a = 0.25, b = 0.25, cp = 0.10, n = 500, verbose = FALSE)
    expect_equal(crrRes, pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10, n = 500, verbose = FALSE))
    expect_s3_class(crrRes, c("pwrss", "z", "med", "defunct"))
    expect_named(crrRes, c("parms", "test", "ncp", "power", "n"))
    expect_equal(crrRes[["parms"]],
                 list(a = 0.25, b = 0.25, cp = 0.1, sdx = 1, sdm = 1, sdy = 1, r2m.x = 0.0625, r2y.mx = 0.0725, n = 500,
                      power = NULL, alpha = 0.05, alternative = "two.sided", mc = TRUE, nsims = 1000, ndraws = 1000,
                      verbose = FALSE))
    expect_equal(crrRes[c("test", "ncp", "power", "n")],
                 list(test = c("z", "joint", "monte.carlo"), ncp = NULL, power = setNames(c(0.980643926, 0.979153469, 0.982074489,
                      0.999805534, 1), mthNme), n = setNames(rep(500, 5), mthNme)))

    crrAsc <- capture.output(pwrss.z.mediation(a = 0.25, b = 0.25, cp = 0.10, n = 500))
    expect_equal(crrAsc, c(" Indirect Effect in a Mediation Model",
                           " ====================================",
                           "      method non-centrality   n     power",
                           "       Sobel       4.027198 500 0.9806439",
                           "      Aroian       3.996535 500 0.9791535",
                           "     Goodman       4.058577 500 0.9820745",
                           "       Joint             NA 500 0.9998055",
                           " Monte Carlo             NA 500 1.0000000",
                           " ------------------------------------",
                           " Type 1 error rate: 0.05", ""))

    crrRes <- pwrss.z.mediation(a = 0.25, b = 0.25, cp = 0.10, power = 0.90, verbose = FALSE)
    expect_s3_class(crrRes, c("pwrss", "z", "med", "defunct"))
    expect_named(crrRes, c("parms", "test", "ncp", "power", "n"))
    expect_equal(crrRes, pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10, power = 0.90, verbose = FALSE))
    expect_equal(crrRes[["parms"]],
                 list(a = 0.25, b = 0.25, cp = 0.1, sdx = 1, sdm = 1, sdy = 1, r2m.x = 0.0625, r2y.mx = 0.0725, n = NULL,
                      power = 0.90, alpha = 0.05, alternative = "two.sided", mc = TRUE, nsims = 1000, ndraws = 1000,
                      verbose = FALSE))
    expect_equal(crrRes[c("test", "ncp", "power", "n")],
                 list(test = c("z", "joint", "monte.carlo"), ncp = NULL, power = setNames(c(0.900055537, 0.900469752,
                      0.900836948, NA, NA), mthNme), n = setNames(c(324, 332, 317, NA, NA), mthNme)))

    crrRes <- pwrss.z.mediation(a = 0.25, b = 0.25, cp = 0.10, n = 500, alternative = "less", verbose = FALSE)
    expect_s3_class(crrRes, c("pwrss", "z", "med", "defunct"))
    expect_named(crrRes, c("parms", "test", "ncp", "power", "n"))
    expect_equal(crrRes[["parms"]],
                 list(a = 0.25, b = 0.25, cp = 0.1, sdx = 1, sdm = 1, sdy = 1, r2m.x = 0.0625, r2y.mx = 0.0725, n = 500,
                      power = NULL, alpha = 0.05, alternative = "one.sided", mc = TRUE, nsims = 1000, ndraws = 1000,
                      verbose = FALSE))
    expect_equal(crrRes[c("test", "ncp", "power", "n")],
                 list(test = c("z", "joint", "monte.carlo"), ncp = NULL, power = setNames(c(0.991398589, 0.990655608, 0.992104778,
                      0.999946621, 1), mthNme), n = setNames(rep(500, 5), mthNme)))
})
