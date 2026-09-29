test_that("regression.logistic.R works", {
    # power.z.logistic (= pwrss.z.logistic / pwrss.z.logreg)
    crrRes <- power.z.logistic(base.prob = 0.15, prob = 0.20, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.2, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41666667, mean = 2.7836765, sd = 0.976071559, null.mean = 0, null.sd = 1,
                      vcf = 1, z.alpha = 1.959964 * c(-1, 1), power = 0.800639587, n = 511))
    expect_equal(crrRes, suppressMessages(pwrss.z.logistic(p0 = 0.15, p1 = 0.20, alpha = 0.05, power = 0.80, distribution = "normal", verbose = FALSE)))
    expect_equal(crrRes, suppressMessages(pwrss.z.logreg(p0 = 0.15, p1 = 0.20, alpha = 0.05, power = 0.80, distribution = "normal", verbose = FALSE)))

    crrRes <- power.z.logistic(base.prob = 0.15, prob = 0.20, alpha = 0.05, power = 0.80, distribution = "normal",
                               method = "demidenko", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.2, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41666667, mean = 2.802677943, sd = 1, null.mean = 0, null.sd = 1, vcf = 0,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800306736, n = 518))

    crrRes <- power.z.logistic(base.prob = 0.15, prob = 0.20, alpha = 0.05, power = 0.80, distribution = "normal",
                               method = "hsieh", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.2, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "hsieh", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41666667, mean = 2.80158522, sd = 1, null.mean = 0, null.sd = 1, vcf = NA,
                      z.alpha = 1.959964, power = 0.8, n = 508))

    crrRes <- power.z.logistic(base.prob = 0.15, prob = 0.20, alpha = 0.05, power = 0.80, distribution = "normal",
                               alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.2, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "one.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41666667, mean = 2.4690014, sd = 0.976071559, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.64485363, power = 0.800763585, n = 402))
    expect_equal(crrRes, pwrss.z.logistic(p0 = 0.15, p1 = 0.20, alpha = 0.05, power = 0.80, distribution = "normal",
                                          alternative = "less", verbose = FALSE))

    crrRes <- power.z.logistic(base.prob = 0.15, prob = 0.20, alpha = 0.05, n = 511, distribution = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.2, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 511, power = NULL, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41666667, mean = 2.7836765, sd = 0.976071559, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800639587, n = 511))

    crrRes <- power.z.logistic(base.prob = 0.15, alpha = 0.05, n = 511, power = 0.800, distribution = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 511, power = 0.8, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416231856, mean = 2.7812787, sd = 0.97611469, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.799942043, n = 511))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.783678356, sd = 0.976071526, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800640122, n = 511))

    crrRes <- power.z.logistic(beta0 = -1.734601, beta1 = 0.3483067, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = NULL, odds.ratio = NULL, beta0 = -1.734601, beta1 = 0.3483067,
                          req.sign = "+", n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05,
                          alternative = "two.sided", method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41666667, mean = 2.7836766, sd = 0.976071558, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800639615, n = 511))

    crrRes <- power.z.logistic(base.prob = 0.15, beta1 = 0.3483067, alpha = 0.05, power = 0.80,
                               distribution = list(dist = "normal", mean = 10, sd = 2), verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = 0.3483067,
                          req.sign = "+", n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05,
                          alternative = "two.sided", method = "demidenko(vc)",
                          distribution = list(dist = "normal", mean = 10, sd = 2), ceil.n = TRUE, verbose = 0,
                          utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41666667, mean = 2.7241170, sd = 0.9028072, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.801340544, n = 134))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, power = 0.80, distribution = "bernoulli", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "bernoulli", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.7957685, sd = 0.992726219, null.mean = 0, null.sd = 1, vcf = 0.85,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.80008686, n = 1816))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, power = 0.80, distribution = "bernoulli",
                               method = "demidenko", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko", distribution = "bernoulli", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.801919840, sd = 1, null.mean = 0, null.sd = 1, vcf = 0,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800094628, n = 1824))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, power = 0.80, distribution = "binomial",
                               method = "hsieh", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "hsieh", distribution = "binomial", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.8015852, sd = 1, null.mean = 0, null.sd = 1, vcf = NA,
                      z.alpha = 1.959964, power = 0.8, n = 1811))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, n = 1811, distribution = "binomial",
                               method = "hsieh", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 1811, power = NULL, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "hsieh", distribution = "binomial", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.80179607, sd = 1, null.mean = 0, null.sd = 1, vcf = NA,
                      z.alpha = 1.959964, power = 0.80005902, n = 1811))

    crrRes <- power.z.logistic(base.prob = 0.15, req.sign = "+", alpha = 0.05, n = 1811, power = 0.8,
                               distribution = "binomial", method = "hsieh", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 1811, power = 0.8, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "hsieh", distribution = "binomial", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416553726, mean = 2.8015852, sd = 1, null.mean = 0, null.sd = 1, vcf = NA,
                      z.alpha = 1.959964, power = 0.8, n = 1811))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, power = 0.80,
                               distribution = list(dist = "bernoulli", prob = 0.30), verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = list(dist = "bernoulli", prob = 0.3), ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.8293272, sd = 1.03295386, null.mean = 0, null.sd = 1, vcf = 0.85,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.80000377, n = 2114))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, r.squared.pred = 0.3345431, alpha = 0.05,
                               power = 0.80, distribution = list(dist = "bernoulli", prob = 0.25077), verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0.3345431, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = list(dist = "bernoulli", prob = 0.25077), ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.83846121, sd = 1.0435957, null.mean = 0, null.sd = 1,
                      vcf = 0.85, z.alpha = 1.959964 * c(-1, 1), power = 0.80005174, n = 3532))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, power = 0.80,
                               distribution = "poisson", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "poisson", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.84692934, sd = 1.0514093, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800555, n = 370))

    crrRes <- power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667, alpha = 0.05, n = 370,
                               distribution = "poisson", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.15, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 370, power = NULL, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "poisson", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.84692934, sd = 1.0514093, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800555, n = 370))

    crrRes <- power.z.logistic(base.prob = 0.20, odds.ratio = 1.416667, alpha = 0.05, power = 0.80,
                               distribution = "uniform", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.20, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "uniform", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.79958919, sd = 0.99751063, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800028722, n = 4402))

    crrRes <- power.z.logistic(base.prob = 0.20, odds.ratio = 1.416667, alpha = 0.05, n = 4402,
                               distribution = "uniform", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.20, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 4402, power = NULL, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "uniform", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.79958919, sd = 0.99751063, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800028722, n = 4402))

    crrRes <- power.z.logistic(base.prob = 0.20, alpha = 0.05, power = 0.8, n = 4402, distribution = "uniform", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.20, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 4402, power = 0.8, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "uniform", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41664954, mean = 2.7994858, sd = 0.9975108, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.79999967, n = 4402))

    crrRes <- power.z.logistic(base.prob = 0.30, odds.ratio = 1.416667, alpha = 0.05, power = 0.80,
                               distribution = "exponential", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.30, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "exponential", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.77189915, sd = 0.964603926, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800030653, n = 289))

    crrRes <- power.z.logistic(base.prob = 0.30, odds.ratio = 1.416667, alpha = 0.05, n = 289,
                               distribution = "exponential", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.30, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 289, power = NULL, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "exponential", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.77189915, sd = 0.964603926, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.800030653, n = 289))

    crrRes <- power.z.logistic(base.prob = 0.30, alpha = 0.05, power = 0.8, n = 289, distribution = "exponential", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.30, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 289, power = 0.8, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "exponential", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.41664968, mean = 2.77181467, sd = 0.9646098, null.mean = 0, null.sd = 1, vcf = 1,
                      z.alpha = 1.959964 * c(-1, 1), power = 0.8000047, n = 289))

    crrRes <- power.z.logistic(base.prob = 0.25, odds.ratio = 1.416667, alpha = 0.05, power = 0.80,
                               distribution = "lognormal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.25, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "lognormal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.60678920, sd = 0.75787988, null.mean = 0, null.sd = 1,
                      vcf = 0.75, z.alpha = 1.959964 * c(-1, 1), power = 0.80329972, n = 119))

    crrRes <- power.z.logistic(base.prob = 0.25, odds.ratio = 1.416667, alpha = 0.05, n = 119,
                               distribution = "lognormal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.25, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 119, power = NULL, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "lognormal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.60678920, sd = 0.75787988, null.mean = 0, null.sd = 1,
                      vcf = 0.75, z.alpha = 1.959964 * c(-1, 1), power = 0.80329972, n = 119))

    crrRes <- power.z.logistic(base.prob = 0.25, odds.ratio = 1.416667, alpha = 0.05, n = 119,
                               distribution = "lognormal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.25, odds.ratio = 1.416667, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 119, power = NULL, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "lognormal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.416667, mean = 2.60678920, sd = 0.75787988, null.mean = 0, null.sd = 1,
                      vcf = 0.75, z.alpha = 1.959964 * c(-1, 1), power = 0.80329972, n = 119))

    # example 29.4 from the GPower manual (first example, normal distribution)
    crrRes <- power.z.logistic(base.prob = 0.5, odds.ratio = 1.5, r.squared.pred = 0, alpha = 0.05, power = 0.95,
                               method = "hsieh", distribution = list(dist = "normal", mean = 0, sd = 1), verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = NULL, base.prob = 0.5, odds.ratio = 1.5, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.95, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "hsieh", distribution = list(dist = "normal", mean = 0, sd = 1), ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 1.5, mean = 3.60481761, sd = 1, null.mean = 0, null.sd = 1, vcf = NA,
                      z.alpha = 1.959964, power = 0.95, n = 317))
    expect_equal(power.z.logistic(base.prob = 0.5, odds.ratio = 1.5, r.squared.pred = 0, alpha = 0.05, power = 0.95,
                                  method = "demidenko(vc)", distribution = list(dist = "normal", mean = 0, sd = 1), verbose = 0)$n, 337)
    expect_equal(power.z.logistic(base.prob = 0.5, odds.ratio = 1.5, r.squared.pred = 0, alpha = 0.05, power = 0.95,
                                  method = "demidenko", distribution = list(dist = "normal", mean = 0, sd = 1), verbose = 0)$n, 355)
    # results are identical: z crit = 1.959964, power = 0.950486, n = hsieh: 317, demidenko: 337 (vc) / 355

    # example 29.4 from the GPower manual (second example, binomial distribution)
    crrRes <- power.z.logistic(base.prob = 0.1, prob = 0.05, r.squared.pred = 0, alpha = 0.05, power = 0.95,
                               method = "hsieh", distribution = list(dist = "binomial", prob = 0.5, size = 1), verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "logistic"))
    expect_named(crrRes, c("parms", "test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(prob = 0.05, base.prob = 0.1, odds.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.95, r.squared.predictor = 0, alpha = 0.05, alternative = "two.sided",
                          method = "hsieh", distribution = list(dist = "binomial", prob = 0.5, size = 1), ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "odds.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", odds.ratio = 0.473684211, mean = 3.60481761, sd = 1, null.mean = 0, null.sd = 1, vcf = NA,
                      z.alpha = 1.959964, power = 0.95, n = 1437))
    expect_equal(power.z.logistic(base.prob = 0.1, prob = 0.05, r.squared.pred = 0, alpha = 0.05, power = 0.95, method = "demidenko(vc)",
                                  distribution = list(dist = "binomial", prob = 0.5, size = 1), verbose = 0)$n, 1437)
    expect_equal(power.z.logistic(base.prob = 0.1, prob = 0.05, r.squared.pred = 0, alpha = 0.05, power = 0.95, method = "demidenko",
                                  distribution = list(dist = "binomial", prob = 0.5, size = 1), verbose = 0)$n, 1498)
    # results are identical: z crit = 1.959964, power = 0.950068, n = hsieh: 1437, demidenko: 1437 (vc) / 1498

    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.20, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0),
                 "`prob` can not have the same value as `base.prob`.")
    expect_error(power.z.logistic(beta0 = 0.15, beta1 = 0.15, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0),
                 "`beta1` can not have the same value as `beta0`.")
    expect_error(power.z.logistic(base.prob = 0.15, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0),
                 paste("Specify `base.prob` & `prob`\n  or `base.prob` & `odds.ratio`\n  or `base.prob` & `beta1`\n  or",
                       "`beta0` & `beta1`\n  or `base.prob` & `n` & `power` \\(the latter calculates `odds.ratio` as effect size\\)."))
    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.40, alpha = 0.05, distribution = "normal", verbose = 0),
                 "Exactly two of the parameters `\\(base.prob, odds.ratio\\)`, `n`, or `power` must be given, one has to be NULL.")
    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.40, alpha = 0.05, power = 0.80, n = 50, distribution = "normal", verbose = 0),
                 "Exactly two of the parameters `\\(base.prob, odds.ratio\\)`, `n`, or `power` must be given, one has to be NULL.")
    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.40, alpha = 0.05, power = 0.80,
                                  distribution = list(dist = "normal", mean = 0, sd = 1, err = 1), verbose = 0),
                 "Unknown input type for `distribution`")
    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.40, alpha = 0.05, power = 0.80,
                                  distribution = list(dist = "normal", mean = 0, sdev = 1), verbose = 0),
                 "Unknown input type for `distribution`")
    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.40, alpha = 0.05, power = 0.80, distribution = NA, verbose = 0),
                 "Unknown input type for `distribution`")
    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.40, alpha = 0.05, power = 0.80,
                                  distribution = list(dist = "binomial", size = 2, prob = 0.5), method = "hsieh", verbose = 0),
                 "Hsieh et al. \\(1998\\) is valid only for a binary covariate or a continuous covariate following normal distribution.")
    expect_error(power.z.logistic(base.prob = 0.20, prob = 0.40, alpha = 0.05, power = 0.80, distribution = "poisson",
                                  method = "hsieh", verbose = 0),
                 "Hsieh et al. \\(1998\\) is valid only for a binary covariate or a continuous covariate following normal distribution.")
    expect_error(power.z.logistic(base.prob = 0.1, alpha = 0.01, power = 0.99, n = 10, verbose = 2),
                 "Design is not feasible. Try `req.sign` = \"-\".")
    expect_error(power.z.logistic(base.prob = 0.1, alpha = 0.01, power = 0.99, n = 2, method = "hsieh", verbose = 0),
                 "Design is not feasible. Try `req.sign` = \"-\".")
    expect_message(power.z.logistic(base.prob = 0.15, prob = 0.20, odds.ratio = 1.5, alpha = 0.05, power = 0.80,
                                    distribution = "normal", verbose = 0),
                   "Using `base.prob` and `prob`, ignoring any specifications to `odds.ratio`, `beta0`, or `beta1`.")
    expect_message(power.z.logistic(base.prob = 0.15, odds.ratio = 1.5, beta1 = 0.80, alpha = 0.05, power = 0.80,
                                    distribution = "normal", verbose = 0),
                   "Using `base.prob` and `odds.ratio`, ignoring any specifications to `prob`, `beta0`, or `beta1`.")
    expect_message(power.z.logistic(base.prob = 0.15, beta1 = 0.80, beta0 = 0.2, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0),
                   "Using `base.prob` and `beta1`, ignoring any specifications to `prob`, `beta0`, or `odds.ratio`.")
    expect_message(power.z.logistic(prob = 0.2, beta1 = 0.80, beta0 = 0.2, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0),
                   "Using `beta0` and `beta1`, ignoring any specifications to `base.prob`, `prob`, or `odds.ratio`.")
})
