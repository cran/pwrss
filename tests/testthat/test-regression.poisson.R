# power.z.poisson (= pwrss.z.poisson) ----------------------------------------------------------------------------------
test_that("power.z.poisson / pwrss.z.poisson work", {
    crrRes <- power.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, power = 0.80, dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = NULL, rate.ratio = NULL, beta0 = 0.5, beta1 = -0.1, req.sign = "+", n = NULL,
                          power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05,
                          alternative = "two.sided", method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.80250964, sd = 0.99999976,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80025971, n = 474))
    expect_equal(crrRes, pwrss.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0))
    expect_equal(crrRes, pwrss.z.poisreg(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, power = 0.80, distribution = "normal", verbose = 0))

    crrRes <- power.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, alternative = "one.sided", power = 0.80,
                              dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = NULL, rate.ratio = NULL, beta0 = 0.5, beta1 = -0.1, req.sign = "+", n = NULL,
                          power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "one.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.4893946, sd = 0.99999976,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = -1.64485362695147, power = 0.800816471, n = 374))
    expect_equal(crrRes, pwrss.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, alternative = "less", power = 0.80,
                                         distribution = "normal", verbose = 0))
    expect_equal(crrRes, pwrss.z.poisreg(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, alternative = "less", power = 0.80,
                                         distribution = "normal", verbose = 0))

    crrRes <- power.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, method = "demidenko", power = 0.80,
                              dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = NULL, rate.ratio = NULL, beta0 = 0.5, beta1 = -0.1, req.sign = "+", n = NULL,
                          power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.80250964, sd = 1,
                      null.mean = 0, null.sd = 1, vcf = 0, z.alpha = c(-1.959964, 1.959964), power = 0.80025966, n = 474))

    crrRes <- power.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, method = "signorini", power = 0.80,
                              dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = NULL, rate.ratio = NULL, beta0 = 0.5, beta1 = -0.1, req.sign = "+", n = NULL,
                          power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "two.sided",
                          method = "signorini", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.80, sd = 0.997503122,
                      null.mean = 0, null.sd = 1, vcf = NA, z.alpha = c(-1.959964, 1.959964), power = 0.800145761, n = 784))
#    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
#                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.80251031, sd = 1.00000024,
#                      null.mean = 0, null.sd = 1, vcf = NA, z.alpha = c(-1.959964, 1.959964), power = 0.80025979, n = 474))

    crrRes <- power.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, n = 474, dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = NULL, rate.ratio = NULL, beta0 = 0.5, beta1 = -0.1, req.sign = "+", n = 474,
                          power = NULL, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.80250964, sd = 0.99999976,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80025971, n = 474))

    crrRes <- power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80,
                              dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = exp(0.50), rate.ratio = exp(-0.10), beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05,
                          alternative = "two.sided", method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.80250964, sd = 0.99999976,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80025971, n = 474))

    crrRes <- power.z.poisson(base.rate = exp(0.50), req.sign = "-", alpha = 0.05, power = 0.80, n = 474,
                              dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = exp(0.50), rate.ratio = NULL, beta0 = NULL, beta1 = NULL, req.sign = "-", n = 474,
                          power = 0.8, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "normal", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = 0.9048634, mean = -2.80170092, sd = 0.99999976,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.8000334, n = 474))

    crrRes <- power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80,
                              dist = list(dist = "normal", mean = 10, sd = 2), verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = exp(0.50), rate.ratio = exp(-0.10), beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05,
                          alternative = "two.sided", method = "demidenko(vc)",
                          distribution = list(dist = "normal", mean = 10, sd = 2), ceil.n = TRUE, verbose = 0,
                          utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.8055045, sd = 0.999998665,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.8010967, n = 318))

    crrRes <- power.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, power = 0.80, dist = "bernoulli", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = NULL, rate.ratio = NULL, beta0 = 0.5, beta1 = -0.1, req.sign = "+", n = NULL,
                          power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "bernoulli", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.800627935, sd = 0.9987513,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80002719, n = 2003))

    crrRes <- power.z.poisson(beta0 = 0.50, beta1 = -0.10, alpha = 0.05, n = 2003, dist = "bernoulli", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = NULL, rate.ratio = NULL, beta0 = 0.5, beta1 = -0.1, req.sign = "+", n = 2003,
                          power = NULL, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "two.sided",
                          method = "demidenko(vc)", distribution = "bernoulli", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.800627935, sd = 0.9987513,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80002719, n = 2003))

    crrRes <- power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80,
                              dist = "bernoulli", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = exp(0.50), rate.ratio = exp(-0.10), beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05,
                          alternative = "two.sided", method = "demidenko(vc)", distribution = "bernoulli", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.800627935, sd = 0.9987513,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80002719, n = 2003))

    crrRes <- power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80,
                              dist = list(dist = "bernoulli", prob = 0.30), verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = exp(0.50), rate.ratio = exp(-0.10), beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = NULL, power = 0.80, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05,
                          alternative = "two.sided", method = "demidenko(vc)",
                          distribution = list(dist = "bernoulli", prob = 0.3), ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.50), rate.ratio = exp(-0.10), mean = -2.78435929, sd = 0.979183926,
                      null.mean = 0, null.sd = 1, vcf = 1, z.alpha = c(-1.959964, 1.959964), power = 0.80008449, n = 2404))

    # example 30.3 ("binomial") from the GPower manual
    crrRes <- power.z.poisson(base.rate = 0.85, rate.ratio = 1.3, alpha = 0.05, power = 0.95, alternative = "one.sided",
                              method = "signorini", dist = "binomial", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = 0.85, rate.ratio = 1.3, beta0 = NULL, beta1 = NULL, req.sign = "+", n = NULL,
                          power = 0.95, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05, alternative = "one.sided",
                          method = "signorini", distribution = "binomial", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = 0.85, rate.ratio = 1.3, mean = 3.193011903, sd = 0.940539943,
                      null.mean = 0, null.sd = 1, vcf = NA, z.alpha = 1.64485363, power = 0.95012135, n = 697))
    expect_equal(power.z.poisson(base.rate = 0.85, rate.ratio = 1.3, alpha = 0.05, power = 0.95, alternative = "one.sided",
                                 method = "demidenko(vc)", dist = "binomial", verbose = 0)$n, 649)
    expect_equal(power.z.poisson(base.rate = 0.85, rate.ratio = 1.3, alpha = 0.05, power = 0.95, alternative = "one.sided",
                                 method = "demidenko", dist = "binomial", verbose = 0)$n, 655)
    # results are identical to GPower

    # example 30.3 ("normal") from the GPower manual
    crrRes <- power.z.poisson(base.rate = exp(0.5), rate.ratio = exp(-0.1), alpha = 0.05, n = 200, alternative = "two.sided",
                              method = "signorini", dist = "normal", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "z", "poisson"))
    expect_named(crrRes,
                 c("parms", "test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n"))
    expect_identical(crrRes[["parms"]],
                     list(base.rate = exp(0.5), rate.ratio = exp(-0.1), beta0 = NULL, beta1 = NULL, req.sign = "+",
                          n = 200, power = NULL, r.squared.predictor = 0, mean.exposure = 1, alpha = 0.05,
                          alternative = "two.sided", method = "signorini", distribution = "normal", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "base.rate", "rate.ratio", "mean", "sd", "null.mean", "null.sd", "vcf", "z.alpha", "power", "n")],
                 list(test = "z", base.rate = exp(0.5), rate.ratio = exp(-0.1), mean = -1.414213562, sd = 0.997503122,
                      null.mean = 0, null.sd = 1, vcf = NA, z.alpha = 1.95996398 * c(-1, 1), power = 0.292508332, n = 200))
    # results are NOT identical to GPower: power is expected to be 0.444593


    # example 30.3 from the GPower manual, Table 2
    pois_pwr <- function(b1, mth, dst) {
      power.z.poisson(beta0 = 0.5, beta1 = b1, alpha = 0.05, n = 200, alternative = "two.sided", method = mth, dist = dst, verbose = 0)$power
    }
    expect_equal(vapply(seq(2, 4) * -0.05, pois_pwr, numeric(1), "demidenko",     "normal"),
                 c(0.444590667, 0.782076179, 0.956208934))
    expect_equal(vapply(seq(2, 4) * -0.05, pois_pwr, numeric(1), "demidenko(vc)", "normal"),
                 c(0.444590654, 0.782076336, 0.956209146))
    expect_equal(vapply(seq(2, 4) * -0.05, pois_pwr, numeric(1), "signorini",     "normal"),
                 c(0.29250833, 0.5644723, 0.8098093))
#                c(0.443, 0.779, 0.954)) # values in Table 2
    expect_equal(vapply(c(1, 2, 4) * -0.2, pois_pwr, numeric(1), "demidenko",     "uniform"),
                 c(0.169362111, 0.474561989, 0.928484710))
    expect_equal(vapply(c(1, 2, 4) * -0.2, pois_pwr, numeric(1), "demidenko(vc)", "uniform"),
                 c(0.169105785, 0.474456655, 0.931624379))
#    expect_equal(vapply(c(1, 2, 4) * -0.2, pois_pwr, numeric(1), "signorini",     "uniform"),
#                 c(0.169863789, 0.477672343, 0.932597708))
    expect_equal(vapply(seq(1, 3) * -0.05, pois_pwr, numeric(1), "demidenko",     "lognormal"),
                 c(0.320220319, 0.695340894, 0.889820796))
    expect_equal(vapply(seq(1, 3) * -0.05, pois_pwr, numeric(1), "demidenko(vc)", "lognormal"),
                 c(0.290596024, 0.746178854, 0.955110662))
#    expect_equal(vapply(seq(1, 3) * -0.05, pois_pwr, numeric(1), "signorini",     "lognormal"),
#                 c(0.478044368, 0.873283908, 0.977554656))
    expect_equal(vapply(seq(1, 2) * -0.20, pois_pwr, numeric(1), "demidenko",     list(dist = "poisson", lambda = 0.5)),
                 c(0.602872827, 0.972089832))
    expect_equal(vapply(seq(1, 2) * -0.20, pois_pwr, numeric(1), "demidenko(vc)", list(dist = "poisson", lambda = 0.5)),
                 c(0.613395917, 0.990250577))
#    expect_equal(vapply(seq(1, 2) * -0.20, pois_pwr, numeric(1), "signorini",     list(dist = "poisson", lambda = 0.5)),
#                 c(0.672691558, 0.988326785))
    expect_equal(vapply(seq(1, 2) * -0.20, pois_pwr, numeric(1), "demidenko",     list(dist = "bernoulli", prob = 0.2)),
                 c(0.267784010, 0.691931262))
    expect_equal(vapply(seq(1, 2) * -0.20, pois_pwr, numeric(1), "demidenko(vc)", list(dist = "bernoulli", prob = 0.2)),
                 c(0.254388220, 0.716318213))
    expect_equal(vapply(seq(1, 2) * -0.20, pois_pwr, numeric(1), "signorini",     list(dist = "bernoulli", prob = 0.2)),
                 c(0.320887143, 0.788431959))
    expect_equal(vapply(c(-0.4, -0.7),     pois_pwr, numeric(1), "demidenko",     list(dist = "exponential", rate = 3)),
                 c(0.518688552, 0.871497202))
    expect_equal(vapply(c(-0.4, -0.7),     pois_pwr, numeric(1), "demidenko(vc)", list(dist = "exponential", rate = 3)),
                 c(0.521140107, 0.918937579))
#   expect_equal(vapply(c(-0.4, -0.7),     pois_pwr, numeric(1), "signorini",     list(dist = "exponential", rate = 3)),
#                 c(0.609332289, 0.933746771))
    # results are not identical to GPower for Signorini; for Demidenko (both with and without VC), the values are correct
    # cf. the GPower 3.1 manual (June 1, 2023; 30.2 Options, p. 76):
    # Results of Monte-Carlo simulations indicate that the accuracy of the procedure based on the work of Demidenko
    # (2007) is comparable to that of the enumeration procedure for N > 200, whereas errors of the procedure proposed
    # by Signorini (1991) can be quite large. We thus recommend to use the procedure based on Demidenko (2007) as the
    # standard procedure.

    # table from PASS (p. 870-5)
    pois_pwr <- function(n, rr) {
      pwrss::power.z.poisson(base.rate = 1, rate.ratio = rr, alpha = 0.05, n = n, alternative = "two.sided", method = "signorini",
                             mean.exposure = 1, dist = list(dist = "normal", mean = 3.2, sd = 2.1), verbose = 0)$power
    }
    expect_equal(vapply(seq(5, 50, 5), pois_pwr, numeric(1), 1.3),
                 c(0.116036529, 0.360429468, 0.612368678, 0.795997955, 0.904027378,
                   0.958760537, 0.983554409, 0.993840309, 0.997813065, 0.999258670))
    expect_equal(vapply(seq(5, 50, 5), pois_pwr, numeric(1), 1.5),
                 c(0.448902366, 0.953538538, 0.998922421, 0.999988721, 0.999999931, 1, 1, 1, 1, 1))
    # table from PASS (p. 870-8)
    pois_n <- function(power) {
      pwrss::power.z.poisson(base.rate = 0.85, rate.ratio = 1.3, alpha = 0.05, power = power, alternative = "one.sided",
                             method = "signorini", mean.exposure = 1, dist = "binomial", verbose = 0)$n
    }
    expect_equal(vapply(c(0.80, 0.90, 0.95), pois_n, numeric(1)),
                 c(406, 556, 697))

    expect_error(power.z.poisson(beta0 = 0.50, alpha = 0.05, power = 0.80, verbose = 0),
                 paste("Specify `base.rate` & `rate.ratio`\n  or `beta0` & `beta1`\n  or `base.rate` & `n` & `power`",
                       "\\(the latter calculates `rate.ratio` as effect size\\)."))
    expect_message(power.z.poisson(beta0 = 0.50, beta1 = -0.10, base.rate = exp(0.50), alpha = 0.05, power = 0.80, verbose = 0),
                   "Using `beta0` and `beta1`, ignoring any specifications to `base.rate` or `rate.ratio`.")
    expect_message(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), beta0 = 0.50, alpha = 0.05, power = 0.80, verbose = 0),
                   "Using `base.rate` and `rate.ratio`, ignoring any specifications to `beta0` or `beta1`.")
    expect_message(power.z.poisson(base.rate = exp(0.50), beta0 = 0.50, alpha = 0.05, power = 0.80, n = 474, verbose = 0),
                   "Calculating the effect size \\(`rate.ratio`\\), ignoring any specifications to `rate.ratio`, `beta0` or `beta1`.")
    expect_error(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(0.50), alpha = 0.05, power = 0.80, verbose = 0),
                 "`beta0` / `base.rate` can not have the same value as `beta1` / `rate.ratio`.")
    expect_error(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, verbose = 0),
                 "Exactly two of the parameters `\\(base.rate, rate.ratio\\)`, `n`, or `power` must be given, one has to be NULL.")
    expect_error(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80, n = 200, verbose = 0),
                 "Exactly two of the parameters `\\(base.rate, rate.ratio\\)`, `n`, or `power` must be given, one has to be NULL.")
    expect_error(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80,
                                  distribution = list(dist = "normal", mean = 0, sd = 1, err = 1), verbose = 0),
                 "Unknown input type for `distribution`")
    expect_error(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80,
                                  distribution = list(dist = "normal", mean = 0, sdev = 1), verbose = 0),
                 "Unknown input type for `distribution`")
    expect_error(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80, distribution = NA, verbose = 0),
                 "Unknown input type for `distribution`")
    expect_error(power.z.poisson(base.rate = exp(0.50), rate.ratio = exp(-0.10), alpha = 0.05, power = 0.80,
                                 method = "signorini", dist = "uniform", verbose = 0),
                 "Distribution type is not supported by the Signorini procedure.")
})
