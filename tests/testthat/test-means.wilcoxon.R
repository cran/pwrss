# power.np.wilcoxon (= pwrss.np.2groups) -------------------------------------------------------------------------------
test_that("power.np.wilcoxon / pwrss.np.2groups work", {
    crrRes <- power.np.wilcoxon(d = 0.25, power = 0.80, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.25, n = c(n1 = 265, n2 = 265), power = 0.801445824, t.alpha = 1.964680946 * c(-1, 1),
                      ncp = 2.81211864, null.ncp = 0, df = 504.112719))
    expect_equal(crrRes, pwrss.np.2groups(mu1 = 0.25, sd1 = 1, power = 0.80, verbose = FALSE))

    crrRes <- power.np.wilcoxon(d = 0.25, n2 = 265, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 265, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.25, n = c(n1 = 265, n2 = 265), power = 0.801445824, t.alpha = 1.964680946 * c(-1, 1),
                      ncp = 2.81211864, null.ncp = 0, df = 504.112719))
    expect_equal(crrRes, pwrss.np.2groups(mu1 = 0.25, sd1 = 1, n2 = 265, verbose = FALSE))

    crrRes <- power.np.wilcoxon(n2 = 265, power = 0.80, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 265, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.249538994, n = c(n1 = 265, n2 = 265), power = 0.8, t.alpha = 1.964680946 * c(-1, 1),
                      ncp = 2.806933, null.ncp = 0, df = 504.112719))

    crrRes <- power.np.wilcoxon(d = 0.20, power = 0.80, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = c(n1 = 412, n2 = 412), power = 0.8000253, t.alpha = 1.9629911 * c(-1, 1),
                      ncp = 2.8051061, null.ncp = 0, df = 784.862039))

    crrRes <- power.np.wilcoxon(d = 0.20, n2 = 412, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 412, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = c(n1 = 412, n2 = 412), power = 0.8000253, t.alpha = 1.9629911 * c(-1, 1),
                      ncp = 2.8051061, null.ncp = 0, df = 784.862039))

    crrRes <- power.np.wilcoxon(n2 = 412, power = 0.80, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 412, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.199993547, n = c(n1 = 412, n2 = 412), power = 0.8, t.alpha = 1.9629911 * c(-1, 1),
                      ncp = 2.80501562, null.ncp = 0, df = 784.862039))

    crrRes <- power.np.wilcoxon(req.sign = "-", n2 = 412, power = 0.80, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "-", n.ratio = 1, n2 = 412, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.199993547, n = c(n1 = 412, n2 = 412), power = 0.8, t.alpha = 1.9629911 * c(-1, 1),
                      ncp = -2.80501562, null.ncp = 0, df = 784.862039))

    crrRes <- power.np.wilcoxon(d = 0.20, power = 0.80, alternative = "two.sided", design = "independent",
                                method = "noether", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "z"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "noether", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df")],
                 list(test = "z", d = 0.20, n = c(n1 = 414, n2 = 414), power = 0.800274145, z.alpha = 1.959964 * c(-1, 1),
                      mean = 2.80256143, sd = 1, null.mean = 0, null.sd = 1, df = Inf))

    crrRes <- power.np.wilcoxon(d = 0.20, n2 = 414, alternative = "two.sided", design = "independent",
                                method = "noether", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "z"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 414, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "noether", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df")],
                 list(test = "z", d = 0.20, n = c(n1 = 414, n2 = 414), power = 0.800274145, z.alpha = 1.959964 * c(-1, 1),
                      mean = 2.80256143, sd = 1, null.mean = 0, null.sd = 1, df = Inf))

    crrRes <- power.np.wilcoxon(n2 = 414, power = 0.80, alternative = "two.sided", design = "independent",
                                method = "noether", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "z"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 414, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "normal",
                          method = "noether", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df")],
                 list(test = "z", d = 0.199929622, n = c(n1 = 414, n2 = 414), power = 0.8, z.alpha = 1.959964 * c(-1, 1),
                      mean = 2.80158179, sd = 1, null.mean = 0, null.sd = 1, df = Inf))

    crrRes <- power.np.wilcoxon(d = 0.25, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.25, n = c(n1 = 208, n2 = 208), power = 0.800182327, t.alpha = 1.64871793,
                      ncp = 2.491393743, null.ncp = 0, df = 395.25074))

    crrRes <- power.np.wilcoxon(d = 0.25, n2 = 208, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 208, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.25, n = c(n1 = 208, n2 = 208), power = 0.800182327, t.alpha = 1.64871793,
                      ncp = 2.491393743, null.ncp = 0, df = 395.25074))

    crrRes <- power.np.wilcoxon(n2 = 208, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 208, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.249934519, n = c(n1 = 208, n2 = 208), power = 0.8, t.alpha = 1.64871793,
                      ncp = 2.49074119, null.ncp = 0, df = 395.25074))

    crrRes <- power.np.wilcoxon(d = 0.20, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = c(n1 = 325, n2 = 325), power = 0.800613556, t.alpha = 1.64732018,
                      ncp = 2.491393743, null.ncp = 0, df = 618.704278))
    expect_equal(crrRes, pwrss.np.2groups(mu1 = 0.20, sd1 = 1, power = 0.8, alternative = "greater", verbose = FALSE))

    crrRes <- power.np.wilcoxon(d = 0.10, margin = -0.05, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.10, null.d = 0, margin = -0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.10, n = c(n1 = 576, n2 = 576), power = 0.8001976, t.alpha = 0.816102757,
                      ncp = 1.658371917, null.ncp = -0.829185959, df = 1098.078967))

    crrRes <- power.np.wilcoxon(d = 0.20, margin = -0.05, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = -0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = c(n1 = 208, n2 = 208), power = 0.800795363, t.alpha = 1.1486718,
                      ncp = 1.993115, null.ncp = -0.49827875, df = 395.25074))

    crrRes <- power.np.wilcoxon(d = 0.10, margin = 0.05, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.10, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.10, n = c(n1 = 5184, n2 = 5184), power = 0.800002944, t.alpha = 4.13322535,
                      ncp = 4.9751158, null.ncp =  2.487557876, df = 9898.7107))

    crrRes <- power.np.wilcoxon(d = 0.10, margin = 0.05, n2 = 5184, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.10, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 5184, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.10, n = c(n1 = 5184, n2 = 5184), power = 0.800002944, t.alpha = 4.13322535,
                      ncp = 4.9751158, null.ncp =  2.487557876, df = 9898.7107))

    crrRes <- power.np.wilcoxon(margin = 0.05, n2 = 5184, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 5184, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.099999789, n = c(n1 = 5184, n2 = 5184), power = 0.80000, t.alpha = 4.13322535,
                      ncp = 4.97510523, null.ncp = 2.48755788, df = 9898.71070))

    crrRes <- power.np.wilcoxon(d = 0.20, margin = 0.05, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = c(n1 = 578, n2 = 578), power = 0.800541149, t.alpha = 2.4783293,
                      ncp = 3.3224971, null.ncp = 0.83062427, df = 1101.8987))

    crrRes <- power.np.wilcoxon(d = 0.20, margin = 0.05, n2 = 578, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 578, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = c(n1 = 578, n2 = 578), power = 0.800541149, t.alpha = 2.4783293,
                      ncp = 3.3224971, null.ncp = 0.83062427, df = 1101.8987))

    crrRes <- power.np.wilcoxon(margin = 0.05, n2 = 578, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 578, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.199883388, n = c(n1 = 578, n2 = 578), power = 0.8, t.alpha = 2.47832927,
                      ncp = 3.32055987, null.ncp = 0.830624270, df = 1101.89869))

    crrRes <- power.np.wilcoxon(d = 0, margin = c(-0.05, 0.05), power = 0.80, alternative = "two.one.sided",
                                design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n.ratio = 1, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          distribution = "normal", method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0, n = c(n1 = 7175, n2 = 7175), power = 0.80000988, t.alpha = 1.2816415 * c(-1, 1),
                      ncp = 0, null.ncp = 2.92652104 * c(-1, 1), df = 13701.2406))
    expect_equal(crrRes, pwrss.np.2groups(mu1 = 0, sd1 = 1, margin = 0.05, power = 0.8, alternative = "equivalent", verbose = FALSE))

    crrRes <- power.np.wilcoxon(d = 0, margin = c(-0.05, 0.05), n2 = 7175, alternative = "two.one.sided",
                                design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n.ratio = 1, n2 = 7175,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          distribution = "normal", method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0, n = c(n1 = 7175, n2 = 7175), power = 0.80000988, t.alpha = 1.2816415 * c(-1, 1),
                      ncp = 0, null.ncp = 2.92652104 * c(-1, 1), df = 13701.2406))

    crrRes <- power.np.wilcoxon(d = 0, margin = c(-0.05, 0.05), power = 0.80, alternative = "two.one.sided",
                                design = "independent", method = "noether", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "z"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n.ratio = 1, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          distribution = "normal", method = "noether", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df")],
                 list(test = "z", d = 0, n = c(n1 = 7178, n2 = 7178), power = 0.80004138, z.alpha = 1.28166946 * c(-1, 1),
                      mean = 0, sd = 1, null.mean = 2.926523083 * c(-1, 1), null.sd = 1, df = Inf))

    crrRes <- power.np.wilcoxon(d = 0, margin = c(-0.05, 0.05), n2 = 7178, alternative = "two.one.sided",
                                design = "independent", method = "noether", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "z"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n.ratio = 1, n2 = 7178,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          distribution = "normal", method = "noether", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "z.alpha", "mean", "sd", "null.mean", "null.sd", "df")],
                 list(test = "z", d = 0, n = c(n1 = 7178, n2 = 7178), power = 0.80004138, z.alpha = 1.28166946 * c(-1, 1),
                      mean = 0, sd = 1, null.mean = 2.926523083 * c(-1, 1), null.sd = 1, df = Inf))

    crrRes <- power.np.wilcoxon(d = 0.05, margin = c(-0.10, -0.05), power = 0.80, alternative = "two.one.sided",
                                design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.10, -0.05), req.sign = "+", n.ratio = 1, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          distribution = "normal", method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.05, n = c(n1 = 1644, n2 = 1644), power = 0.80001254, t.alpha = c(-4.7655851, 0.5592075),
                      ncp = 1.400849902, null.ncp = c(-2.801699804, -1.400849902), df = 3137.8087))

    crrRes <- power.np.wilcoxon(d = 0.05, margin = c(-0.10, -0.05), n2 = 1644, alternative = "two.one.sided",
                                design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.10, -0.05), req.sign = "+", n.ratio = 1, n2 = 1644,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          distribution = "normal", method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.05, n = c(n1 = 1644, n2 = 1644), power = 0.80001254, t.alpha = c(-4.7655851, 0.5592075),
                      ncp = 1.400849902, null.ncp = c(-2.801699804, -1.400849902), df = 3137.8087))

    crrRes <- suppressWarnings(power.np.wilcoxon(margin = c(-0.05, 0.05), req.sign = "0", n2 = 7175, power = 0.80,
                                                 alternative = "two.one.sided", design = "independent", verbose = 0))
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "0", n.ratio = 1, n2 = 7175,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          distribution = "normal", method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0, n = c(n1 = 7175, n2 = 7175), power = 0.80000988, t.alpha = 1.2816415 * c(-1, 1),
                      ncp = 0, null.ncp = 2.92652104 * c(-1, 1), df = 13701.2406))

    crrRes <- power.np.wilcoxon(d = -0.25, power = 0.80, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.25, n = 134, power = 0.80138478, t.alpha = 1.97882545 * c(-1, 1),
                      ncp = -2.8279915, null.ncp = 0, df = 126.9605742))

    crrRes <- power.np.wilcoxon(d = -0.25, n2 = 134, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 134, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.25, n = 134, power = 0.80138478, t.alpha = 1.97882545 * c(-1, 1),
                      ncp = -2.8279915, null.ncp = 0, df = 126.9605742))

    crrRes <- power.np.wilcoxon(n2 = 134, power = 0.8, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 134, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.249558443, n = 134, power = 0.8, t.alpha = 1.97882545 * c(-1, 1),
                      ncp = 2.82299662, null.ncp = 0, df = 126.9605742))

    crrRes <- power.np.wilcoxon(d = -0.20, power = 0.80, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.20, n = 208, power = 0.800946514, t.alpha = 1.97204047 * c(-1, 1),
                      ncp = -2.818690256, null.ncp = 0, df = 197.62537))

    crrRes <- power.np.wilcoxon(d = -0.20, n2 = 208, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 208, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.20, n = 208, power = 0.800946514, t.alpha = 1.97204047 * c(-1, 1),
                      ncp = -2.818690256, null.ncp = 0, df = 197.62537))

    crrRes <- power.np.wilcoxon(n2 = 208, power = 0.80, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 208, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.199758584, n = 208, power = 0.8, t.alpha = 1.97204047 * c(-1, 1),
                      ncp = 2.81528788, null.ncp = 0, df = 197.62537))

    crrRes <- power.np.wilcoxon(d = -0.25, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.25, n = 106, power = 0.80326579, t.alpha = -1.66019985,
                      ncp = -2.51523537, null.ncp = 0, df = 100.2225438))

    crrRes <- power.np.wilcoxon(d = -0.25, n2 = 106, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 106, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.25, n = 106, power = 0.80326579, t.alpha = -1.66019985,
                      ncp = -2.51523537, null.ncp = 0, df = 100.2225438))

    crrRes <- power.np.wilcoxon(n2 = 106, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 106, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.2488268, n = 106, power = 0.8, t.alpha = 1.66019985,
                      ncp = 2.503431843, null.ncp = 0, df = 100.2225438))

    crrRes <- power.np.wilcoxon(d = -0.20, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.20, n = 164, power = 0.80153238, t.alpha = -1.654704869,
                      ncp = -2.5028661, null.ncp = 0, df = 155.608464))

    crrRes <- power.np.wilcoxon(d = -0.20, n2 = 164, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 164, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.20, n = 164, power = 0.80153238, t.alpha = -1.654704869,
                      ncp = -2.5028661, null.ncp = 0, df = 155.608464))

    crrRes <- power.np.wilcoxon(d = -0.10, margin = 0.05, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.10, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.10, n = 289, power = 0.801088074, t.alpha = -0.81596483,
                      ncp = -1.66124854, null.ncp = 0.83062427, df = 274.97467))

    crrRes <- power.np.wilcoxon(d = -0.10, margin = 0.05, n2 = 289, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.10, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 289, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.10, n = 289, power = 0.801088074, t.alpha = -0.81596483,
                      ncp = -1.66124854, null.ncp = 0.83062427, df = 274.97467))

    crrRes <- power.np.wilcoxon(d = 0.20, margin = 0.05, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = 291, power = 0.80004359, t.alpha = 2.4897522,
                      ncp = 3.33397379, null.ncp = 0.833493447, df = 276.88453))

    crrRes <- power.np.wilcoxon(d = 0.20, margin = 0.05, n2 = 291, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 291, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = 291, power = 0.80004359, t.alpha = 2.4897522,
                      ncp = 3.33397379, null.ncp = 0.833493447, df = 276.88453))

    crrRes <- power.np.wilcoxon(margin = 0.05, n2 = 291, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 291, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.199990608, n = 291, power = 0.8, t.alpha = 2.4897522,
                      ncp = 3.33381722, null.ncp = 0.833493447, df = 276.88453))

    crrRes <- power.np.wilcoxon(d = -0.10, margin = -0.05, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.10, null.d = 0, margin = -0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.10, n = 2599, power = 0.80004255, t.alpha = -4.13902135,
                      ncp = -4.98182916, null.ncp = -2.49091458, df = 2480.86218))

    crrRes <- power.np.wilcoxon(d = -0.10, margin = -0.05, n2 = 2599, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.10, null.d = 0, margin = -0.05, req.sign = "+", n.ratio = 1, n2 = 2599, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.10, n = 2599, power = 0.80004255, t.alpha = -4.13902135,
                      ncp = -4.98182916, null.ncp = -2.49091458, df = 2480.86218))

    crrRes <- power.np.wilcoxon(margin = 0.05, n2 = 2599, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05, req.sign = "+", n.ratio = 1, n2 = 2599, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.099996944, n = 2599, power = 0.8000, t.alpha = 4.13902135,
                      ncp = 4.981676920, null.ncp = 2.49091458, df = 2480.86218))

    crrRes <- power.np.wilcoxon(d = 0.20, margin = -0.05, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = -0.05, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = 105, power = 0.80239028, t.alpha = 1.1525527,
                      ncp = 2.00267435, null.ncp = -0.500668588, df = 99.2676141))

    crrRes <- power.np.wilcoxon(d = 0.20, margin = -0.05, n2 = 105, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = -0.05, req.sign = "+", n.ratio = 1, n2 = 105, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.20, n = 105, power = 0.80239028, t.alpha = 1.1525527,
                      ncp = 2.00267435, null.ncp = -0.500668588, df = 99.2676141))

    crrRes <- power.np.wilcoxon(margin = -0.05, n2 = 105, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = -0.05, req.sign = "+", n.ratio = 1, n2 = 105, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.199141404, n = 105, power = 0.8, t.alpha = 1.1525527,
                      ncp = 1.994076911, null.ncp = -0.500668588, df = 99.2676141))

    crrRes <- power.np.wilcoxon(d = 0, margin = c(-0.05, 0.05), power = 0.80, alternative = "two.one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n.ratio = 1, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "paired",
                          distribution = "normal", method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0, n = 3589, power = 0.800132137, t.alpha = 1.28217542 * c(-1, 1),
                      ncp = 0, null.ncp = 2.9271328 * c(-1, 1), df = 3426.24254))

    crrRes <- power.np.wilcoxon(d = 0, margin = c(-0.05, 0.05), n2 = 3589, alternative = "two.one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n.ratio = 1, n2 = 3589,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "paired",
                          distribution = "normal", method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0, n = 3589, power = 0.800132137, t.alpha = 1.28217542 * c(-1, 1),
                      ncp = 0, null.ncp = 2.9271328 * c(-1, 1), df = 3426.24254))

    crrRes <- power.np.wilcoxon(d = -0.25, power = 0.80, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.25, n = 134, power = 0.80138478, t.alpha = 1.97882545 * c(-1, 1),
                      ncp = -2.8279915, null.ncp = 0, df = 126.9605742))

    crrRes <- power.np.wilcoxon(d = -0.25, n2 = 134, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 134, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.25, n = 134, power = 0.80138478, t.alpha = 1.97882545 * c(-1, 1),
                      ncp = -2.8279915, null.ncp = 0, df = 126.9605742))

    crrRes <- power.np.wilcoxon(n2 = 134, power = 0.80, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 134, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.249558443, n = 134, power = 0.8, t.alpha = 1.97882545 * c(-1, 1),
                      ncp = 2.82299662, null.ncp = 0, df = 126.9605742))

    crrRes <- power.np.wilcoxon(d = -0.20, power = 0.80, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.20, n = 208, power = 0.800946514, t.alpha = 1.97204047 * c(-1, 1),
                      ncp = -2.818690256, null.ncp = 0, df = 197.62537))

    crrRes <- power.np.wilcoxon(d = -0.25, power = 0.80, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.25, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.25, n = 106, power = 0.80326579, t.alpha = -1.66019985,
                      ncp = -2.51523537, null.ncp = 0, df = 100.2225438))

    crrRes <- power.np.wilcoxon(d = -0.20, power = 0.80, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = NULL, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.20, n = 164, power = 0.80153238, t.alpha = -1.654704869,
                      ncp = -2.5028661, null.ncp = 0, df = 155.608464))

    crrRes <- power.np.wilcoxon(d = -0.20, n2 = 164, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 164, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = -0.20, n = 164, power = 0.80153238, t.alpha = -1.654704869,
                      ncp = -2.5028661, null.ncp = 0, df = 155.608464))

    crrRes <- power.np.wilcoxon(n2 = 164, power = 0.80, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 164, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.199559679, n = 164, power = 0.8, t.alpha = 1.654704869,
                      ncp = 2.497355746, null.ncp = 0, df = 155.608464))

    crrRes <- power.np.wilcoxon(d = 0.1, n2 = 649, design = "one.sample", alternative = "one.sided", distribution = "normal",
                                verbose = 0) # example 22.3 from the GPower manual
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.1, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 649, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", distribution = "normal",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.1, n = 649, power = 0.800078016, t.alpha = 1.64732, ncp = 2.489476548, null.ncp = 0,
                      df = 618.74935))
    # uses a different calculation (Lehmann without continuity correction in GPower, Guenther in pwrss)
    # results are (roughly) identical: t / z crit ~ 1.65, power ~ 0.8001

    crrRes <- power.np.wilcoxon(d = 0.8, n2 = 11, design = "one.sample", alternative = "one.sided",
                                distribution = "laplace", verbose = 0) # example 22.3 from the GPower manual
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.8, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 11, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", distribution = "laplace",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.8, n = 11, power = 0.92759865, t.alpha = 1.7493444, ncp = 3.2496154, null.ncp = 0, df = 15.5))
    # uses a different calculation (Lehmann without continuity correction in GPower, Guenther in pwrss)
    # results are quite off: t / z crit = 1.6448536 (GPower) vs. 1.7493444 (pwrss), power = 0.8308341 (GPower) vs. 0.92759865 (pwrss)

    crrRes <- power.np.wilcoxon(d = 1.13842, n2 = 5, design = "paired", alternative = "two.sided", distribution = "laplace",
                                verbose = 0) # example 23.3 from the GPower manual
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 1.13842, null.d = 0, margin = 0, req.sign = "+", n.ratio = 1, n2 = 5, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "paired", distribution = "laplace",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 1.13842, n = 5, power = 0.75156987, t.alpha = 2.40200201 * c(-1, 1), ncp = 3.1176916,
                      null.ncp = 0, df = 6.5))
    # uses a different calculation (Lehmann without continuity correction in GPower, Guenther in pwrss)
    # results are quite off: t / z crit = 1.9599640 (GPower) vs. 2.40200201 (pwrss), power = 0.8534958 (GPower) vs. 0.75156987 (pwrss)

    crrRes <- power.np.wilcoxon(d = 0.375, n.ratio = 2, n2 = 67, design = "independent", alternative = "two.sided",
                                distribution = "laplace", verbose = 0) # example 24.3.1 from the GPower manual
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.375, null.d = 0, margin = 0, req.sign = "+", n.ratio = 2, n2 = 67, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "independent", distribution = "laplace",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.375, n = c(n1 = 134, n2 = 67), power = 0.86426579, t.alpha = 1.96791632 * c(-1, 1),
                      ncp = 3.0695073, null.ncp = 0, df = 299.5))
    # uses a different calculation (Lehmann without continuity correction in GPower, Guenther in pwrss)
    # results are a bit off: t / z crit = 1.9599640 (GPower) vs. 1.96791632 (pwrss), power = 0.8472103 (GPower) vs. 0.86426579 (pwrss)

    crrRes <- power.np.wilcoxon(d = 0.2911838, n.ratio = 2, n2 = 67, design = "independent", alternative = "one.sided",
                                distribution = "logistic", verbose = 0) # example 24.3.2 from the GPower manual, using d instead of p1
    expect_s3_class(crrRes, c("pwrss", "np", "wilcoxon", "t"))
    expect_named(crrRes, c("parms", "test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2911838, null.d = 0, margin = 0, req.sign = "+", n.ratio = 2, n2 = 67, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "independent", distribution = "logistic",
                          method = "guenther", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "d", "n", "power", "t.alpha", "ncp", "null.ncp", "df")],
                 list(test = "t", d = 0.2911838, n = c(n1 = 134, n2 = 67), power = 0.65052962, t.alpha = 1.65185978,
                      ncp = 2.03792219, null.ncp = 0, df = 218.421165))
    # uses a different calculation (Lehmann without continuity correction in GPower, Guenther in pwrss)
    # results are slightly off: t / z crit = 1.6448536 (GPower) vs. 1.65185978 (pwrss), power = 0.6461276 (GPower) vs. 0.65052962 (pwrss)

    expect_equal(fmt_test_wilcoxon("independent"),
                 "Wilcoxon Rank-Sum Test (Independent Samples) \n(Wilcoxon-Mann-Whitney or Mann-Whitney U Test)")
    expect_equal(fmt_test_wilcoxon("paired"),
                 "Wilcoxon Signed-Rank Test (Paired Samples)")
    expect_equal(fmt_test_wilcoxon("one.sample"),
                 "Wilcoxon Signed-Rank Test (One Sample)")

    expect_error(power.np.wilcoxon(d = 0.10, margin = -11, power = 0.80, alternative = "one.sided", verbose = 0),
                 "Possibly incorrect value for `margin` \\(should be within -10 ... 10\\).")
    expect_error(power.np.wilcoxon(d = -0.20, power = 0.80, alternative = "two.sided", design = "paired", method = "noether", verbose = 0),
                 "Specify `method` = \"guenther\" to request Wilcoxon signed-rank test for matched pairs.")
    expect_error(power.np.wilcoxon(d = 0, alpha = 1e-4, power = 0.99, alternative = "two.sided"),
                 "Design is not feasible.")
    expect_error(power.np.wilcoxon(n2 = 4, alpha = 1e-4, power = 0.99, alternative = "two.sided"),
                 "Design is not feasible.")
    expect_error(power.np.wilcoxon(d = 0.2, n2 = 2, alpha = 1e-4, alternative = "two.sided"),
                 "Degrees of freedom can not be smaller than 3.")
    expect_error(pwrss.np.2means(), "This function is no longer available. Please use `power.np.wilcoxon\\(\\)`.")
    expect_warning(power.np.wilcoxon(margin = c(-0.05, 0.05), req.sign = "0", n2 = 20, power = 0.80,
                                     alternative = "two.one.sided", design = "independent", verbose = 0),
                   "The target power rate cannot be achieved within the null bounds.")
})
