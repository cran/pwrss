# power.t.student (= pwrss.t.mean / pwrss.t.2means) --------------------------------------------------------------------
test_that("power.t.student / pwrss.t.mean / pwrss.t.2means work", {
    crrRes <- power.t.student(d = 0.20, power = 0.80, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2,  null.d = 0, margin = 0,  req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 786, ncp = 2.80713377, null.ncp = 0, t.alpha = 1.96298672 * c(-1, 1),
                      d = 0.20, power = 0.800593128, n = c(n1 = 394, n2 = 394), n.total = 788))

    crrRes <- power.t.student(d = 0.20, n2 = 394, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2,  null.d = 0, margin = 0, req.sign = "+", n2 = 394, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 786, ncp = 2.80713377, null.ncp = 0, t.alpha = 1.96298672 * c(-1, 1),
                      d = 0.20, power = 0.800593128, n = c(n1 = 394, n2 = 394), n.total = 788))

    crrRes <- power.t.student(n2 = 394, power = 0.80, alternative = "two.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,  null.d = 0, margin = 0, req.sign = "+", n2 = 394, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 786, ncp = 2.80501064, null.ncp = 0, t.alpha = 1.96298672 * c(-1, 1),
                      d = 0.199848734, power = 0.8, n = c(n1 = 394, n2 = 394), n.total = 788))

    crrRes <- power.t.student(d = 0.20, power = 0.80, alternative = "two.sided", design = "independent",
                              claim.basis = "smd.ci", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2,  null.d = 0, margin = 0, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", claim.basis = "smd.ci",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 790, ncp = 2.80722243, null.ncp = 0, t.alpha = 1.96297139 * c(-1, 1),
                      d = 0.20, power = 0.80062273, n = c(n1 = 396, n2 = 396), n.total = 792))

    crrRes <- power.t.student(d = 0.20, n2 = 396, alternative = "two.sided", design = "independent",
                              claim.basis = "smd.ci", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2,   null.d = 0, margin = 0, req.sign = "+", n2 = 396, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "independent", claim.basis = "smd.ci",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 790, ncp = 2.80722243, null.ncp = 0, t.alpha = 1.96297139 * c(-1, 1),
                      d = 0.20, power = 0.80062273, n = c(n1 = 396, n2 = 396), n.total = 792))

    crrRes <- power.t.student(n2 = 396, power = 0.80, alternative = "two.sided", design = "independent",
                              claim.basis = "smd.ci", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,   null.d = 0, margin = 0, req.sign = "+", n2 = 396, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "independent", claim.basis = "smd.ci",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 790, ncp = 2.80499325, null.ncp = 0, t.alpha = 1.96297139 * c(-1, 1),
                      d = 0.199840388, power = 0.8, n = c(n1 = 396, n2 = 396), n.total = 792))

    crrRes <- power.t.student(d = 0.20, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2,   null.d = 0, margin = 0, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 618, ncp = 2.48997992, null.ncp = 0, t.alpha = 1.647323,
                      d = 0.20, power = 0.8002178, n = c(n1 = 310, n2 = 310), n.total = 620))

    crrRes <- power.t.student(d = 0.20, n2 = 310, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2,   null.d = 0, margin = 0, req.sign = "+", n2 = 310, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 618, ncp = 2.48997992, null.ncp = 0, t.alpha = 1.647323,
                      d = 0.20, power = 0.8002178, n = c(n1 = 310, n2 = 310), n.total = 620))

    crrRes <- power.t.student(n2 = 310, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,   null.d = 0, margin = 0, req.sign = "+", n2 = 310, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 618, ncp = 2.48920085, null.ncp = 0, t.alpha = 1.647323,
                      d = 0.199937424, power = 0.8, n = c(n1 = 310, n2 = 310), n.total = 620))

    crrRes <- power.t.student(d = 0.20, margin = -0.05, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2, null.d = 0, margin = -0.05, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 396, ncp = 1.994993734, null.ncp = -0.49874843, t.alpha = 1.14819672,
                      d = 0.20, power = 0.80145163, n = c(n1 = 199, n2 = 199), n.total = 398))

    crrRes <- power.t.student(d = 0.20, margin = -0.05, n2 = 199, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2, null.d = 0, margin = -0.05, req.sign = "+", n2 = 199, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 396, ncp = 1.994993734, null.ncp = -0.49874843, t.alpha = 1.14819672,
                      d = 0.20, power = 0.80145163, n = c(n1 = 199, n2 = 199), n.total = 398))

    crrRes <- power.t.student(margin = -0.05, n2 = 199, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = -0.05, req.sign = "+", n2 = 199, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 396, ncp = 1.989792938, null.ncp = -0.49874843, t.alpha = 1.14819672,
                      d = 0.199478615, power = 0.8, n = c(n1 = 199, n2 = 199), n.total = 398))

    crrRes <- power.t.student(d = 0.20, margin = 0.05, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2, null.d = 0, margin = 0.05, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 1102, ncp = 3.32264955, null.ncp = 0.830662386, t.alpha = 2.4783672,
                      d = 0.20, power = 0.800573124, n = c(n1 = 552, n2 = 552), n.total = 1104))

    crrRes <- power.t.student(d = 0.20, margin = 0.05, n2 = 552, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.2, null.d = 0, margin = 0.05, req.sign = "+", n2 = 552, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 1102, ncp = 3.32264955, null.ncp = 0.830662386, t.alpha = 2.4783672,
                      d = 0.20, power = 0.800573124, n = c(n1 = 552, n2 = 552), n.total = 1104))

    crrRes <- power.t.student(margin = 0.05, n2 = 552, power = 0.80, alternative = "one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05, req.sign = "+", n2 = 552, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 1102, ncp = 3.3205978, null.ncp = 0.830662386, t.alpha = 2.4783672,
                      d = 0.1998765, power = 0.8, n = c(n1 = 552, n2 = 552), n.total = 1104))

    crrRes <- power.t.student(margin = -0.05, req.sign = "-", n2 = 552, power = 0.80, alternative = "one.sided",
                              design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = -0.05, req.sign = "-", n2 = 552, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "independent", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 1102, ncp = -3.3205978, null.ncp = -0.830662386, t.alpha = -2.4783672,
                      d = -0.1998765, power = 0.8, n = c(n1 = 552, n2 = 552), n.total = 1104))

    crrRes <- power.t.student(d = 0, margin = c(-0.05, 0.05), power = 0.80, alternative = "two.one.sided",
                              design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n2 = NULL, n.ratio = 1,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 13702, ncp = 0, null.ncp = 2.92660213 * c(-1, 1), t.alpha = 1.2817226 * c(-1, 1),
                      d = 0, power = 0.800038332, n = c(n1 = 6852, n2 = 6852), n.total = 13704))

    crrRes <- power.t.student(d = 0, margin = c(-0.05, 0.05), n2 = 6852, alternative = "two.one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n2 = 6852, n.ratio = 1,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 13702, ncp = 0, null.ncp = 2.92660213 * c(-1, 1), t.alpha = 1.2817226 * c(-1, 1),
                      d = 0, power = 0.800038332, n = c(n1 = 6852, n2 = 6852), n.total = 13704))

    crrRes <- power.t.student(d = 0.05, margin = c(-0.05, -0.10), power = 0.80, alternative = "two.one.sided",
                              design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.05, -0.10), req.sign = "+", n2 = NULL, n.ratio = 1,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 3138, ncp = 1.400892573, null.ncp = c(-1.400892573, -2.801785145),
                      t.alpha = c(-4.76567037, 0.55916477), d = 0.05, power = 0.800036434, n = c(n1 = 1570, n2 = 1570),
                      n.total = 3140))

    crrRes <- power.t.student(d = 0.05, margin = c(-0.05, -0.10), n2 = 1570, alternative = "two.one.sided", design = "independent", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.05, -0.10), req.sign = "+", n2 = 1570, n.ratio = 1,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 3138, ncp = 1.400892573, null.ncp = c(-1.400892573, -2.801785145),
                      t.alpha = c(-4.76567037, 0.55916477), d = 0.05, power = 0.800036434, n = c(n1 = 1570, n2 = 1570),
                      n.total = 3140))

    crrRes <- suppressWarnings(power.t.student(margin = c(-0.05, 0.05), req.sign = "0", n2 = 6852, power = 0.80,
                                               alternative = "two.one.sided", design = "independent", verbose = 0))
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "0", n2 = 6852, n.ratio = 1,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "independent",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 13702, ncp = 0, null.ncp = 2.92660213 * c(-1, 1), t.alpha = 1.2817226 * c(-1, 1),
                      d = 0, power = 0.800038332, n = c(n1 = 6852, n2 = 6852), n.total = 13704))

    crrRes <- power.t.student(d = -0.20, power = 0.80, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,   null.d = 0, margin = 0, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 198, ncp = -2.8213472, null.ncp = 0, t.alpha = 1.97201748 * c(-1, 1),
                      d = -0.20, power = 0.801691024, n = 199, n.total = 199))
    expect_equal(crrRes, pwrss.t.2means(mu1 = -0.20, sd1 = 1, power = 0.80, alternative = "not equal", paired = TRUE, verbose = FALSE))

    crrRes <- power.t.student(d = -0.20, n2 = 199, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,   null.d = 0, margin = 0, req.sign = "+", n2 = 199, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 198, ncp = -2.8213472, null.ncp = 0, t.alpha = 1.97201748 * c(-1, 1),
                      d = -0.20, power = 0.801691024, n = 199, n.total = 199))

    crrRes <- power.t.student(n2 = 199, power = 0.80, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,   null.d = 0, margin = 0, req.sign = "+", n2 = 199, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 198, ncp = 2.81526177, null.ncp = 0, t.alpha = 1.97201748 * c(-1, 1),
                      d = 0.199568616, power = 0.8, n = 199, n.total = 199))

    crrRes <- power.t.student(d = -0.20, power = 0.80, alternative = "two.sided", claim.basis = "smd.ci", design = "paired",
                              verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,   null.d = 0, margin = 0, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", claim.basis = "smd.ci",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 202, ncp = -2.82148619, null.ncp = 0, t.alpha = 1.97177738 * c(-1, 1),
                      d = -0.20, power = 0.801805282, n = 203, n.total = 203))

    crrRes <- power.t.student(d = -0.20, n2 = 203, alternative = "two.sided", claim.basis = "smd.ci", design = "paired",
                              verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,   null.d = 0, margin = 0, req.sign = "+", n2 = 203, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "paired", claim.basis = "smd.ci",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 202, ncp = -2.82148619, null.ncp = 0, t.alpha = 1.97177738 * c(-1, 1),
                      d = -0.20, power = 0.801805282, n = 203, n.total = 203))

    crrRes <- power.t.student(n2 = 203, power = 0.80, alternative = "two.sided", claim.basis = "smd.ci", design = "paired",
                              verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,   null.d = 0, margin = 0, req.sign = "+", n2 = 203, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "paired", claim.basis = "smd.ci",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 202, ncp = 2.814989096, null.ncp = 0, t.alpha = 1.97177738 * c(-1, 1),
                      d = 0.19953028, power = 0.8, n = 203, n.total = 203))

    crrRes <- power.t.student(d = -0.20, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,   null.d = 0, margin = 0, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 155, ncp = -2.4979992, null.ncp = 0, t.alpha = -1.654743774,
                      d = -0.20, power = 0.80016732, n = 156, n.total = 156))

    crrRes <- power.t.student(d = -0.20, n2 = 156, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,   null.d = 0, margin = 0, req.sign = "+", n2 = 156, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 155, ncp = -2.4979992, null.ncp = 0, t.alpha = -1.654743774,
                      d = -0.20, power = 0.80016732, n = 156, n.total = 156))

    crrRes <- power.t.student(n2 = 156, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,   null.d = 0, margin = 0, req.sign = "+", n2 = 156, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 155, ncp = 2.497398749, null.ncp = 0, t.alpha = 1.654743774,
                      d = 0.199951925, power = 0.8, n = 156, n.total = 156))

    crrRes <- power.t.student(d = 0.20, margin = 0.05, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 277, ncp = 3.3346664, null.ncp = 0.8336666, t.alpha = 2.489922,
                      d = 0.20, power = 0.80018923, n = 278, n.total = 278))

    crrRes <- power.t.student(d = 0.20, margin = 0.05, n2 = 278, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", n2 = 278, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 277, ncp = 3.3346664, null.ncp = 0.8336666, t.alpha = 2.489922,
                      d = 0.20, power = 0.80018923, n = 278, n.total = 278))

    crrRes <- power.t.student(margin = 0.05, n2 = 278, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05, req.sign = "+", n2 = 278, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 277, ncp = 3.3339865, null.ncp = 0.8336666, t.alpha = 2.489922,
                      d = 0.199959223, power = 0.8, n = 278, n.total = 278))

    crrRes <- power.t.student(d = 0.20, margin = -0.05, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = -0.05, req.sign = "+", n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 99, ncp = 2, null.ncp = -0.5, t.alpha = 1.1532522,
                      d = 0.20, power = 0.801453942, n = 100, n.total = 100))

    crrRes <- power.t.student(d = 0.20, margin = -0.05, n2 = 100, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = -0.05, req.sign = "+", n2 = 100, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 99, ncp = 2, null.ncp = -0.5, t.alpha = 1.1532522,
                      d = 0.20, power = 0.801453942, n = 100, n.total = 100))

    crrRes <- power.t.student(margin = -0.05, n2 = 100, power = 0.80, alternative = "one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = -0.05, req.sign = "+", n2 = 100, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 99, ncp = 1.99477777, null.ncp = -0.5, t.alpha = 1.1532522,
                      d = 0.199477777, power = 0.8, n = 100, n.total = 100))

    crrRes <- power.t.student(d = 0, margin = c(-0.05, 0.05), power = 0.80, alternative = "two.one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n2 = NULL, n.ratio = 1,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "paired",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 3426, ncp = 0, null.ncp = 2.92702921 * c(-1, 1), t.alpha = 1.282071859 * c(-1, 1),
                      d = 0, power = 0.800095812, n = 3427, n.total = 3427))

    crrRes <- power.t.student(d = 0, margin = c(-0.05, 0.05), n2 = 3427, alternative = "two.one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", n2 = 3427, n.ratio = 1,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "paired",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 3426, ncp = 0, null.ncp = 2.92702921 * c(-1, 1), t.alpha = 1.282071859 * c(-1, 1),
                      d = 0, power = 0.800095812, n = 3427, n.total = 3427))

    crrRes <- power.t.student(d = 0.05, margin = c(-0.05, -0.10), power = 0.80, alternative = "two.one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.05, -0.10), req.sign = "+", n2 = NULL, n.ratio = 1,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "paired",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 785, ncp = 1.401784577, null.ncp = c(-1.401784577, -2.803569154),
                      t.alpha = c(-4.7792355, 0.5585519), d = 0.05, power = 0.800477147, n = 786, n.total = 786))

    crrRes <- power.t.student(d = 0.05, margin = c(-0.05, -0.10), n2 = 786, alternative = "two.one.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.05, -0.10), req.sign = "+", n2 = 786, n.ratio = 1,
                          power = NULL, alpha = 0.05, alternative = "two.one.sided", design = "paired",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 785, ncp = 1.401784577, null.ncp = c(-1.401784577, -2.803569154),
                      t.alpha = c(-4.7792355, 0.5585519), d = 0.05, power = 0.800477147, n = 786, n.total = 786))

    crrRes <- power.t.student(d = -0.20, power = 0.80, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,  null.d = 0, margin = 0, req.sign = "+",  n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "one.sample", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 198, ncp = -2.8213472, null.ncp = 0, t.alpha = 1.97201748 * c(-1, 1),
                      d = -0.20, power = 0.801691024, n = 199, n.total = 199))
    expect_equal(crrRes, pwrss.t.mean(mu = -0.20, sd = 1, power = 0.8, alternative = "not equal", verbose = FALSE))

    crrRes <- power.t.student(d = -0.20, n2 = 199, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,  null.d = 0, margin = 0, req.sign = "+",  n2 = 199, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "one.sample", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 198, ncp = -2.8213472, null.ncp = 0, t.alpha = 1.97201748 * c(-1, 1),
                      d = -0.20, power = 0.801691024, n = 199, n.total = 199))
    expect_equal(crrRes, pwrss.t.mean(mu = -0.20, sd = 1, n = 199, alternative = "not equal", verbose = FALSE))

    crrRes <- power.t.student(n2 = 199, power = 0.80, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,  null.d = 0, margin = 0, req.sign = "+",  n2 = 199, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "two.sided", design = "one.sample", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 198, ncp = 2.81526177, null.ncp = 0, t.alpha = 1.97201748 * c(-1, 1),
                      d = 0.199568616, power = 0.8, n = 199, n.total = 199))

    crrRes <- power.t.student(d = -0.20, power = 0.80, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,  null.d = 0, margin = 0, req.sign = "+",  n2 = NULL, n.ratio = 1, power = 0.80,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 155, ncp = -2.4979992, null.ncp = 0, t.alpha = -1.654743774,
                      d = -0.20, power = 0.80016732, n = 156, n.total = 156))
    expect_equal(crrRes, pwrss.t.mean(mu = -0.20, sd = 1, power = 0.8, alternative = "less", verbose = FALSE))

    crrRes <- power.t.student(d = -0.20, n2 = 156, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = -0.20,  null.d = 0, margin = 0, req.sign = "+",  n2 = 156, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 155, ncp = -2.4979992, null.ncp = 0, t.alpha = -1.654743774,
                      d = -0.20, power = 0.80016732, n = 156, n.total = 156))
    expect_equal(crrRes, pwrss.t.mean(mu = -0.20, sd = 1, n = 156, alternative = "less", verbose = FALSE))

    crrRes <- power.t.student(n2 = 156, power = 0.80, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,  null.d = 0, margin = 0, req.sign = "+",  n2 = 156, n.ratio = 1, power = 0.80, alpha = 0.05,
                      alternative = "one.sided", design = "one.sample", claim.basis = "md.pval", ceil.n = TRUE,
                      verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 155, ncp = 2.497398749, null.ncp = 0, t.alpha = 1.654743774,
                      d = 0.199951925, power = 0.8, n = 156, n.total = 156))

    crrRes <- power.t.student(d = 0.20, margin = c(-0.1, 0.1), power = 0.8, alternative = "two.one.sided",
                              design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = c(-0.1, 0.1), req.sign = "+", n2 = NULL, n.ratio = 1,
                          power = 0.80, alpha = 0.05, alternative = "two.one.sided", design = "one.sample",
                          claim.basis = "md.pval", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 796, ncp = 5.6462377, null.ncp = 2.82311884 * c(-1, 1), t.alpha = 4.79868895 * c(-1, 1),
                      d = 0.20, power = 0.800381536, n = 797, n.total = 797))
    expect_equal(crrRes, pwrss.t.mean(mu = 0.20, margin = 0.1, sd = 1, power = 0.8, alternative = "equivalent", verbose = FALSE))

    # example 19.3 from the GPower manual
    crrRes <- power.t.student(d = 0.421637021, n2 = 50, alternative = "two.sided", design = "paired", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.421637021,  null.d = 0, margin = 0, req.sign = "+",  n2 = 50, n.ratio = 1, power = NULL,
                          alpha = 0.05, alternative = "two.sided", design = "paired", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 49, ncp = 2.981424, null.ncp = 0, t.alpha = 2.00957524 * c(-1, 1),
                      d = 0.421637021, power = 0.8321145, n = 50, n.total = 50))
    # the results are identical: ncp ~ 2.981424, t.crit ~ 2.009575, power ~ 0.832114, n = 50

    # example 20.3 from the GPower manual
    crrRes <- power.t.student(d = 0.625, power = 0.95, alternative = "one.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.625,  null.d = 0, margin = 0, req.sign = "+",  n2 = NULL, n.ratio = 1, power = 0.95,
                          alpha = 0.05, alternative = "one.sided", design = "one.sample", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 29, ncp = 3.423266, null.ncp = 0, t.alpha = 1.699127027, d = 0.625,
                      power = 0.95514436, n = 30, n.total = 30))
    # the results are identical: ncp ~ 3.423266, t.crit ~ 1.699127, power ~ 0.955144, n = 30

    # example 20.3 from the GPower manual
    crrRes <- power.t.student(d = 0.1, power = 0.90, alpha = 0.01, alternative = "two.sided", design = "one.sample", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "student"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.1,  null.d = 0, margin = 0, req.sign = "+",  n2 = NULL, n.ratio = 1, power = 0.90,
                          alpha = 0.01, alternative = "two.sided", design = "one.sample", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 1491, ncp = 3.86264158, null.ncp = 0, t.alpha = 2.57913076 * c(-1, 1),
                      d = 0.10, power = 0.90016873, n = 1492, n.total = 1492))
    # the results are identical: ncp ~ 3.862642, t.crit ~ 2.579131, power ~ 0.900169, n = 1492

    expect_error(power.t.student(d = 1e-6, power = 0.99, alpha = 1e-6, alternative = "two.sided", design = "independent"),
                 "Design is not feasible.")
    expect_error(power.t.student(power = 0.99, n2 = 5, alpha = 1e-6, alternative = "two.sided", design = "independent"),
                 "Design is not feasible.")
    expect_error(power.t.student(d = 0.2, n2 = 2, alpha = 1e-6, alternative = "two.sided", design = "independent"),
                 "Degrees of freedom can not be smaller than 3.")
    expect_warning(power.t.student(margin = c(-0.05, 0.05), req.sign = "0", n2 = 20, power = 0.80,
                                   alternative = "two.one.sided", design = "independent", verbose = 0),
                   "The target power rate cannot be achieved within the null bounds.")
    expect_warning(pwrss.t.mean(mu = 0.20, margin = 0.1, sd = 1, power = 0.8, alternative = "not equal", verbose = FALSE),
                   "Margin is forced to be 0 for the 'two.sided' test.")
})

# power.t.welch --------------------------------------------------------------------------------------------------------
test_that("power.t.welch works", {
    crrRes <- power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 2, power = 0.80, alternative = "two.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 629.77715, ncp = 2.81109027, null.ncp = 0, t.alpha = 1.96373795 * c(-1, 1),
                      d = 0.20, power = 0.801457968, n = c(n1 = 474, n2 = 237), n.total = 711))

    crrRes <- power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 2, n2 = 237, alternative = "two.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = 237,
                          power = NULL, alpha = 0.05, alternative = "two.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 629.77715, ncp = 2.81109027, null.ncp = 0, t.alpha = 1.96373795 * c(-1, 1),
                      d = 0.20, power = 0.801457968, n = c(n1 = 474, n2 = 237), n.total = 711))

    crrRes <- power.t.welch(n.ratio = 2, var.ratio = 2, n2 = 237, power = 0.80, alternative = "two.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = 237,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 629.77715, ncp = 2.805863, null.ncp = 0, t.alpha = 1.96373795 * c(-1, 1),
                      d = 0.199628097, power = 0.8, n = c(n1 = 474, n2 = 237), n.total = 711))

    crrRes <- power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 2, power = 0.80, alternative = "two.sided", claim.basis = "smd.ci", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", claim.basis = "smd.ci", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 632.44382, ncp = 2.80857462, null.ncp = 0, t.alpha = 1.963722 * c(-1, 1),
                      d = 0.20, power = 0.800762167, n = c(n1 = 476, n2 = 238), n.total = 714))

    crrRes <- power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 2, n2 = 238, alternative = "two.sided", claim.basis = "smd.ci", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = 238,
                          power = NULL, alpha = 0.05, alternative = "two.sided", claim.basis = "smd.ci", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 632.44382, ncp = 2.80857462, null.ncp = 0, t.alpha = 1.963722 * c(-1, 1),
                      d = 0.20, power = 0.800762167, n = c(n1 = 476, n2 = 238), n.total = 714))

    crrRes <- power.t.welch(n.ratio = 2, var.ratio = 2, n2 = 238, power = 0.80, alternative = "two.sided", claim.basis = "smd.ci", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = 238,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", claim.basis = "smd.ci", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 632.44382, ncp = 2.805844912, null.ncp = 0, t.alpha = 1.963722 * c(-1, 1),
                      d = 0.199804448, power = 0.8, n = c(n1 = 476, n2 = 238), n.total = 714))

    crrRes <- power.t.welch(req.sign = "-", n.ratio = 2, var.ratio = 2, n2 = 238, power = 0.80, alternative = "two.sided",
                            claim.basis = "smd.ci", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,  null.d = 0, margin = 0, req.sign = "-",  var.ratio = 2, n.ratio = 2, n2 = 238,
                          power = 0.80, alpha = 0.05, alternative = "two.sided", claim.basis = "smd.ci", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 632.44382, ncp = -2.805844912, null.ncp = 0, t.alpha = 1.963722 * c(-1, 1),
                      d = -0.199804448, power = 0.8, n = c(n1 = 476, n2 = 238), n.total = 714))

    crrRes <- power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 2, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 493.776978, ncp = 2.4904277, null.ncp = 0, t.alpha = 1.64794541,
                      d = 0.20, power = 0.80015083, n = c(n1 = 372, n2 = 186), n.total = 558))

    crrRes <- power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 2, n2 = 186, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = 186,
                          power = NULL, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 493.776978, ncp = 2.4904277, null.ncp = 0, t.alpha = 1.64794541,
                      d = 0.20, power = 0.80015083, n = c(n1 = 372, n2 = 186), n.total = 558))

    crrRes <- power.t.welch(n.ratio = 2, var.ratio = 2, n2 = 186, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL,  null.d = 0, margin = 0, req.sign = "+",  var.ratio = 2, n.ratio = 2, n2 = 186,
                          power = 0.80, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 493.776978, ncp = 2.4898881, null.ncp = 0, t.alpha = 1.64794541,
                      d = 0.199956665, power = 0.8, n = c(n1 = 372, n2 = 186), n.total = 558))

    crrRes <- power.t.welch(d = 0.20, margin = -0.05, n.ratio = 2, var.ratio = 2, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = -0.05, req.sign = "+", var.ratio = 2, n.ratio = 2, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 315.109859, ncp = 1.99221018, null.ncp = -0.49805255, t.alpha = 1.149433434,
                      d = 0.20, power = 0.8003316, n = c(n1 = 238, n2 = 119), n.total = 357))

    crrRes <- power.t.welch(d = 0.20, margin = -0.05, n.ratio = 2, var.ratio = 2, n2 = 119, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = -0.05, req.sign = "+", var.ratio = 2, n.ratio = 2, n2 = 119,
                          power = NULL, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 315.109859, ncp = 1.99221018, null.ncp = -0.49805255, t.alpha = 1.149433434,
                      d = 0.20, power = 0.8003316, n = c(n1 = 238, n2 = 119), n.total = 357))

    crrRes <- power.t.welch(margin = -0.05, n.ratio = 2, var.ratio = 2, n2 = 119, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = -0.05, req.sign = "+", var.ratio = 2, n.ratio = 2, n2 = 119,
                          power = 0.80, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 315.109859, ncp = 1.991023943, null.ncp = -0.49805255, t.alpha = 1.149433434,
                      d = 0.19988091, power = 0.8, n = c(n1 = 238, n2 = 119), n.total = 357))

    crrRes <- power.t.welch(d = 0.20, margin = 0.05, n.ratio = 2, var.ratio = 2, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", var.ratio = 2, n.ratio = 2, n2 = NULL,
                          power = 0.80, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 880.444, ncp = 3.32198134, null.ncp = 0.830495336, t.alpha = 2.478918093,
                      d = 0.20, power = 0.80018977, n = c(n1 = 662, n2 = 331), n.total = 993))

    crrRes <- power.t.welch(d = 0.20, margin = 0.05, n.ratio = 2, var.ratio = 2, n2 = 331, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.20, null.d = 0, margin = 0.05, req.sign = "+", var.ratio = 2, n.ratio = 2, n2 = 331,
                          power = NULL, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 880.444, ncp = 3.32198134, null.ncp = 0.830495336, t.alpha = 2.478918093,
                      d = 0.20, power = 0.80018977, n = c(n1 = 662, n2 = 331), n.total = 993))

    crrRes <- power.t.welch(margin = 0.05, n.ratio = 2, var.ratio = 2, n2 = 331, power = 0.80, alternative = "one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05, req.sign = "+", var.ratio = 2, n.ratio = 2, n2 = 331,
                          power = 0.80, alpha = 0.05, alternative = "one.sided", claim.basis = "md.pval", ceil.n = TRUE,
                          verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 880.444, ncp = 3.3213021, null.ncp = 0.830495336, t.alpha = 2.478918093,
                      d = 0.19995911, power = 0.8, n = c(n1 = 662, n2 = 331), n.total = 993))

    crrRes <- power.t.welch(d = 0, margin = c(-0.05, 0.05), n.ratio = 2, var.ratio = 2, power = 0.80, alternative = "two.one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", var.ratio = 2, n.ratio = 2,
                          n2 = NULL, power = 0.80, alpha = 0.05, alternative = "two.one.sided", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 10960.44441, ncp = 0, null.ncp = 2.926554671 * c(-1, 1), t.alpha = 1.28166865 * c(-1, 1),
                      d = 0, power = 0.800013981, n = c(n1 = 8222, n2 = 4111), n.total = 12333))

    crrRes <- power.t.welch(d = 0, margin = c(-0.05, 0.05), n.ratio = 2, var.ratio = 2, n2 = 4111, alternative = "two.one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "+", var.ratio = 2, n.ratio = 2,
                          n2 = 4111, power = NULL, alpha = 0.05, alternative = "two.one.sided", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 10960.44441, ncp = 0, null.ncp = 2.926554671 * c(-1, 1), t.alpha = 1.28166865 * c(-1, 1),
                      d = 0, power = 0.800013981, n = c(n1 = 8222, n2 = 4111), n.total = 12333))

    crrRes <- power.t.welch(d = 0.05, margin = c(-0.10, -0.05), n.ratio = 2, var.ratio = 2, power = 0.80,
                            alternative = "two.one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.10, -0.05), req.sign = "+", var.ratio = 2, n.ratio = 2,
                          n2 = NULL, power = 0.80, alpha = 0.05, alternative = "two.one.sided", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 2509.77762, ncp = 1.400942178, null.ncp = c(-2.801884357, -1.400942178),
                      t.alpha = c(-4.766751724, 0.559138526), d = 0.05, power = 0.80005932, n = c(n1 = 1884, n2 = 942),
                      n.total = 2826))

    crrRes <- power.t.welch(d = 0.05, margin = c(-0.10, -0.05), n.ratio = 2, var.ratio = 2, n2 = 942,
                            alternative = "two.one.sided", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = 0.05, null.d = 0, margin = c(-0.10, -0.05), req.sign = "+", var.ratio = 2, n.ratio = 2,
                          n2 = 942, power = NULL, alpha = 0.05, alternative = "two.one.sided", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 2509.77762, ncp = 1.400942178, null.ncp = c(-2.801884357, -1.400942178),
                      t.alpha = c(-4.766751724, 0.559138526), d = 0.05, power = 0.80005932, n = c(n1 = 1884, n2 = 942),
                      n.total = 2826))

    crrRes <- suppressWarnings(power.t.welch(margin = c(-0.05, 0.05), req.sign = "0", var.ratio = 2, n.ratio = 2, n2 = 4111,
                                             power = 0.80, alternative = "two.one.sided", verbose = 0))
    expect_s3_class(crrRes, c("pwrss", "t", "welch"))
    expect_named(crrRes, c("parms", "test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(d = NULL, null.d = 0, margin = 0.05 * c(-1, 1), req.sign = "0", var.ratio = 2, n.ratio = 2,
                          n2 = 4111, power = 0.80, alpha = 0.05, alternative = "two.one.sided", claim.basis = "md.pval",
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "df", "ncp", "null.ncp", "t.alpha", "d", "power", "n", "n.total")],
                 list(test = "t", df = 10960.44441, ncp = 0, null.ncp = 2.926554671 * c(-1, 1), t.alpha = 1.28166865 * c(-1, 1),
                      d = 0, power = 0.800013981, n = c(n1 = 8222, n2 = 4111), n.total = 12333))

    expect_equal(power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 1, power = 0.80, alternative = "two.sided", verbose = 0),
                 pwrss.t.2means(mu1 = 0.20, sd1 = 1, kappa = 2, power = 0.8, alternative = "not equal", verbose = FALSE))

    expect_equal(power.t.welch(d = 0.20, n.ratio = 2, var.ratio = 1, power = 0.80, alternative = "one.sided", verbose = 0),
                 pwrss.t.2means(mu1 = 0.20, sd1 = 1, kappa = 2, power = 0.8, alternative = "less", verbose = FALSE))

    expect_equal(power.t.welch(d = 0.20, margin = c(-0.1, 0.1), n.ratio = 2, var.ratio = 1, power = 0.80, alternative = "two.one.sided", verbose = 0),
                 pwrss.t.2means(mu1 = 0.20, sd1 = 1, margin = 0.10, kappa = 2, power = 0.8, alternative = "equivalent", verbose = FALSE))

    expect_equal(power.t.welch(d = 0.20, n.ratio = 1 / 2, var.ratio = 1, n2 = 400, alternative = "two.sided", verbose = 0),
                 pwrss.t.2means(mu1 = 0.20, sd1 = 1, kappa = 1 / 2, n2 = 400, alternative = "not equal", verbose = FALSE))

    expect_equal(power.t.welch(d = 0.20, n.ratio = 1 / 2, var.ratio = 1, n2 = 400, alternative = "one.sided", verbose = 0),
                 pwrss.t.2means(mu1 = 0.20, sd1 = 1, kappa = 1 / 2, n2 = 400, alternative = "less", verbose = FALSE))

    expect_equal(power.t.welch(d = 0.20, margin = c(-0.1, 0.1), n.ratio = 2, var.ratio = 1, n2 = 800, alternative = "two.one.sided", verbose = 0),
                 pwrss.t.2means(mu1 = 0.20, sd1 = 1, margin = 0.1, kappa = 2, n2 = 800, alternative = "equivalent", verbose = FALSE))

    # expect_error(power.t.welch(margin = c(-0.05, 0.05), n2 = 7517, power = 0.80, alternative = "two.one.sided"),
    #              "Determining the effect size is not possible if `alternative` is \"two.one.sided\".")
    expect_error(power.t.welch(d = 1e-6, power = 0.99, alpha = 1e-6, alternative = "two.sided"),
                 "Design is not feasible.")
    expect_error(power.t.welch(n2 = 5, power = 0.99, alpha = 1e-6, alternative = "two.sided"),
                 "Design is not feasible.")
    expect_error(power.t.welch(d = 0.2, n2 = 2, alpha = 1e-6, alternative = "two.sided"),
                 "Degrees of freedom can not be smaller than 3.")
    expect_warning(power.t.welch(margin = c(-0.05, 0.05), req.sign = "0", n2 = 20, power = 0.80, alternative = "two.one.sided", verbose = 0),
                   "The target power rate cannot be achieved within the null bounds.")
    expect_warning(pwrss.t.2means(mu1 = 0.20, sd1 = 1, power = 0.8, welch.df = FALSE, alternative = "not equal", verbose = FALSE),
                   "Forcing welch.df = TRUE.")
    expect_warning(pwrss.t.2means(mu1 = -0.20, sd1 = 1, margin = 0.1, power = 0.80, alternative = "not equal", verbose = FALSE),
                   "Margin is forced to be 0 for the 'two.sided' test.")
})

# pwrss.z.mean, pwrss.z.2means (not longer supported) ------------------------------------------------------------------
test_that("pwrss.z.mean / pwrss.z.2means return error", {
    expect_error(pwrss.z.mean(),   "This function is no longer available. Please use `power.t.student\\(\\)`.")
    expect_error(pwrss.z.2means(), "This function is no longer available. Please use `power.t.student\\(\\)` or `power.t.welch\\(\\)`.")
})
