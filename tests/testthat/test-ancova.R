# power.f.ancova (= pwrss.f.ancova) ------------------------------------------------------------------------------------
test_that("power.f.ancova / pwrss.f.ancova work", {
    crrRes <- power.f.ancova(eta.squared = 0.059, factor.levels = 2, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.059, null.eta.squared = 0, factor.levels = 2, target.effect = NULL, k.covariates = 0,
                          n.total = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.059, df1 = 1, df2 = 126, ncp = 8.02550478, null.ncp = 0,
                      f.alpha = 3.91632464, power = 0.8027032, n.total = 128))
    expect_identical(crrRes, pwrss.f.ancova(eta2 = 0.059, n.levels = 2, power = 0.80, alpha = 0.05, verbose = 0))
    expect_identical(crrRes, pwrss.f.ancova(f2 = 0.059 / (1 - 0.059), n.levels = 2, power = 0.80, alpha = 0.05, verbose = 0))

    crrRes <- power.f.ancova(eta.squared = 0.059, factor.levels = 4, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.059, null.eta.squared = 0, factor.levels = 4, target.effect = NULL, k.covariates = 0,
                          n.total = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(4)", eta.squared = 0.059, df1 = 3, df2 = 176, ncp = 11.2858661, null.ncp = 0,
                      f.alpha = 2.65593888, power = 0.805367138, n.total = 180))

    crrRes <- power.f.ancova(eta.squared = 0.059, factor.levels = 4, alpha = 0.05, n.total = 180, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.059, null.eta.squared = 0, factor.levels = 4, target.effect = NULL, k.covariates = 0,
                          n.total = 180, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(4)", eta.squared = 0.059, df1 = 3, df2 = 176, ncp = 11.2858661, null.ncp = 0,
                      f.alpha = 2.65593888, power = 0.805367138, n.total = 180))

    crrRes <- power.f.ancova(factor.levels = 4, alpha = 0.05, n.total = 180, power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = 4, target.effect = NULL, k.covariates = 0,
                          n.total = 180, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(4)", eta.squared = 0.058319032, df1 = 3, df2 = 176, ncp = 11.14753955, null.ncp = 0,
                      f.alpha = 2.65593888, power = 0.8, n.total = 180))

    crrRes <- power.f.ancova(eta.squared = 0.030, factor.levels = c(2, 2), power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.030, null.eta.squared = 0, factor.levels = c(2, 2), target.effect = NULL, k.covariates = 0,
                          n.total = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2):B(2)", eta.squared = 0.030, df1 = 1, df2 = 252, ncp = 7.91752577, null.ncp = 0,
                      f.alpha = 3.87862445, power = 0.800416655, n.total = 256))

    crrRes <- power.f.ancova(eta.squared = 0.030, factor.levels = c(2, 2), n.total = 256, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.030, null.eta.squared = 0, factor.levels = c(2, 2), target.effect = NULL, k.covariates = 0,
                          n.total = 256, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2):B(2)", eta.squared = 0.030, df1 = 1, df2 = 252, ncp = 7.91752577, null.ncp = 0,
                      f.alpha = 3.87862445, power = 0.800416655, n.total = 256))

    crrRes <- power.f.ancova(factor.levels = c(2, 2), n.total = 256, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = c(2, 2), target.effect = NULL, k.covariates = 0,
                          n.total = 256, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2):B(2)", eta.squared = 0.029969086, df1 = 1, df2 = 252, ncp = 7.909115,
                      null.ncp = 0, f.alpha = 3.87862445, power = 0.8, n.total = 256))

    crrRes <- power.f.ancova(eta.squared = 0.048, factor.levels = 2, k.covariates = 1, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.048, null.eta.squared = 0, factor.levels = 2, target.effect = NULL, k.covariates = 1,
                          n.total = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.048, df1 = 1, df2 = 155, ncp = 7.96638655, null.ncp = 0,
                      f.alpha = 3.90215432, power = 0.8009416, n.total = 158))

    crrRes <- power.f.ancova(eta.squared = 0.048, factor.levels = 2, k.covariates = 1, n.total = 158, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.048, null.eta.squared = 0, factor.levels = 2, target.effect = NULL, k.covariates = 1,
                          n.total = 158, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.048, df1 = 1, df2 = 155, ncp = 7.96638655, null.ncp = 0,
                      f.alpha = 3.90215432, power = 0.8009416, n.total = 158))

    crrRes <- power.f.ancova(factor.levels = 2, k.covariates = 1, n.total = 158, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = 2, target.effect = NULL, k.covariates = 1,
                          n.total = 158, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.047890306, df1 = 1, df2 = 155, ncp = 7.9472653,
                      null.ncp = 0, f.alpha = 3.90215432, power = 0.8, n.total = 158))

    crrRes <- power.f.ancova(eta.squared = 0.020, factor.levels = c(2, 2), k.covariates = 1, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.020, null.eta.squared = 0, factor.levels = c(2, 2), target.effect = NULL, k.covariates = 1,
                          n.total = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2):B(2)", eta.squared = 0.020, df1 = 1, df2 = 383, ncp = 7.9183673, null.ncp = 0,
                      f.alpha = 3.8658527, power = 0.80148462, n.total = 388))

    crrRes <- power.f.ancova(eta.squared = 0.020, factor.levels = c(2, 2, 2), k.covariates = 1, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.020, null.eta.squared = 0, factor.levels = c(2, 2, 2), target.effect = NULL, k.covariates = 1,
                          n.total = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2):B(2):C(2)", eta.squared = 0.020, df1 = 1, df2 = 383, ncp = 8, null.ncp = 0,
                      f.alpha = 3.86585275, power = 0.8054821, n.total = 392))

    crrRes <- power.f.ancova(eta.squared = 0.020, factor.levels = c(2, 2, 2), k.covariates = 1, n.total = 392, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.020, null.eta.squared = 0, factor.levels = c(2, 2, 2), target.effect = NULL, k.covariates = 1,
                          n.total = 392, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2):B(2):C(2)", eta.squared = 0.020, df1 = 1, df2 = 383, ncp = 8, null.ncp = 0,
                      f.alpha = 3.86585275, power = 0.8054821, n.total = 392))

    crrRes <- power.f.ancova(factor.levels = c(2, 2, 2), k.covariates = 1, n.total = 392, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = c(2, 2, 2), target.effect = NULL, k.covariates = 1,
                          n.total = 392, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(2):B(2):C(2)", eta.squared = 0.0197265246, df1 = 1, df2 = 383,
                      ncp = 7.8884085, null.ncp = 0, f.alpha = 3.86585275, power = 0.8, n.total = 392))

    etaSqd <- f.to.etasq(0.25, verbose = 0)$eta.squared
    crrRes <- power.f.ancova(eta.squared = etaSqd, factor.levels = 10, alpha = 0.05, power = 0.95, verbose = 0) # example 10.4 from GPower
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = etaSqd, null.eta.squared = 0, factor.levels = 10, target.effect = NULL,
                          k.covariates = 0, n.total = NULL, power = 0.95, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A(10)", eta.squared = 0.05882353, df1 = 9, df2 = 380, ncp = 24.375, null.ncp = 0,
                      f.alpha = 1.90453773, power = 0.95236341, n.total = 390))
    # results identical: power ~ 0.952363, n = 390, ncp = 24.375

    etaSqd <- f.to.etasq(0.7066856, verbose = 0)$eta.squared
    crrRes <- power.f.ancova(eta.squared = etaSqd, factor.levels = c(3, 3, 4), target.effect = "A", alpha = 0.05,
                             n.total = 108, verbose = 0) # example 11.3.1 from GPower, main effect
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = etaSqd, null.eta.squared = 0, factor.levels = c(3, 3, 4), target.effect = "A",
                          k.covariates = 0, n.total = 108, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A from A(3):B(3):C(4)", eta.squared = 0.333068578, df1 = 2, df2 = 72,
                      ncp = 53.935690, null.ncp = 0, f.alpha = 3.12390745, power = 0.9999994, n.total = 108))
    # results identical: power ~ 0.99999, n = 108, ncp ~ 53.935690, f.alpha ~ 3.123907

    etaSqd <- f.to.etasq(0.2450722, verbose = 0)$eta.squared
    crrRes <- power.f.ancova(eta.squared = etaSqd, factor.levels = c(3, 3, 4), target.effect = "A:B", alpha = 0.05,
                             n.total = 108, verbose = 0) # example 11.3.1 from GPower, 2-way interact.
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = etaSqd, null.eta.squared = 0, factor.levels = c(3, 3, 4), target.effect = "A:B",
                          k.covariates = 0, n.total = 108, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A:B from A(3):B(3):C(4)", eta.squared = 0.056657511, df1 = 4, df2 = 72,
                      ncp = 6.48652139, null.ncp = 0, f.alpha = 2.49891858, power = 0.47563458, n.total = 108))
    # results identical: power ~ 0.475635, n = 108, ncp ~ 6.486521, f.alpha ~ 2.498919

    etaSqd <- f.to.etasq(0.3288016, verbose = 0)$eta.squared
    crrRes <- power.f.ancova(eta.squared = etaSqd, factor.levels = c(3, 3, 4),
                             target.effect = "A:B:C", alpha = 0.05, n.total = 108, verbose = 0) # example 11.3.1 from GPower, 3-way interact.
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = etaSqd, null.eta.squared = 0, factor.levels = c(3, 3, 4), target.effect = "A:B:C",
                      k.covariates = 0, n.total = 108, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A:B:C from A(3):B(3):C(4)", eta.squared = 0.097562917, df1 = 12, df2 = 72,
                      ncp = 11.675933, null.ncp = 0, f.alpha = 1.889242097, power = 0.51344242, n.total = 108))
    # results identical: power ~ 0.513442, n = 108, ncp ~ 11.675933, f.alpha ~ 3.123907

    etaSqd <- f.to.etasq(0.1, verbose = 0)$eta.squared
    crrRes <- power.f.ancova(eta.squared = etaSqd, factor.levels = c(3, 2, 5), target.effect = "A:C", alpha = 0.05,
                             n.total = 2283, verbose = 0) # example 11.3.2 from GPower
    expect_s3_class(crrRes, c("pwrss", "f", "ancova"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = etaSqd, null.eta.squared = 0, factor.levels = c(3, 2, 5), target.effect = "A:C",
                          k.covariates = 0, n.total = 2283, power = NULL, alpha = 0.05, ceil.n = TRUE, verbose = 0,
                          utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.total")],
                 list(test = "F", effect = "A:C from A(3):B(2):C(5)", eta.squared = 0.00990099, df1 = 8, df2 = 2253,
                      ncp = 22.83, null.ncp = 0, f.alpha = 1.94250745, power = 0.950077975, n.total = 2283))
    # -> GPower does ceil.n for n.total not per group -> when using power = 0.95 as parameter, GPower returns a N = 2283,
    #    pwrss returns 2310 (i.e., `ceiling(2283 / 30) * 30` = 2310)
    # if n.total = 2283 (instead power = 0.95), results identical: power ~ 0.950078, n = 2283, ncp ~ 22.83, f.alpha ~ 1.942507

    expect_error(power.f.ancova(eta.squared = 0.059, factor.levels = 2, alpha = 0.05, verbose = 0),
                 "Exactly two of the parameters `eta.squared`, `n.total`, or `power` must be given, one has to be NULL.")
    expect_error(power.f.ancova(eta.squared = 0.059, factor.levels = 2, power = 0.80, alpha = 0.05, n.total = 1000, verbose = 0),
                 "Exactly two of the parameters `eta.squared`, `n.total`, or `power` must be given, one has to be NULL.")
    expect_error(power.f.ancova(eta.squared = 0.1, factor.levels = c(2, 2), target.effect = "A:B:C", k.covariates = 0, alpha = 0.01, power = 0.80),
                 paste("Invalid specification of `target.effect`. It must be either a single letter \"A\", \"B\" or \"C\"",
                       "\\(depending on the length of `factor.levels`\\), assessing a main effect, or a combination of",
                       "these letters separated by \":\", e.g., \"A:B\", assessing an interaction."))
    expect_error(power.f.ancova(eta.squared = 1e-9, factor.levels = 8, k.covariates = 0, alpha = 0.01, power = 0.80),
                 "Design is not feasible.")
    expect_error(power.f.ancova(factor.levels = 8, k.covariates = 0, alpha = 0.01, power = 0.80, n.total = 9),
                 "Design is not feasible.")
    expect_error(power.f.ancova(eta.squared = 0.059, factor.levels = rep(2, 4), power = 0.80, alpha = 0.05, verbose = 0),
                 "More than three-way ANOVA or ANCOVA is not allowed at the moment.")
    expect_error(pwrss.f.ancova(eta2 = 0.059, f2 = 0.059 / (1 - 0.059), n.levels = 2, power = 0.80, alpha = 0.05, verbose = 0),
                 "Effect size conflict for the alternative. Specify only either `eta2` or `f2`.")
})

# power.f.ancova.keppel ------------------------------------------------------------------------------------------------
test_that("power.f.ancova.keppel works", {
    crrRes <-  power.f.ancova.keppel(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), p.vector = rep(0.50, 2),
                                     k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "keppel"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), n.vector = NULL, p.vector = rep(0.50, 2),
                          factor.levels = NULL, r.squared = 0.50, k.covariates = 1, power = 0.80, alpha = 0.05,
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.111111111, f = 0.353553391, df1 = 1, df2 = 63, ncp = 8.25,
                      null.ncp = 0, f.alpha = 3.993364924, power = 0.807379456, n.vector = c(33, 33), n.total = 66))

    crrRes <-  power.f.ancova.keppel(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), n.vector = rep(33, 2), p.vector = rep(0.50, 2),
                                     k.covariates = 1, r.squared = 0.50, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "keppel"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), n.vector = rep(33, 2), p.vector = rep(0.50, 2),
                          factor.levels = NULL, r.squared = 0.50, k.covariates = 1, power = NULL, alpha = 0.05,
                          ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.111111111, f = 0.353553391, df1 = 1, df2 = 63, ncp = 8.25,
                      null.ncp = 0, f.alpha = 3.993364924, power = 0.807379456, n.vector = rep(33, 2), n.total = 66))

    crrRes <-  power.f.ancova.keppel(mu.vector = c(0.50, 0.00, 0.40, 0.20), sd.vector = rep(1, 4), n.vector = NULL, p.vector = rep(0.25, 4),
                                     factor.levels = 4, k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "keppel"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.50, 0.00, 0.40, 0.20), sd.vector = rep(1, 4), n.vector = NULL,
                          p.vector = rep(0.25, 4), factor.levels = 4, r.squared = 0.50, k.covariates = 1, power = 0.80,
                          alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(4)", eta.squared = 0.068684517, f = 0.27156951, df1 = 3, df2 = 147, ncp = 11.21,
                      null.ncp = 0, f.alpha = 2.6661488, power = 0.800524082, n.vector = rep(38, 4), n.total = 152))

    crrRes <-  power.f.ancova.keppel(mu.vector = c(0.50, 0.00, 0.40, 0.20), sd.vector = rep(1, 4), n.vector = rep(38, 4),
                                     p.vector = rep(0.25, 4), factor.levels = 4, k.covariates = 1, r.squared = 0.50,
                                     alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "keppel"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.50, 0.00, 0.40, 0.20), sd.vector = rep(1, 4), n.vector = rep(38, 4),
                          p.vector = rep(0.25, 4), factor.levels = 4, r.squared = 0.50, k.covariates = 1, power = NULL,
                          alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(4)", eta.squared = 0.068684517, f = 0.27156951, df1 = 3, df2 = 147, ncp = 11.21,
                      null.ncp = 0, f.alpha = 2.6661488, power = 0.800524082, n.vector = rep(38, 4), n.total = 152))

    expect_error(power.f.ancova.keppel(mu.vector = NULL, k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05, verbose = 0),
                 "Provide a vector of means \\(`mu.vector`\\) with its length equal to number of groups.")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00, 0.20, 0.30), sd.vector = rep(1, 4), p.vector = rep(0.50, 4),
                                       factor.levels = c(2, 2), k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05),
                 "Factorial designs are not allowed in Keppel's approach.")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00, 0.20, 0.30), sd.vector = rep(1, 4), p.vector = rep(0.50, 4),
                                       factor.levels = 6, k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05),
                 "Length of the vector of means \\(`mu.vector`\\) does not match number of levels.")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), p.vector = rep(0.50, 2),
                                       k.covariates = 1, r.squared = 1.01, power = 0.80, alpha = 0.05, verbose = 0),
                 "R-squared \\(explanatory power of covariates\\) takes a value between 0 and 1.")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), p.vector = rep(0.50, 2),
                                       k.covariates = 0, r.squared = 0.50, power = 0.80, alpha = 0.05, verbose = 0),
                 "Explanatory power of covariates is expected to be non-zero when number of covariates is non-zero.")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), p.vector = rep(0.50, 2),
                                       k.covariates = 1, r.squared = 0.50, alpha = 0.05, verbose = 0),
                 "Exactly one of the parameters `n.vector` or `power` must be given, one has to be NULL.")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), n.vector = rep(50, 2), p.vector = rep(0.50, 2),
                                       k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05, verbose = 0),
                 "Exactly one of the parameters `n.vector` or `power` must be given, one has to be NULL.")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00, 0.20, 0.30), sd.vector = rep(1, 4),
                                       factor.levels = 4, k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05),
                 "`p.vector` can not be NULL when sample size is requested")
    expect_error(power.f.ancova.keppel(mu.vector = c(0.50, 0.00, 0.20, 0.30), sd.vector = rep(1, 4), p.vector = rep(0.24, 4),
                                       factor.levels = 4, k.covariates = 1, r.squared = 0.50, power = 0.80, alpha = 0.05),
                 "The elements of the `p.vector` should sum to 1")
})

# factorial.contrasts --------------------------------------------------------------------------------------------------
test_that("factorial.contrasts works", {
    crrRes <- factorial.contrasts(factor.levels = 3, coding = "deviation", verbose = 0)
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, 3)
    expect_equal(crrRes$factor.data, data.frame(A = gl(3, 1)))
    expect_equal(crrRes$model.matrix, model.matrix(~ A, data.frame(A = as.factor(seq(3))), contrasts = list(A = contr.sum(3, 3))))
    expect_equal(crrRes$contrast.matrix, matrix(c(2 / 3, -1 / 3, -1 / 3, 2 / 3, -1 / 3, -1 / 3), nrow = 2,
                                                dimnames = list(c("A1", "A2"), c("A1", "A2", "A3"))))
    expect_equal(capture.output(factorial.contrasts(factor.levels = 3, coding = "deviation", verbose = 1)),
                 c("       A1     A2     A3", "A1  0.667 -0.333 -0.333", "A2 -0.333  0.667 -0.333"))

    crrRes <- factorial.contrasts(factor.levels = 3, coding = "deviation", intercept = TRUE, verbose = 0)
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, 3)
    expect_equal(crrRes$factor.data, data.frame(A = gl(3, 1)))
    expect_equal(crrRes$model.matrix, model.matrix(~ A, data.frame(A = as.factor(seq(3))), contrasts = list(A = contr.sum(3, 3))))
    expect_equal(crrRes$contrast.matrix, matrix(c(1 / 3, 2 / 3, -1 / 3, 1 / 3, -1 / 3, 2 / 3, 1 / 3, -1 / 3, -1 / 3), nrow = 3,
                                                dimnames = list(c("(Intercept)", "A1", "A2"), c("A1", "A2", "A3"))))
    expect_equal(capture.output(factorial.contrasts(factor.levels = 3, coding = "deviation", verbose = 1)),
                 c("       A1     A2     A3", "A1  0.667 -0.333 -0.333", "A2 -0.333  0.667 -0.333"))

    crrRes <- factorial.contrasts(factor.levels = 3, coding = "treatment", verbose = 0)
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, 3)
    expect_equal(crrRes$factor.data, data.frame(A = gl(3, 1)))
    expect_equal(crrRes$model.matrix, model.matrix(~ A, data.frame(A = as.factor(seq(3))), contrasts = list(A = contr.treatment(3, 3))))
    expect_equal(crrRes$contrast.matrix, matrix(c(1, 0, 0, 1, -1, -1), nrow = 2, dimnames = list(c("A1", "A2"), c("A1", "A2", "A3"))))
    expect_equal(capture.output(factorial.contrasts(factor.levels = 3, coding = "treatment", verbose = 1)),
                 c("   A1 A2 A3", "A1  1  0 -1", "A2  0  1 -1"))

    crrRes <- factorial.contrasts(factor.levels = 3, base = 1, coding = "treatment", verbose = 0)
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, 3)
    expect_equal(crrRes$factor.data, data.frame(A = gl(3, 1)))
    expect_equal(crrRes$model.matrix, model.matrix(~ A, data.frame(A = as.factor(seq(3))), contrasts = list(A = contr.treatment(3, 1))))
    expect_equal(crrRes$contrast.matrix, matrix(c(-1, -1, 1, 0, 0, 1), nrow = 2, dimnames = list(c("A2", "A3"), c("A1", "A2", "A3"))))
    expect_equal(capture.output(factorial.contrasts(factor.levels = 3, coding = "treatment", verbose = 1)),
                 c("   A1 A2 A3", "A1  1  0 -1", "A2  0  1 -1"))

    crrRes <- factorial.contrasts(factor.levels = 3, coding = "helmert", verbose = 0)
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, 3)
    expect_equal(crrRes$factor.data, data.frame(A = gl(3, 1)))
    expect_equal(crrRes$model.matrix, model.matrix(~ A, data.frame(A = as.factor(seq(3))), contrasts = list(A = contr.helmert(3, 3))))
    expect_equal(crrRes$contrast.matrix, matrix(c(-1 / 2, -1 / 6, 1 / 2, -1 / 6, 0, 1 / 3), nrow = 2,
                                                dimnames = list(c("A1", "A2"), c("A1", "A2", "A3"))))
    expect_equal(capture.output(factorial.contrasts(factor.levels = 3, coding = "helmert", verbose = 1)),
                 c("       A1     A2    A3", "A1 -0.500  0.500 0.000", "A2 -0.167 -0.167 0.333"))

    crrRes <- factorial.contrasts(factor.levels = 3, coding = "poly", verbose = 0)
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, 3)
    expect_equal(crrRes$factor.data, data.frame(A = gl(3, 1)))
    expect_equal(crrRes$model.matrix, model.matrix(~ A, data.frame(A = as.factor(seq(3))), contrasts = list(A = contr.poly)))
    expect_equal(crrRes$contrast.matrix, t(matrix(as.vector(contr.poly(3)), ncol = 2, dimnames = list(c("A1", "A2", "A3"), c("A.L", "A.Q")))))
    expect_equal(capture.output(factorial.contrasts(factor.levels = 3, coding = "poly", verbose = 1)),
                 c("        A1     A2    A3", "A.L -0.707  0.000 0.707", "A.Q  0.408 -0.816 0.408"))

    crrRes <- suppressMessages(factorial.contrasts(factor.levels = c(4, 3), coding.scheme = "deviation", intercept = FALSE, verbose = 0))
    crrAsc <- suppressMessages(capture.output(factorial.contrasts(factor.levels = c(4, 3), coding.scheme = "deviation",
                                                                  intercept = FALSE, verbose = 1)))
    crrMdM <- model.matrix(~ A * B, data.frame(A = gl(4, 3), B = gl(3, 1, 3 * 4)), contrasts = list(A = contr.sum, B = contr.sum))
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, c(4, 3))
    expect_equal(crrRes$factor.data, data.frame(A = gl(4, 3), B = gl(3, 1, 3 * 4)))
    expect_equal(crrRes$model.matrix, crrMdM)
    expect_equal(unname(crrRes$contrast.matrix), unname(solve(crrMdM)[-1, ]))
    expect_equal(dimnames(crrRes$contrast.matrix),
                 list(c("A1", "A2", "A3", "B1", "B2", "A1:B1", "A2:B1", "A3:B1", "A1:B2", "A2:B2", "A3:B2"),
                      c("A1:B1", "A1:B2", "A1:B3", "A2:B1", "A2:B2", "A2:B3", "A3:B1", "A3:B2", "A3:B3", "A4:B1", "A4:B2", "A4:B3")))
    expect_equal(crrAsc, c("       A1:B1  A1:B2  A1:B3  A2:B1  A2:B2  A2:B3  A3:B1  A3:B2  A3:B3  A4:B1",
                           "A1     0.250  0.250  0.250 -0.083 -0.083 -0.083 -0.083 -0.083 -0.083 -0.083",
                           "A2    -0.083 -0.083 -0.083  0.250  0.250  0.250 -0.083 -0.083 -0.083 -0.083",
                           "A3    -0.083 -0.083 -0.083 -0.083 -0.083 -0.083  0.250  0.250  0.250 -0.083",
                           "B1     0.167 -0.083 -0.083  0.167 -0.083 -0.083  0.167 -0.083 -0.083  0.167",
                           "B2    -0.083  0.167 -0.083 -0.083  0.167 -0.083 -0.083  0.167 -0.083 -0.083",
                           "A1:B1  0.500 -0.250 -0.250 -0.167  0.083  0.083 -0.167  0.083  0.083 -0.167",
                           "A2:B1 -0.167  0.083  0.083  0.500 -0.250 -0.250 -0.167  0.083  0.083 -0.167",
                           "A3:B1 -0.167  0.083  0.083 -0.167  0.083  0.083  0.500 -0.250 -0.250 -0.167",
                           "A1:B2 -0.250  0.500 -0.250  0.083 -0.167  0.083  0.083 -0.167  0.083  0.083",
                           "A2:B2  0.083 -0.167  0.083 -0.250  0.500 -0.250  0.083 -0.167  0.083  0.083",
                           "A3:B2  0.083 -0.167  0.083  0.083 -0.167  0.083 -0.250  0.500 -0.250  0.083",
                           "       A4:B2  A4:B3", "A1    -0.083 -0.083", "A2    -0.083 -0.083", "A3    -0.083 -0.083",
                           "B1    -0.083 -0.083", "B2     0.167 -0.083", "A1:B1  0.083  0.083", "A2:B1  0.083  0.083",
                           "A3:B1  0.083  0.083", "A1:B2 -0.167  0.083", "A2:B2 -0.167  0.083", "A3:B2 -0.167  0.083"))

    crrRes <- suppressMessages(factorial.contrasts(factor.levels = c(4, 3), coding.scheme = "poly", intercept = FALSE, verbose = 0))
    crrAsc <- suppressMessages(capture.output(factorial.contrasts(factor.levels = c(4, 3), coding.scheme = "poly",
                                                                  intercept = FALSE, verbose = 1)))
    crrMdM <- model.matrix(~ A * B, data.frame(A = gl(4, 3), B = gl(3, 1, 3 * 4)), contrasts = list(A = contr.poly, B = contr.poly))
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, c(4, 3))
    expect_equal(crrRes$factor.data, data.frame(A = gl(4, 3), B = gl(3, 1, 3 * 4)))
    expect_equal(crrRes$model.matrix, crrMdM)
    expect_equal(unname(crrRes$contrast.matrix), unname(solve(crrMdM)[-1, ]))
    expect_equal(dimnames(crrRes$contrast.matrix),
                 list(c("A.L", "A.Q", "A.C", "B.L", "B.Q", "A.L:B.L", "A.Q:B.L", "A.C:B.L", "A.L:B.Q", "A.Q:B.Q", "A.C:B.Q"),
                      c("A1:B1", "A1:B2", "A1:B3", "A2:B1", "A2:B2", "A2:B3", "A3:B1", "A3:B2", "A3:B3", "A4:B1", "A4:B2", "A4:B3")))
    expect_equal(crrAsc, c("         A1:B1  A1:B2  A1:B3  A2:B1  A2:B2  A2:B3  A3:B1  A3:B2  A3:B3  A4:B1",
                           "A.L     -0.224 -0.224 -0.224 -0.075 -0.075 -0.075  0.075  0.075  0.075  0.224",
                           "A.Q      0.167  0.167  0.167 -0.167 -0.167 -0.167 -0.167 -0.167 -0.167  0.167",
                           "A.C     -0.075 -0.075 -0.075  0.224  0.224  0.224 -0.224 -0.224 -0.224  0.075",
                           "B.L     -0.177  0.000  0.177 -0.177  0.000  0.177 -0.177  0.000  0.177 -0.177",
                           "B.Q      0.102 -0.204  0.102  0.102 -0.204  0.102  0.102 -0.204  0.102  0.102",
                           "A.L:B.L  0.474  0.000 -0.474  0.158  0.000 -0.158 -0.158  0.000  0.158 -0.474",
                           "A.Q:B.L -0.354  0.000  0.354  0.354  0.000 -0.354  0.354  0.000 -0.354 -0.354",
                           "A.C:B.L  0.158  0.000 -0.158 -0.474  0.000  0.474  0.474  0.000 -0.474 -0.158",
                           "A.L:B.Q -0.274  0.548 -0.274 -0.091  0.183 -0.091  0.091 -0.183  0.091  0.274",
                           "A.Q:B.Q  0.204 -0.408  0.204 -0.204  0.408 -0.204 -0.204  0.408 -0.204  0.204",
                           "A.C:B.Q -0.091  0.183 -0.091  0.274 -0.548  0.274 -0.274  0.548 -0.274  0.091",
                           "         A4:B2 A4:B3", "A.L      0.224 0.224", "A.Q      0.167 0.167", "A.C      0.075 0.075",
                           "B.L      0.000 0.177", "B.Q     -0.204 0.102", "A.L:B.L  0.000 0.474", "A.Q:B.L  0.000 0.354",
                           "A.C:B.L  0.000 0.158", "A.L:B.Q -0.548 0.274", "A.Q:B.Q -0.408 0.204", "A.C:B.Q -0.183 0.091"))

    crrRes <- suppressMessages(factorial.contrasts(factor.levels = c(4, 3), coding.scheme = c("deviation", "poly"), intercept = FALSE, verbose = 0))
    crrAsc <- suppressMessages(capture.output(factorial.contrasts(factor.levels = c(4, 3), coding.scheme = c("deviation", "poly"),
                                                                  intercept = FALSE, verbose = 1)))
    crrMdM <- model.matrix(~ A * B, data.frame(A = gl(4, 3), B = gl(3, 1, 3 * 4)), contrasts = list(A = contr.sum, B = contr.poly))
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, c(4, 3))
    expect_equal(crrRes$factor.data, data.frame(A = gl(4, 3), B = gl(3, 1, 3 * 4)))
    expect_equal(crrRes$model.matrix, crrMdM)
    expect_equal(unname(crrRes$contrast.matrix), unname(solve(crrMdM)[-1, ]))
    expect_equal(dimnames(crrRes$contrast.matrix),
                 list(c("A1", "A2", "A3", "B.L", "B.Q", "A1:B.L", "A2:B.L", "A3:B.L", "A1:B.Q", "A2:B.Q", "A3:B.Q"),
                      c("A1:B1", "A1:B2", "A1:B3", "A2:B1", "A2:B2", "A2:B3", "A3:B1", "A3:B2", "A3:B3", "A4:B1", "A4:B2", "A4:B3")))
    expect_equal(crrAsc, c("        A1:B1  A1:B2  A1:B3  A2:B1  A2:B2  A2:B3  A3:B1  A3:B2  A3:B3  A4:B1",
                           "A1      0.250  0.250  0.250 -0.083 -0.083 -0.083 -0.083 -0.083 -0.083 -0.083",
                           "A2     -0.083 -0.083 -0.083  0.250  0.250  0.250 -0.083 -0.083 -0.083 -0.083",
                           "A3     -0.083 -0.083 -0.083 -0.083 -0.083 -0.083  0.250  0.250  0.250 -0.083",
                           "B.L    -0.177  0.000  0.177 -0.177  0.000  0.177 -0.177  0.000  0.177 -0.177",
                           "B.Q     0.102 -0.204  0.102  0.102 -0.204  0.102  0.102 -0.204  0.102  0.102",
                           "A1:B.L -0.530  0.000  0.530  0.177  0.000 -0.177  0.177  0.000 -0.177  0.177",
                           "A2:B.L  0.177  0.000 -0.177 -0.530  0.000  0.530  0.177  0.000 -0.177  0.177",
                           "A3:B.L  0.177  0.000 -0.177  0.177  0.000 -0.177 -0.530  0.000  0.530  0.177",
                           "A1:B.Q  0.306 -0.612  0.306 -0.102  0.204 -0.102 -0.102  0.204 -0.102 -0.102",
                           "A2:B.Q -0.102  0.204 -0.102  0.306 -0.612  0.306 -0.102  0.204 -0.102 -0.102",
                           "A3:B.Q -0.102  0.204 -0.102 -0.102  0.204 -0.102  0.306 -0.612  0.306 -0.102",
                           "        A4:B2  A4:B3", "A1     -0.083 -0.083", "A2     -0.083 -0.083", "A3     -0.083 -0.083",
                           "B.L     0.000  0.177", "B.Q    -0.204  0.102", "A1:B.L  0.000 -0.177", "A2:B.L  0.000 -0.177",
                           "A3:B.L  0.000 -0.177", "A1:B.Q  0.204 -0.102", "A2:B.Q  0.204 -0.102", "A3:B.Q  0.204 -0.102"))

    crrRes <- suppressMessages(factorial.contrasts(factor.levels = c(3, 2, 2), coding.scheme = "poly", intercept = FALSE, verbose = 0))
    crrAsc <- suppressMessages(capture.output(factorial.contrasts(factor.levels = c(3, 2, 2), coding.scheme = "poly",
                                                                  intercept = FALSE, verbose = 1)))
    crrMdM <- model.matrix(~ A * B * C, data.frame(A = gl(3, 4), B = gl(2, 2, 12), C = gl(2, 1, 12)),
                           contrasts = list(A = contr.poly, B = contr.poly, C = contr.poly))
    expect_type(crrRes, "list")
    expect_named(crrRes, c("factor.levels", "factor.data", "model.matrix", "contrast.matrix"))
    expect_equal(crrRes$factor.levels, c(3, 2, 2))
    expect_equal(crrRes$factor.data, data.frame(A = gl(3, 4), B = gl(2, 2, 12), C = gl(2, 1, 12)))
    expect_equal(crrRes$model.matrix, crrMdM)
    expect_equal(unname(crrRes$contrast.matrix), unname(solve(crrMdM)[-1, ]))
    expect_equal(dimnames(crrRes$contrast.matrix),
                 list(c("A.L", "A.Q", "B.L", "C.L", "A.L:B.L", "A.Q:B.L", "A.L:C.L", "A.Q:C.L", "B.L:C.L", "A.L:B.L:C.L", "A.Q:B.L:C.L"),
                      c("A1:B1:C1", "A1:B1:C2", "A1:B2:C1", "A1:B2:C2", "A2:B1:C1", "A2:B1:C2",
                        "A2:B2:C1", "A2:B2:C2", "A3:B1:C1", "A3:B1:C2", "A3:B2:C1", "A3:B2:C2")))
    expect_equal(crrAsc, c("            A1:B1:C1 A1:B1:C2 A1:B2:C1 A1:B2:C2 A2:B1:C1 A2:B1:C2 A2:B2:C1",
                           "A.L           -0.177   -0.177   -0.177   -0.177    0.000    0.000    0.000",
                           "A.Q            0.102    0.102    0.102    0.102   -0.204   -0.204   -0.204",
                           "B.L           -0.118   -0.118    0.118    0.118   -0.118   -0.118    0.118",
                           "C.L           -0.118    0.118   -0.118    0.118   -0.118    0.118   -0.118",
                           "A.L:B.L        0.250    0.250   -0.250   -0.250    0.000    0.000    0.000",
                           "A.Q:B.L       -0.144   -0.144    0.144    0.144    0.289    0.289   -0.289",
                           "A.L:C.L        0.250   -0.250    0.250   -0.250    0.000    0.000    0.000",
                           "A.Q:C.L       -0.144    0.144   -0.144    0.144    0.289   -0.289    0.289",
                           "B.L:C.L        0.167   -0.167   -0.167    0.167    0.167   -0.167   -0.167",
                           "A.L:B.L:C.L   -0.354    0.354    0.354   -0.354    0.000    0.000    0.000",
                           "A.Q:B.L:C.L    0.204   -0.204   -0.204    0.204   -0.408    0.408    0.408",
                           "            A2:B2:C2 A3:B1:C1 A3:B1:C2 A3:B2:C1 A3:B2:C2",
                           "A.L            0.000    0.177    0.177    0.177    0.177",
                           "A.Q           -0.204    0.102    0.102    0.102    0.102",
                           "B.L            0.118   -0.118   -0.118    0.118    0.118",
                           "C.L            0.118   -0.118    0.118   -0.118    0.118",
                           "A.L:B.L        0.000   -0.250   -0.250    0.250    0.250",
                           "A.Q:B.L       -0.289   -0.144   -0.144    0.144    0.144",
                           "A.L:C.L        0.000   -0.250    0.250   -0.250    0.250",
                           "A.Q:C.L       -0.289   -0.144    0.144   -0.144    0.144",
                           "B.L:C.L        0.167    0.167   -0.167   -0.167    0.167",
                           "A.L:B.L:C.L    0.000    0.354   -0.354   -0.354    0.354",
                           "A.Q:B.L:C.L   -0.408    0.204   -0.204   -0.204    0.204"))

    expect_equal(factorial.contrasts(factor.levels = c(3, 3), coding = rep("deviation", 2), verbose = -1),
                 factorial.contrasts(factor.levels = c(3, 3), coding = rep("deviation", 2), base = 3, verbose = -1))
    expect_equal(factorial.contrasts(factor.levels = c(3, 3), coding = rep("deviation", 2), verbose = -1),
                 factorial.contrasts(factor.levels = c(3, 3), coding = rep("deviation", 2), base = rep(3, 3), verbose = -1))
    expect_message(factorial.contrasts(factor.levels = 3, coding = c("deviation", "poly", "helmert", "sum"), verbose = 0),
                   "Provide as many coding schemes as number of factors. Using the first 1.")
    crrMsg <- capture_messages(factorial.contrasts(factor.levels = c(2, 2, 2), coding = "deviation", verbose = 0))
    expect_equal(crrMsg, c("Assuming the same coding scheme applies to the other factor(s)\n",
                           paste0("Elements of `mu.vector`, `sd.vector`, `n.vector` or `p.vector` should follow this specific order:\n",
                                  "A1:B1:C1  A1:B1:C2  A1:B2:C1  A1:B2:C2  A2:B1:C1  A2:B1:C2  A2:B2:C1  A2:B2:C2  \n\n")))
    expect_error(factorial.contrasts(factor.levels = 3, coding = "invalid", verbose = 0),
                 "Contrast type \"invalid\" not supported at the moment.")
    expect_error(factorial.contrasts(factor.levels = rep(2, 4), coding = rep("poly", 4), verbose = 0),
                 "This version supports only up to three-way interactions.")
})

# power.f.ancova.shieh / power.t.contrasts / power.t.contrast ----------------------------------------------------------
test_that("power.f.ancova.shieh / power.t.contrasts / power.t.contrast work", {
    crrRes <- power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), n.vector = rep(150, 2), r.squared = 0.50,
                                   k.covariates = 1, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.2, 0.0), sd.vector = rep(1, 2), n.vector = rep(150, 2), p.vector = NULL,
                          factor.levels = NULL, r.squared = 0.5, k.covariates = 1, contrast.matrix = NULL, power = NULL,
                          alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.019543331, f = 0.141183873, df1 = 1, df2 = 297, ncp = 5.97986577,
                      null.ncp = 0, f.alpha = 3.87295916, power = 0.683501582, n.vector = c(150, 150), n.total = 300))

    crrRes <- power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.50, 2), r.squared = 0.50,
                                   k.covariates = 1, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.2, 0.0), sd.vector = rep(1, 2), n.vector = NULL, p.vector = rep(0.5, 2),
                          factor.levels = NULL, r.squared = 0.5, k.covariates = 1, contrast.matrix = NULL, power = 0.80,
                          alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.01955905041, f = 0.141241773, df1 = 1, df2 = 393, ncp = 7.89989848,
                      null.ncp = 0, f.alpha = 3.86522919, power = 0.800619012, n.vector = c(198, 198), n.total = 396))

    crrRes <- power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.50, 2), r.squared = 0.50,
                                   k.covariates = 2, power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.2, 0.0), sd.vector = rep(1, 2), n.vector = NULL, p.vector = rep(0.5, 2),
                          factor.levels = NULL, r.squared = 0.5, k.covariates = 2, contrast.matrix = NULL, power = 0.80,
                          alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(2)", eta.squared = 0.0195109904, f = 0.141064681, df1 = 1, df2 = 394, ncp = 7.91989924,
                      null.ncp = 0, f.alpha = 3.86516859, power = 0.801596275, n.vector = c(199, 199), n.total = 398))

    crrRes <- power.f.ancova.shieh(mu.vector = c(0.20, 0.25, 0.15, 0.05), sd.vector = rep(1, 4), p.vector = rep(0.25, 4),
                                   factor.levels = c(2, 2), r.squared = 0.50, k.covariates = 1, alpha = 0.05,
                                   power = 0.80, verbose = -1)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.20, 0.25, 0.15, 0.05), sd.vector = rep(1, 4), n.vector = NULL,
                          p.vector = rep(0.25, 4), factor.levels = c(2, 2), r.squared = 0.5, k.covariates = 1,
                          contrast.matrix = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = -1, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(2):B(2)", eta.squared = 0.0028036103, f = 0.05302351, df1 = 1, df2 = 2791, ncp = 7.86093347,
                      null.ncp = 0, f.alpha = 3.84479279, power = 0.800332660, n.vector = rep(699, 4), n.total = 2796))

    crrRes <- power.f.ancova.shieh(mu.vector = c(0.20, 0.25, 0.30, 0.15, 0.05, 0.10, 0.00, 0.05), sd.vector = rep(1, 8),
                                   p.vector = rep(1 / 8, 8), factor.levels = c(2, 2, 2), r.squared = 0.50, k.covariates = 1,
                                   power = 0.80, alpha = 0.05, verbose = -1)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.20, 0.25, 0.30, 0.15, 0.05, 0.10, 0.00, 0.05), sd.vector = rep(1, 8),
                          n.vector = NULL, p.vector = rep(1 / 8, 8), factor.levels = c(2, 2, 2), r.squared = 0.5,
                          k.covariates = 1, contrast.matrix = NULL, power = 0.80, alpha = 0.05, ceil.n = TRUE,
                          verbose = -1, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(2):B(2):C(2)", eta.squared = 0.0012482409, f = 0.035352524, df1 = 1, df2 = 6279,
                      ncp = 7.85874841, null.ncp = 0, f.alpha = 3.84294023, power = 0.800373625, n.vector = rep(786, 8), n.total = 6288))

    mtxCnt <- factorial.contrasts(factor.levels = 3, coding = "treatment", verbose = 0)$contrast.matrix
    crrRes <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), p.vector = rep(1 / 3, 3),
                                   contrast.matrix = mtxCnt, r.squared = 0.50, k.covariates = 1, power = 0.80,
                                   alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL,
                          p.vector = rep(1 / 3, 3), factor.levels = NULL, r.squared = 0.5, k.covariates = 1,
                          contrast.matrix = mtxCnt, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(3)", eta.squared = 0.0077115847, f = 0.0881562, df1 = 2, df2 = 1241, ncp = 9.67553677,
                      null.ncp = 0, f.alpha = 3.00297552, power = 0.800762614, n.vector = rep(415, 3), n.total = 1245))
    crrCnt <- power.t.contrasts(crrRes, adjust.alpha = "fdr", verbose = 0)
    expect_s3_class(crrCnt, c("pwrss", "t", "contrasts"))
    expect_named(crrCnt, c("parms", "test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = rep(415, 3),
                          p.vector = rep(1 / 3, 3), r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt,
                          power = NULL, alpha = 0.05, adjust.alpha = "fdr", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power")],
                 list(test = "t", contrast = seq(2), comparison = c("A1 <=> A3", "A2 <=> A3"), psi = c(-0.05, 0.10),
                      d = c(-0.070710678, 0.141421356), ncp = c(-1.018167467, 2.036334933), df = rep(1241, 2),
                      t.alpha = rep(2.24412588, 2), n.total = rep(1245, 2), power = c(0.110986043, 0.417964282)))

    mtxCnt <- factorial.contrasts(factor.levels = 3, coding = "helmert", verbose = 0)$contrast.matrix
    crrRes <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), p.vector = rep(1 / 3, 3),
                                   contrast.matrix = mtxCnt, r.squared = 0.50, k.covariates = 1, alpha = 0.05,
                                   power = 0.80, verbose = 0)
    crrCnt <- power.t.contrasts(crrRes, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL,
                          p.vector = rep(1 / 3, 3), factor.levels = NULL, r.squared = 0.5, k.covariates = 1,
                          contrast.matrix = mtxCnt, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(3)", eta.squared = 0.0077115847, f = 0.0881562, df1 = 2, df2 = 1241, ncp = 9.67553677,
                      null.ncp = 0, f.alpha = 3.00297552, power = 0.800762614, n.vector = rep(415, 3), n.total = 1245))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrasts"))
    expect_named(crrCnt, c("parms", "test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = rep(415, 3),
                          p.vector = rep(1 / 3, 3), r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt,
                          power = NULL, alpha = 0.05, adjust.alpha = "none", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power")],
                 list(test = "t", contrast = seq(2), comparison = c("A2 <=> A1", "A3 <=> A1 A2"), psi = c(0.075, -0.0083333333),
                      d = c(0.106066017, -0.011785113), ncp = c(3.0545024, -0.5878393), df = rep(1241, 2),
                      t.alpha = rep(1.9618774, 2), n.total = rep(1245, 2), power = c(0.86262155, 0.09036882)))
    crrCnt <- do.call(power.t.contrast,
                      c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates", "alpha")],
                        list(contrast.vector = mtxCnt[1, ], power = 0.8, verbose = 0)))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL, p.vector = rep(1 / 3, 3),
                          contrast.vector = mtxCnt[1, ], r.squared = 0.5, k.covariates = 1, power = 0.80, alpha = 0.05,
                          tukey.kramer = FALSE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = 0.075, d = 0.106066017, df = 1046, t.alpha = c(-1.9622345, 1.9622345), ncp = 2.8049032,
                      null.ncp = 0, power = 0.800208158, n.vector = rep(350, 3), n.total = 1050))

    mtxCnt <- factorial.contrasts(factor.levels = 3, coding = "poly", verbose = 0)$contrast.matrix
    crrRes <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), p.vector = rep(1 / 3, 3),
                                   contrast.matrix = mtxCnt, r.squared = 0.50, k.covariates = 1, alpha = 0.05,
                                   power = 0.80, verbose = 0)
    crrCnt <- power.t.contrasts(crrRes, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL, p.vector = rep(1 / 3, 3),
                          factor.levels = NULL, r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt, power = 0.80,
                          alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(3)", eta.squared = 0.0077115847, f = 0.0881562, df1 = 2, df2 = 1241, ncp = 9.67553677,
                      null.ncp = 0, f.alpha = 3.00297552, power = 0.800762614, n.vector = rep(415, 3), n.total = 1245))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrasts"))
    expect_named(crrCnt, c("parms", "test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = rep(415, 3),
                          p.vector = rep(1 / 3, 3), r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt,
                          power = NULL, alpha = 0.05, adjust.alpha = "none", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power")],
                 list(test = "t", contrast = seq(2), comparison = c("A3 <=> A1", "A1 A3 <=> A2"), psi = c(0.035355339, -0.102062073),
                      d = c(0.05, -0.144337567), ncp = c(1.01816747, -2.9391963), df = rep(1241, 2),
                      t.alpha = rep(1.9618774, 2), n.total = rep(1245, 2), power = c(0.174400616, 0.835704937)))
    expect_equal(crrCnt, do.call(power.t.contrasts,
                                 c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates", "alpha",
                                                       "contrast.matrix")], list(n.vector = crrRes$n.vector, verbose = 0))))

    crrCnt <- do.call(power.t.contrasts,
                      c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates", "alpha")],
                        list(contrast.matrix = mtxCnt[1, ], power = 0.80, verbose = 0))) # linear trend as vector
    expect_s3_class(crrCnt, c("pwrss", "t", "contrasts"))
    expect_named(crrCnt, c("parms", "test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL,
                          p.vector = rep(1 / 3, 3), r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt[1, ],
                          power = 0.80, alpha = 0.05, adjust.alpha = "none", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power")],
                 list(test = "t", contrast = 1, comparison = "A3 <=> A1", psi = 0.035355339, d = 0.05, ncp = 2.80208252,
                      df = 9419, t.alpha = 1.96021588, n.total = 9423, power = 0.80006019))

    crrCnt <- do.call(power.t.contrast,
                      c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates", "alpha")],
                        list(contrast.vector = mtxCnt[1, ], power = 0.8, verbose = 0)))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL,
                          p.vector = rep(1 / 3, 3), contrast.vector = mtxCnt[1, ], r.squared = 0.5, k.covariates = 1,
                          power = 0.80, alpha = 0.05, tukey.kramer = FALSE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = 0.035355339, d = 0.05, df = 9419, t.alpha = c(-1.96021588, 1.96021588), ncp = 2.80208252,
                      null.ncp = 0, power = 0.80006019, n.vector = rep(3141, 3), n.total = 9423))

    mtxCnt <- factorial.contrasts(factor.levels = 3, coding = "poly", verbose = 0)$contrast.matrix
    crrRes <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), p.vector = rep(1 / 3, 3),
                                   contrast.matrix = mtxCnt, r.squared = 0.50, k.covariates = 2, alpha = 0.05,
                                   power = 0.80, verbose = 0)
    crrCnt <- power.t.contrasts(crrRes, adjust.alpha = "tukey", verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL,
                          p.vector = rep(1 / 3, 3), factor.levels = NULL, r.squared = 0.5, k.covariates = 2,
                          contrast.matrix = mtxCnt, power = 0.80, alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(3)", eta.squared = 0.0077054286, f = 0.08812073, df1 = 2, df2 = 1240, ncp = 9.66775275,
                      null.ncp = 0, f.alpha = 3.0029814, power = 0.80041063, n.vector = rep(415, 3), n.total = 1245))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrasts"))
    expect_named(crrCnt, c("parms", "test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = rep(415, 3),
                          p.vector = rep(1 / 3, 3), r.squared = 0.5, k.covariates = 2, contrast.matrix = mtxCnt,
                          power = NULL, alpha = 0.05, adjust.alpha = "tukey", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power")],
                 list(test = "t", contrast = seq(2), comparison = c("A3 <=> A1", "A1 A3 <=> A2"), psi = c(0.035355339, -0.102062073),
                      d = c(0.05, -0.144337567), ncp = c(1.017741531, -2.937966734), df = rep(1240, 2), t.alpha = rep(2.34652362, 2),
                          n.total = rep(1245, 2), power = c(0.09267207, 0.72283065)))
    expect_equal(crrCnt, do.call(power.t.contrasts,
                                 c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates",
                                                       "alpha", "contrast.matrix")],
                                   list(n.vector = crrRes$n.vector, adjust.alpha = "tukey", verbose = 0))))

    crrCnt <- do.call(power.t.contrast,
                      c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates", "alpha")],
                        list(contrast.vector = crrRes$parms$contrast.matrix[2, ], power = 0.8, verbose = 0)))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL, p.vector = rep(1 / 3, 3),
                          contrast.vector = mtxCnt[2, ], r.squared = 0.5, k.covariates = 2, power = 0.80, alpha = 0.05,
                          tukey.kramer = FALSE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = -0.102062073, d = -0.144337567, df = 1132, t.alpha = c(-1.96206183, 1.96206183),
                      ncp = -2.8074314, null.ncp = 0, power = 0.80096821, n.vector = rep(379, 3), n.total = 1137))

    # custom contrasts
    mtxCnt <- rbind(c(A1 = 1, A2 = -0.50, A3 = -0.50), c(A1 = 0.50, A2 = 0.50, A3 = -1))
    crrRes <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), p.vector = rep(1 / 3, 3),
                                   contrast.matrix = mtxCnt, r.squared = 0.50, k.covariates = 1, alpha = 0.05,
                                   power = 0.80, verbose = 0)
    crrCnt <- power.t.contrasts(crrRes, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "ancova", "shieh"))
    expect_named(crrRes,
                 c("parms", "test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power",
                   "n.vector", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL, p.vector = rep(1 / 3, 3),
                          factor.levels = NULL, r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt, power = 0.8,
                          alpha = 0.05, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "eta.squared", "f", "df1", "df2", "ncp", "null.ncp", "f.alpha", "power", "n.vector", "n.total")],
                 list(test = "F", effect = "A(3)", eta.squared = 0.0077115847, f = 0.0881562, df1 = 2, df2 = 1241, ncp = 9.67553677,
                      null.ncp = 0, f.alpha = 3.00297552, power = 0.800762614, n.vector = rep(415, 3), n.total = 1245))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrasts"))
    expect_named(crrCnt, c("parms", "test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = rep(415, 3),
                          p.vector = rep(1 / 3, 3), r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt,
                          power = NULL, alpha = 0.05, adjust.alpha = "none", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power")],
                 list(test = "t", contrast = seq(2), comparison = c("A1 <=> A2 A3", "A1 A2 <=> A3"), psi = c(-0.1, 0.025),
                      d = c(-0.141421356, 0.035355339), ncp = c(-2.35135704, 0.58783926), df = rep(1241, 2),
                      t.alpha = rep(1.9618774, 2), n.total = rep(1245, 2), power = c(0.651581934, 0.090368815)))
    expect_equal(crrCnt, do.call(power.t.contrasts,
                                 c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates",
                                                       "alpha", "contrast.matrix")],
                                   list(n.vector = crrRes$n.vector, verbose = 0))))

    crrCnt <- do.call(power.t.contrasts,
                      c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates", "alpha", "contrast.matrix")],
                        list(power = crrRes$power, verbose = 0)))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrasts"))
    expect_named(crrCnt, c("parms", "test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL, p.vector = rep(1 / 3, 3),
                          r.squared = 0.5, k.covariates = 1, contrast.matrix = mtxCnt, power = crrRes$power, alpha = 0.05,
                          adjust.alpha = "none", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "contrast", "comparison", "psi", "d", "ncp", "df", "t.alpha", "n.total", "power")],
                 list(test = "t", contrast = seq(2), comparison = c("A1 <=> A2 A3", "A1 A2 <=> A3"), psi = c(-0.1, 0.025),
                      d = c(-0.141421356, 0.035355339), ncp = c(-2.806340906, 2.804411203), df = c(1769, 28310),
                      t.alpha = c(1.961305911, 1.960047784), n.total = c(1773, 28314), power = c(0.800904361, 0.800764602)))

    crrCnt <- do.call(power.t.contrast,
                      c(crrRes[["parms"]][c("mu.vector", "sd.vector", "p.vector", "r.squared", "k.covariates", "alpha")],
                        list(contrast.vector = crrRes$parms$contrast.matrix[1, ], power = 0.8, verbose = 0)))
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = NULL, p.vector = rep(1 / 3, 3),
                          contrast.vector = mtxCnt[1, ], r.squared = 0.5, k.covariates = 1, power = 0.80, alpha = 0.05,
                          tukey.kramer = FALSE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = -0.1, d = -0.141421356, df = 1766, t.alpha = c(-1.96130819, 1.96130819), ncp = -2.80396433,
                      null.ncp = 0, power = 0.8002398, n.vector = rep(590, 3), n.total = 1770))

    crrCnt <- power.t.contrast(mu.vector = c(0.20, 0.10), sd.vector = rep(1, 2), n.vector = rep(1500, 2), contrast.vector = c(1, -1),
                               r.squared = 0.50, k.covariates = 1, alpha = 0.05, tukey.kramer = TRUE, verbose = 0)
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.20, 0.10), sd.vector = rep(1, 2), n.vector = rep(1500, 2), p.vector = NULL,
                          contrast.vector = c(1, -1), r.squared = 0.5, k.covariates = 1, power = NULL, alpha = 0.05,
                          tukey.kramer = TRUE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = 0.1, d = 0.141421356, df = 2997, t.alpha = c(-1.96075583, 1.96075583), ncp = 3.87233747,
                      null.ncp = 0, power = 0.972006163, n.vector = rep(1500, 2), n.total = 3000))

    crrCnt <- power.t.contrast(mu.vector = c(0.20, 0.10), sd.vector = rep(1, 2), n.vector = rep(1500, 2), contrast.vector = c(1, -1),
                               r.squared = 0.50, k.covariates = 1, alpha = 0.05, tukey.kramer = FALSE, verbose = 0)
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.20, 0.10), sd.vector = rep(1, 2), n.vector = rep(1500, 2), p.vector = NULL,
                          contrast.vector = c(1, -1), r.squared = 0.5, k.covariates = 1, power = NULL, alpha = 0.05,
                          tukey.kramer = FALSE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = 0.1, d = 0.141421356, df = 2997, t.alpha = c(-1.96075583, 1.96075583), ncp = 3.87233747,
                      null.ncp = 0, power = 0.972006163, n.vector = rep(1500, 2), n.total = 3000))

    crrCnt <- power.t.contrast(mu.vector = c(0.20, 0.10), sd.vector = rep(1, 2), p.vector = rep(0.5, 2), contrast.vector = c(1, -1),
                               r.squared = 0.50, k.covariates = 1, power = 0.972006163, alpha = 0.05, tukey.kramer = FALSE, verbose = 0)
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(0.20, 0.10), sd.vector = rep(1, 2), n.vector = NULL, p.vector = rep(0.5, 2),
                          contrast.vector = c(1, -1), r.squared = 0.5, k.covariates = 1, power = 0.972006163, alpha = 0.05,
                          tukey.kramer = FALSE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = 0.1, d = 0.141421356, df = 2997, t.alpha = c(-1.96075583, 1.96075583), ncp = 3.87233747,
                      null.ncp = 0, power = 0.972006163, n.vector = rep(1500, 2), n.total = 3000))

    crrCnt <- power.t.contrast(mu.vector = c(1.5, 2, 3, 4), sd.vector = rep(2, 4), n.vector = rep(5, 4), contrast.vector = c(-3, -1, 1, 3),
                               k.covariates = 0, alpha = 0.05, verbose = 0) # example 11.3.3 from GPower
    expect_s3_class(crrCnt, c("pwrss", "t", "contrast"))
    expect_named(crrCnt, c("parms", "test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total"))
    expect_identical(crrCnt[["parms"]],
                     list(mu.vector = c(1.5, 2, 3, 4), sd.vector = rep(2, 4), n.vector = rep(5, 4), p.vector = NULL,
                          contrast.vector = c(-3, -1, 1, 3), r.squared = 0, k.covariates = 0, power = NULL, alpha = 0.05,
                          tukey.kramer = FALSE, ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrCnt[c("test", "psi", "d", "df", "t.alpha", "ncp", "null.ncp", "power", "n.vector", "n.total")],
                 list(test = "t", psi = 8.5, d = 4.25, df = 16, t.alpha = 2.1199053 * c(-1, 1), ncp = 2.125,
                      null.ncp = 0, power = 0.514737, n.vector = rep(5, 4), n.total = 20))
    # results are nearly identical, but numerically different since GPower uses the F-distribition
    # power ~ 0.514736, crit. F / t = 4.493998 - sqrt(4.493998) = 2.1199053, ncp = 4.515617 / sqrt(4.515617) ~ 0.2125

    mtxCnt <- c(A1 = 1, A2 = -0.50, A3 = -0.50)
    expect_equal(power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = rep(400, 3), alpha = 0.05,
                                      contrast.matrix = mtxCnt, r.squared = 0.50, k.covariates = 1, verbose = 0)[-1], # [-1] excludes parms
                 power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), sd.vector = rep(1, 3), n.vector = rep(400, 3), alpha = 0.05,
                                      contrast.matrix = t(as.matrix(mtxCnt)), r.squared = 0.50, k.covariates = 1, verbose = 0)[-1])
    expect_equal(power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5, 2), contrast.vector = c(1, -1),
                                  r.squared = 0.50, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0)[-1],
                 power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5, 2),
                                  contrast.vector = matrix(c(1, -1), ncol = 2), r.squared = 0.50, k.covariates = 1,
                                  power = 0.8, alpha = 0.05, verbose = 0)[-1])

    expect_error(power.f.ancova.shieh(mu.vector = NULL, r.squared = 0.50, k.covariates = 1, alpha = 0.05, verbose = 0),
                 "Provide a vector of means \\(`mu.vector`\\) with its length equal to the number of groups.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), n.vector = rep(150, 2), factor.levels = 3,
                                      r.squared = 0.50, k.covariates = 1, alpha = 0.05, verbose = 0),
                 "Provide a vector of means \\(`mu.vector`\\) with its length equal to the the product of `factor.levels`.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), n.vector = rep(150, 2),
                                      r.squared = -0.01, k.covariates = 1, alpha = 0.05, verbose = 0),
                 "R-squared \\(explanatory power of covariates\\) takes a value between 0 and 1.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), n.vector = rep(150, 2),
                                      r.squared = 0.50, k.covariates = 1, alpha = -0.01, verbose = 0),
                 "Type 1 error rate \\(alpha\\) takes a value between 0 and 1.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), n.vector = rep(150, 2), contrast.matrix = data.frame,
                                      r.squared = 0.50, k.covariates = 1, alpha = 0.05, verbose = 0),
                 "Contrast coefficients are not provided in the form of a matrix.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), n.vector = rep(150, 2), contrast.matrix = 1,
                                      r.squared = 0.50, k.covariates = 1, alpha = 0.05, verbose = 0),
                 "The number of columns in the contrast matrix should match number of groups.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), n.vector = rep(150, 2), contrast.matrix = diag(2),
                                      r.squared = 0.50, k.covariates = 1, alpha = 0.05, verbose = 0),
                 "The number of rows in the contrast matrix should be less than or equal to number of groups minus one.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2),
                                      r.squared = 0.50, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                 "`p.vector` can not be NULL when sample size is requested.")
    expect_error(power.f.ancova.shieh(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5 - 1e-5, 2),
                                      r.squared = 0.50, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                 "The elements of the `p.vector` should sum to 1.")
    expect_error(power.f.ancova.shieh(mu.vector = runif(16, 0, 1), sd.vector = rep(1, 16), p.vector = rep(1 / 16, 16),
                                      factor.levels = c(2, 2, 2, 2), contrast.matrix = matrix(runif(240), ncol = 16),
                                      r.squared = 0.50, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                 "More than three-way ANCOVA is not allowed at the moment.")
    expect_error(power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5, 2), contrast.vector = c(1, -1),
                                  r.squared = -0.01, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                 "R-squared \\(explanatory power of covariates\\) takes a value between 0 and 1.")
    expect_error(power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5, 2), contrast.vector = data.frame,
                                  r.squared = 0.50, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                 "`contrast.vector` must be either a vector or a matrix.")
    expect_error(power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5, 2), contrast.vector = c(0.5, 0.5, 1),
                                  r.squared = 0.50, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                 "The number of columns / elements in the contrast vector should match number of groups.")
    expect_error(power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5, 2), contrast.vector = diag(2),
                                  r.squared = 0.50, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                 "The number of rows in the contrast matrix should be one.")
    expect_error(power.t.contrast(mu.vector = c(0.001, 0, 0.0005), sd.vector = rep(1, 3), p.vector = rep(1 / 3, 3), contrast.vector = c(1, 0, -1),
                                  r.squared = 0, k.covariates = 1e8, power = 0.99, alpha = 1e-9, verbose = 0),
                 "Design not feasible.")
    expect_warning(power.t.contrast(mu.vector = c(0.20, 0.20), sd.vector = rep(1, 2), p.vector = rep(0.5, 2), contrast.vector = c(1, -1),
                                    r.squared = 0.5, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                   "Using infinity \\(maximum integer number as defined in R\\) for `n.total` because `psi` = 0.")
    expect_error(power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), contrast.vector = c(1, -1),
                                  r.squared = 0.5, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                   "The `p.vector` can not be NULL when the sample size is requested.")
    expect_error(power.t.contrast(mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5 - 1e-5, 2), contrast.vector = c(1, -1),
                                  r.squared = 0.5, k.covariates = 1, power = 0.8, alpha = 0.05, verbose = 0),
                   "The elements of the `p.vector` should sum to 1.")
    expect_error(power.t.contrasts(data.frame),
                 "This function only works with an object of type `pwrss`, `f`, `ancova`, and `shieh`.")
})

test_that("ANCOVA helper functions work", {
    expect_equal(calc.df2(100, 12, 2), 86)
    expect_equal(calc.df2(100, 12, 0), 88)
    expect_equal(calc.df2(100,  0, 2), 98)
    expect_equal(fmt_test_ancova(1, 0), "One-Way Analysis of Variance (F-Test)")
    expect_equal(fmt_test_ancova(2, 1), "Two-Way Analysis of Covariance (F-Test)")
    expect_equal(fmt_test_ancova(3, 0), "Three-Way Analysis of Variance (F-Test)")
    expect_warning(check_var.ratio(c(0.2, 1, 2), rep(20, 3)),
                   "Interpretation of results may no longer be valid when variances differ beyond sampling error.")
    expect_null(check_var.ratio(c(1.1, 1, 0.9), rep(20, 3)))
})
