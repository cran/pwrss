# power.f.mixed.anova (= pwrss.f.rmanova) ------------------------------------------------------------------------------
test_that("power.f.mixed.anova / pwrss.f.rmanova work", {
    crrRes <- power.f.mixed.anova(eta.squared = 0.022, factor.levels = c(1, 2), rho.within = 0.50, effect = "within",
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.022, null.eta.squared = 0, factor.levels = c(1, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = NULL,
                          power = 0.80, alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(1)", df1 = 1, df2 = 89, ncp = 8.09815951, null.ncp = 0,
                      f.alpha = 3.94808435, eta.squared = 0.022, power = 0.803705778, n.total = 90))
    expect_equal(crrRes, pwrss.f.rmanova(eta2 = 0.022, n.levels = 1, n.rm = 2, corr.rm = 0.50, type = "within",
                                         power = 0.80, alpha = 0.05, verbose = 0))
    expect_equal(crrRes, pwrss.f.rmanova(f2 = 0.0224948875, n.levels = 1, n.rm = 2, corr.rm = 0.50, type = "within",
                                         power = 0.80, alpha = 0.05, verbose = 0))

    crrRes <- power.f.mixed.anova(eta.squared = 0.022, factor.levels = c(1, 2), rho.within = 0.50, effect = "within",
                                  n.total = 120, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.022, null.eta.squared = 0, factor.levels = c(1, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 120, power = NULL,
                          alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(1)", df1 = 1, df2 = 119, ncp = 10.797546, null.ncp = 0,
                      f.alpha = 3.9207955, eta.squared = 0.022, power = 0.9030908, n.total = 120))
    expect_equal(crrRes, pwrss.f.rmanova(eta2 = 0.022, n.levels = 1, n.rm = 2, corr.rm = 0.50, type = "within",
                                         n = 120, alpha = 0.05, verbose = 0))
    expect_equal(crrRes, pwrss.f.rmanova(f2 = 0.0224948875, n.levels = 1, n.rm = 2, corr.rm = 0.50, type = "within",
                                         n = 120, alpha = 0.05, verbose = 0))

    crrRes <- power.f.mixed.anova(factor.levels = c(1, 2), rho.within = 0.50, effect = "within", n.total = 90,
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = c(1, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 90, power = 0.80,
                          alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(1)", df1 = 1, df2 = 89, ncp = 8.02174, null.ncp = 0,
                      f.alpha = 3.94808435, eta.squared = 0.021796919, power = 0.8, n.total = 90))

    crrRes <- suppressWarnings(power.f.mixed.anova(eta.squared = 0.08255, factor.levels = c(1, 2), rho.within = NA,
                                                   effect = "within", power = 0.80, alpha = 0.05, verbose = 0))
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.08255, null.eta.squared = 0, factor.levels = c(1, 2),
                          factor.type = c("between", "within"), rho.within = NA, epsilon = 1, n.total = NULL,
                          power = 0.80, alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(1)", df1 = 1, df2 = 89, ncp = 8.097989, null.ncp = 0,
                      f.alpha = 3.94808435, eta.squared = 0.08255, power = 0.803697575, n.total = 90))

    crrRes <- suppressWarnings(power.f.mixed.anova(eta.squared = 0.08255, factor.levels = c(1, 2), rho.within = NA,
                                                   effect = "within", n.total = 90, alpha = 0.05, verbose = 0))
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.08255, null.eta.squared = 0, factor.levels = c(1, 2),
                          factor.type = c("between", "within"), rho.within = NA, epsilon = 1, n.total = 90, power = NULL,
                          alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(1)", df1 = 1, df2 = 89, ncp = 8.097989, null.ncp = 0,
                      f.alpha = 3.94808435, eta.squared = 0.08255, power = 0.803697575, n.total = 90))

    crrRes <- power.f.mixed.anova(eta.squared = 0.059, factor.levels = c(2, 1), effect = "between", alpha = 0.05,
                                  power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.059, null.eta.squared = 0, factor.levels = c(2, 1),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = NULL,
                          power = 0.80, alpha = 0.05, effect = "between", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2)|W(1)", df1 = 1, df2 = 126, ncp = 8.0255048, null.ncp = 0,
                      f.alpha = 3.9163246, eta.squared = 0.059, power = 0.8027032, n.total = 128))

    crrRes <- power.f.mixed.anova(eta.squared = 0.059, factor.levels = c(2, 1), effect = "between", alpha = 0.05,
                                  n.total = 128, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.059, null.eta.squared = 0, factor.levels = c(2, 1),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 128,
                          power = NULL, alpha = 0.05, effect = "between", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2)|W(1)", df1 = 1, df2 = 126, ncp = 8.0255048, null.ncp = 0,
                      f.alpha = 3.9163246, eta.squared = 0.059, power = 0.8027032, n.total = 128))

    crrRes <- power.f.mixed.anova(factor.levels = c(2, 1), effect = "between", alpha = 0.05, n.total = 128, power = 0.80, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = c(2, 1),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 128,
                          power = 0.80, alpha = 0.05, effect = "between", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2)|W(1)", df1 = 1, df2 = 126, ncp = 7.9702414, null.ncp = 0,
                      f.alpha = 3.9163246, eta.squared = 0.058617543, power = 0.8, n.total = 128))

    crrRes <- power.f.mixed.anova(eta.squared = 0.038, factor.levels = c(2, 2),  rho.within = 0.50, effect = "between",
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.038, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = NULL,
                          power = 0.80, alpha = 0.05, effect = "between", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2)|W(2)", df1 = 1, df2 = 150, ncp = 8.005544, null.ncp = 0,
                      f.alpha = 3.9042019, eta.squared = 0.038, power = 0.802695047, n.total = 152))

    crrRes <- power.f.mixed.anova(eta.squared = 0.038, factor.levels = c(2, 2), rho.within = 0.50, effect = "between",
                                  n.total = 152, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.038, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 152,
                          power = NULL, alpha = 0.05, effect = "between", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2)|W(2)", df1 = 1, df2 = 150, ncp = 8.005544, null.ncp = 0,
                      f.alpha = 3.9042019, eta.squared = 0.038, power = 0.802695047, n.total = 152))

    crrRes <- power.f.mixed.anova(factor.levels = c(2, 2), rho.within = 0.50, effect = "between", n.total = 152,
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 152,
                          power = 0.80, alpha = 0.05, effect = "between", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2)|W(2)", df1 = 1, df2 = 150, ncp = 7.9505852, null.ncp = 0,
                      f.alpha = 3.9042019, eta.squared = 0.037748974, power = 0.8, n.total = 152))

    crrRes <- power.f.mixed.anova(eta.squared = 0.01, factor.levels = c(2, 2), rho.within = 0.50, effect = "within",
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.01, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = NULL,
                          power = 0.80, alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(2)", df1 = 1, df2 = 196, ncp = 8.00, null.ncp = 0,
                      f.alpha = 3.889341, eta.squared = 0.01, power = 0.8036086, n.total = 198))

    crrRes <- power.f.mixed.anova(eta.squared = 0.01, factor.levels = c(2, 2), rho.within = 0.50, effect = "within",
                                  n.total = 198, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.01, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 198,
                          power = NULL, alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(2)", df1 = 1, df2 = 196, ncp = 8.00, null.ncp = 0,
                      f.alpha = 3.889341, eta.squared = 0.01, power = 0.8036086, n.total = 198))

    crrRes <- power.f.mixed.anova(factor.levels = c(2, 2), rho.within = 0.50, effect = "within", n.total = 198,
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 198,
                          power = 0.80, alpha = 0.05, effect = "within", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "W(2)|B(2)", df1 = 1, df2 = 196, ncp = 7.92649015, null.ncp = 0,
                      f.alpha = 3.889341, eta.squared = 0.0099090232, power = 0.8, n.total = 198))

    crrRes <- power.f.mixed.anova(eta.squared = 0.01, factor.levels = c(2, 2), rho.within = 0.50, effect = "interaction",
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.01, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = NULL,
                          power = 0.80, alpha = 0.05, effect = "interaction", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2):W(2)", df1 = 1, df2 = 196, ncp = 8.00, null.ncp = 0,
                      f.alpha = 3.889341, eta.squared = 0.01, power = 0.8036086, n.total = 198))

    crrRes <- power.f.mixed.anova(eta.squared = 0.01, factor.levels = c(2, 2), rho.within = 0.50, effect = "interaction",
                                  n.total = 198, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = 0.01, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 198,
                          power = NULL, alpha = 0.05, effect = "interaction", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2):W(2)", df1 = 1, df2 = 196, ncp = 8.00, null.ncp = 0,
                      f.alpha = 3.889341, eta.squared = 0.01, power = 0.8036086, n.total = 198))

    crrRes <- power.f.mixed.anova(factor.levels = c(2, 2), rho.within = 0.50, effect = "interaction", n.total = 198,
                                  power = 0.80, alpha = 0.05, verbose = 0)
    expect_s3_class(crrRes, c("pwrss", "f", "anova_mixed"))
    expect_named(crrRes, c("parms", "test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total"))
    expect_identical(crrRes[["parms"]],
                     list(eta.squared = NULL, null.eta.squared = 0, factor.levels = c(2, 2),
                          factor.type = c("between", "within"), rho.within = 0.5, epsilon = 1, n.total = 198,
                          power = 0.80, alpha = 0.05, effect = "interaction", ceil.n = TRUE, verbose = 0, utf = FALSE))
    expect_equal(crrRes[c("test", "effect", "df1", "df2", "ncp", "null.ncp", "f.alpha", "eta.squared", "power", "n.total")],
                 list(test = "F", effect = "B(2):W(2)", df1 = 1, df2 = 196, ncp = 7.92649015, null.ncp = 0,
                      f.alpha = 3.889341, eta.squared = 0.0099090232, power = 0.8, n.total = 198))

    expect_error(power.f.mixed.anova(eta.squared = 0.022, factor.levels = c(2, 2), factor.type = c("within", "within"),
                                     rho.within = 0.50, effect = "within", power = 0.80, alpha = 0.05, verbose = 0),
                 paste("The `factor.type` argument must be specified as either c\\('between', 'within'\\) or",
                       "c\\('within', 'between'\\), indicating the order in which the corresponding values in",
                       "`factor.levels` are interpreted - specifically, which factor is treated as between-subjects",
                       "and which as within-subjects."))
    expect_error(power.f.mixed.anova(eta.squared = 0.022, factor.levels = c(2, 2), factor.type = c("within", "wrong"),
                                     rho.within = 0.50, effect = "within", power = 0.80, alpha = 0.05, verbose = 0),
                 paste("The `factor.type` argument must be specified as either c\\('between', 'within'\\) or",
                       "c\\('within', 'between'\\), indicating the order in which the corresponding values in",
                       "`factor.levels` are interpreted - specifically, which factor is treated as between-subjects",
                       "and which as within-subjects."))
    expect_error(power.f.mixed.anova(eta.squared = 0.022, factor.levels = c(2, 2, 2), factor.type = c("within", "between"),
                                     rho.within = 0.50, effect = "within", power = 0.80, alpha = 0.05, verbose = 0),
                 "Exactly two factors are allowed in this procedure.")
    expect_error(power.f.mixed.anova(eta.squared = 0.022, factor.levels = c(2, 2), rho.within = 0.50, effect = "within",
                                     epsilon = 0.1, power = 0.80, alpha = 0.05, verbose = 0),
                 "Incorrect value for the non-sphericity correction factor \\(`epsilon`\\).")
    expect_error(power.f.mixed.anova(factor.levels = c(1, 2), rho.within = NA, effect = "within", n.total = 90,
                                     power = 0.80, alpha = 0.05, verbose = 0),
                 "When the effect size shall be estimated, `rho.within` needs to be defined.")
    expect_error(power.f.mixed.anova(eta.squared = 1e-12, factor.levels = c(4, 4), rho.within = 0.50, effect = "within",
                                     power = 0.99, alpha = 0.05, verbose = 0),
                 "Design is not feasible.")
    expect_error(power.f.mixed.anova(factor.levels = c(4, 4), rho.within = 0.50, effect = "within", power = 0.80,
                                     n.total = 2, alpha = 0.05, verbose = 0),
                 "Design is not feasible.")
    expect_error(power.f.mixed.anova(eta.squared = 0.022, factor.levels = c(1, 1), rho.within = 0.50, effect = "within",
                                     power = 0.80, alpha = 0.05, verbose = 0),
                 "Design is not feasible.")
    expect_error(pwrss.f.rmanova(eta2 = 0.022, f2 = 0.0224948875, n.levels = 1, n.rm = 2, type = "within", n = 120, verbose = 0),
                 "Effect size conflict for the alternative. Specify only either `eta2` or `f2`.")
})
