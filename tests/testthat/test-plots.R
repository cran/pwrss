test_that("plots work", {

    # uses check.snap4plot in helper.R =================================================================================
    # to create snaps for the first time, use `crtSnp <- TRUE` before running `devtools::test()`

    # binomial-tests ---------------------------------------------------------------------------------------------------
    check.snap4plot("binom_plot_1", power.binom.test, list(size = 200, prob = 0.6, null.prob = 0.5,         alternative = "one.sided"))
    check.snap4plot("binom_plot_2", power.binom.test, list(size = 200, prob = 0.4, null.prob = 0.5,         alternative = "one.sided"))
    check.snap4plot("binom_plot_3", power.binom.test, list(size = 200, prob = 0.4, null.prob = 0.5,         alternative = "two.sided"))
    check.snap4plot("binom_plot_4", power.binom.test, list(size = 200, prob = 0.5, null.prob = c(0.4, 0.6), alternative = "two.one.sided"))
    check.snap4plot("binom_plot_5", power.binom.test, list(size = 200, prob = 0.5, null.prob = c(0.3, 0.4), alternative = "two.one.sided"))

    # chisq-tests ------------------------------------------------------------------------------------------------------
    check.snap4plot("chisq_plot_1", power.chisq.test, list(ncp = 20, df =    1))
    check.snap4plot("chisq_plot_2", power.chisq.test, list(ncp = 20, df =    2))
    check.snap4plot("chisq_plot_3", power.chisq.test, list(ncp = 20, df =    4))
    check.snap4plot("chisq_plot_4", power.chisq.test, list(ncp = 20, df =   10))
    check.snap4plot("chisq_plot_5", power.chisq.test, list(ncp = 20, df =  100))
    check.snap4plot("chisq_plot_6", power.chisq.test, list(ncp = 20, df = 1000))
    check.snap4plot("chisq_plot_7", power.chisq.test, list(ncp = 40, df =  100))

    # F-tests ----------------------------------------------------------------------------------------------------------
    check.snap4plot("F_plot_01", power.f.test, list(ncp = 1, df1 = 1, df2 =   10))
    check.snap4plot("F_plot_02", power.f.test, list(ncp = 1, df1 = 1, df2 =  100))
    check.snap4plot("F_plot_03", power.f.test, list(ncp = 1, df1 = 1, df2 = 1000))
    check.snap4plot("F_plot_04", power.f.test, list(ncp = 2, df1 = 1, df2 =   10))
    check.snap4plot("F_plot_05", power.f.test, list(ncp = 2, df1 = 1, df2 =  100))
    check.snap4plot("F_plot_06", power.f.test, list(ncp = 2, df1 = 1, df2 = 1000))
    check.snap4plot("F_plot_07", power.f.test, list(ncp = 4, df1 = 1, df2 =   10))
    check.snap4plot("F_plot_08", power.f.test, list(ncp = 4, df1 = 1, df2 =  100))
    check.snap4plot("F_plot_09", power.f.test, list(ncp = 4, df1 = 1, df2 = 1000))
    check.snap4plot("F_plot_10", power.f.test, list(ncp = 8, df1 = 1, df2 =   10))
    check.snap4plot("F_plot_11", power.f.test, list(ncp = 8, df1 = 1, df2 =  100))
    check.snap4plot("F_plot_12", power.f.test, list(ncp = 8, df1 = 1, df2 = 1000))
    check.snap4plot("F_plot_13", power.f.test, list(ncp = 4, df1 = 2, df2 =   10))
    check.snap4plot("F_plot_14", power.f.test, list(ncp = 4, df1 = 2, df2 =  100))
    check.snap4plot("F_plot_15", power.f.test, list(ncp = 4, df1 = 2, df2 = 1000))
    check.snap4plot("F_plot_16", power.f.test, list(ncp = 8, df1 = 4, df2 =   10))
    check.snap4plot("F_plot_17", power.f.test, list(ncp = 8, df1 = 4, df2 =  100))
    check.snap4plot("F_plot_18", power.f.test, list(ncp = 8, df1 = 4, df2 = 1000))

    # Lambda-prime distribution ----------------------------------------------------------------------------------------
    check.snap4plot("lp_plot_1", power.lp.test, list(ncp =  1.96, df = 100))
    check.snap4plot("lp_plot_2", power.lp.test, list(ncp =  1.96, df = 100, alternative = "one.sided"))
    check.snap4plot("lp_plot_3", power.lp.test, list(ncp =  0, null.ncp = c(-2, 2), df = 100, alternative = "two.one.sided"))
    check.snap4plot("lp_plot_4", power.lp.test, list(ncp =  2, null.ncp = c(-1, 1), df = 100, alternative = "two.one.sided"))
    check.snap4plot("lp_plot_5", power.lp.test, list(ncp = -2, null.ncp = c(-1, 1), df = 100, alternative = "two.one.sided"))
    check.snap4plot("lp_plot_6", power.lp.test, list(ncp = -1.96, df = 100))

    # t-tests ----------------------------------------------------------------------------------------------------------
    check.snap4plot("t_plot_1", power.t.test, list(ncp =  1.96,                   df = 100, alternative = "two.sided"))
    check.snap4plot("t_plot_2", power.t.test, list(ncp = -1.96,                   df = 100, alternative = "two.sided"))
    check.snap4plot("t_plot_3", power.t.test, list(ncp =  1.96,                   df = 100, alternative = "one.sided"))
    check.snap4plot("t_plot_4", power.t.test, list(ncp =  0, null.ncp = c(-2, 2), df = 100, alternative = "two.one.sided"))
    check.snap4plot("t_plot_5", power.t.test, list(ncp =  2, null.ncp = c(-1, 1), df = 100, alternative = "two.one.sided"))
    check.snap4plot("t_plot_6", power.t.test, list(ncp = -2, null.ncp = c(-1, 1), df = 100, alternative = "two.one.sided"))

    # z-tests ----------------------------------------------------------------------------------------------------------
    check.snap4plot("z_plot_1", power.z.test, list(mean =  1.96,                    alternative = "two.sided"))
    check.snap4plot("z_plot_2", power.z.test, list(mean = -1.96,                    alternative = "two.sided"))
    check.snap4plot("z_plot_3", power.z.test, list(mean =  1.96,                    alternative = "one.sided"))
    check.snap4plot("z_plot_4", power.z.test, list(mean =  0, null.mean = c(-2, 2), alternative = "two.one.sided"))
    check.snap4plot("z_plot_5", power.z.test, list(mean =  2, null.mean = c(-1, 1), alternative = "two.one.sided"))

    # plots.R (taking a results-obj and using it for plotting) ---------------------------------------------------------
    # binom: proportions.onetwo (only exact.oneprop) -------------------------------------------------------------------
    check.snap4plot("plot_oneprop.binom", power.exact.oneprop, list(power = 0.8, prob = 0.45, null.prob = 0.50, alternative = "one.sided"))
    # chisq: chisq.gof -------------------------------------------------------------------------------------------------
    mtxW <- probs.to.w(cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078), c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101)), verbose = 0)
    check.snap4plot("plot_gof", power.chisq.gof, list(power = 0.8, w = mtxW$w, df = mtxW$df))
    # F: ancova, keppel, shieh, mixed.anova, (f.)regression ------------------------------------------------------------
    check.snap4plot("plot_ancova",        power.f.ancova,     list(power = 0.8, eta.squared = 0.059, factor.levels = 2))
    check.snap4plot("plot_ancova.keppel", power.f.ancova.keppel,
                    list(power = 0.8, mu.vector = c(0.50, 0.00), sd.vector = rep(1, 2), p.vector = rep(0.50, 2),
                         k.covariates = 1, r.squared = 0.50))
    check.snap4plot("plot_ancova.shieh",  power.f.ancova.shieh,
                    list(power = 0.8, mu.vector = c(0.20, 0), sd.vector = rep(1, 2), p.vector = rep(0.5, 2),
                         r.squared = 0.50, k.covariates = 1))
    check.snap4plot("plot_mixed.anova",   power.f.mixed.anova,
                    list(power = 0.8, eta.squared = 0.022, factor.levels = c(1, 2), rho.within = 0.50, effect = "within"))
    check.snap4plot("plot_regression.f",  power.f.regression, list(power = 0.8, r.squared = 0.15, k.total = 3))

    # t: student, welch, wilcoxon, (t.)regression ----------------------------------------------------------------------
    check.snap4plot("plot_student",      power.t.student,    list(d = 0.20, n2 = 100, alternative = "two.sided", design = "independent"))
    check.snap4plot("plot_welch",        power.t.welch,      list(d = 0.20, n2 = 100, n.ratio = 2, var.ratio = 2, alternative = "two.sided"))
    check.snap4plot("plot_wilcoxon",     power.np.wilcoxon,  list(d = 0.25, n2 = 100, alternative = "two.sided", design = "independent"))
    check.snap4plot("plot_regression.t", power.t.regression, list(power = 0.8, beta = 0.20, k.total = 5, r.squared = 0.30))

    # z: proportions.onetwo, correlations (steiger, twocors, onecor), logistic, poisson, mediation ---------------------
    check.snap4plot("plot_oneprop.z",  power.z.oneprop,  list(power = 0.8, prob = 0.45, null.prob = 0.50, alternative = "one.sided"))
    check.snap4plot("plot_twoprops.z", power.z.twoprops, list(power = 0.8, prob1 = 0.65, prob2 = 0.60, alternative = "one.sided"))
    check.snap4plot("plot_twocors.steiger", power.z.twocors.steiger,
                    list(power = 0.8, rho12 = 0.35, rho13 = 0.45, rho23 = 0.05, alternative = "two.sided", common.index = TRUE))
    check.snap4plot("plot_twocors",   power.z.twocors,   list(power = 0.8, rho1 = 0.20, rho2 = 0.30, alternative = "two.sided"))
    check.snap4plot("plot_onecor",    power.z.onecor,    list(power = 0.8, rho = 0.20, alternative = "two.sided"))
    check.snap4plot("plot_logistic",  power.z.logistic,  list(power = 0.8, base.prob = 0.15, prob = 0.20, distribution = "normal"))
    check.snap4plot("plot_mediation", power.z.mediation, list(power = 0.8, beta.a = 0.25, beta.b = 0.25, beta.cp = 0.10))
    check.snap4plot("plot_poisson",   power.z.poisson,   list(power = 0.8, beta0 = 0.50, beta1 = -0.10, dist = "normal"))

    # errors -----------------------------------------------------------------------------------------------------------
    expect_error(plot(pwrss.z.mediation(power = 0.8, a = 0.25, b = 0.25, cp = 0.10, verbose = FALSE)),
                 "Plotting is no longer available for this type of object.")
    expect_error(plot(power.exact.fisher(power = 0.8, prob1 = 0.60, prob2 = 0.40, verbose = 0)),
                 "Plotting is not available for exact tests.")
    expect_error(plot.pwrss(NULL), "Not an object of the type 'pwrss'.")
    expect_error(.plot.binom.t1t2(size = 5, prob = 0.4, alternative = "two.sided"),
                 "Number of trials should be greater than 10 for plotting.")
    expect_error(.plot.chisq.t1t2(ncp = 40, df = Inf),
                 "`df` must be numeric, finite, have a value of at least 1 and have a length of 1.")
    expect_error(.plot.f.t1t2(ncp = 1, df1 = 0, df2 = 10),
                 "`df1` must be numeric, finite, have a value of at least 1 and have a length of 1.")
    expect_error(.plot.f.t1t2(ncp = 1, df1 = 1, df2 = Inf),
                 "`df2` must be numeric, finite, have a value of at least 1 and have a length of 1.")
    expect_error(.plot.lp.t1t2(ncp = 1, df = NA),
                 "`df` must be numeric, have a value of at least 1 and have a length of 1.")
    expect_error(.plot.t.t1t2(ncp = 2, df = 0),
                 "`df` must be numeric, have a value of at least 1 and have a length of 1.")
})
