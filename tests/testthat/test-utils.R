# inflate.sample -------------------------------------------------------------------------------------------------------
test_that("inflate.sample works", {
    expect_equal(vapply(2 ^ seq(10), inflate.sample, numeric(1), rate = 0.05, verbose = 0), c(3, 5,  9, 17, 34, 68, 135, 270, 539, 1078))
    expect_equal(vapply(2 ^ seq(10), inflate.sample, numeric(1), rate = 0.20, verbose = 0), c(3, 5, 10, 20, 40, 80, 160, 320, 640, 1280))
    expect_equal(vapply(2 ^ seq(10), inflate.sample, numeric(1), rate = 0.20, ceil.n = FALSE, verbose = 0),
                 c(2.5, 5, 10, 20, 40, 80, 160, 320, 640, 1280))
    expect_equal(capture.output(invisible(vapply(2 ^ seq(10), inflate.sample, numeric(1), rate = 0.20))), "35102040801603206401280")
})

# etasq.to.f -----------------------------------------------------------------------------------------------------------
test_that("etasq.to.f works", {
    expect_equal(etasq.to.f(0.009900990, verbose = 0), list(f.squared = 0.0100, f = 0.10, eta.squared = 0.009900990))
    expect_equal(etasq.to.f(0.058823530, verbose = 0), list(f.squared = 0.0625, f = 0.25, eta.squared = 0.058823530))
    expect_equal(etasq.to.f(0.137931034, verbose = 0), list(f.squared = 0.1600, f = 0.40, eta.squared = 0.137931034))
    expect_equal(capture.output(etasq.to.f(0.009900990)), c("  f.squared           f eta.squared ", " 0.01000000  0.10000000  0.00990099 "))
})

# f.to.etasq -----------------------------------------------------------------------------------------------------------
test_that("f.to.etasq works", {
    expect_equal(f.to.etasq(0.10, verbose = 0), list(eta.squared = 0.009900990, f.squared = 0.0100, f = 0.10))
    expect_equal(f.to.etasq(0.25, verbose = 0), list(eta.squared = 0.058823530, f.squared = 0.0625, f = 0.25))
    expect_equal(f.to.etasq(0.40, verbose = 0), list(eta.squared = 0.137931034, f.squared = 0.1600, f = 0.40))
    expect_equal(capture.output(f.to.etasq(0.10)), c("eta.squared   f.squared           f ", " 0.00990099  0.01000000  0.10000000 "))
})

# rsq.to.f -------------------------------------------------------------------------------------------------------------
test_that("rsq.to.f works", {
    expect_equal(rsq.to.f(0.009900990),
                 list(f.squared = 0.0100, f = 0.10, r.squared.full = 0.009900990, r.squared.reduced = 0))
    expect_equal(rsq.to.f(0.058823530),
                 list(f.squared = 0.0625, f = 0.25, r.squared.full = 0.058823530, r.squared.reduced = 0))
    expect_equal(rsq.to.f(0.137931034),
                 list(f.squared = 0.1600, f = 0.40, r.squared.full = 0.137931034, r.squared.reduced = 0))
    expect_error(rsq.to.f(0.1, 0.2), "Expecting `r.squared.full` > `r.squared.reduced`.")
    expect_equal(capture.output(rsq.to.f(0.137931034, verbose = 1)),
                 c("        f.squared                 f    r.squared.full r.squared.reduced ",
                   "         0.160000          0.400000          0.137931          0.000000 "))
})

# f.to.rsq -------------------------------------------------------------------------------------------------------------
test_that("f.to.rsq works", {
    expect_equal(f.to.rsq(0.10),
                 list(f.squared = 0.0100, f = 0.10, r.squared.full = 0.009900990, r.squared.reduced = 0))
    expect_equal(f.to.rsq(0.25),
                 list(f.squared = 0.0625, f = 0.25, r.squared.full = 0.058823530, r.squared.reduced = 0))
    expect_equal(f.to.rsq(0.40),
                 list(f.squared = 0.1600, f = 0.40, r.squared.full = 0.137931034, r.squared.reduced = 0))
    expect_equal(f.to.rsq(0.40, r.squared.full = 0.00862069, verbose = 0),
                 list(f.squared = 0.1600, f = 0.40, r.squared.full = 0.00862069, r.squared.reduced = -0.15))
    expect_equal(capture.output(f.to.rsq(0.40, verbose = 1)),
                 c("        f.squared                 f    r.squared.full r.squared.reduced ",
                   "         0.160000          0.400000          0.137931          0.000000 "))
})

# cor.to.z -------------------------------------------------------------------------------------------------------------
test_that("cor.to.z works", {
    expect_equal(cor.to.z(1.000,  verbose = 0), list(z =  Inf,         rho =  1.000))
    expect_equal(cor.to.z(-1.000, verbose = 0), list(z = -Inf,         rho = -1.000))
    expect_equal(cor.to.z(0.999,  verbose = 0), list(z =  3.800201170, rho =  0.999))
    expect_equal(cor.to.z(-0.999, verbose = 0), list(z = -3.800201170, rho = -0.999))
    expect_equal(cor.to.z(0.990,  verbose = 0), list(z =  2.646652410, rho =  0.990))
    expect_equal(cor.to.z(-0.990, verbose = 0), list(z = -2.646652410, rho = -0.990))
    expect_equal(cor.to.z(0.950,  verbose = 0), list(z =  1.831780820, rho =  0.950))
    expect_equal(cor.to.z(-0.950, verbose = 0), list(z = -1.831780820, rho = -0.950))
    expect_equal(cor.to.z(0.900,  verbose = 0), list(z =  1.472219490, rho =  0.900))
    expect_equal(cor.to.z(-0.900, verbose = 0), list(z = -1.472219490, rho = -0.900))
    expect_equal(cor.to.z(0.500,  verbose = 0), list(z =  0.549306144, rho =  0.500))
    expect_equal(cor.to.z(-0.500, verbose = 0), list(z = -0.549306144, rho = -0.500))
    expect_equal(cor.to.z(0.200,  verbose = 0), list(z =  0.202732554, rho =  0.200))
    expect_equal(cor.to.z(-0.200, verbose = 0), list(z = -0.202732554, rho = -0.200))
    expect_equal(cor.to.z(0.100,  verbose = 0), list(z =  0.100335348, rho =  0.100))
    expect_equal(cor.to.z(-0.100, verbose = 0), list(z = -0.100335348, rho = -0.100))
    expect_equal(cor.to.z(0.000,  verbose = 0), list(z =  0.000000000, rho =  0.000))
    expect_equal(capture.output(cor.to.z(1.000)), c("  z rho ", "Inf   1 "))
})

# z.to.cor -------------------------------------------------------------------------------------------------------------
test_that("z.to.cor works", {
    expect_equal(z.to.cor(0.1,  verbose = 0), list(rho =  0.099667995, z =  0.1))
    expect_equal(z.to.cor(-0.1, verbose = 0), list(rho = -0.099667995, z = -0.1))
    expect_equal(z.to.cor(1.0,  verbose = 0), list(rho =  0.761594160, z =  1.0))
    expect_equal(z.to.cor(-1.0, verbose = 0), list(rho = -0.761594160, z = -1.0))
    expect_equal(z.to.cor(2.0,  verbose = 0), list(rho =  0.964027580, z =  2.0))
    expect_equal(z.to.cor(-2.0, verbose = 0), list(rho = -0.964027580, z = -2.0))
    expect_equal(z.to.cor(4.0,  verbose = 0), list(rho =  0.999329300, z =  4.0))
    expect_equal(z.to.cor(-4.0, verbose = 0), list(rho = -0.999329300, z = -4.0))
    expect_equal(capture.output(z.to.cor(0.1)), c("       rho          z ", "0.09966799 0.10000000 "))
    expect_error(z.to.cor(Inf, verbose = 0), "All elements of `z` need to be valid numeric values \\(numeric, and finite\\)")
})

# cors.to.q ------------------------------------------------------------------------------------------------------------
test_that("cors.to.q works", {
    expect_equal(cors.to.q(rho2 = 0.571202682, rho1 = 0.5, verbose = 0),
                 list(q = -0.1, delta = -0.071202682, rho1 = 0.5, rho2 = 0.571202682))
    expect_equal(cors.to.q(rho2 = 0.690706810, rho1 = 0.5, verbose = 0),
                 list(q = -0.3, delta = -0.190706810, rho1 = 0.5, rho2 = 0.690706810))
    expect_equal(cors.to.q(rho2 = 0.781536455, rho1 = 0.5, verbose = 0),
                 list(q = -0.5, delta = -0.281536455, rho1 = 0.5, rho2 = 0.781536455))
    expect_equal(capture.output(cors.to.q(rho2 = 0.571202682, rho1 = 0.5)),
                 c("          q       delta        rho1        rho2 ", "-0.10000000 -0.07120268  0.50000000  0.57120268 "))
})

# q.to.cors ------------------------------------------------------------------------------------------------------------
test_that("q.to.cors works", {
    expect_equal(q.to.cors(q = 0.10, rho1 = 0.5, verbose = 0),
                 list(q =  0.1, delta = -0.071202682, rho1 = 0.5, rho2 = 0.571202682))
    expect_equal(q.to.cors(q = 0.30, rho1 = 0.5, verbose = 0),
                 list(q =  0.3, delta = -0.190706810, rho1 = 0.5, rho2 = 0.690706810))
    expect_equal(q.to.cors(q = 0.50, rho1 = 0.5, verbose = 0),
                 list(q =  0.5, delta = -0.281536455, rho1 = 0.5, rho2 = 0.781536455))
    expect_equal(q.to.cors(q = 0.10, rho2 = 0.5, verbose = 0),
                 list(q =  0.1, delta = -0.078671512, rho1 = 0.421328488, rho2 = 0.5))
    expect_equal(q.to.cors(q = 0.30, rho2 = 0.5, verbose = 0),
                 list(q =  0.3, delta = -0.255733683, rho1 = 0.244266317, rho2 = 0.5))
    expect_equal(q.to.cors(q = 0.50, rho2 = 0.5, verbose = 0),
                 list(q =  0.5, delta = -0.450733773, rho1 = 0.049266227, rho2 = 0.5))
    expect_equal(capture.output(q.to.cors(q = 0.10, rho1 = 0.5)),
                 c("          q       delta        rho1        rho2 ", " 0.10000000 -0.07120268  0.50000000  0.57120268 "))
    expect_equal(capture.output(q.to.cors(q = 0.10, rho2 = 0.5)),
                 c("          q       delta        rho1        rho2 ", " 0.10000000 -0.07867151  0.42132849  0.50000000 "))
    expect_error(q.to.cors(q = 0.10, verbose = 0), "Both `rho1` and `rho2` can not be NULL.")
    expect_error(q.to.cors(q = 0.10, rho1 = 0.5, rho2 = 0.3, verbose = 0), "Exactly one of the `rho1` or `rho2` should be NULL.")
})

# d.to.cles ------------------------------------------------------------------------------------------------------------
test_that("d.to.cles works", {
    expect_equal(d.to.cles(0.2,                        verbose = 0), list(cles = 0.556231458, d = 0.2))
    expect_equal(d.to.cles(0.5,                        verbose = 0), list(cles = 0.638163195, d = 0.5))
    expect_equal(d.to.cles(0.8,                        verbose = 0), list(cles = 0.714196178, d = 0.8))
    expect_equal(d.to.cles(0.2, design = "paired",     verbose = 0), list(cles = 0.579259710, d = 0.2))
    expect_equal(d.to.cles(0.5, design = "paired",     verbose = 0), list(cles = 0.691462461, d = 0.5))
    expect_equal(d.to.cles(0.8, design = "paired",     verbose = 0), list(cles = 0.788144600, d = 0.8))
    expect_equal(d.to.cles(0.2, design = "one.sample", verbose = 0), list(cles = 0.579259710, d = 0.2))
    expect_equal(d.to.cles(0.5, design = "one.sample", verbose = 0), list(cles = 0.691462461, d = 0.5))
    expect_equal(d.to.cles(0.8, design = "one.sample", verbose = 0), list(cles = 0.788144600, d = 0.8))
    expect_equal(capture.output(d.to.cles(0.2, verbose = 1)), c("     cles         d ", "0.5562315 0.2000000 "))
    expect_equal(capture.output(d.to.cles(0.5, verbose = 1)), c("     cles         d ", "0.6381632 0.5000000 "))
    expect_equal(capture.output(d.to.cles(0.8, verbose = 1)), c("     cles         d ", "0.7141962 0.8000000 "))
})

# cles.to.d ------------------------------------------------------------------------------------------------------------
test_that("cles.to.d works", {
    expect_equal(cles.to.d(0.556231458,                        verbose = 0), list(d = 0.2, cles = 0.556231458))
    expect_equal(cles.to.d(0.638163195,                        verbose = 0), list(d = 0.5, cles = 0.638163195))
    expect_equal(cles.to.d(0.714196178,                        verbose = 0), list(d = 0.8, cles = 0.714196178))
    expect_equal(cles.to.d(0.579259710, design = "paired",     verbose = 0), list(d = 0.2, cles = 0.579259710))
    expect_equal(cles.to.d(0.691462461, design = "paired",     verbose = 0), list(d = 0.5, cles = 0.691462461))
    expect_equal(cles.to.d(0.788144600, design = "paired",     verbose = 0), list(d = 0.8, cles = 0.788144600))
    expect_equal(cles.to.d(0.579259710, design = "one.sample", verbose = 0), list(d = 0.2, cles = 0.579259710))
    expect_equal(cles.to.d(0.691462461, design = "one.sample", verbose = 0), list(d = 0.5, cles = 0.691462461))
    expect_equal(cles.to.d(0.788144600, design = "one.sample", verbose = 0), list(d = 0.8, cles = 0.788144600))
    expect_equal(capture.output(cles.to.d(0.5562315, verbose = 1)), c("        d      cles ", "0.2000002 0.5562315 "))
    expect_equal(capture.output(cles.to.d(0.6381632, verbose = 1)), c("        d      cles ", "0.5000000 0.6381632 "))
    expect_equal(capture.output(cles.to.d(0.7141962, verbose = 1)), c("        d      cles ", "0.8000001 0.7141962 "))
})

# means.to.d -----------------------------------------------------------------------------------------------------------
test_that("means.to.d", {
    expect_equal(suppressWarnings(means.to.d(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 15, n2 = 30, verbose = 0)),
                 list(parms =
                      list(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 15, n.ratio = 1, n2 = 30, paired = FALSE, rho.paired = 0.5, verbose = 0),
                      d = 0.223606800, pooled.sd = 11.1803400, var.ratio = 1 / 9, n1 = 30, n2 = 30))
    expect_warning(means.to.d(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 15, n2 = 30, verbose = 0),
                   "Interpretation of Cohen's d may no longer be valid when variances differ beyond sampling error\\.")
    expect_equal(means.to.d(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 6, n2 = 30, verbose = 0),
                 list(parms =
                      list(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 6, n.ratio = 1, n2 = 30, paired = FALSE, rho.paired = 0.5, verbose = 0),
                      d = 0.452678730, pooled.sd =  5.52268051, var.ratio = 0.694444444, n1 = 30, n2 = 30))
    expect_equal(suppressWarnings(means.to.d(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 15, n2 = 30, paired = TRUE, verbose = 0)),
                 list(parms =
                      list(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 15, n.ratio = 1, n2 = 30, paired = TRUE, rho.paired = 0.5, verbose = 0),
                      d = 0.188982237, pooled.sd = 13.2287566, var.ratio = 1 / 9, n1 = 30, n2 = 30))
    expect_warning(means.to.d(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 15, n2 = 30, paired = TRUE, verbose = 0),
                   "Interpretation of Cohen's d may no longer be valid when variances differ beyond sampling error\\.")
    expect_equal(means.to.d(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 6, n2 = 30, paired = TRUE, verbose = 0),
                 list(parms =
                      list(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 6, n.ratio = 1, n2 = 30, paired = TRUE, rho.paired = 0.5, verbose = 0),
                      d = 0.449013255, pooled.sd =  5.56776436, var.ratio = 0.694444444, n1 = 30, n2 = 30))
    expect_equal(capture.output(suppressWarnings(means.to.d(mu1 = 20, mu2 = 17.5, sd1 = 5, sd2 = 15, n2 = 30))),
                 c("        d ", "0.2236068 "))
})

# probs.to.h -------------------------------------------------------------------------------------------------------
test_that("probs.to.h works", {
    expect_equal(probs.to.h(prob1 = 0.56,   prob2 = 0.50, verbose = 0), list(h = 0.120289882, prob1 = 0.56,   prob2 = 0.50))
    expect_equal(probs.to.h(prob1 = 0.60,   prob2 = 0.40, verbose = 0), list(h = 0.402715842, prob1 = 0.60,   prob2 = 0.40))
    expect_equal(probs.to.h(prob1 = 0.25,   prob2 = 0.15, verbose = 0), list(h = 0.251798721, prob1 = 0.25,   prob2 = 0.15))
    expect_equal(probs.to.h(prob1 = 8 / 15, prob2 = 0.40, verbose = 0), list(h = 0.268074069, prob1 = 8 / 15, prob2 = 0.40))
    expect_equal(capture.output(probs.to.h(prob1 = 0.56,   prob2 = 0.50)), c("        h     prob1     prob2 ", "0.1202899 0.5600000 0.5000000 "))
    expect_equal(capture.output(probs.to.h(prob1 = 0.60,   prob2 = 0.40)), c("        h     prob1     prob2 ", "0.4027158 0.6000000 0.4000000 "))
    expect_equal(capture.output(probs.to.h(prob1 = 0.25,   prob2 = 0.15)), c("        h     prob1     prob2 ", "0.2517987 0.2500000 0.1500000 "))
    expect_equal(capture.output(probs.to.h(prob1 = 8 / 15, prob2 = 0.40)), c("        h     prob1     prob2 ", "0.2680741 0.5333333 0.4000000 "))
})

# joint.probs.2x2 ------------------------------------------------------------------------------------------------------
test_that("joint.probs.2x2 works", {
    expect_equal(joint.probs.2x2(prob1 = 0.51, prob2 = 0.49, rho = 0.4141414, verbose = 0),
                 list(parms = list(prob1 = 0.51, prob2 = 0.49, rho = 0.4141414, verbose = 0),
                      rho.min = -1, rho.max = 0.9607843,
                      prob11 = 0.353393936, prob10 = 0.156606064, prob01 = 0.136606064, prob00 = 0.353393936))
    expect_equal(joint.probs.2x2(prob1 = 0.55, prob2 = 0.45, rho = 0.4141414, verbose = 0),
                 list(parms = list(prob1 = 0.55, prob2 = 0.45, rho = 0.4141414, verbose = 0),
                      rho.min = -1, rho.max = 0.818181818,
                      prob11 = 0.35, prob10 = 0.200000004, prob01 = 0.100000003, prob00 = 0.35))
    expect_equal(joint.probs.2x2(prob1 = 0.60, prob2 = 0.40, rho = 0.4141414, verbose = 0),
                 list(parms = list(prob1 = 0.60, prob2 = 0.40, rho = 0.4141414, verbose = 0),
                      rho.min = -1, rho.max = 2 / 3,
                      prob11 = 0.339393936, prob10 = 0.260606064, prob01 = 0.060606064, prob00 = 0.339393936))
    expect_equal(joint.probs.2x2(prob1 = 2 / 3, prob2 = 1 / 3, rho = 0.4141414, verbose = 0),
                 list(parms = list(prob1 = 2 / 3, prob2 = 1 / 3, rho = 0.4141414, verbose = 0),
                      rho.min = -1, rho.max = 0.5,
                      prob11 = 0.314253644, prob10 = 0.352413022, prob01 = 0.019079689, prob00 = 0.314253644))
    expect_equal(joint.probs.2x2(prob1 = 0.90, prob2 = 0.10, rho = 0.1, verbose = 0),
                 list(parms = list(prob1 = 0.90, prob2 = 0.10, rho = 0.1, verbose = 0),
                      rho.min = -1, rho.max = 1 / 9,
                      prob11 = 0.0990, prob10 = 0.8010, prob01 = 0.0010, prob00 = 0.0990))
    expect_equal(capture.output(joint.probs.2x2(prob1 = 0.51, prob2 = 0.49, rho = 0.4141414)),
                 c("   rho.min    rho.max     prob11     prob10     prob01     prob00 ",
                   "-1.0000000  0.9607843  0.3533939  0.1566061  0.1366061  0.3533939 "))
    expect_error(joint.probs.2x2(prob1 = 0.51, prob2 = 0.49, rho = 2, verbose = 0),
                 "Argument `rho` does not have a valid correlation value \\(must be length 1, numeric, >= -1, and <= 1\\)")
    expect_error(joint.probs.2x2(prob1 = 0.51, prob2 = 0.49, rho = 1, verbose = 0),
                 "Combination of `prob1`, `prob2` and `rho` is not feasible.\n`rho` should be between -1 and 0.961")
})

# marginal.probs.2x2 ---------------------------------------------------------------------------------------------------
test_that("marginal.probs.2x2 works", {
    expect_equal(marginal.probs.2x2(prob11 = 0.3534, prob10 = 0.1566, prob01 = 0.1366, prob00 = 0.3534, verbose = 0),
                 list(parms = list(prob11 = 0.3534, prob10 = 0.1566, prob01 = 0.1366, prob00 = 0.3534, verbose = 0),
                      prob1 = 0.51, prob2 = 0.49, rho = 0.414165666))
    expect_equal(marginal.probs.2x2(prob11 = 0.3500, prob10 = 0.2000, prob01 = 0.1000, prob00 = 0.3500, verbose = 0),
                 list(parms = list(prob11 = 0.3500, prob10 = 0.2000, prob01 = 0.1000, prob00 = 0.3500, verbose = 0),
                      prob1 = 0.55, prob2 = 0.45, rho = 0.414141414))
    expect_equal(marginal.probs.2x2(prob11 = 0.3394, prob10 = 0.2606, prob01 = 0.0606, prob00 = 0.3394, verbose = 0),
                 list(parms = list(prob11 = 0.3394, prob10 = 0.2606, prob01 = 0.0606, prob00 = 0.3394, verbose = 0),
                      prob1 = 0.60, prob2 = 0.40, rho = 0.41416667))
    expect_equal(marginal.probs.2x2(prob11 = 0.3143, prob10 = 0.3524, prob01 = 0.0190, prob00 = 0.3143, verbose = 0),
                 list(parms = list(prob11 = 0.3143, prob10 = 0.3524, prob01 = 0.0190, prob00 = 0.3143, verbose = 0),
                      prob1 = round(2 / 3, 4), prob2 = round(1 / 3, 4), rho = 0.41442073))
    expect_equal(marginal.probs.2x2(prob11 = 0.0990, prob10 = 0.8010, prob01 = 0.0010, prob00 = 0.0990, verbose = 0),
                 list(parms = list(prob11 = 0.0990, prob10 = 0.8010, prob01 = 0.0010, prob00 = 0.0990, verbose = 0),
                      prob1 = 0.90, prob2 = 0.10, rho = 0.1))
    expect_equal(capture.output(marginal.probs.2x2(prob11 = 0.3534, prob10 = 0.1566, prob01 = 0.1366, prob00 = 0.3534)),
                 c("    prob1     prob2       rho ", "0.5100000 0.4900000 0.4141657 "))
    expect_error(marginal.probs.2x2(prob11 = 0.35, prob10 = 0.21, prob01 = 0.10, prob00 = 0.35, verbose = 0),
                 "Joint probabilities must sum to 1.")
    expect_warning(marginal.probs.2x2(prob11 = 0.35, prob10 = 0.65, prob01 = 0.0, prob00 = 0.0, verbose = 0),
                   "Undefined correlation: division by zero in denominator.")
})

# probs.to.w -----------------------------------------------------------------------------------------------------------
test_that("probs.to.w works", {
    expect_equal(probs.to.w(c(0.28, 0.72), rep(0.5, 2), verbose = 0),
                 list(w = 0.44, df = 1, prob.matrix = c(0.28, 0.72), null.prob.matrix = rep(0.5, 2)))
    expect_equal(probs.to.w(c(0.28, 0.72), verbose = 0),
                 list(w = 0.44, df = 1, prob.matrix = c(0.28, 0.72), null.prob.matrix = rep(0.5, 2)))
    crrPrb <- rbind(c(0.056, 0.132), c(0.944, 0.868))
    expect_equal(probs.to.w(crrPrb, verbose = 0),
                 list(w = 0.130213368, df = 1, prob.matrix = crrPrb, null.prob.matrix = outer(rowSums(crrPrb), colSums(crrPrb)) / sum(crrPrb)))
    crrPrb <- cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078), c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101))
    expect_equal(probs.to.w(crrPrb, verbose = 0),
                 list(w = 0.030220075, df = 4, prob.matrix = crrPrb, null.prob.matrix = outer(rowSums(crrPrb), colSums(crrPrb)) / sum(crrPrb)))
    expect_equal(capture.output(probs.to.w(c(0.28, 0.72), rep(0.5, 2))), c("   w   df ", "0.44 1.00 "))
    expect_equal(capture.output(probs.to.w(rbind(c(0.056, 0.132), c(0.944, 0.868)))), c("        w        df ", "0.1302134 1.0000000 "))
    expect_equal(capture.output(probs.to.w(cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078), c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101)))),
                 c("         w         df ", "0.03022008 4.00000000 "))
    expect_error(probs.to.w(c(0, 1.1)), "Matrix elements outside of \\[0, 1\\] range.")
    expect_error(probs.to.w(c(0.2, 0.8), c(0, 1.1)), "Matrix elements outside of \\[0, 1\\] range.")
    expect_error(probs.to.w(c(0.1, 0.2, 0.3, 0.4), c(0.5, 0.5)), "Length of `prob.matrix` and `null.prob.matrix` should match.")
    expect_error(probs.to.w(c(0.2, 0.3, 0.4, 0.5)), "Cell probabilities should sum to 1.")
    expect_error(probs.to.w(matrix(c(0.2, 0.3, 0.4, 0.5), nrow = 2), matrix(c(0.2, 0.3, 0.4, 0.5), nrow = 1)),
                 "Dimensions for `prob.matrix` and `null.prob.matrix` do not match.")
    expect_error(probs.to.w(data.frame(A = c(0.2, 0.3, 0.4, 0.5))), "`prob.matrix` must be either a vector or a matrix.")
})

# prob.limits.paired ---------------------------------------------------------------------------------------------------
test_that("prob.limits.paired works", {
    expect_equal(prob.limits.paired(prob1 = 0.40, prob2 = NULL), c(prob.min = 0.1429, prob.max = 0.7272))
    expect_equal(prob.limits.paired(prob1 = 0.50, prob2 = NULL), c(prob.min = 0.2000, prob.max = 0.7999))
    expect_equal(prob.limits.paired(prob1 = NULL, prob2 = 0.50), c(prob.min = 0.2000, prob.max = 0.7999))
    expect_equal(prob.limits.paired(prob1 = 0.10, prob2 = NULL), c(prob.min = 0.0271, prob.max = 0.3076))
    expect_equal(prob.limits.paired(prob1 = 0.01, prob2 = NULL), c(prob.min = 0.0026, prob.max = 0.0388))
    expect_equal(prob.limits.paired(prob1 = 0.10, prob2 = NULL, rho = 1.0),   c(prob.min = 0.1000, prob.max = 0.1000))
    expect_equal(prob.limits.paired(prob1 = 0.10, prob2 = NULL, rho = 0.8),   c(prob.min = 0.0664, prob.max = 0.1479))
    expect_equal(prob.limits.paired(prob1 = 0.10, prob2 = NULL, rho = 0.2),   c(prob.min = 0.0045, prob.max = 0.7352))
    expect_equal(prob.limits.paired(prob1 = 0.10, prob2 = NULL, rho = 0.1),   c(prob.min = 0.0012, prob.max = 0.9174))
    expect_equal(prob.limits.paired(prob1 = 0.10, prob2 = NULL, rho = 0.01),  c(prob.min = 0.0001, prob.max = 0.9991))
    expect_equal(prob.limits.paired(prob1 = 0.10, prob2 = NULL, rho = 1e-12), c(prob.min = 0.0001, prob.max = 0.9999))
    expect_error(prob.limits.paired(prob1 = NULL, prob2 = NULL),
                 "Exactly one of the parameters `prob1` or `prob2` must be given, one has to be NULL.")
    expect_error(prob.limits.paired(prob1 = 1e-6, prob2 = NULL),
                 "No feasible values found. Check that rho is achievable given the fixed probability.")
})

# means.to.etasq -------------------------------------------------------------------------------------------------------
test_that("means.to.etasq works", {
    expect_equal(means.to.etasq(mu.vector = c(0.50, 0), sd.vector = rep(1, 2), n.vector = rep(33, 2), k.cov = 1,
                                r.squared = 0.50, verbose = 0),
                 list(f = 0.353553391, eta.squared = 1 / 9, df1 = 1, df2 = 63, ncp = 8.25))
    expect_equal(means.to.etasq(mu.vector = c(0.50, 0), sd.vector = rep(1, 2), n.vector = rep(33, 2), verbose = 0),
                 list(f = 0.25, eta.squared = 0.05882353, df1 = 1, df2 = 64, ncp = 4.125))
    expect_equal(means.to.etasq(mu.vector = seq(-1, 1, 0.5), sd.vector = rep(1, 5), n.vector = rep(20, 5), k.cov = 1,
                                r.squared = 0.50, verbose = 0),
                 list(f = 1, eta.squared = 0.5, df1 = 4, df2 = 94, ncp = 100))
    expect_equal(means.to.etasq(mu.vector = seq(-1, 1, 0.5), sd.vector = rep(1, 5), n.vector = rep(20, 5), verbose = 0),
                 list(f = sqrt(0.5), eta.squared = 1 / 3, df1 = 4, df2 = 95, ncp = 50))
    expect_output(means.to.etasq(mu.vector = seq(-1, 1, 0.5), sd.vector = rep(1, 5), n.vector = rep(20, 5), verbose = 1),
                  "f eta.squared         df1         df2         ncp \n  0.7071068   0.3333333   4.0000000  95.0000000  50.0000000")
    expect_error(means.to.etasq(mu.vector = NULL, sd.vector = rep(1, 5), n.vector = rep(20, 5), k.cov = 0, verbose = 0),
                 "Provide a vector of means \\(`mu.vector`\\) with its length equal to number of groups.")
    expect_error(means.to.etasq(mu.vector = seq(-1, 1, 0.5), sd.vector = rep(1, 5), n.vector = rep(20, 5),
                                factor.levels = c(2, 2), verbose = 0),
                 "Factorial designs are not allowed.")
    expect_error(means.to.etasq(mu.vector = seq(-1, 1, 0.5), sd.vector = rep(1, 5), n.vector = rep(20, 5),
                                factor.levels = 4, verbose = 0),
                 "Length of the vector of means \\(`mu.vector`\\) does not match number of levels.")
    expect_error(means.to.etasq(mu.vector = seq(-1, 1, 0.5), sd.vector = rep(1, 5), n.vector = rep(20, 5), k.cov = 1,
                                r.squared = -1e-12, verbose = 0),
                 "R-squared \\(explanatory power of covariates\\) takes a value between 0 and 1.")
    expect_error(means.to.etasq(mu.vector = seq(-1, 1, 0.5), sd.vector = rep(1, 5), n.vector = rep(20, 5), k.cov = 0,
                                r.squared = 0.50, verbose = 0),
                 "Explanatory power of covariates is expected to be non-zero when number of covariates is non-zero.")
})
