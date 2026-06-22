## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      comment = "#>",
                      collapse = TRUE,
                      message = FALSE,
                      warning = FALSE,
                      fig.height = 5,
                      fig.width = 7,
                      fig.align = "center",
                      out.width = "100%")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(pwrss)

## ----message = FALSE, eval = FALSE, warning = FALSE---------------------------
# install.packages("pwrss")
# library(pwrss)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
data(mtcars)

model <- lm(mpg ~ hp + wt, data = mtcars)
summary(model)

power.t.test(ncp = -3.519, # t-value for hp variable
             df = 29,      # residual degrees of freedom
             alpha = 0.05, # type 1 error rate
             alternative = "two.sided",
             plot = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.t.test(ncp = NULL,     # calculate this
             req.sign = "-", # direction of the test statistics
             df = 29,        # residual degrees of freedom
             alpha = 0.05,   # type 1 error rate
             power = 0.80,
             alternative = "two.sided",
             plot = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# define parameters
d <- 0.50           # effect size under alternative
null.d <- 0         # effect size under null
margin <- 0.10      # smallest meaningful diff between d and null.d

p <- 0.50           # proportion of subjects in the intervention school
n <- 120            # total sample size
g <- 1              # number of covariates
r.squared <- 0.50   # explanatory power of covariates
rtx.squared <- 0.10 # squared point-biserial cor between trt dummy and outcome

# calculate the standard error
std.error <- sqrt((1 - r.squared) / (p * (1 - p) * n * (1 - rtx.squared)))

# calculate non-centrality parameters
ncp <- (d - null.d) / std.error
null.ncp <- margin / std.error

# calculate power
power.t.test(ncp = ncp,
             null.ncp = null.ncp,
             df = n - g - 2,
             alpha = 0.05,
             alternative = "one.sided",
             plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
data(warpbreaks)

model <- glm(breaks ~ wool + tension, data = warpbreaks,
             family = poisson(link = "log"))
summary(model)

power.z.test(mean = -3.994, # z-value for wool B
             alpha = 0.05,  # type 1 error rate
             alternative = "two.sided",
             plot = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.z.test(mean = NULL,   # calculate this
             alpha = 0.05,  # type 1 error rate
             power = 0.80,
             alternative = "two.sided",
             plot = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# define parameters
rs <- 0.30   # spearman rho rank cor under alternative
null.rs <- 0 # spearman rho rank cor under null
n <- 100     # sample size

# apply Fisher's Z transformation
z.rs <- cor.to.z(rs)$z
z.null.rs <- cor.to.z(null.rs)$z

# calculate the standard error
z.std.error <- sqrt(1.06 / (n - 3))

# calculate the non-centrality parameter
ncp <- (z.rs - z.null.rs) / z.std.error

# calculate power
power.z.test(mean = ncp,
             alpha = 0.05,
             alternative = "two.sided",
             plot = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
data(trees)

# model <- aov(Girth ~ Height, data = trees)
model <- lm(Girth ~ Height, data = trees)
summary(model)

power.f.test(ncp = 10.71,  # non-centrality under alternative
             df1 = 1,      # numerator degrees of freedom
             df2 = 29,     # denominator degrees of freedom
             alpha = 0.05, # type 1 error rate
             plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.f.test(ncp = NULL,   # calculate this
             df1 = 1,      # numerator degrees of freedom
             df2 = 29,     # denominator degrees of freedom
             power = 0.80,
             alpha = 0.05, # type 1 error rate
             plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
data(HairEyeColor)

table <- margin.table(HairEyeColor, c(1, 2))
print(table)

chisq.test(table)

power.chisq.test(ncp = 138.29, # X-squared
                 df = 9,       # degrees of freedom
                 alpha = 0.05, # type 1 error rate
                 plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.chisq.test(ncp = NULL,    # calculate this
                 df = 9,        # degrees of freedom
                 power = 0.80,
                 alpha = 0.05,  # type 1 error rate
                 plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
data(infert)

fit.reduced <- glm(case ~ age + parity + spontaneous,
                   data = infert,
                   family = binomial)

fit.full <- glm(case ~ age + parity + spontaneous + induced,
                data = infert,
                family = binomial)

anova(fit.reduced, fit.full)

power.chisq.test(ncp = 18.463,
                 df = 1,
                 alpha = 0.05,
                 plot = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.chisq.test(ncp = NULL,    # calculate this
                 df = 1,        # degrees of freedom
                 power = 0.80,
                 alpha = 0.05,  # type 1 error rate
                 plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
data(faithful)

long <- faithful$eruptions > 3
n.success <- sum(long)
n.total <- nrow(faithful)

binom.test(n.success, n.total, p = 0.50)

power.binom.test(size = n.total,             # number of eruptions
                 prob = n.success / n.total, # prob. of occurrence under alt.
                 null.prob = 0.50,           # prob. of occurrence under null
                 alpha = 0.05,
                 alternative = "one.sided",
                 plot = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.binom.test(size = NULL,                 # number of tosses needed
                 prob = 0.50,                 # prob. of head under alt.
                 null.prob = c(0.495, 0.505), # equivalence margins
                 power = 0.80,
                 alpha = 0.05,                # type 1 error rate
                 alternative = "two.one.sided",
                 plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.binom.test(size = NULL,                 # number of replications needed
                 prob = 0.05,                 # prob. of falsely rejecting null
                 null.prob = c(0.045, 0.055), # equivalence margins
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.one.sided",
                 plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.binom.test(size = NULL,                 # number of replications needed
                 prob = 0.80,                 # prob. of correctly rejecting null
                 null.prob = c(0.795, 0.805), # equivalence margins
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.one.sided",
                 plot = FALSE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = FALSE, warning = FALSE----
power.t.student(d = 0.20, power = 0.80) |>
  plot()

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = TRUE----
power.t.student(d = -0.20,
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                design = "independent")

# account for attrition
inflate.sample(n = 310, rate = 0)    # control
inflate.sample(n = 310, rate = 0.05) # treatment

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = TRUE----
power.t.welch(d = 0.20,
              power = 0.80,
              n.ratio = 2,
              var.ratio = 1.5,
              alpha = 0.05,
              alternative = "two.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = TRUE----
power.np.wilcoxon(d = 0.20,
                  power = 0.80,
                  n.ratio = 1,
                  alpha = 0.05,
                  alternative = "two.sided",
                  design = "independent",
                  distribution = "normal")

inflate.sample(n = 412, rate = 0.05) # treatment

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = 0.20,
                margin = -0.05,
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                design = "independent")

# consider 0.05 attrition rate
inflate.sample(n = 398, rate = 0.05)

# robust parametric
power.t.welch(d = 0.20,
              margin = -0.05,
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alpha = 0.05,
              alternative = "one.sided")

# non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = -0.05,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "one.sided",
                  design = "independent")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = 0.20,
                margin = 0.05,
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                design = "independent")

# consider 0.05 attrition rate
inflate.sample(n = 1104, rate = 0.05)

# robust parametric
power.t.welch(d = 0.20,
              margin = 0.05,
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alpha = 0.05,
              alternative = "one.sided")

# non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = 0.05,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "one.sided",
                  design = "independent")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = 0,
                margin = c(-0.10, 0.10),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.one.sided",
                design = "independent")

# consider 0.05 attrition rate
inflate.sample(n = 3428, rate = 0.05)

# robust parametric
power.t.welch(d = 0,
              margin = c(-0.10, 0.10),
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alpha = 0.05,
              alternative = "two.one.sided")

# non-parametric
power.np.wilcoxon(d = 0,
                  margin = c(-0.10, 0.10),
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.one.sided",
                  design = "independent")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = 0.20,
                margin = c(-0.05, 0.05),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.one.sided",
                design = "independent")

# consider 0.05 attrition rate
inflate.sample(n = 1400, rate = 0.05)

# robust parametric
power.t.welch(d = 0.20,
              margin = c(-0.05, 0.05),
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alpha = 0.05,
              alternative = "two.one.sided")

# non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = c(-0.05, 0.05),
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.one.sided",
                  design = "independent")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = TRUE----
power.t.student(d = 0.20,
                power = 0.80,
                alpha = 0.05,
                alternative = "two.sided",
                design = "paired")

# consider an attrition rate of 0.05
inflate.sample(n = 199, rate = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = TRUE----
power.np.wilcoxon(d = 0.20,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.sided",
                  design = "paired",
                  distribution = "normal")

# consider an attrition rate of 0.05
inflate.sample(n = 208, rate = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = 0,
                margin = -0.10,
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                design = "paired")

# consider 0.05 attrition rate
inflate.sample(n = 619, rate = 0.05)

# non-parametric
power.np.wilcoxon(d = 0,
                  margin = -0.10,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "one.sided",
                  design = "paired")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = 0.20,
                margin = 0.10,
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                design = "paired")

# consider 0.05 attrition rate
inflate.sample(n = 627, rate = 0.05)

# non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = 0.10,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "one.sided",
                  design = "paired")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = 0,
                margin = c(-0.10, 0.10),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.one.sided",
                design = "paired")

# consider 0.05 attrition rate
inflate.sample(n = 858, rate = 0.05)

# non-parametric
power.np.wilcoxon(d = 0,
                  margin = c(-0.10, 0.10),
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.one.sided",
                  design = "paired")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric (an example report is provided below)
power.t.student(d = -0.20,
                margin = c(-0.10, 0.10),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.one.sided",
                design = "paired")

# consider 0.05 attrition rate
inflate.sample(n = 797, rate = 0.05)

# non-parametric
power.np.wilcoxon(d = -0.20,
                  margin = c(-0.05, 0.05),
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.one.sided",
                  design = "paired")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = TRUE----
power.t.student(d = 0.40,
                power = 0.80,
                alpha = 0.05,
                alternative = "two.sided",
                design = "one.sample")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = TRUE----
power.np.wilcoxon(d = -0.20,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.sided",
                  design = "one.sample",
                  distribution = "normal")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# parametric
power.t.student(d = 0.40,
                margin = 0.10,
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                design = "one.sample")

# non-parametric
power.np.wilcoxon(d = 0.40,
                  margin = 0.10,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "one.sided",
                  design = "one.sample")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
d <- (7.5 - 7.5) / 1.2
margin <- c((6.5 - 7.5) / 1.2, (8.5 - 7.5) / 1.2)

# parametric
power.t.student(d = d,
                margin = margin,
                power = 0.80,
                alpha = 0.05,
                alternative = "two.one.sided",
                design = "one.sample")

# non-parametric
power.np.wilcoxon(d = d,
                  margin = margin,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.one.sided",
                  design = "one.sample")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
d <- (9 - 7.5) / 1.2
margin <- c((6.5 - 7.5) / 1.2, (8.5 - 7.5) / 1.2)

# parametric
power.t.student(d = d,
                margin = margin,
                power = 0.80,
                alpha = 0.05,
                alternative = "two.one.sided",
                design = "one.sample")

# non-parametric
power.np.wilcoxon(d = d,
                  margin = margin,
                  power = 0.80,
                  alpha = 0.05,
                  alternative = "two.one.sided",
                  design = "one.sample")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.oneprop(prob = 0.80, # probability of success under alternative
                null.prob = 0.90, # probability of success under null
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided")


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.oneprop(prob = 0.80, # probability of success under alternative
                null.prob = 0.90, # probability of success under null
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                arcsine = TRUE)


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.oneprop(prob = 0.80, # probability of success under alternative
                null.prob = 0.90, # probability of success under null
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                correct = TRUE)


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.oneprop(prob = 0.80, # probability of success under alternative
                null.prob = 0.90, # probability of success under null
                power = 0.80,
                alpha = 0.05,
                alternative = "one.sided",
                std.error = "alternative")


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.exact.oneprop(prob = 0.002, # probability of success under alternative
                    null.prob = 0, # probability of success under null
                    power = 0.80,
                    alpha = 0.05,
                    alternative = "one.sided")


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# z-test approach
power.z.twoprops(prob1 = 0.60, prob2 = 0.50,
                 power = 0.80, arcsine = TRUE)

# find Cohen's h
probs.to.h(prob1 = 0.60, prob2 = 0.50)

# t-test approach
power.t.student(d = 0.2013579, power = 0.80)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.15,
                 prob2 = 0.10,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.15,
                 prob2 = 0.10,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "one.sided",
                 arcsine = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.15,
                 prob2 = 0.10,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "one.sided",
                 correct = TRUE)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.15,
                 prob2 = 0.10,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "one.sided",
                 std.error = "unpooled")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.exact.twoprops(prob1 = 0.60,
                     prob2 = 0.50,
                     power = 0.80,
                     alpha = 0.05,
                     alternative = "one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.01,
                 prob2 = 0.02,
                 margin = 0.01,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "one.sided")

# consider 5% attrition rate
inflate.sample(n = 914, rate = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.45,
                 prob2 = 0.50,
                 margin = -0.01,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "one.sided")

# consider 5% attrition rate
inflate.sample(n = 3852, rate = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.10,
                 prob2 = 0.10,
                 margin = c(-0.02, 0.02),
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.13,
                 prob2 = 0.10,
                 margin = c(-0.01, 0.01),
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twoprops(prob1 = 0.50,
                 prob2 = 0.40,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.sided",
                 paired = TRUE,
                 rho.paired = 0.50)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.exact.twoprops(prob1 = 0.50,
                     prob2 = 0.40,
                     power = 0.80,
                     alpha = 0.05,
                     alternative = "two.sided",
                     paired = TRUE,
                     rho.paired = 0.50)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.onecor(rho = 0.20,
               null.rho = 0.10,
               power = 0.80,
               alpha = 0.05,
               alternative = "one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# z-test approach
power.z.twocors(rho1 = 0.20, rho2 = 0.10, power = .80)

# find Cohen's q
cors.to.q(rho1 = 0.20, rho2 = 0.10)

# t-test approach
power.t.student(d = 0.1023972, power = 0.80)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.twocors(rho1 = 0.10,
                rho2 = 0,
                power = .80,
                alpha = 0.05,
                alternative = "one.sided")

# calculate Cohen's h
cors.to.q(rho1 = 0.10, rho2 = 0)

# t-test approximation
power.t.student(d = 0.1003353,
                power = .80,
                alpha = 0.05,
                alternative = "one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# example data for one common index
# compare cor(V1, V2) to cor(V1, V3)

# subject    V1       V2      V3
# <int>    <dbl>    <dbl>    <dbl>
#   1       1.2      2.3      0.8
#   2      -0.0      1.1      0.7
#   3       1.9     -0.4     -2.3
#   4       0.7      1.3      0.4
#   5       2.1     -0.1      0.8
#   ...     ...      ...      ...
#   1000   -0.5      2.7     -1.7

# V1: socio-economic status (common)
# V2: pre-test
# V3: post-test

power.z.twocors.steiger(rho12 = 0.50,
                        rho13 = 0.40,
                        rho23 = 0.70,
                        power = 0.80,
                        alpha = 0.05,
                        alternative = "one.sided",
                        common.index = TRUE)

# calculate the effect size
cors.to.q(rho1 = 0.40, rho2 = 0.50)

# adjust the sample for 5% attrition rate
inflate.sample(n = 286, rate = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# example data for no common index
# compare cor(V1, V2) to cor(V3, V4)

# subject    V1       V2       V3       V4
# <int>    <dbl>    <dbl>    <dbl>    <dbl>
#   1       1.2      2.3      0.8      1.2
#   2      -0.0      1.1      0.7      0.9
#   3       1.9     -0.4     -2.3     -0.1
#   4       0.7      1.3      0.4     -0.3
#   5       2.1     -0.1      0.8      2.7
#   ...     ...      ...      ...      ...
#   1000   -0.5      2.7     -1.7      0.8

# V1: pre-test reading
# V2: pre-test math
# V3: post-test reading
# V4: post-test math

power.z.twocors.steiger(rho12 = 0.50, # cor(V1, V2)
                        rho13 = 0.70,
                        rho23 = 0.40,
                        rho14 = 0.40,
                        rho24 = 0.70,
                        rho34 = 0.60, # cor(V3, V4)
                        power = 0.80,
                        alpha = 0.05,
                        alternative = "one.sided",
                        common.index = FALSE)

# calculate the effect size
cors.to.q(rho1 = 0.60, rho2 = 0.50)

# adjust the sample for 5% attrition rate
inflate.sample(n = 317, rate = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.regression(r.squared = 0.10,
                   k.total = 3, # number of total predictors
                   power = 0.80,
                   alpha = 0.05)


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.regression(r.squared = 0.10,
                   margin = 0.05,
                   k.total = 3,
                   power = 0.80,
                   alpha = 0.05)


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.regression(r.squared.change = 0.10,
                   k.total = 5, # number of total predictors
                   k.tested = 2, # number of tested predictors
                   power = 0.80,
                   alpha = 0.05)


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.regression(r.squared.change = 0.10,
                   margin = 0.05,
                   k.total = 5, # number of total predictors
                   k.tested = 2, # number of tested predictors
                   power = 0.80,
                   alpha = 0.05)


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.t.regression(beta = 0.20,
                   k.total = 3,
                   r.squared = 0.30,
                   power = .80,
                   alpha = 0.05,
                   alternative = "two.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.t.regression(beta = 0.60,
                   sd.outcome = 12,
                   sd.predictor = 4,
                   k.total = 3,
                   r.squared = 0.30,
                   power = .80,
                   alpha = 0.05,
                   alternative = "two.sided")


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
p <- 0.50
sd.predictor <- sqrt(p * (1 - p))

power.t.regression(beta = 0.20,
                   sd.predictor = sd.predictor,
                   k.total = 3,
                   r.squared = 0.30,
                   power = .80,
                   alpha = 0.05,
                   alternative = "two.sided")


## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
p <- 0.50
sd.predictor <- sqrt(p * (1 - p))

power.t.regression(beta = 0.20,
                   null.beta = 0.10,
                   margin = -0.05,
                   sd.predictor = sd.predictor,
                   k.total = 3,
                   r.squared = 0.30,
                   power = .80,
                   alpha = 0.05,
                   alternative = "one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
p <- 0.50
sd.predictor <- sqrt(p * (1 - p))

power.t.regression(beta = 0.20,
                   null.beta = 0.10,
                   margin = 0.05,
                   sd.predictor = sd.predictor,
                   k.total = 3,
                   r.squared = 0.30,
                   power = .80,
                   alpha = 0.05,
                   alternative = "one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
p <- 0.50
sd.predictor <- sqrt(p * (1 - p))

power.t.regression(beta = 0.20,
                   null.beta = 0.20,
                   margin = c(-0.05, 0.05),
                   sd.predictor = sd.predictor,
                   k.total = 3,
                   r.squared = 0.30,
                   power = .80,
                   alpha = 0.05,
                   alternative = "two.one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
p <- 0.50
sd.predictor <- sqrt(p * (1 - p))

power.t.regression(beta = 0.20,
                   margin = c(-0.05, 0.05),
                   sd.predictor = sd.predictor,
                   k.total = 3,
                   r.squared = 0.30,
                   power = .80,
                   alpha = 0.05,
                   alternative = "two.one.sided")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.logistic(prob = 0.10,
                 base.prob = 0.15,
                 r.squared.predictor = 0.20,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.sided",
                 distribution = "normal")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.logistic(odds.ratio = 0.6296,
                 base.prob = 0.15,
                 r.squared.predictor = 0.20,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.sided",
                 distribution = "normal")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.logistic(beta1 = -0.4626,
                 base.prob = 0.15,
                 r.squared.predictor = 0.20,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.sided",
                 distribution = "normal")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
distribution <- list(dist = "normal", mean = 20, sd = 8)

power.z.logistic(beta1 = -0.4626,
                 base.prob = 0.15,
                 r.squared.predictor = 0.20,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.sided",
                 distribution = distribution)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.logistic(beta1 = -0.4626,
                 base.prob = 0.15,
                 r.squared.predictor = 0.20,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.sided",
                 distribution = "bernoulli")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
distribution <- list(dist = "bernoulli", prob = 0.40)

power.z.logistic(beta1 = -0.4626,
                 base.prob = 0.15,
                 r.squared.predictor = 0.20,
                 power = 0.80,
                 alpha = 0.05,
                 alternative = "two.sided",
                 distribution = distribution)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.poisson(beta0 = 0.50,
                beta1 = -0.10,
                power = 0.80,
                alpha = 0.05,
                alternative = "two.sided",
                distribution = "normal")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.sided",
                distribution = "normal")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
distribution <- list(dist = "normal", mean = 20, sd = 8)

power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.sided",
                distribution = distribution)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.sided",
                distribution = "bernoulli")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
distribution <- list(dist = "bernoulli", prob = 0.40)

power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                power = 0.80,
                alpha = 0.05,
                alternative = "two.sided",
                distribution = distribution)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# mediation model with base R-squared values
power.z.mediation(beta.a = 0.25,
                  beta.b = 0.25,
                  power = 0.80,
                  alpha = 0.05,
                  method = "sobel")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
# base R-squared values are 0 (zero)
# do not specify 'cp'
power.z.mediation(beta.a = 0.25,
                  beta.b = 0.25,
                  r.squared.mediator = 0,
                  r.squared.outcome = 0,
                  power = 0.80,
                  alpha = 0.05,
                  method = "sobel")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
p <- 0.50 # proportion of subjects in one of the groups
sd.predictor <- sqrt(p * (1 - p))

power.z.mediation(beta.a = 0.25,
                  beta.b = 0.25,
                  sd.predictor = sd.predictor,
                  power = 0.80,
                  alpha = 0.05,
                  method = "sobel")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# binary X
p <- 0.50 # proportion of subjects in one of the groups
sd.predictor <- sqrt(p * (1 - p))

# set seed to ensure the same results throughout different runs of examples.Rmd
set.seed(1)

power.z.mediation(beta.a = 0.25,
                  beta.b = 0.25,
                  sd.predictor = sd.predictor,
                  n = 300,
                  alpha = 0.05,
                  method = "joint")

power.z.mediation(beta.a = 0.25,
                  beta.b = 0.25,
                  sd.predictor = sd.predictor,
                  n = 300,
                  alpha = 0.05,
                  method = "monte.carlo")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# continuous X
power.z.mediation(beta.a = 0.25,
                  beta.b = 0.25,
                  r.squared.mediator = 0.50,
                  r.squared.outcome = 0.50,
                  power = 0.80,
                  alpha = 0.05,
                  method = "sobel")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# binary X
p <- 0.50 # proportion of subjects in one of the groups
sd.predictor <- sqrt(p * (1 - p))

power.z.mediation(beta.a = 0.25,
                  beta.b = 0.25,
                  sd.predictor = sd.predictor,
                  r.squared.outcome = 0.50,
                  power = 0.80,
                  alpha = 0.05,
                  method = "sobel")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova(eta.squared = 0.059,
               factor.levels = 2,
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova(eta.squared = 0.03,
               factor.levels = c(2, 2),
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova(eta.squared = 0.02,
               factor.levels = c(2, 2, 3),
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova(eta.squared = 0.02,
               null.eta.squared = 0.01,
               factor.levels = c(2, 2, 3),
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova.keppel(mu.vector = c(0.50, 0), # vector of means
                      sd.vector = c(1, 1), # vector of standard deviations
                      p.vector = c(0.50, 0.50), # sample allocation rates
                      power = 0.80,
                      alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.mixed.anova(eta.squared = 0.059,
                    factor.levels = c(2, 1), # c("between", "within")
                    power = 0.80,
                    alpha = 0.05,
                    effect = "between")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.mixed.anova(eta.squared = 0.022,
                    factor.levels = c(1, 2), # c("between", "within")
                    power = 0.80,
                    alpha = 0.05,
                    rho.within = 0.50,
                    effect = "within")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.mixed.anova(eta.squared = 0.038,
                    factor.levels = c(2, 2), # c("between", "within")
                    power = 0.80,
                    alpha = 0.05,
                    rho.within = 0.50,
                    effect = "between")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.mixed.anova(eta.squared = 0.022,
                    factor.levels = c(2, 2), # c("between", "within")
                    power = 0.80,
                    alpha = 0.05,
                    rho.within = 0.50,
                    effect = "within")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.mixed.anova(eta.squared = 0.01,
                    factor.levels = c(2, 2), # c("between", "within")
                    power = 0.80,
                    alpha = 0.05,
                    rho.within = 0.50,
                    effect = "interaction")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.mixed.anova(eta.squared = 0.05,
                    factor.levels = c(2, 2), # c("between", "within")
                    power = 0.80,
                    alpha = 0.05,
                    rho.within = NA,
                    effect = "between")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.f.mixed.anova(eta.squared = 0.05,
                    null.eta.squared = 0.01,
                    factor.levels = c(2, 2), # c("between", "within")
                    power = 0.80,
                    alpha = 0.05,
                    rho.within = NA,
                    effect = "between")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova(eta.squared = 0.048,
               factor.levels = 2,
               k.covariates = 1,
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova(eta.squared = 0.02,
               factor.levels = c(2, 2),
               k.covariates = 1,
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova(eta.squared = 0.01,
               factor.levels = c(2, 2, 3),
               k.covariates = 1,
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.f.ancova(eta.squared = 0.048,
               null.eta.squared = 0.01,
               factor.levels = 2,
               k.covariates = 1,
               power = 0.80,
               alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova.keppel(mu.vector = c(0.318, 0), # vector of adjusted means
                      sd.vector = c(1, 1), # vector of unadjusted standard deviations
                      p.vector = c(0.50, 0.50), # sample allocation rates
                      r.squared = 0.50,
                      k.covariates = 1,
                      power = 0.80,
                      alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova.shieh(mu.vector = c(0.318, 0), # vector of adjusted means
                      sd.vector = c(1, 1), # vector of unadjusted standard deviations
                      p.vector = c(0.50, 0.50), # sample allocation rates
                      r.squared = 0.50,
                      k.covariates = 1,
                      power = 0.80,
                      alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
power.f.ancova.shieh(mu.vector = c(0.30, 0.09, 0.05, 0.245), # vector of adjusted means
                     sd.vector = c(1, 1, 1, 1), # vector of unadjusted standard deviations
                     p.vector = c(0.25, 0.25, 0.25, 0.25), # sample allocation rates
                     factor.levels = c(2, 2),
                     r.squared = 0.50,
                     k.covariates = 1,
                     power = 0.80,
                     alpha = 0.05)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
contr.obj <- factorial.contrasts(factor.levels = 3,
                                 coding = "treatment")

# get contrast matrix from the contrast object
contrast.matrix <- contr.obj$contrast.matrix

# calculate sample size given design characteristics
design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                               sd.vector = c(1, 1, 1), # unadjusted standard deviations
                               p.vector = c(1 / 3, 1 / 3, 1 / 3), # sample allocation rate
                               contrast.matrix = contrast.matrix,
                               r.squared = 0.50,
                               k.covariates = 1,
                               alpha = 0.05,
                               power = 0.80)

# power of planned contrasts, adjusted for alpha level
power.t.contrasts(design, adjust.alpha = "fdr")

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
contr.obj <- factorial.contrasts(factor.levels = 3,
                                 coding = "helmert")

# get contrast matrix from the contrast object
contrast.matrix <- contr.obj$contrast.matrix

# calculate sample size given design characteristics
design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                               sd.vector = c(1, 1, 1), # unadjusted standard deviations
                               p.vector = c(1 / 3, 1 / 3, 1 / 3), # allocation, should sum to 1
                               contrast.matrix = contrast.matrix,
                               r.squared = 0.50, # proportion of variance explained by covariates
                               k.covariates = 1, # number of covariates
                               alpha = 0.05,
                               power = 0.80)

# power of planned contrasts
power.t.contrasts(design)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
contr.obj <- factorial.contrasts(factor.levels = 3,
                                 coding = "poly")

# get contrast matrix from the contrast object
contrast.matrix <- contr.obj$contrast.matrix

# calculate sample size given design characteristics
design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                               sd.vector = c(1, 1, 1), # unadjusted standard deviations
                               p.vector = c(1 / 3, 1 / 3, 1 / 3), # allocation, should sum to 1
                               contrast.matrix = contrast.matrix,
                               r.squared = 0.50, # proportion of variance explained by covariates
                               k.covariates = 1, # number of covariates
                               alpha = 0.05,
                               power = 0.80)

# power of the planned contrasts
power.t.contrasts(design)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
# custom contrasts
contrast.matrix <- rbind(
  cbind(A1 = 1, A2 = -0.50, A3 = -0.50),
  cbind(A1 = 0.50, A2 = 0.50, A3 = -1)
)
# labels are not required for custom contrasts,
# but they make it easier to understand power.t.contrasts() output

# calculate sample size given design characteristics
design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                               sd.vector = c(1, 1, 1), # unadjusted standard deviations
                               p.vector = c(1 / 3, 1 / 3, 1 / 3), # allocation, should sum to 1
                               contrast.matrix = contrast.matrix,
                               r.squared = 0.50, # proportion of variance explained by covariates
                               k.covariates = 1, # number of covariates
                               alpha = 0.05,
                               power = 0.80)

# power of the planned contrasts
power.t.contrasts(design)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
prob.matrix <- c(0.28, 0.72)

probs.to.w(prob.matrix = prob.matrix)
# Degrees of freedom is k - 1 for Cohen's w.

power.chisq.gof(w = 0.44,
                df = 1,
                alpha = 0.05,
                power = 0.80)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
prob.matrix <- rbind(c(0.056, 0.132),
                  c(0.944, 0.868))
colnames(prob.matrix) <- c("Girl", "Boy")
rownames(prob.matrix) <- c("ADHD", "No ADHD")
prob.matrix

probs.to.w(prob.matrix = prob.matrix)
# Degrees of freedom is 1 for Phi coefficient.

power.chisq.gof(w = 0.1302134,
                df = 1,
                alpha = 0.05,
                power = 0.80)

## ----message = FALSE, fig.width = 7, fig.height = 5, results = TRUE-----------
prob.matrix <- cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078),
                  c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101))
rownames(prob.matrix) <- c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
colnames(prob.matrix) <- c("Female", "Male")
prob.matrix

probs.to.w(prob.matrix = prob.matrix)
# Degrees of freedom is (nrow - 1) * (ncol - 1) for Cramer's V.

power.chisq.gof(w = 0.03022008,
                df = 4,
                alpha = 0.05,
                power = 0.80)

