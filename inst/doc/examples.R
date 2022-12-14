## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(pwrss)

## ----  message = FALSE, eval = FALSE------------------------------------------
#  install.packages("pwrss")
#  library(pwrss)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = FALSE---------
power.t.test(ncp = 1.96, df = 99, alpha = 0.05, alternative = "equivalent")
power.z.test(ncp = 1.96, alpha = 0.05, alternative = "equivalent")
power.f.test(ncp = 2, df1 = 2, df2 = 98, alpha = 0.05)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = FALSE---------
# comparing two means
design1 <- pwrss.t.2means(mu1 = 0.20, margin = -0.05, paired = TRUE,
                         power = 0.80, alternative = "non-inferior")
plot(design1)

# ANCOVA design
design2 <- pwrss.f.ancova(eta2 = 0.10, n.levels = c(2,3),
                         power = .80)
plot(design2)

# indirect effect in mediation analysis
design3 <- pwrss.z.med(a = 0.10, b = 0.20, cp = 0.10,
                         power = .80)
plot(design3)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# it is sufficient to provide standardized difference between two groups for 'mu1'
# e.g. Cohen's d = 0.50 for mu1, because mu2 = 0 by default
pwrss.t.2means(mu1 = 0.50, power = .80)
pwrss.t.2means(mu1 = 0.50, power = .80, paired  = TRUE)

# it is sufficient to provide pooled standard deviation for sd1
# because sd2 = sd1 by default
pwrss.t.2means(mu1 = 10, mu2 = 5, sd1 = 10, power = .80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------

# one-way ANOVA
pwrss.f.ancova(n.levels = 2, n.covariates = 0,
               power = .80, eta2 = 0.10)
# one-way ANCOVA
pwrss.f.ancova(n.levels = 2, n.covariates = 1,
               power = .80, eta2 = 0.08)

# two-way ANOVA
pwrss.f.ancova(n.levels = c(2,3), n.covariates = 0,
               power = .80, eta2 = 0.10)
# two-way ANCOVA
pwrss.f.ancova(n.levels = c(2,3), n.covariates = 1,
               power = .80, eta2 = 0.08)

# three-way ANOVA
pwrss.f.ancova(n.levels = c(2,3,2), n.covariates = 0,
               power = .80, eta2 = 0.10)
# three-way ANCOVA
pwrss.f.ancova(n.levels = c(2,3,2), n.covariates = 1,
               power = .80, eta2 = 0.08)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------

pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                 repmeasures.r = 0.50, type = "between", power = 0.80)

pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                repmeasures.r = 0.50, epsilon = 1,
                type = "within", power = 0.80)

pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                repmeasures.r = 0.50, epsilon = 1,
                type = "interaction", power = 0.80)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------

# regression coefficient for a continuous predictor
pwrss.t.reg(beta = 0.20, r2 = 0.041, k = 1,
            power = 0.80)
pwrss.z.reg(beta = 0.20, r2 = 0.041, 
            power = 0.80)

# regression coefficient for a binary predictor
p <- 0.50 # proportion of subjects in one group
pwrss.t.reg(beta = 0.20, r2 = 0.041, k = 1,
            sdx = sqrt(p*(1-p)), power = 0.80)

# R-squared against zero
pwrss.f.reg(r2 = 0.041, k = 1,
            power = 0.80)

# R-squared difference against zero in hierarchical regression
pwrss.f.reg(r2 = 0.041, k = 5, m = 3,
            power = 0.80)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------

# indirect effect in mediation analysis
# X (cont.), M (cont.) , Z (cont.)
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            n = 200)

# X (binary), M (cont.) , Z (cont.)
p <- 0.50 # proportion of subjects in one group
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            sdx = sqrt(p*(1-p)), n = 200)

# covariate adjusted mediator and outcome model
# X (cont.), M (cont.) , Z (cont.)
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            r2m.x = 0.50, r2y.mx = 0.50, n = 200)

# covariate adjusted mediator and outcome model
# X (binary), M (cont.) , Z (cont.)
p <- 0.50 # proportion of subjects in one group
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            sdx = sqrt(p*(1-p)), 
            r2m.x = 0.50, r2y.mx = 0.50, n = 200)


