## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(pwrss)

## ----  message = FALSE, eval = FALSE, warning = FALSE-------------------------
#  install.packages("pwrss")
#  library(pwrss)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.t.test(ncp = 1.96, df = 99, alpha = 0.05,
             alternative = "equivalent", plot = TRUE)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.z.test(ncp = 1.96, alpha = 0.05, 
             alternative = "not equal", plot = TRUE)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.chisq.test(ncp = 15, df = 20,
                 alpha = 0.05, plot = TRUE)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.f.test(ncp = 3, df1 = 2, df2 = 98,
             alpha = 0.05, plot = TRUE)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
power.t.test(ncp = c(0.50, 1.00, 1.50, 2.00, 2.50), plot = FALSE,
             df = 99, alpha = 0.05, alternative = "not equal")

power.z.test(alpha = c(0.001, 0.010, 0.025, 0.050), plot = FALSE,
             ncp = 1.96, alternative = "greater")

power.chisq.test(df = c(80, 90, 100, 120, 150, 200), plot = FALSE, 
                 ncp = 2, alpha = 0.05)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = FALSE, warning = FALSE----
# comparing two means
design1 <- pwrss.t.2means(mu1 = 0.20, margin = -0.05, paired = TRUE,
                          power = 0.80, alpha = 0.05,
                          alternative = "non-inferior")
plot(design1)

# ANCOVA design
design2 <- pwrss.f.ancova(eta2 = 0.10, n.levels = c(2,3),
                          power = .80, alpha = 0.05)
plot(design2)

# indirect effect in mediation analysis
design3 <- pwrss.z.med(a = 0.10, b = 0.20, cp = 0.10,
                       power = .80, alpha = 0.05)
plot(design3)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, kappa = 1, 
               n2 = 50, alpha = 0.05,
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, kappa = 1, 
               power = .80, alpha = 0.05,
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 28, sd1 = 10.198, kappa = 1,
               power = .80, alpha = 0.05,
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 0.196, kappa = 1,
               power = .80, alpha = 0.05, 
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, 
               paired = TRUE, paired.r = 0.50,
               power = .80, alpha = 0.05,
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 28, sd1 = 10.583, 
               paired = TRUE, paired.r = 0.50,
               power = .80, alpha = 0.05,
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 0.1883, paired = TRUE, paired.r = 0.50,
               power = .80, alpha = 0.05, 
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.np.2groups(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, kappa = 1, 
               power = .80, alpha = 0.05,
               alternative = "not equal")



## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.np.2groups(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, 
               paired = TRUE, paired.r = 0.50,
               power = .80, alpha = 0.05,
               alternative = "not equal")


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, 
               margin = -1, power = 0.80,
               alternative = "non-inferior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, 
               margin = 1, power = 0.80,
               alternative = "superior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.t.2means(mu1 = 30, mu2 = 30, sd1 = 12, sd2 = 8, 
               margin = 1, power = 0.80,
               alternative = "equivalent")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.np.2groups(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, 
               margin = -1, power = 0.80,
               alternative = "non-inferior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.np.2groups(mu1 = 30, mu2 = 28, sd1 = 12, sd2 = 8, 
               margin = 1, power = 0.80,
               alternative = "superior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE, warning = FALSE----
pwrss.np.2groups(mu1 = 30, mu2 = 30, sd1 = 12, sd2 = 8, 
               margin = 1, power = 0.80,
               alternative = "equivalent")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.reg(r2 = 0.30, k = 3, power = 0.80, alpha = 0.05)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.reg(r2 = 0.15, k = 5, m = 2, power = 0.80, alpha = 0.05)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.t.reg(beta1 = 0.20, k = 3, r2 = 0.30, 
            power = .80, alpha = 0.05, alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.t.reg(beta1 = 0.60, sdy = 12, sdx = 4, k = 3, r2 = 0.30, 
            power = .80, alpha = 0.05, alternative = "not equal")


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
p <- 0.50
pwrss.t.reg(beta1 = 0.20, k = 3, r2 = 0.30, sdx = sqrt(p*(1-p)),
            power = .80, alpha = 0.05, alternative = "not equal")


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
p <- 0.50
pwrss.t.reg(beta1 = 0.20, beta0 = 0.10, margin = -0.05, 
            k = 3, r2 = 0.30, sdx = sqrt(p*(1-p)),
            power = .80, alpha = 0.05, alternative = "non-inferior")


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
p <- 0.50
pwrss.t.reg(beta1 = 0.20, beta0 = 0.10, margin = 0.05, 
            k = 3, r2 = 0.30, sdx = sqrt(p*(1-p)),
            power = .80, alpha = 0.05, alternative = "superior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
p <- 0.50
pwrss.t.reg(beta1 = 0.20, beta0 = 0.20, margin = 0.05, 
            k = 3, r2 = 0.30, sdx = sqrt(p*(1-p)),
            power = .80, alpha = 0.05, alternative = "equivalent")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.logreg(p0 = 0.15, p1 = 0.10, r2.other.x = 0.20,
               power = 0.80, alpha = 0.05, 
               dist = "normal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.logreg(p0 = 0.15, odds.ratio = 0.6296, r2.other.x = 0.20,
               alpha = 0.05, power = 0.80,
               dist = "normal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626, r2.other.x = 0.20,
               alpha = 0.05, power = 0.80,
               dist = "normal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
dist.x <- list(dist = "normal", mean = 25, sd = 8)

pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626, r2.other.x = 0.20,
               alpha = 0.05, power = 0.80,
               dist = dist.x)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626, r2.other.x = 0.20,
               alpha = 0.05, power = 0.80,
               dist = "bernoulli")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
dist.x <- list(dist = "bernoulli", prob = 0.40)

pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626, r2.other.x = 0.20,
               alpha = 0.05, power = 0.80,
               dist = dist.x)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.poisreg(beta0 = 0.50, beta1 = -0.10,
                power = 0.80, alpha = 0.05, 
                dist = "normal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                power = 0.80, alpha = 0.05, 
                dist = "normal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
dist.x <- list(dist = "normal", mean = 25, sd = 8)

pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = dist.x)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.poisreg(beta0 = 0.50, beta1 = -0.10,
                alpha = 0.05, power = 0.80,
                dist = "bernoulli")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = "bernoulli")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
dist.x <- list(dist = "bernoulli", prob = 0.40)

pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = dist.x)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# mediation model with base R-squared values
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            power = 0.80, alpha = 0.05)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# base R-squared values are 0 (zero)
# do not specify 'cp' 
pwrss.z.med(a = 0.25, b = 0.25, 
            r2m = 0, r2y = 0,
            power = 0.80, alpha = 0.05)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
p <- 0.50 # proportion of subjects in one of the groups
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10, 
            sdx = sqrt(p*(1-p)), 
            power = 0.80, alpha = 0.05)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# binary X
p <- 0.50 # proportion of subjects in one of the groups
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10, 
            sdx = sqrt(p*(1-p)), 
            n = 300, alpha = 0.05)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# continuous X
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            r2m = 0.50, r2y = 0.50, 
            power = 0.80, alpha = 0.05)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# binary X
p <- 0.50 # proportion of subjects in one of the groups
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            sdx = sqrt(p*(1-p)), r2y = 0.50, 
            power = 0.80, alpha = 0.05)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.ancova(eta2 = 0.059, n.levels = 2,
               power = .80, alpha = 0.05)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.ancova(eta2 = 0.048, n.levels = 2, n.cov = 1,
               alpha = 0.05, power = .80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.ancova(eta2 = 0.03, n.levels = c(2,2),
               alpha = 0.05, power = 0.80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.ancova(eta2 = 0.02, n.levels = c(2,2), n.cov = 1,
               alpha = 0.05, power = .80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.ancova(eta2 = 0.02, n.levels = c(2,2,3),
               alpha = 0.05, power = 0.80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.ancova(eta2 = 0.01, n.levels = c(2,2,3), n.cov = 1,
               alpha = 0.05, power = .80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.rmanova(eta2 = 0.059,  n.levels = 2, n.rm = 1,
                power = 0.80, alpha = 0.05,
                type = "between")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.rmanova(eta2 = 0.022,  n.levels = 1, n.rm = 2,
                power = 0.80, alpha = 0.05,
                corr.rm = 0.50, type = "within")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.rmanova(eta2 = 0.038,  n.levels = 2, n.rm = 2,
                power = 0.80, alpha = 0.05, 
                corr.rm = 0.50, type = "between")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.rmanova(eta2 = 0.022,  n.levels = 2, n.rm = 2,
                power = 0.80, alpha = 0.05, 
                corr.rm = 0.50, type = "within")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.f.rmanova(eta2 = 0.01,  n.levels = 2, n.rm = 2,
                power = 0.80, alpha = 0.05, 
                corr.rm = 0.50, type = "interaction")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
prob.mat <- c(0.28, 0.72) 
pwrss.chisq.gofit(p1 = c(0.28, 0.72), 
                  p0 = c(0.50, 0.50),
                  alpha = 0.05, power = 0.80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.chisq.gofit(w = 0.44, df = 1,
                  alpha = 0.05, power = 0.80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
prob.mat <- rbind(c(0.056, 0.132),
                  c(0.944, 0.868))
colnames(prob.mat) <- c("Girl", "Boy")
rownames(prob.mat) <- c("ADHD", "No ADHD")
prob.mat
pwrss.chisq.gofit(p1 = prob.mat,
                  alpha = 0.05, power = 0.80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.chisq.gofit(w = 0.1302134, df = 1,
                  alpha = 0.05, power = 0.80)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
prob.mat <- cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078),
                  c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101))
rownames(prob.mat) <- c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
colnames(prob.mat) <- c("Female", "Male")
prob.mat
pwrss.chisq.gofit(p1 = prob.mat,
                  alpha = 0.05, power = 0.80)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.chisq.gofit(w = 0.03022008, df = 4,
                  alpha = 0.05, power = 0.80)


## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.corr(r = 0.20, r0 = 0.10,
             power = 0.80, alpha = 0.05, 
             alternative = "greater")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.corr(r = 0.20, r0 = 0,
             power = 0.80, alpha = 0.05, 
             alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2corrs(r1 = 0.30, r2 = 0.20,
               power = .80, alpha = 0.05, 
               alternative = "greater")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2corrs(r1 = 0.30, r2 = 0.20,
               power = .80, alpha = 0.05, 
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# normal approximation
pwrss.z.prop(p = 0.45, p0 = 0.50,
             alpha = 0.05, power = 0.80,
             alternative = "less",
             arcsin.trans = FALSE)

# arcsine transformation
pwrss.z.prop(p = 0.45, p0 = 0.50,
             alpha = 0.05, power = 0.80,
             alternative = "less",
             arcsin.trans = TRUE)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.prop(p = 0.45, p0 = 0.50,
             alpha = 0.05, power = 0.80,
             alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.prop(p = 0.45, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "non-inferior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.prop(p = 0.55, p0 = 0.50, margin = -0.01,
             alpha = 0.05, power = 0.80,
             alternative = "non-inferior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.prop(p = 0.45, p0 = 0.50, margin = -0.01,
             alpha = 0.05, power = 0.80,
             alternative = "superior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.prop(p = 0.55, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "superior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.prop(p = 0.50, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "equivalent")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
# normal approximation
pwrss.z.2props(p1 = 0.45, p2 = 0.50,
               alpha = 0.05, power = 0.80,
               alternative = "less",
               arcsin.trans = FALSE)
# arcsine transformation
pwrss.z.2props(p1 = 0.45, p2 = 0.50,
               alpha = 0.05, power = 0.80,
               alternative = "less",
               arcsin.trans = TRUE)

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2props(p1 = 0.45, p2 = 0.50,
               alpha = 0.05, power = 0.80,
               alternative = "not equal")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2props(p1 = 0.45, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "non-inferior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2props(p1 = 0.55, p2 = 0.50,  margin = -0.01,
               alpha = 0.05, power = 0.80,
               alternative = "non-inferior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2props(p1 = 0.45, p2 = 0.50, margin = -0.01,
               alpha = 0.05, power = 0.80,
               alternative = "superior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2props(p1 = 0.55, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "superior")

## ---- message = FALSE, fig.width = 7, fig.height = 5, results = TRUE----------
pwrss.z.2props(p1 = 0.50, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "equivalent")

