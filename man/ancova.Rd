\name{power.f.ancova}
\alias{power.f.ancova}
\alias{pwrss.f.ancova}

\title{Power Analysis for One-, Two-, Three-Way ANOVA/ANCOVA Using Effect Size (F-Test)}

\description{
Calculates power or sample size for one-way, two-way, or three-way ANOVA/ANCOVA. Set \code{k.cov = 0} for ANOVA, and \code{k.cov > 0} for ANCOVA. Note that in the latter, the effect size (\code{eta.squared} should be obtained from the relevant ANCOVA model, which is already adjusted for the explanatory power of covariates (thus, an additional R-squared argument is not required as an input).

Note that R has a partial matching feature which allows you to specify shortened versions of arguments, such as \code{k} or \code{k.cov} instead of \code{k.covariates}.

Formulas are validated using G*Power and tables in PASS documentation.
}

\usage{

power.f.ancova(eta.squared,
               null.eta.squared = 0,
               factor.levels = 2,
               k.covariates = 0,
               n.total = NULL,
               power = NULL,
               alpha = 0.05,
               ceiling = TRUE,
               verbose = TRUE,
               pretty = FALSE)

}

\arguments{
  \item{eta.squared}{(partial) eta-squared for the alternative.}
  \item{null.eta.squared}{(partial) eta-squared for the null.}
  \item{factor.levels}{integer; number of levels or groups in each factor. For example, for two factors each having two levels or groups use e.g. c(2, 2), for three factors each having two levels (groups) use e.g. c(2, 2, 2).}
  \item{k.covariates}{integer; number of covariates in the ANCOVA model.}
  \item{n.total}{integer; total sample size}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{ceiling}{logical; if \code{FALSE} sample size in each cell is not rounded up.}
  \item{verbose}{logical; if \code{FALSE} no output is printed on the console.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (F-Test).}
  \item{df1}{numerator degrees of freedom.}
  \item{df2}{denominator degrees of freedom.}
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null.}
  \item{f.alpha}{critical value.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n.total}{total sample size.}
}


\examples{
#############################################
#              one-way ANOVA                #
#############################################

# Cohen's d = 0.50 between treatment and control
# translating into Eta-squared = 0.059

# estimate sample size using ANOVA approach
power.f.ancova(eta.squared = 0.059,
               factor.levels = 2,
               alpha = 0.05, power = .80)

# estimate sample size using regression approach(F-Test)
power.f.regression(r.squared = 0.059,
                   k.total = 1,
                   alpha = 0.05, power = 0.80)

# estimate sample size using regression approach (T-Test)
p <- 0.50 # proportion of sample in treatment (allocation rate)
power.t.regression(beta = 0.50, r.squared = 0,
                   k.total = 1,
                   sd.predictor = sqrt(p*(1-p)),
                   alpha = 0.05, power = 0.80)

# estimate sample size using t test approach
power.t.student(d = 0.50, alpha = 0.05, power = 0.80)

#############################################
#              two-way ANOVA                #
#############################################

# a researcher is expecting a partial Eta-squared = 0.03
# for interaction of treatment (Factor A) with
# gender consisting of two levels (Factor B)

power.f.ancova(eta.squared = 0.03,
               factor.levels = c(2,2),
               alpha = 0.05, power = 0.80)

# estimate sample size using regression approach (F test)
# one dummy for treatment, one dummy for gender, and their interaction (k = 3)
# partial Eta-squared is equivalent to the increase in R-squared by adding
# only the interaction term (m = 1)
power.f.regression(r.squared = 0.03,
                   k.total = 3, k.test = 1,
                   alpha = 0.05, power = 0.80)

#############################################
#              one-way ANCOVA               #
#############################################

# a researcher is expecting an adjusted difference of
# Cohen's d = 0.45 between treatment and control after
# controllling for the pretest (k.cov = 1)
# translating into Eta-squared = 0.048

power.f.ancova(eta.squared = 0.048,
               factor.levels = 2,
               k.covariates = 1,
               alpha = 0.05, power = .80)

#############################################
#              two-way ANCOVA               #
#############################################

# a researcher is expecting an adjusted partial Eta-squared = 0.02
# for interaction of treatment (Factor A) with
# gender consisting of two levels (Factor B)

power.f.ancova(eta.squared = 0.02,
               factor.levels = c(2,2),
               k.covariates = 1,
               alpha = 0.05, power = .80)
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
