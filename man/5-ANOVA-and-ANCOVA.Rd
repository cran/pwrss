\name{pwrss.f.ancova}
\alias{pwrss.f.anova}
\alias{pwrss.f.ancova}


\title{Analysis of (Co)Variance (F test)}

\description{
Calculates statistical power or minimum required sample size for one-way, two-way, or three-way ANOVA/ANCOVA. Set \code{n.covariates = 0} for ANOVA, and \code{n.covariates > 0} for ANCOVA. Note that in each case, the effect size (partial) (\code{eta2} or \code{f2}) should be obtained from the relevant model.

Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.
}

\usage{
pwrss.f.ancova(eta2 = 0.01, f2 = eta2 / (1 - eta2),
               n.way = length(n.levels),
               n.levels = 2, n.covariates = 0, alpha = 0.05,
               n = NULL, power = NULL, verbose = TRUE)

}

\arguments{
  \item{eta2}{expected Eta-squared}
  \item{f2}{expected Cohen's f2 (an alternative to \code{eta2} specification). f2 = eta2 / (1 - eta2)}
  \item{n.way}{1 for one-way, 2 for two-way, 3 for three-way ANOVA or ANCOVA. The default takes its value from the length of \code{n.levels}}
  \item{n.levels}{number of levels (groups) in each factor. For example, for two factors each having two levels (groups) use e.g. c(2,2), for three factors each having two levels (groups) use e.g. c(2,2,2)}
  \item{n.covariates}{number of covariates in the ANCOVA model}
  \item{n}{total sample size}
  \item{alpha}{probability of type I error}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{verbose}{if \code{FALSE} no output is printed on the console}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (F test)}
  \item{df1}{numerator degrees of freedom}
  \item{df2}{denominator degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{total sample size}
}


\examples{
#############################################
#              one-way ANOVA                #
#############################################

# a researcher is expecting a difference of
# Cohen's d = 0.50 between treatment and control
# translating into Eta-squared = 0.059

# estimate sample size using ANOVA approach
pwrss.f.ancova(eta2 = 0.059, n.levels = 2,
               alpha = 0.05, power = .80)

# estimate sample size using regression approach(F test)
pwrss.f.reg(r2 = 0.059, k = 1,
            alpha = 0.05, power = 0.80)

# estimate sample size using regression approach (t test)
p <- 0.50 # proportion of sample in treatment
pwrss.t.reg(beta1 = 0.50, r2 = 0,
            k = 1, sdx = sqrt(p*(1-p)),
            alpha = 0.05, power = 0.80)

# estimate sample size using t test approach
pwrss.t.2means(mu1 = 0.50,
               alpha = 0.05, power = 0.80)

#############################################
#              two-way ANOVA                #
#############################################

# a researcher is expecting a partial Eta-squared = 0.03
# for interaction of treatment (Factor A) with
# gender consisting of two levels (Factor B)

pwrss.f.ancova(eta2 = 0.03, n.levels = c(2,2),
               alpha = 0.05, power = 0.80)

# estimate sample size using regression approach (F test)
# one dummy for treatment, one dummy for gender, and their interaction (k = 3)
# partial Eta-squared is equivalent to the increase in R-squared by adding
# only the interaction term (m = 1)
pwrss.f.reg(r2 = 0.03, k = 3, m = 1,
            alpha = 0.05, power = 0.80)

#############################################
#              one-way ANCOVA               #
#############################################

# a researcher is expecting an adjusted difference of
# Cohen's d = 0.45 between treatment and control after
# controllling for the pretest (n.cov = 1)
# translating into Eta-squared = 0.048

pwrss.f.ancova(eta2 = 0.048, n.levels = 2, n.cov = 1,
               alpha = 0.05, power = .80)

#############################################
#              two-way ANCOVA               #
#############################################

# a researcher is expecting an adjusted partial Eta-squared = 0.02
# for interaction of treatment (Factor A) with
# gender consisting of two levels (Factor B)

pwrss.f.ancova(eta2 = 0.02, n.levels = c(2,2), n.cov = 1,
               alpha = 0.05, power = .80)
}

\references{
Bulus, M., & Polat, C. (in press). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi. \url{https://osf.io/ua5fc/download/}

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
