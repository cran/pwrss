\name{pwrss.f.rmanova}
\alias{pwrss.f.rmanova}


\title{Repeated Measures Analysis of Variance (F test)}

\description{
Calculates statistical power or minimum required sample size for one-way Repeated Measures Analysis of Variance (RM-ANOVA).

Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.
}

\usage{
pwrss.f.rmanova(eta2 = 0.10, f2 = eta2/(1 - eta2),
                corr.rm = 0.50, n.levels = 2, n.rm = 2,
                epsilon = 1, alpha = 0.05,
                type = c("between","within","interaction"),
                n = NULL, power = NULL, verbose = TRUE)

}

\arguments{
  \item{eta2}{expected (partial) Eta-squared}
  \item{f2}{expected Cohen's f-squared (an alternative to \code{eta2} specification). \code{f2 = eta2 / (1 - eta2)}}
  \item{corr.rm}{expected correlation between repeated measures. For example, for pretest/posttest designs, this is the correlation between pretest and posttest scores regardless of group membership. The default is 0.50}
  \item{n.levels}{number of levels (groups). For example, for randomized controlled trials with two arms (treatment/control) it takes a value of 2}
  \item{n.rm}{number of measurements. For example, for pretest/posttest designs it takes a value of 2. When there is a follow-up test it takes a value of 3}
  \item{epsilon}{non-sperhicity correction factor, default is 1 (means no violation of sphericity). Lower bound for this argument is \code{epsilon} = 1 / (\code{n.rm} - 1)}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{type}{the effect to be tested: "between", "within", or "interaction". The type of the effect depends on the hypothesis test. If the interest is in the group effect after controlling for the time effect use "between"; if the interest is the time effect after controlling for the group membership use "within"; if the interest is in the group x time interaction use "interaction"}
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
######################################################
# pretest-posttest design with treatment group only  #
######################################################

# a researcher is expecting a difference of Cohen's d = 0.30
# between posttest and pretest score translating into
# Eta-squared = 0.022
pwrss.f.rmanova(eta2 = 0.022,  n.levels = 1, n.rm = 2,
                corr.rm = 0.50, type = "within",
                alpha = 0.05, power = 0.80)

# paired t-test approach
pwrss.t.2means(mu1 = 0.30, mu2 = 0,
               sd1 = 1, sd2 = 1,
               paired = TRUE, paired.r = 0.50,
               alpha = 0.05, power = 0.80)

##########################################################
# posttest only design with treatment and control groups #
##########################################################

# a researcher is expecting a difference of Cohen's d = 0.50
# on the posttest score between treatment and control groups
# translating into Eta-squared = 0.059
pwrss.f.rmanova(eta2 = 0.059,  n.levels = 2, n.rm = 1,
                type = "between",
                alpha = 0.05, power = 0.80)

# independent t-test approach
pwrss.t.2means(mu1 = 0.50, mu2 = 0,
               sd1 = 1, sd2 = 1,
               alpha = 0.05, power = 0.80)

#############################################################
# pretest-posttest design with treatment and control groups #
#############################################################

# a researcher is expecting a difference of Cohen's d = 0.40
# on the posttest score between treatment and control groups
# after controlling for the pretest translating into
# partial Eta-squared = 0.038
pwrss.f.rmanova(eta2 = 0.038,  n.levels = 2, n.rm = 2,
                corr.rm = 0.50, type = "between",
                alpha = 0.05, power = 0.80)

# regression approach
p <- 0.50 # proportion of subjects in treatment group
pwrss.t.reg(beta1 = 0.40, r2 = 0.25, k = 2,
            sdx = sqrt(p*(1-p)),
            alpha = 0.05, power = 0.80)

# a researcher is expecting an interaction effect
# (between groups and time) of Eta-squared = 0.01
pwrss.f.rmanova(eta2 = 0.01,  n.levels = 2, n.rm = 2,
                corr.rm = 0.50, type = "interaction",
                alpha = 0.05, power = 0.80)
}

\references{
Bulus, M., & Polat, C. (in press). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi. \url{https://osf.io/ua5fc/download/}
}
