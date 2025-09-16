\name{power.f.mixed.anova}
\alias{power.f.mixed.anova}
\alias{pwrss.f.rmanova}


\title{Power Analysis for Mixed-Effects Analysis of Variance (F-Test)}

\description{
  Calculates power or sample size for mixed-effects ANOVA design with two factors (between and within). When there is only one group observed over time, this design is often referred to as repeated-measures ANOVA.

  Formulas are validated using G*Power and tables in PASS documentation.

  NOTE: The \code{pwrss.f.rmanova()} function is deprecated and will no longer be supported, but it will remain available as a wrapper for \code{power.f.mixed.anova()} during the transition period.
}

\usage{
 power.f.mixed.anova(eta.squared, null.eta.squared = 0,
                    factor.levels = c(2, 2),
                    factor.type = c("between", "within"),
                    rho.within = 0.50, epsilon = 1,
                    n.total = NULL, power = NULL, alpha = 0.05,
                    effect = c("between", "within", "interaction"),
                    ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{eta.squared}{(partial) eta-squared for the alternative.}
  \item{null.eta.squared}{(partial) eta-squared for the null.}
  \item{rho.within}{Correlation between repeated measures. For example, for pretest/post-test designs, this is the correlation between pretest and post-test scores regardless of group membership. The default is 0.50. If \code{eta.squared} is already adjusted for this correlation specify 'rho.within = NA'.}
  \item{factor.levels}{vector; integer; length of two representing the number of levels for groups and measures. For example, in randomized controlled trials with two arms (treatment and control) where pre-test, post-test, and follow-up test are administered, this would be represented as c(2, 3).}
  \item{factor.type}{vector; character; length of two indicating the order of between-subject and within-subject factors. By default, the first value represents the between-subject factor and the second value represents the within-subject factor. This argument is rarely needed, except when unsure which element in 'factor.levels' represent between-subject or within-subject factors. Therefore, specify the 'factor.levels' accordingly.}
  \item{epsilon}{non-sphericity correction factor, default is 1 (means no violation of sphericity). Lower bound for this argument is \code{epsilon} = 1 / (\code{factor.levels[2]} - 1).}
  \item{n.total}{integer; total sample size.}
 \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{effect}{character; the effect of interest: "between", "within", or "interaction".}
  \item{ceiling}{logical;  \code{TRUE} by default. If \code{FALSE} sample size in each group is NOT rounded up.}
  \item{verbose}{logical; \code{TRUE} by default. If \code{FALSE} no output is printed on the console.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (F-Test).}
  \item{df1}{numerator degrees of freedom.}
  \item{df2}{denominator degrees of freedom.}
  \item{ncp}{non-centrality parameter under alternative.}
  \item{null.ncp}{non-centrality parameter under null.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n.total}{total sample size.}
}


\examples{

######################################################
# pretest-post-test design with treatment group only  #
######################################################

# a researcher is expecting a difference of Cohen's d = 0.30
# between post-test and pretest score translating into
# Eta-squared = 0.022

# adjust effect size for correlation with 'rho.within'
power.f.mixed.anova(eta.squared = 0.022,
                    factor.levels = c(1, 2), # 1 between 2 within
                    rho.within = 0.50,
                    effect = "within",
                    alpha = 0.05, power = 0.80)

# if effect size is already adjusted for correlation
# use 'rho.within = NA'
power.f.mixed.anova(eta.squared = 0.08255,
                    factor.levels = c(1, 2), # 1 between 2 within
                    rho.within = NA,
                    effect = "within",
                    alpha = 0.05, power = 0.80)

##########################################################
# post-test only design with treatment and control groups #
##########################################################

# a researcher is expecting a difference of Cohen's d = 0.50
# on the post-test score between treatment and control groups
# translating into Eta-squared = 0.059
power.f.mixed.anova(eta.squared = 0.059,
                    factor.levels = c(2, 1),  # 2 between 1 within
                    effect = "between",
                    alpha = 0.05, power = 0.80)


#############################################################
# pretest-post-test design with treatment and control groups #
#############################################################

# a researcher is expecting a difference of Cohen's d = 0.40
# on the post-test score between treatment and control groups
# after controlling for the pretest translating into
# partial Eta-squared = 0.038
power.f.mixed.anova(eta.squared = 0.038,
                    factor.levels = c(2, 2),  # 2 between 2 within
                    rho.within = 0.50,
                    effect = "between",
                    alpha = 0.05, power = 0.80)

# a researcher is expecting an interaction effect
# (between groups and time) of Eta-squared = 0.01
power.f.mixed.anova(eta.squared = 0.01,
                    factor.levels = c(2, 2),  # 2 between 2 within
                    rho.within = 0.50,
                    effect = "interaction",
                    alpha = 0.05, power = 0.80)

# a researcher is expecting an interaction effect
# (between groups and time) of Eta-squared = 0.01
power.f.mixed.anova(eta.squared = 0.01,
                    factor.levels = c(2, 2),  # 2 between 2 within
                    rho.within = 0.50,
                    effect = "within",
                    alpha = 0.05, power = 0.80)
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}
}
