\name{power.f.ancova.keppel}
\alias{power.f.ancova.keppel}

\title{Power Analysis for One-Way ANOVA/ANCOVA Using Means and Standard Deviations (F test)}

\description{
Calculates power or sample size for one-way ANOVA/ANCOVA. Set \code{k.cov = 0} for one-way ANOVA (without any pretest or covariate adjustment). Set \code{k.cov > 0} in combination with  \code{r2 > 0} for one-way ANCOVA (with pretest or covariate adjustment).

Note that R has a partial matching feature which allows you to specify shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of \code{k.covariates}.

Formulas are validated using PASS documentation.
}

\usage{
power.f.ancova.keppel(mu.vector, sd.vector,
                      n.vector = NULL, p.vector = NULL,
                      factor.levels = length(mu.vector),
                      r.squared = 0, k.covariates = 0,
                      power = NULL, alpha = 0.05,
                      ceiling = TRUE, verbose = TRUE,
                      pretty = FALSE)
}

\arguments{
  \item{mu.vector}{vector of adjusted means (or estimated marginal means) for each level of a factor.}
  \item{sd.vector}{vector of unadjusted standard deviations for each level of a factor.}
  \item{n.vector}{vector of sample sizes for each level of a factor.}
  \item{p.vector}{vector of proportion of total sample size in each level of a factor. These proportions should sum to one.}
  \item{factor.levels}{integer; number of levels or groups in each factor. For example, for two factors each having two levels or groups use e.g. c(2, 2), for three factors each having two levels or groups use e.g. c(2, 2, 2)}
  \item{r.squared}{explanatory power of covariates (R-squared) in the ANCOVA model. The default is \code{r.squared = 0}, which means an ANOVA model would be of interest.}
  \item{k.covariates}{integer; number of covariates in the ANCOVA model. The default is \code{k.cov = 0}, which means an ANOVA model would be of interest.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
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
# required sample size to detect a mean difference of
# Cohen's d = 0.50 for a one-way two-group design
power.f.ancova.keppel(mu.vector = c(0.50, 0), # marginal means
                      sd.vector = c(1, 1), # unadjusted standard deviations
                      n.vector = NULL, # sample size (will be calculated)
                      p.vector = c(0.50, 0.50), # balanced allocation
                      k.cov = 1, # number of covariates
                      r.squared = 0.50, # explanatory power of covariates
                      alpha = 0.05, # Type 1 error rate
                      power = .80)

# effect size approach
power.f.ancova(eta.squared = 0.111, # effect size that is already adjusted for covariates
               factor.levels = 2, # one-way ANCOVA with two levels (groups)
               k.covariates = 1, # number of covariates
               alpha = 0.05, # Type 1 error rate
               power = .80)

# regression approach
p <- 0.50
power.t.regression(beta = 0.50,
                   sd.predictor = sqrt(p * (1 - p)),
                   sd.outcome = 1,
                   k.total = 1,
                   r.squared = 0.50,
                   n = NULL, power = 0.80)
}

\references{
Keppel, G., & Wickens, T. D. (2004). Design and analysis: A researcher's handbook (4th ed.). Pearson.
}
