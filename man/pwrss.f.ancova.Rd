\name{pwrss.f.ancova}
\alias{pwrss.f.anova}
\alias{pwrss.f.ancova}


\title{Analysis of Variance (ANOVA) or Analysis of Covariance (ANCOVA) (F test)}

\description{
Calculates statistical power or minimum required sample size for One-way/Two-way/Three-way Analysis of Variance (ANOVA) and Analysis of Covariance (ANCOVA). Set \code{n.covariates = 0} for ANOVA, and \code{n.covariates > 0} for ANCOVA. Note that in each case, the effect size (\code{eta2} or \code{f2}) should be obtained from the relevant model.
}

\usage{
pwrss.f.ancova(eta2 = 0.01, f2 = eta2 / (1 - eta2),
               n.way = length(n.levels),
               n.levels = 2, n.covariates = 0, alpha = 0.05,
               n = NULL, power = NULL, verbose = TRUE)

}

\arguments{
  \item{eta2}{Expected Eta-squared}
  \item{f2}{Cohen's f2 (an alternative to \code{eta2} specification). f2 = eta2 / (1 - eta2)}
  \item{n.way}{1 for one-way, 2 for two-way, 3 for three-way ANOVA or ANCOVA. The default takes its value from the length of \code{n.levels}}
  \item{n.levels}{number of levels (groups) in each factor. For two factors each having two levels (groups) use e.g. c(2,2), for three factors each having two levels (groups) use e.g. c(2,2,2)}
  \item{n.covariates}{number of covariates in the ANCOVA model}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{df1}{numerator degrees of freedom}
  \item{df2}{denominator degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}


\examples{
# one-way ANOVA
pwrss.f.ancova(n.levels = 2,
               n = 26, eta2 = 0.10)
pwrss.f.ancova(n.levels = 2,
               power = .80, eta2 = 0.10)

# two-way ANOVA
pwrss.f.ancova(n.levels = c(2,3),
               n = 26, eta2 = 0.10)
pwrss.f.ancova(n.levels = c(2,3),
               power = .80, eta2 = 0.10)

# three-way ANOVA
pwrss.f.ancova(n.levels = c(2,3,2),
               n = 26, eta2 = 0.10)
pwrss.f.ancova(n.levels = c(2,3,2),
               power = .80, eta2 = 0.10)

# one-way ANCOVA
pwrss.f.ancova(n.levels = 2, n.covariates = 1,
               n = 26, eta2 = 0.08)
pwrss.f.ancova(n.levels = 2, n.covariates = 1,
               power = .80, eta2 = 0.08)

# two-way ANCOVA
pwrss.f.ancova(n.levels = c(2,3), n.covariates = 1,
               n = 26, eta2 = 0.10)
pwrss.f.ancova(n.levels = c(2,3), n.covariates = 1,
               power = .80, eta2 = 0.08)

# three-way ANCOVA
pwrss.f.ancova(n.levels = c(2,3,2), n.covariates = 1,
               n = 26, eta2 = 0.10)
pwrss.f.ancova(n.levels = c(2,3,2), n.covariates = 1,
               power = .80, eta2 = 0.08)
}

\references{
Add references.
}
