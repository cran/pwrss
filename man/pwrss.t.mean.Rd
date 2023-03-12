\name{One Mean (z and t Test)}
\alias{pwrss.z.mean}
\alias{pwrss.t.mean}

\title{A Mean against a Constant (z or t Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a mean against a constant mean.

Formulas are validated using \url{http://powerandsamplesize.com/}, and tables in PASS documentation.
}

\usage{
pwrss.t.mean(mu, sd = 1, mu0 = 0, margin = 0, alpha = 0.05,
             alternative = c("not equal","greater","less",
                             "equivalent","non-inferior","superior"),
             n = NULL, power = NULL, verbose = TRUE)

pwrss.z.mean(mu, sd = 1, mu0 = 0, margin = 0, alpha = 0.05,
             alternative = c("not equal","greater","less",
                             "equivalent","non-inferior","superior"),
             n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{mu}{expected mean}
  \item{sd}{expected standard devation}
  \item{mu0}{constant to be tested against (a mean)}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin: boundry of \code{mu - mu0} that is practically insignificant)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{df}{degrees of freedom}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}

\examples{
# Example: A researcher is expecting a score of 23
# on Beck depression inventory (BDI) which is
# 0.50 standard devation above the threshold value 20
# (assume standard deviation of BDI scores is 6).

# to find that a score of 23 is greater than the threshold 20
pwrss.t.mean(mu = 23, mu0 = 20, sd = 6,
             alpha = 0.05, power = 0.80,
             alternative = "greater")
# standardized formulation
pwrss.t.mean(mu = 0.50, mu0 = 0, sd = 1,
             alpha = 0.05, power = 0.80,
             alternative = "greater")
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Under review. \url{https://osf.io/ua5fc}
}
