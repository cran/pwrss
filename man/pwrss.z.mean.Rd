\name{pwrss.z.mean}
\alias{pwrss.z.mean}

\title{A Mean against a Constant (One Sample z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a mean against a constant mean.
}

\usage{
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
  \item{margin}{non-inferority, superiority, or equivalence margin (margin = p - p0)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}

\examples{
pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             n = 155, alternative = "greater")
pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             power = 0.80, alternative = "greater")

pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             n = 155, alternative = "not equal")
pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             power = 0.80, alternative = "not equal")

pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             n = 155, margin = -0.05, alternative = "non-inferior")
pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             power = 0.80, margin = -0.05, alternative = "non-inferior")

pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             n = 155, margin = 0.05, alternative = "superior")
pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             power = 0.80, margin = 0.05, alternative = "superior")

pwrss.z.mean(mu = 0.40, mu0 = 0.20, sd = 1,
             power = 0.80, margin = 0.05, alternative = "equivalent")

}

\references{
Add references.
}
