\name{pwrss.z.corr}
\alias{pwrss.z.corr}


\title{A Correlation against 0 (zero) (One Sample z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a correlation against 0 (zero).
}

\usage{
pwrss.z.corr(r = 0.50, r0 = 0, alpha = 0.05,
             alternative = c("not equal","greater","less"),
             n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{r}{expected correlation}
  \item{r0}{constant proportion to be compared}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", or "less"}
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
pwrss.z.corr(r = .20, n = 194)
pwrss.z.corr(r = .20, power = .80, alternative = "not equal")
}

\references{
Add references.
}
