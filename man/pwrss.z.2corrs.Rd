\name{pwrss.z.2corrs}
\alias{pwrss.z.2corrs}


\title{Difference between Two Correlations (Independent Samples z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two independent correlations.
}

\usage{
pwrss.z.2corrs(r1 = 0.50, r2 = 0.30,
               alpha = 0.05, kappa = 1,
               alternative = c("not equal","greater","less"),
               n2 = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{r1}{expected correlation in the first group}
  \item{r2}{expected correlation in the second group}
  \item{n2}{sample size in the second group}
  \item{kappa}{n1/n2}
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
pwrss.z.2corrs(r1 = .20, r2 = 0.30, n = 194)
pwrss.z.2corrs(r1 = .20, r2 = 0.30, power = .80)
}

\references{
Add references.
}
