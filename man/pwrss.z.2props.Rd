\name{pwrss.z.2props}
\alias{pwrss.z.2props}

\title{Difference between Two Proportions (Independent Samples z Test)}

\description{
  Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two proportions.
}

\usage{
pwrss.z.2props(p1, p2, margin = 0, arcsin.trans = TRUE, kappa = 1, alpha = 0.05,
               alternative = c("not equal","greater","less",
                               "equivalent","non-inferior","superior"),
               n2 = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{p1}{expected proportion in the first group}
  \item{p2}{expected proportion in the second group}
  \item{arcsin.trans}{if \code{TRUE} uses arcsin transformation (default), if \code{FALSE} uses normal approximation}
  \item{kappa}{n1/n2}
  \item{n2}{sample size in the second group}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin = p1 - p2)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size for the first and second group}
}


\examples{
pwrss.z.2props(p1 = 0.65, p2 = 0.85,
               n2 = 55, alternative = "less")
pwrss.z.2props(p1 = 0.65, p2 = 0.85,
               power = 0.80, alternative = "less")

pwrss.z.2props(p1 = 0.65, p2 = 0.85,
               n2 = 55, alternative = "not equal")
pwrss.z.2props(p1 = 0.65, p2 = 0.85,
               power = 0.80, alternative = "not equal")

pwrss.z.2props(p1 = 0.65, p2 = 0.85, n2 = 55,
               margin = -.10, alternative = "non-inferior")
pwrss.z.2props(p1 = 0.65, p2 = 0.85, power = 0.80,
               margin = -.10, alternative = "non-inferior")

pwrss.z.2props(p1 = 0.65, p2 = 0.85, n2 = 55,
               margin = .10, alternative = "superior")
pwrss.z.2props(p1 = 0.65, p2 = 0.85, power = 0.80,
               margin = .10, alternative = "superior")

pwrss.z.2props(p1 = 0.65, p2 = 0.85, power = 0.80,
               margin = .10, alternative = "equivalent")
}

\references{
Add references.
}
