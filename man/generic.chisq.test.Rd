\name{power.chisq}
\alias{power.chisq.test}
\alias{power.chisq}

\title{Statistical Power for the Generic Chi-square Test}

\description{
Calculates power for the generic chi-square test with (optional) Type 1 and Type 2 error plots.
}

\usage{
power.chisq.test(ncp, null.ncp = 0, df, alpha = 0.05,
                 plot = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null.}
  \item{df}{integer; degrees of freedom. For example, for the test of independence df = (nrow - 1)*(ncol - 1).}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{plot}{logical; \code{FALSE} switches off Type 1 and Type 2 error plot. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{power}{statistical power \eqn{(1-\beta)}.}
}

\examples{
# power is defined as the probability of observing Chi-square-statistics
# greater than the critical  value
power.chisq.test(ncp = 20, df = 100, alpha = 0.05)
}
