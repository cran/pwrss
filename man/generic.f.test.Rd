\name{power.f}
\alias{power.f.test}
\alias{power.f}

\title{Statistical Power for the Generic F-Test}

\description{
Calculates power for the generic F-Test with (optional) Type 1 and Type 2 error plots.
}

\usage{
power.f.test(ncp, null.ncp = 0, df1, df2, alpha = 0.05,
             plot = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{df1}{integer; numerator degrees of freedom.}
  \item{df2}{integer; denominator degrees of freedom.}
  \item{plot}{logical; \code{FALSE} switches off Type 1 and Type 2 error plot. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{df1}{numerator degrees of freedom.}
  \item{df2}{denominator degrees of freedom.}
  \item{ncp}{non-centrality parameter under alternative.}
  \item{ncp.null}{non-centrality parameter under null.}
  \item{f.alpha}{critical value(s).}
  \item{power}{statistical power \eqn{(1-\beta)}.}
}

\examples{
# power is defined as the probability of observing F-statistics
# greater than the critical value
power.f.test(ncp = 1, df1 = 4, df2 = 100, alpha = 0.05)
}
