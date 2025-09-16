\name{power.t}
\alias{power.t.test}
\alias{power.t}

\title{Statistical Power for the Generic T-Test}

\description{
Calculates power for the generic T-Test with (optional) Type 1 and Type 2 error plots.
}

\usage{
power.t.test(ncp, null.ncp = 0, df, alpha = 0.05,
             alternative = c("two.sided", "one.sided", "two.one.sided"),
             plot = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null. When alternative = "two.one.sided", the function expects two values in the form c(lower, upper). If a single value is provided, it is interpreted as the absolute bound and automatically expanded to c(-value, +value).}
  \item{df}{degrees of freedom.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; direction or type of the hypothesis test: "one.sided", "two.sided", or "two.one.sided". "two.one.sided" is used for equivalence and minimal effect testing.}
  \item{plot}{logical; \code{FALSE} switches off Type 1 and Type 2 error plot. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
\item{df}{degrees of freedom.}
\item{ncp}{non-centrality parameter under alternative.}
\item{ncp.null}{non-centrality parameter under null.}
\item{t.alpha}{critical value(s).}
\item{power}{statistical power \eqn{(1-\beta)}.}
}

\examples{
# two-sided
# power defined as the probability of observing a test statistic
# greater than the positive critical value OR
# less than the negative critical value
power.t.test(ncp = 1.96, df = 100, alpha = 0.05,
             alternative = "two.sided")

# one-sided
# power is defined as the probability of observing a test statistic
# greater than the critical value
power.t.test(ncp = 1.96, df = 100, alpha = 0.05,
             alternative = "one.sided")

# equivalence
# power is defined as the probability of observing a test statistic
# greater than the upper critical value (for the lower bound) AND
# less than the lower critical value (for the upper bound)
power.t.test(ncp = 0, df = 100,
             null.ncp = c(-2, 2), alpha = 0.05,
             alternative = "two.one.sided")

# minimal effect testing
# power is defined as the probability of observing a test statistic
# greater than the upper critical value (for the upper bound) OR
# less than the lower critical value (for the lower bound).
power.t.test(ncp = 2, df = 100,
             null.ncp = c(-1, 1), alpha = 0.05,
             alternative = "two.one.sided")
}
