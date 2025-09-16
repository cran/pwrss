\name{power.z}
\alias{power.z.test}
\alias{power.z}

\title{Statistical Power for the Generic Z-Test}

\description{
Calculates power for the generic Z-Test with (optional) Type 1 and Type 2 error plots.
}

\usage{
power.z.test(mean = NULL, sd = 1, null.mean = 0, null.sd = 1,
             alpha = 0.05, alternative = c("two.sided",
                                           "one.sided", "two.one.sided"),
             plot = TRUE, verbose = TRUE, pretty = FALSE, ...)
}

\arguments{
  \item{mean}{mean of the alternative.}
  \item{sd}{standard deviation of the alternative. Do not change this value except when some sort of variance correction is applied (e.g. as in logistic and Poisson regressions).}
  \item{null.mean}{mean of the null. When alternative = "two.one.sided", the function expects two values in the form c(lower, upper). If a single value is provided, it is interpreted as the absolute bound and automatically expanded to c(-value, +value).}
  \item{null.sd}{standard deviation of the null. Do not change this value except when some sort of correction is applied.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; direction or type of the hypothesis test: "one.sided", "two.sided", or "two.one.sided". "two.one.sided" is used for equivalence and minimal effect testing.}
  \item{plot}{logical; \code{FALSE} switches off Type 1 and Type 2 error plot. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{...}{legacy inputs will be mapped to their corresponding arguments (silent). e.g. \code{ncp}}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
\item{mean}{mean of the alternative distribution.}
\item{sd}{standard deviation of the alternative distribution.}
\item{null.mean}{mean of the null distribution.}
\item{null.sd}{standard deviation of the null distribution.}
\item{z.alpha}{critical value(s).}
\item{power}{statistical power \eqn{(1-\beta)}.}
}

\examples{
# two-sided
# power defined as the probability of observing z-statistics
# greater than the positive critical t value OR
# less than the negative critical t value
power.z.test(mean = 1.96, alpha = 0.05,
             alternative = "two.sided")

# one-sided
# power is defined as the probability of observing z-statistics
# greater than the critical t value
power.z.test(mean = 1.96, alpha = 0.05,
             alternative = "one.sided")

# equivalence
# power is defined as the probability of observing a test statistic
# greater than the upper critical value (for the lower bound) AND
# less than the lower critical value (for the upper bound)
power.z.test(mean = 0, null.mean = c(-2, 2), alpha = 0.05,
             alternative = "two.one.sided")

# minimal effect testing
# power is defined as the probability of observing a test statistic
# greater than the upper critical value (for the upper bound) OR
# less than the lower critical value (for the lower bound).
power.z.test(mean = 2, null.mean = c(-1, 1), alpha = 0.05,
             alternative = "two.one.sided")
}
