\name{power.binom}
\alias{power.binom.test}
\alias{power.binom}

\title{Power Analysis for the Generic Binomial Test}

\description{
  Calculates power for the generic binomial test with (optional) Type 1 and Type 2 error plots.
}

\usage{
power.binom.test(size, prob, null.prob = 0.5, alpha = 0.05,
                 alternative = c("two.sided", "one.sided", "two.one.sided"),
                 plot = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{size}{number of trials (zero or more).}
  \item{prob}{probability of success on each trial under alternative.}
  \item{null.prob}{probability of success on each trial under null.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{direction or type of the hypothesis test: "two.sided", "one.sided", or "two.one.sided". For non-inferiority or superiority tests, add or subtract the margin from the null hypothesis value and use alternative = "one.sided".}
  \item{plot}{logical; \code{FALSE} switches off Type 1 and Type 2 error plot. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{size}{number of trials (zero or more).}
  \item{prob}{probability of success on each trial under alternative.}
  \item{null.prob}{probability of success on each trial under null.}
  \item{binom.alpha}{critical value(s).}
  \item{power}{statistical power \eqn{(1-\beta)}.}
}

\examples{
# one-sided
power.binom.test(size = 200, prob = 0.6, null.prob = 0.5,
                 alpha = 0.05, alternative = "one.sided")

# two-sided
power.binom.test(size = 200, prob = 0.4, null.prob = 0.5,
                 alpha = 0.05, alternative = "two.sided")

# equivalence
power.binom.test(size = 200, prob = 0.5, null.prob = c(0.4, 0.6),
                 alpha = 0.05, alternative = "two.one.sided")
}

