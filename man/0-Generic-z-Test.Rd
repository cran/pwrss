\name{power.z.test}
\alias{power.z.test}

\title{Statistical Power for the Generic z Test}

\description{
Calculates statistical power for the generic z test with (optional) Type I and Type II error plots.
Unlike other more specific functions \code{power.z.test()} function allows multiple values for one parameter at a time (only when \code{plot = FALSE}).
}

\usage{
power.z.test(ncp, alpha = 0.05,
             alternative = c("not equal", "greater", "less",
                             "non-inferior", "superior", "equivalent"),
             plot = TRUE, plot.main = NULL, plot.sub = NULL,
             verbose = TRUE)
}

\arguments{
  \item{ncp}{non-centrality parameter (lambda)}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior". The same non-centrality parameters will produce the same power rates for "greater", "less", "non-inferior", and "superior" tests. Different labels have been used for consistency. However, it should be noted that the non-centrality parameter should conform to the specific test type}
  \item{plot}{if \code{TRUE} plots Type I and Type II error}
  \item{plot.main}{plot title}
  \item{plot.sub}{plot subtitle}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{power}{statistical power \eqn{(1-\beta)}}
}

\examples{
# power defined as the probability of observing z-statistics
# greater than the positive critical t value OR
# less than the negative critical t value
power.z.test(ncp = 1.96, alpha = 0.05,
             alternative = "not equal")

# power is defined as the probability of observing z-statistics
# greater than the critical t value
power.z.test(ncp = 1.96, alpha = 0.05,
             alternative = "greater")

# power is defined as the probability of observing z-statistics
# greater than the critical t value where the non-centrality parameter
# for the alternative distribution is adjusted for the non-inferiority margin
power.z.test(ncp = 1.98, alpha = 0.05,
             alternative = "non-inferior")

# power is defined as the probability of observing z-statistics
# greater than the critical t value where the non-centrality parameter
# for the alternative distribution is adjusted for the superiority margin
power.z.test(ncp = 1.94, alpha = 0.05,
             alternative = "superior")

# power is defined as the probability of observing z-statistics
# less than the positive critical t value AND
# greater than the negative critical t value
# the non-centrality parameter is for the null distribution
# and is derived from the equivalence margins (lower and upper)
power.z.test(ncp = 1.96, alpha = 0.05,
             alternative = "equivalent")
# or, define lower and upper bound with rbind()
power.z.test(ncp = rbind(-1.96, 1.96), alpha = 0.05,
             alternative = "equivalent")
}
