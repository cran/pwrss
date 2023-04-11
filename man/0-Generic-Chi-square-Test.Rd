\name{power.chisq.test}
\alias{power.chisq.test}

\title{Statistical Power for the Generic Chi-square Test}

\description{
Calculates statistical power for the generic chi-square test with (optional) Type I and Type II error plots.
Unlike other more specific functions \code{power.chisq.test()} function allows multiple values for one parameter at a time (only when \code{plot = FALSE}).
}

\usage{
power.chisq.test(ncp, df, alpha = 0.05, plot = TRUE,
                 plot.main = NULL, plot.sub = NULL,
                 verbose = TRUE)
}

\arguments{
  \item{ncp}{non-centrality parameter (lambda)}
  \item{df}{degrees of freedom. For example, for the test of homogeneity or independence df = (nrow - 1)*(ncol - 1)}
  \item{alpha}{probability of type I error}
  \item{plot}{if \code{TRUE} plots Type I and Type II error}
  \item{plot.main}{plot title}
  \item{plot.sub}{plot subtitle}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{power}{statistical power \eqn{(1-\beta)}}
}

\examples{
# power is defined as the probability of observing Chi-square-statistics
# greater than the critical Chi-square value
power.chisq.test(ncp = 20, df = 100, alpha = 0.05)

# power of multiple Chi-square-statistics
power.chisq.test(ncp = c(5, 10, 15, 20), plot = FALSE,
                 df = 100, alpha = 0.05)

}
