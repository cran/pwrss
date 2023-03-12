\name{Generic F Test}
\alias{power.f.test}

\title{Statistical Power for the Generic F Test with (optional) Type I and Type II Error Plots}

\description{
Calculates statistical power for the generic F test with (optional) Type I and Type II error plots.
Unlike other more specific functions \code{power.f.test()} function allows specification of multiple values for one parameter at a time (only when \code{plot = FALSE}).
}

\usage{
power.f.test(ncp, df1, df2, alpha = 0.05, plot = TRUE,
             plot.main = NULL, plot.sub = NULL,
             verbose = TRUE)
}

\arguments{
  \item{ncp}{non-centrality parameter (lambda)}
  \item{alpha}{probability of type I error}
  \item{df1}{numerator degrees of freedom for the F test}
  \item{df2}{denominator degrees of freedom for the F test}
  \item{plot}{if \code{TRUE} plots Type I and Type II error plots}
  \item{plot.main}{plot title}
  \item{plot.sub}{plot subtitle}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{power}{statistical power \eqn{(1-\beta)}}
}

\examples{
# power is defined as the probability of observing F-statistics
# greater than the critical F value
power.f.test(ncp = 1, df1 = 4, df2 = 100, alpha = 0.05)

# power of multiple F-statistics
power.f.test(ncp = c(1.0, 1.5, 2.0, 2.5), plot = FALSE,
             df1 = 4, df2 = 100, alpha = 0.05)

}
