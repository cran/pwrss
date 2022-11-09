\name{power.f.test}
\alias{power.f.test}

\title{Statistical Power for the Generic F Test with Type I and Type II Error Plots}

\description{
Calculates statistical power for the generic F test with Type I and Type II error plots.
}

\usage{
power.f.test(ncp, df1, df2, alpha,
             plot = TRUE, plot.main = NULL, plot.sub = NULL)
}

\arguments{
  \item{ncp}{non-centrality parameter (lambda)}
  \item{alpha}{probability of type I error}
  \item{df1}{numerator degrees of freedom for the F test}
  \item{df2}{denominator degrees of freedom for the F test}
  \item{plot}{if \code{TRUE} plots Type I and Type II error plots}
  \item{plot.main}{title}
  \item{plot.sub}{subtitle}
}

\value{
  \item{power}{statistical power \eqn{(1-\beta)}}
}

\examples{
power.f.test(ncp = 1, df1 = 4, df2 = 100, alpha = 0.05)
}
