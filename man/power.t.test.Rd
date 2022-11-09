\name{power.t.test}
\alias{power.t.test}

\title{Statistical Power for the Generic t Test with Type I and Type II Error Plots}

\description{
Calculates statistical power for the generic t test with Type I and Type II error plots.
}

\usage{
power.t.test(ncp, df, alpha, alternative,
             plot = TRUE, plot.main = NULL, plot.sub = NULL)
}

\arguments{
  \item{ncp}{non-centrality parameter (lambda)}
  \item{df}{degrees of freedom for the t test}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior". This information is needed only for power calculation}
  \item{plot}{if \code{TRUE} plots Type I and Type II error plots}
  \item{plot.main}{title}
  \item{plot.sub}{subtitle}
}

\value{
  \item{power}{statistical power \eqn{(1-\beta)}}
}

\examples{
power.t.test(ncp = 2.2, df = 100, alpha = 0.05, alternative = "not equal")
}
