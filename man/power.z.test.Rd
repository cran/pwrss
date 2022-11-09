\name{power.z.test}
\alias{power.z.test}

\title{Statistical Power for the Generic z Test with Type I and Type II Error Plots}

\description{
Calculates statistical power for the generic z test with Type I and Type II error plots.
}

\usage{
power.z.test(ncp, alpha, alternative,
             plot = TRUE, plot.main = NULL, plot.sub = NULL)
}

\arguments{
  \item{ncp}{non-centrality parameter (lambda)}
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
power.z.test(ncp = 2.2, alpha = 0.05, alternative = "not equal")
}
