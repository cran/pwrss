\name{pwrss.z.2means}
\alias{pwrss.z.2means}

\title{Difference between Two Means (Independent Samples z Test)}

\description{
  Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two means.
  For standardized mean difference (Cohen's d) set \code{mu1 = d} and use defaults for \code{mu2}, \code{sd1}, and \code{sd2}.
  If pooled standard deviation (psd) is available set \code{sd1 = psd}.
}

\usage{
pwrss.z.2means(mu1, mu2 = 0, sd1 = 1, sd2 = sd1, margin = 0,
               kappa = 1, alpha = 0.05,
               alternative = c("not equal","greater","less",
                               "equivalent","non-inferior","superior"),
               n2 = NULL, power = NULL, verbose = TRUE)

}

\arguments{
  \item{mu1}{expected mean in the first group}
  \item{mu2}{expected mean in the second group}
  \item{sd1}{pooled standard devation in the first group}
  \item{sd2}{pooled standard devation in the second group}
  \item{kappa}{n1/n2}
  \item{n2}{sample size in the second group}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin = mu1 - mu2) }
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size for the first and second group}
}


\examples{
pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  alternative = "greater")
pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  alternative = "greater")

pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  alternative = "not equal")
pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  alternative = "not equal")


pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  margin = -0.10, alternative = "non-inferior")
pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  margin = -0.10, alternative = "non-inferior")

pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  margin = 0.10, alternative = "superior")
pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  margin = 0.10, alternative = "superior")


pwrss.z.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  margin = 0.10, alternative = "equivalent")
}

\references{
Add references.
}
