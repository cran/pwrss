\name{pwrss.t.2means}
\alias{pwrss.t.2means}

\title{Difference between Two Means (Independent Samples t Test)}

\description{
  Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two means.
  For standardized mean difference (Cohen's d) set \code{mu1 = d} and use defaults for \code{mu2}, \code{sd1}, and \code{sd2}.
  If pooled standard deviation (psd) is available set \code{sd1 = psd}.
}

\usage{
pwrss.t.2means(mu1, mu2 = 0, margin = 0,
               sd1 = ifelse(paired, sqrt(1/(2*(1-paired.r))), 1), sd2 = sd1,
               kappa = 1, paired = FALSE, paired.r = 0.50,
               alpha = 0.05, welch.df = FALSE,
               alternative = c("not equal","greater","less",
                               "equivalent","non-inferior","superior"),
               n2 = NULL, power = NULL, verbose = TRUE)

}

\arguments{
  \item{mu1}{expected mean in the first group}
  \item{mu2}{expected mean in the second group}
  \item{sd1}{standard devation in the first group}
  \item{sd2}{standard devation in the second group}
  \item{kappa}{n1/n2}
  \item{paired}{if \code{TRUE} paired samples t test}
  \item{paired.r}{correlation between the two repeated measures (e.g., pretest and posttest)}
  \item{n2}{sample size in the second group (or for the single group in paired samples)}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{welch.df}{if \code{TRUE} uses Welch degrees of freedom adjustment when groups sizes or variances are not equal (applies to independent samples t test)}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin = mu1 - mu2) }
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{df}{degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}


\examples{
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  alternative = "greater")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  alternative = "greater")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               n2 = 474,  alternative = "greater")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               power = 0.80,  alternative = "greater")

pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  alternative = "not equal")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  alternative = "not equal")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               n2 = 474,  alternative = "not equal")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               power = 0.80,  alternative = "not equal")


pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  margin = -0.10, alternative = "non-inferior")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  margin = -0.10, alternative = "non-inferior")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               n2 = 474,  margin = -0.10, alternative = "non-inferior")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               power = 0.80,  margin = -0.10, alternative = "non-inferior")


pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               n2 = 474,  margin = 0.10, alternative = "superior")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  margin = 0.10, alternative = "superior")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               n2 = 474,  margin = 0.10, alternative = "superior")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               power = 0.80,  margin = 0.10, alternative = "superior")


pwrss.t.2means(mu1 = 0.40, mu2 = 0.20,
               power = 0.80,  margin = 0.10, alternative = "equivalent")
pwrss.t.2means(mu1 = 0.40, mu2 = 0.20, paired = TRUE,
               power = 0.80,  margin = 0.10, alternative = "equivalent")
}

\references{
Add references.
}
