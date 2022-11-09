\name{pwrss.z.prop}
\alias{pwrss.z.prop}


\title{A Proportion against a Constant (One Sample z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a proportion against a constant proportion.
}

\usage{
pwrss.z.prop(p, p0 = 0, margin = 0, arcsin.trans = TRUE, alpha = 0.05,
             alternative = c("not equal","greater","less",
                             "equivalent","non-inferior","superior"),
             n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{p}{expected proportion}
  \item{p0}{constant to be tested against (a proportion)}
  \item{arcsin.trans}{if \code{TRUE} uses arcsin transformation (default), if \code{FALSE} uses normal approximation}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin = p - p0)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}


\examples{
pwrss.z.prop(p = 0.05, p0 = 0.02, n = 190, power = NULL,
             alternative = "greater")
pwrss.z.prop(p = 0.05, p0 = 0.02, power = 0.80,
             alternative = "greater")

pwrss.z.prop(p = 0.05, p0 = 0.02, n = 190, power = NULL,
             alternative = "not equal")
pwrss.z.prop(p = 0.05, p0 = 0.02, power = 0.80,
             alternative = "not equal")

pwrss.z.prop(p = 0.05, p0 = 0.02, n = 190, power = NULL,
             margin = -0.005, alternative = "non-inferior")
pwrss.z.prop(p = 0.05, p0 = 0.02, power = 0.80,
             margin = -0.005, alternative = "non-inferior")

pwrss.z.prop(p = 0.05, p0 = 0.02, n = 190, power = NULL,
             margin = 0.005, alternative = "superior")
pwrss.z.prop(p = 0.05, p0 = 0.02, power = 0.80,
             margin = 0.005, alternative = "superior")

pwrss.z.prop(p = 0.05, p0 = 0.02, power = 0.80,
             margin = 0.005, alternative = "equivalent")
}

\references{
Add references.
}
