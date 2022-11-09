\name{pwrss.f.rmanova}
\alias{pwrss.f.rmanova}


\title{Repeated Measures Analysis of Variance (RM-ANOVA) (F test)}

\description{
Calculates statistical power or minimum required sample size for one-way Repeated Measures Analysis of Variance (RM-ANOVA).
}

\usage{
pwrss.f.rmanova(eta2 = 0.10, f2 = eta2/(1 - eta2),
                repmeasures.r = 0.50, n.levels = 2, n.measurements = 2,
                epsilon = 1, alpha = 0.05,
                type = c("between","within","interaction"),
                n = NULL, power = NULL, verbose = TRUE)

}

\arguments{
  \item{eta2}{Expected Eta-squared}
  \item{f2}{Cohen's f2 (an alternative to \code{eta2} specification). f2 = eta2 / (1 - eta2)}
  \item{repmeasures.r}{correlation between repeated measures}
  \item{n.levels}{number of levels (groups)}
  \item{n.measurements}{number of measurements}
  \item{epsilon}{non-sperhicity correction f, default 1 (means no violation of sphericity). Lower bound for \code{epsilon} = 1 / (\code{n.measurements} - 1)}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{type}{the effect to be tested: "between", "within", or "interaction"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{df1}{numerator degrees of freedom}
  \item{df2}{denominator degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}


\examples{
pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                type = "between", n = 10)
pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                type = "between", power = 0.80)

pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                type = "within", n = 10)
pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                type = "within", power = 0.80)

pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                type = "interaction", n = 10)
pwrss.f.rmanova(eta2 = 0.10,  n.levels = 2, n.measurements = 3,
                type = "interaction", power = 0.80)
}

\references{
Add references.
}
