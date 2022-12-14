\name{pwrss.f.reg}
\alias{pwrss.f.reg}


\title{R-squared Deviation from 0 (zero) or R-squared Difference in Hierarchical Multiple Linear Regression (F Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test R-squared Deviation from 0 (zero) or to test R-squared difference from 0 (zero) in hierarchical regression.
}

\usage{
pwrss.f.reg(r2 = 0.10, f2 = r2 /(1 - r2),
            k = 4, m = k, alpha = 0.05,
            n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{r2}{Expected R-squared}
  \item{f2}{Cohen's f-squared (an alternative to \code{r2} specification). f2 = r2 / (1 - r2)}
  \item{k}{(total) number of predictors}
  \item{m}{(number of predictors in the subset of interest}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
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
# R-squared deviation from 0 (zero)
pwrss.f.reg(r2 = 0.20, k = 5, n = 26)
pwrss.f.reg(r2 = 0.20, k = 5, power = 0.80)

# Model 1 = 3 predictors
# Model 2 = 2 predictors
# k = 5 (total number of predictors)
# m = 2 (predictors whose incrimental R-squared change to be inspected)
pwrss.f.reg(r2 = 0.15, k = 5, m = 2, n = 26)
pwrss.f.reg(r2 = 0.15, k = 5, m = 2, power = 0.80)
}

\references{
Add references.
}
