\name{Linear Regression (F Test)}
\alias{pwrss.f.reg}
\alias{pwrss.f.regression}


\title{Test of R-squared Deviation from 0 (zero) or R-squared Difference in Hierarchical Multiple Linear Regression (F Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test R-squared deviation from 0 (zero) in linear regression or to test R-squared difference from 0 (zero) in hierarchical regression.

Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.
}

\usage{
pwrss.f.reg(r2 = 0.10, f2 = r2 /(1 - r2),
            k = 4, m = k, alpha = 0.05,
            n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{r2}{expected R-squared (or R-squared change)}
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
# EXample 1: A researcher is expecting that
# three variables together explain 15 percent of the variance
# in the outcome (R-squared = 0.15).
pwrss.f.reg(r2 = 0.15, k = 3,
            alpha = 0.05, power = 0.80)

# Example 2: A researcher is expecting that
# adding two more variables will increase R-squared
# from 0.15 (with 3 predictors) to 0.20 (with 3 + 2 predictors)
# k = 5 (total number of predictors)
# m = 2 (predictors whose incremental contribution to R-squared change is of interest)
pwrss.f.reg(r2 = 0.05, k = 5, m = 2,
            alpha = 0.05, power = 0.80)
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Under review. \url{https://osf.io/ua5fc}

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
