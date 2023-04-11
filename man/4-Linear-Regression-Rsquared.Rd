\name{pwrss.f.reg}
\alias{pwrss.f.reg}
\alias{pwrss.f.regression}


\title{Linear Regression: R-squared or R-squared Difference (F Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test R-squared deviation from 0 (zero) in linear regression or to test R-squared difference between two linear regression models. The test of R-squared difference is often used to evaluate incremental contribution of a set of predictors in hierarchical linear regression.

Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.
}

\usage{
pwrss.f.reg(r2 = 0.10, f2 = r2 /(1 - r2),
            k = 1, m = k, alpha = 0.05,
            n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{r2}{expected R-squared (or R-squared change)}
  \item{f2}{expected Cohen's f-squared (an alternative to \code{r2} specification). f2 = r2 / (1 - r2)}
  \item{k}{(total) number of predictors}
  \item{m}{number of predictors in the subset of interest. By default \code{m = k}, which implies that one is interested in the contribution of all predictors, and tests whether R-squared value is different from 0 (zero)}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{verbose}{if \code{FALSE} no output is printed on the console}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (F test)}
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
Bulus, M., & Polat, C. (in press). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi. \url{https://osf.io/ua5fc/download/}

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
