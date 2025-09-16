\name{power.f.regression}
\alias{pwrss.f.reg}
\alias{power.f.reg}
\alias{pwrss.f.regression}
\alias{power.f.regression}


\title{Power Analysis for Linear Regression: R-squared or R-squared Change (F-Test)}

\description{
Calculates power or sample size (only one can be NULL at a time) to test R-squared deviation from 0 (zero) in linear regression or to test R-squared change between two linear regression models. The test of R-squared change is often used to evaluate incremental contribution of a set of predictors in hierarchical linear regression.

Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.

 NOTE: The \code{pwrss.f.reg()} function and its alias \code{pwrss.f.regression} are deprecated, but they will remain available as a wrapper for \code{power.f.regression()} during the transition period.

}

\usage{
power.f.regression(r.squared.change = NULL, margin = 0,
                   k.total, k.tested = k.total,
                   n = NULL, power = NULL, alpha = 0.05,
                   ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{r.squared.change}{R-squared (or R-squared change).}
  \item{margin}{margin - ignorable R-squared (or R-squared change).}
  \item{k.total}{integer; total number of predictors.}
  \item{k.tested}{integer; number of predictors in the subset of interest. By default \code{m.tested = k.total}, which implies that one is interested in the contribution of all predictors, and tests whether R-squared value is different from 0 (zero).}
  \item{n}{integer; sample size.}
 \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
   \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (F-Test).}
  \item{df1}{numerator degrees of freedom.}
  \item{df2}{denominator degrees of freedom.}
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null.}
  \item{f.alpha}{critical value.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n}{sample size.}
}

\examples{
# in the outcome (R-squared = 0.15).
power.f.regression(r.squared = 0.15,
                   k.total = 3, # total number of predictors
                   power = 0.80)

# adding two more variables will increase R-squared
# from 0.15 (with 3 predictors) to 0.25 (with 3 + 2 predictors)
power.f.regression(r.squared.change = 0.10, # R-squared change
                   k.total = 5, # total number of predictors
                   k.tested = 2, # predictors to be tested
                   power = 0.80)
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
