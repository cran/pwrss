\name{power.z.oneprop}
\alias{power.z.oneprop}
\alias{power.exact.oneprop}
\alias{pwrss.z.prop}

\title{Power Analysis for the Test of One Proportion (Normal Approximation and Exact Methods)}

\description{
  Calculates power or sample size (only one can be NULL at a time) for test of a proportion against a constant using normal approximation or exact method.

  Formulas are validated using PASS documentation.

  NOTE: The \code{pwrss.z.prop()} function is deprecated, but it will remain available as a wrapper for the \code{power.z.oneprop()} function during the transition period.

}

\usage{
power.z.oneprop(prob, null.prob = 0.50,
                n = NULL, power = NULL, alpha = 0.05,
                alternative = c("two.sided", "one.sided", "two.one.sided"),
                std.error = c("null", "alternative"),
                arcsine = FALSE, correct = FALSE,
                ceiling = TRUE, verbose = TRUE, pretty = FALSE)

power.exact.oneprop(prob, null.prob = 0.50,
                    n = NULL, power = NULL, alpha = 0.05,
                    alternative = c("two.sided", "one.sided", "two.one.sided"),
                    verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{prob}{probability of success under alternative.}
  \item{null.prob}{probability of success under null.}
  \item{n}{integer; sample size.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{std.error}{character; whether to calculate standard error using "null" or "alternative" value. "null" by default.}
  \item{arcsine}{logical; whether arcsine transformation should be applied. \code{FALSE} by default. Note that when \code{arcsine = TRUE}, any specification to \code{correct} and \code{std.error} will be ignored.}
  \item{correct}{logical; whether Yate's continuity correction should be applied.}
  \item{alternative}{character; the direction or type of the hypothesis test: "two.sided", "one.sided", or "two.one.sided". For non-inferiority or superiority tests, add margin to the null hypothesis value and use \code{alternative = "one.sided"}.}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test ("exact").}
  \item{mean}{mean of the alternative distribution.}
  \item{sd}{standard deviation of the alternative distribution.}
  \item{null.mean}{mean of the null distribution.}
  \item{null.sd}{standard deviation of the null distribution.}
  \item{z.alpha}{critical value(s).}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n}{sample size.}
}


\examples{
# power
power.z.oneprop(prob = 0.45, null.prob = 0.50,
                alpha = 0.05, n = 500,
                alternative = "one.sided")

power.exact.oneprop(prob = 0.45, null.prob = 0.50,
                    alpha = 0.05, n = 500,
                    alternative = "one.sided")



# sample size
power.z.oneprop(prob = 0.45, null.prob = 0.50,
                alpha = 0.05, power = 0.80,
                alternative = "one.sided")

power.exact.oneprop(prob = 0.45, null.prob = 0.50,
                    alpha = 0.05, power = 0.80,
                    alternative = "one.sided")
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}
}
