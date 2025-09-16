\name{power.z.twoprops}
\alias{power.z.twoprops}
\alias{power.z.twoprop}
\alias{power.exact.twoprops}
\alias{power.exact.twoprop}
\alias{power.z.twoprop}
\alias{pwrss.z.2props}
\alias{pwrss.z.2prop}

\title{Power Analysis for Testing Difference Between Two Proportions (Normal Approximation and Exact Methods)}

\description{
  Calculates power or sample size (only one can be NULL at a time) for two proportions using normal approximation method.

  Validated via G*Power and PASS documentation.

  NOTE: The \code{pwrss.z.2props()} function is deprecated, but it will remain available as a wrapper for the \code{power.z.twoprops()} function during the transition period.

}

\usage{
power.z.twoprops(prob1, prob2, margin = 0,
                 n2 = NULL, n.ratio = 1,
                 power = NULL, alpha = 0.05,
                 alternative = c("two.sided", "one.sided", "two.one.sided"),
                 arcsine = FALSE, correct = FALSE,
                 paired = FALSE, rho.paired = 0.50,
                 std.error = c("pooled", "unpooled"),
                 ceiling = TRUE, verbose = TRUE, pretty = FALSE)

power.exact.twoprops(prob1, prob2, n2 = NULL, n.ratio = 1,
                     power = NULL, alpha = 0.05,
                     alternative = c("two.sided", "one.sided"),
                     paired = FALSE, rho.paired = 0.50,
                     method = c("exact", "approximate"),
                     ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{prob1}{probability of success in the first group.}
  \item{prob2}{probability of success in the second group.}
  \item{margin}{ignorable \code{prob1} - \code{prob2} difference. For two one-sided tests provide lower and upper margins in the form of c(lower, upper).}
  \item{n2}{integer; sample size for the second group.}
  \item{n.ratio}{sample size ratio (n1 / n2).}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{paired}{logical; if \code{TRUE} samples are paired. \code{FALSE} by default.}
  \item{rho.paired}{correlation between paired observations.}
  \item{method}{character; whether to use "approximate" or "exact" method. Default is "exact" (only in the \code{power.exact.twoprops()} function).}
  \item{arcsine}{logical; whether arcsine transformation should be applied. Note that this only applies to independent proportions without continuity correction. }
  \item{correct}{logical; whether Yates' continuity correction should be applied to the test statistic. Ignored for the paired test.}
  \item{std.error}{character; whether to calculate standard error using "pooled" or "unpooled" standard deviation. Ignored for the paired test.}
  \item{alternative}{character; direction or type of the hypothesis test: "two.sided", "one.sided", or "two.one.sided".}
  \item{ceiling}{logical; \code{TRUE} rounds up sample size in each group.}
  \item{verbose}{logical; \code{TRUE} prints the output on the console.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
 \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the test, which is "z" or "exact".}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{mean}{mean of the alternative distribution.}
  \item{sd}{standard deviation of the alternative distribution.}
  \item{null.mean}{mean of the null distribution.}
  \item{null.sd}{standard deviation of the null distribution.}
  \item{z.alpha}{critical value(s).}
  \item{n}{sample size in the form of c(n1, n2) (applies to independent proportions).}
  \item{n.total}{total sample size (applies to independent proportions).}
  \item{n.paired}{paired sample size (applies to paired proportions).}
}


\examples{
  # power
  power.z.twoprops(prob1 = 0.65, prob2 = 0.60,
                   alpha = 0.05, n2 = 500,
                   alternative = "one.sided")

  # sample size
  power.z.twoprops(prob1 = 0.65, prob2 = 0.60,
                   alpha = 0.05, power = 0.80,
                   alternative = "one.sided")
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}
}
