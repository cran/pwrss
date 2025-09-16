\name{power.z.onecor}
\alias{power.z.onecor}
\alias{pwrss.z.corr}
\alias{pwrss.z.cor}

\title{Power Analysis for One-Sample Correlation}

\description{
Calculates power or sample size (only one can be NULL at a time) to test a (Pearson) correlation against a constant using Fisher's z transformation.

Formulas are validated using PASS and G*Power.
}

\usage{
power.z.onecor(rho, null.rho = 0,
               n = NULL, power = NULL, alpha = 0.05,
               alternative = c("two.sided", "one.sided"),
               ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{rho}{correlation.}
  \item{null.rho}{correlation when null is true.}
  \item{n}{sample size.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; direction or type of the hypothesis test: "two.sided" or "one.sided".}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (Z-Test)}
  \item{mean}{mean of the alternative distribution.}
  \item{sd}{standard deviation of the alternative distribution.}
  \item{null.mean}{mean of the null distribution.}
  \item{null.sd}{standard deviation of the null distribution.}
  \item{z.alpha}{critical value(s).}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n}{sample size.}
}

\examples{
# expected correlation is 0.20 and it is different from 0
# it could be 0.20 as well as -0.20
power.z.onecor(rho = 0.20,
               power = 0.80,
               alpha = 0.05,
               alternative = "two.sided")

# expected correlation is 0.20 and it is greater than 0.10
power.z.onecor(rho = 0.20, null = 0.10,
               power = 0.80,
               alpha = 0.05,
               alternative = "one.sided")

}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
