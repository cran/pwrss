\name{power.z.twocors}
\alias{power.z.twocors}
\alias{power.z.twocor}
\alias{pwrss.z.2corrs}
\alias{pwrss.z.2corr}

\title{Power Analysis for Independent Correlations}

\description{
Calculates power or sample size (only one can be NULL at a time) to test difference between two independent (Pearson) correlations using Fisher's z transformation.

Formulas are validated using PASS and G*Power.
}

\usage{
power.z.twocors(rho1, rho2,
                n2 = NULL, n.ratio = 1,
                power = NULL, alpha = 0.05,
                alternative = c("two.sided", "one.sided"),
                ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{rho1}{correlation in the first group.}
  \item{rho2}{correlation in the second group.}
  \item{n2}{sample size in the second group. Sample size in the first group can be calculated as \code{n2*kappa}. By default, \code{n1 = n2} because \code{kappa = 1}.}
  \item{n.ratio}{\code{n1/n2} ratio.}
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
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size for the first and second groups, in the form of c(n1, n2).}
}

\examples{
# difference between r1 and r2 is different from zero
# it could be -0.10 as well as 0.10
power.z.twocors(rho1 = .20, rho2 = 0.30,
               alpha = 0.05, power = .80,
               alternative = "two.sided")

# difference between r1 and r2 is greater than zero
power.z.twocors(rho1 = .30, rho2 = 0.20,
               alpha = 0.05, power = .80,
               alternative = "one.sided")
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
