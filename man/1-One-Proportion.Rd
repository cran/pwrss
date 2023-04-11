\name{pwrss.z.prop}
\alias{pwrss.z.prop}

\title{One Proportion against a Constant (z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a proportion against a constant.

Formulas are validated using Monte Carlo simulation, G*Power, \url{http://powerandsamplesize.com/} and tables in PASS documentation.

}

\usage{
pwrss.z.prop(p, p0 = 0, margin = 0, arcsin.trans = FALSE, alpha = 0.05,
             alternative = c("not equal","greater","less",
                             "equivalent","non-inferior","superior"),
             n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{p}{expected proportion}
  \item{p0}{constant to be compared (a proportion)}
  \item{arcsin.trans}{if \code{TRUE} uses Cohen's arcsine transformation, if \code{FALSE} uses normal approximation (default)}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin: boundry of \code{p - p0} that is practically insignificant)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z or t test)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}


\examples{
# Example 1: expecting p - p0 smaller than 0
## one-sided test with normal approximation
pwrss.z.prop(p = 0.45, p0 = 0.50,
             alpha = 0.05, power = 0.80,
             alternative = "less",
             arcsin.trans = FALSE)
## one-sided test with arcsine transformation
pwrss.z.prop(p = 0.45, p0 = 0.50,
             alpha = 0.05, power = 0.80,
             alternative = "less",
             arcsin.trans = TRUE)

# Example 2: expecting p - p0 smaller than 0 or greater than 0
## two-sided test with normal approximation
pwrss.z.prop(p = 0.45, p0 = 0.50,
             alpha = 0.05, power = 0.80,
             alternative = "not equal",
             arcsin.trans = FALSE)
## two-sided test with arcsine transformation
pwrss.z.prop(p = 0.45, p0 = 0.50,
             alpha = 0.05, power = 0.80,
             alternative = "not equal",
             arcsin.trans = TRUE)

# Example 2: expecting p - p0 smaller than 0.01
# when smaller proportion is better
## non-inferiority test with normal approximation
pwrss.z.prop(p = 0.45, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "non-inferior",
             arcsin.trans = FALSE)
## non-inferiority test with arcsine transformation
pwrss.z.prop(p = 0.45, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "non-inferior",
             arcsin.trans = TRUE)

# Example 3: expecting p - p0 greater than -0.01
# when bigger proportion is better
## non-inferiority test with normal approximation
pwrss.z.prop(p = 0.55, p0 = 0.50, margin = -0.01,
             alpha = 0.05, power = 0.80,
             alternative = "non-inferior",
             arcsin.trans = FALSE)
## non-inferiority test with arcsine transformation
pwrss.z.prop(p = 0.55, p0 = 0.50, margin = -0.01,
             alpha = 0.05, power = 0.80,
             alternative = "non-inferior",
             arcsin.trans = TRUE)

# Example 4: expecting p - p0 smaller than -0.01
# when smaller proportion is better
## superiority test with normal approximation
pwrss.z.prop(p = 0.45, p0 = 0.50, margin = -0.01,
             alpha = 0.05, power = 0.80,
             alternative = "superior",
             arcsin.trans = FALSE)
## superiority test with arcsine transformation
pwrss.z.prop(p = 0.45, p0 = 0.50, margin = -0.01,
             alpha = 0.05, power = 0.80,
             alternative = "superior",
             arcsin.trans = TRUE)

# Example 5: expecting p - p0 greater than 0.01
# when bigger proportion is better
## superiority test with normal approximation
pwrss.z.prop(p = 0.55, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "superior",
             arcsin.trans = FALSE)
## superiority test with arcsine transformation
pwrss.z.prop(p = 0.55, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "superior",
             arcsin.trans = TRUE)

# Example 6: expecting p - p0 between -0.01 and 0.01
## equivalence test with normal approximation
pwrss.z.prop(p = 0.50, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "equivalent",
             arcsin.trans = FALSE)
# equivalence test with arcsine transformation
pwrss.z.prop(p = 0.50, p0 = 0.50, margin = 0.01,
             alpha = 0.05, power = 0.80,
             alternative = "equivalent",
             arcsin.trans = TRUE)
}

\references{
Bulus, M., & Polat, C. (in press). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi. \url{https://osf.io/ua5fc/download/}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
