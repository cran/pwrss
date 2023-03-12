\name{Difference between Two Proportions (z Test)}
\alias{pwrss.z.2props}

\title{Difference between Two Proportions (z Test)}

\description{
  Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two proportions.

  Formulas are validated using Monte Carlo simulation, G*Power, \url{http://powerandsamplesize.com/} and tables in PASS documentation.
}

\usage{
pwrss.z.2props(p1, p2, margin = 0, arcsin.trans = FALSE, kappa = 1, alpha = 0.05,
               alternative = c("not equal","greater","less",
                               "equivalent","non-inferior","superior"),
               n2 = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{p1}{expected proportion in the first group}
  \item{p2}{expected proportion in the second group}
  \item{arcsin.trans}{if \code{TRUE} uses arcsin transformation, if \code{FALSE} uses normal approximation (default)}
  \item{kappa}{n1/n2}
  \item{n2}{sample size in the second group}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin: boundry of \code{p1 - p2} that is practically insignificant)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size for the first and second group}
}


\examples{
# Example 1: expecting p1 - p2 smaller than 0
## one-sided test with normal approximation
pwrss.z.2props(p1 = 0.45, p2 = 0.50,
               alpha = 0.05, power = 0.80,
               alternative = "less",
               arcsin.trans = FALSE)
## one-sided test with arcsine transformation
pwrss.z.2props(p1 = 0.45, p2 = 0.50,
               alpha = 0.05, power = 0.80,
               alternative = "less",
               arcsin.trans = TRUE)

# Example 2: expecting p1 - p2 smaller than 0 or greater than 0
## two-sided test with normal approximation
pwrss.z.2props(p1 = 0.45, p2 = 0.50,
               alpha = 0.05, power = 0.80,
               alternative = "not equal",
               arcsin.trans = FALSE)
## two-sided test with arcsine transformation
pwrss.z.2props(p1 = 0.45, p2 = 0.50,
               alpha = 0.05, power = 0.80,
               alternative = "not equal",
               arcsin.trans = TRUE)

# Example 2: expecting p1 - p2 smaller than 0.01
# when smaller proportion is better
## non-inferiority test with normal approximation
pwrss.z.2props(p1 = 0.45, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "non-inferior",
               arcsin.trans = FALSE)
## non-inferiority test with arcsine transformation
pwrss.z.2props(p1 = 0.45, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "non-inferior",
               arcsin.trans = TRUE)

# Example 3: expecting p1 - p2 greater than -0.01
# when bigger proportion is better
## non-inferiority test with normal approximation
pwrss.z.2props(p1 = 0.55, p2 = 0.50, margin = -0.01,
               alpha = 0.05, power = 0.80,
               alternative = "non-inferior",
               arcsin.trans = FALSE)
## non-inferiority test with arcsine transformation
pwrss.z.2props(p1 = 0.55, p2 = 0.50,  margin = -0.01,
               alpha = 0.05, power = 0.80,
               alternative = "non-inferior",
               arcsin.trans = TRUE)

# Example 4: expecting p1 - p2 smaller than -0.01
# when smaller proportion is better
## superiority test with normal approximation
pwrss.z.2props(p1 = 0.45, p2 = 0.50, margin = -0.01,
               alpha = 0.05, power = 0.80,
               alternative = "superior",
               arcsin.trans = FALSE)
## superiority test with arcsine transformation
pwrss.z.2props(p1 = 0.45, p2 = 0.50, margin = -0.01,
               alpha = 0.05, power = 0.80,
               alternative = "superior",
               arcsin.trans = TRUE)

# Example 5: expecting p1 - p2 greater than 0.01
# when bigger proportion is better
## superiority test with normal approximation
pwrss.z.2props(p1 = 0.55, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "superior",
               arcsin.trans = FALSE)
## superiority test with arcsine transformation
pwrss.z.2props(p1 = 0.55, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "superior",
               arcsin.trans = TRUE)

# Example 6: expecting p1 - p2 between -0.01 and 0.01
## equivalence test with normal approximation
pwrss.z.2props(p1 = 0.50, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "equivalent",
               arcsin.trans = FALSE)
# equivalence test with arcsine transformation
pwrss.z.2props(p1 = 0.50, p2 = 0.50, margin = 0.01,
               alpha = 0.05, power = 0.80,
               alternative = "equivalent",
               arcsin.trans = TRUE)
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Under review. \url{https://osf.io/ua5fc}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
