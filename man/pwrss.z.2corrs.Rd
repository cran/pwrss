\name{Difference between Two Correlations (z Test)}
\alias{pwrss.z.2corrs}

\title{Difference between Two Correlations (Independent Samples z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two independent (Pearson) correlations using Fisher's z transformation.

Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.
}

\usage{
pwrss.z.2corrs(r1 = 0.50, r2 = 0.30,
               alpha = 0.05, kappa = 1,
               alternative = c("not equal","greater","less"),
               n2 = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{r1}{expected correlation in the first group}
  \item{r2}{expected correlation in the second group}
  \item{n2}{sample size in the second group}
  \item{kappa}{\code{n1/n2}}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error.}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", or "less"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size for the first and second groups}
}

\examples{
# difference between r1 and r2 is different from zero
# it could be -0.10 as well as 0.10
pwrss.z.2corrs(r1 = .20, r2 = 0.30,
               alpha = 0.05, power = .80,
               alternative = "not equal")

# difference between r1 and r2 is greater than zero
pwrss.z.2corrs(r1 = .30, r2 = 0.20,
               alpha = 0.05, power = .80,
               alternative = "greater")
}

\references{
Bulus, M., & Polat, C. (2022). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Under review. \url{https://osf.io/ua5fc}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
