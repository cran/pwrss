\name{One Correlation (z Test)}
\alias{pwrss.z.corr}


\title{A Correlation against a Constant (One Sample z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a (Pearson) correlation against a constant using Fisher's z transformation.

Formulas are validated using G*Power and tables in PASS documentation.
}

\usage{
pwrss.z.corr(r = 0.50, r0 = 0, alpha = 0.05,
             alternative = c("not equal","greater","less"),
             n = NULL, power = NULL, verbose = TRUE)
}

\arguments{
  \item{r}{expected correlation}
  \item{r0}{constant proportion to be compared}
  \item{n}{sample size}
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
  \item{n}{sample size}
}

\examples{
# expected correlation is 0.20 and it is different from 0
# it could be 0.20 as well as -0.20
pwrss.z.corr(r = 0.20, r0 = 0,
             alpha = 0.05, power = 0.80,
             alternative = "not equal")

# expected correlation is 0.20 and it is greater than 0.10
pwrss.z.corr(r = 0.20, r0 = 0.10,
             alpha = 0.05, power = 0.80,
             alternative = "greater")
}

\references{
Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Under review. \url{https://osf.io/ua5fc}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
