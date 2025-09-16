\name{power.exact.mcnemar}
\alias{power.exact.mcnemar}
\alias{power.exact.twoprops.mcnemar}

\title{Power Analysis for McNemar's Exact Test (Paired Proportions)}

\description{
  Calculates power or sample size for McNemar's test on paired binary outcomes. Approximate and exact methods are available (check references for details).

  Validated via PASS and G*Power.
}

\usage{
  power.exact.mcnemar(prob10, prob01, n.paired = NULL,
                    power = NULL,  alpha = 0.05,
                    alternative = c("two.sided", "one.sided"),
                    method = c("exact", "approximate"),
                    ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{prob10}{(joint) probability of success in case (or after) but failure in matched control (or before). 'prob10' and 'prob01' are known as discordant probs.}
  \item{prob01}{(joint) probability of failure in case (or after) but success in matched control (or before). prob10' and 'prob01' are known as discordant probs.}
  \item{n.paired}{number of pairs, which is sum of cell frequencies in 2 x 2 table (f11 + f10 + f01 + f00), or number of rows in a data frame with matched variables 'case' and 'control' or 'after' and 'before'.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{method}{character; the method used for power calculation. "exact" specifies Fisher's exact test, while "approximate" refers to the z-test based on the normal approximation.}
  \item{alternative}{character; the direction or type of the hypothesis test: "not equal", "greater", "less". Optionally, users can specify 'two.sided' as an alternative to 'not equal' for consistency with other R statistical functions.}
  \item{ceiling}{logical; if \code{TRUE} rounds up sample size in each cell. This procedure assumes symmetry for concordant probs, which are 'p11' and 'p00'). Thus results may differ from other software by a few units. To match results set 'ceiling = FALSE'.}
  \item{verbose}{logical; if \code{FALSE} no output is printed on the console.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the test, which is "exact" or "z".}
  \item{odds.ratio}{odds ratio.}
  \item{mean}{mean of the alternative distribution.}
  \item{sd}{standard deviation of the alternative distribution.}
  \item{null.mean}{mean of the null distribution.}
  \item{null.sd}{standard deviation of the null distribution.}
  \item{z.alpha}{critical value(s).}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n.paired}{paired sample size, which is sum of cell frequencies in 2 x 2 table (f11 + f10 + f01 + f00), or number of rows in a data frame with variables 'case' and 'control' or 'after' and 'before'.}

}

\examples{
# example data for a matched case-control design
# subject  case     control
# <int>    <dbl>    <dbl>
#   1        1        1
#   2        0        1
#   3        1        0
#   4        0        1
#   5        1        1
#   ...     ...      ...
#   100      0        0

# example data for a before-after design
# subject  before   after
# <int>    <dbl>    <dbl>
#   1        1        1
#   2        0        1
#   3        1        0
#   4        0        1
#   5        1        1
#   ...     ...      ...
#   100      0        0

# convert to a 2 x 2 frequency table
freqs <- matrix(c(30, 10, 20, 40), nrow = 2, ncol = 2)
colnames(freqs) <- c("control_1", "control_0")
rownames(freqs) <- c("case_1", "case_0")
freqs

# convert to a 2 x 2 proportion table
props <- freqs / sum(freqs)
props

# discordant pairs (0 and 1, or 1 and 0) in 'props' matrix
# are the sample estimates of prob01 and prob10

# post-hoc exact power
power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10,
                    n.paired = 100, alpha = 0.05,
                    method = "exact", alt = "two.sided")

# required sample size for exact test
# assuming prob01 and prob10 are population parameters
power.exact.mcnemar(prob10 = 0.20, prob01 = 0.10,
                    power = 0.80, alpha = 0.05,
                    method = "exact", alt = "two.sided")

# we may not have 2 x 2 joint probs
# convert marginal probs to joint probs
joint.probs.2x2(prob1 = 0.55, # mean of case group (or after)
                    prob2 = 0.45, # mean of matched control group (or before)
                    # correlation between matched case-control or before-after
                    rho = 0.4141414
)
}

\references{
Bennett, B. M., & Underwood, R. E. (1970). 283. Note: On McNemar's Test for the 2 * 2 Table and Its Power Function. Biometrics, 26(2), 339-343. \doi{10.2307/2529083}

Connor, R. J. (1987). Sample size for testing differences in proportions for the paired-sample design. Biometrics, 43(1), 207-211. \doi{10.2307/2531961}

Duffy, S. W. (1984). Asymptotic and exact power for the McNemar test and its analogue with R controls per case. Biometrics, 40(4) 1005-1015. \doi{10.2307/2531151}

McNemar, Q. (1947). Note on the sampling error of the difference between correlated proportions or percentages. Psychometrika, 12(2), 153-157. \doi{10.1007/BF02295996}

Miettinen, O. S. (1968). The matched pairs design in the case of all-or-none responses. Biometrics, 24(2), 339-352. \doi{10.2307/2528039}
}
