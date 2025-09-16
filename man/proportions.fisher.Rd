\name{power.exact.fisher}
\alias{power.exact.fisher}
\alias{power.exact.twoprops.fisher}

\title{Power Analysis for Fisher's Exact Test (Independent Proportions)}

\description{
  Calculates power or sample size for Fisher's exact test on independent binary outcomes. Approximate and exact methods are available.

  Validated via PASS and G*Power.
}

\usage{
power.exact.fisher(prob1, prob2, n2 = NULL, n.ratio = 1,
                   alpha = 0.05, power = NULL,
                   alternative = c("two.sided", "one.sided"),
                   method = c("exact", "approximate"),
                   ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{prob1}{probability of success in the first group.}
  \item{prob2}{probability of success in the second group.}
  \item{n2}{integer; sample size for the second group.}
  \item{n.ratio}{n1 / n2 ratio.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; direction or type of the hypothesis test: "two.sided" or "one.sided".}
  \item{method}{character; method used for power calculation. "exact" specifies Fisher's exact test, while "approximate" refers to the Z-Test based on the normal approximation.}
  \item{ceiling}{logical; if \code{TRUE} rounds up sample size in each group.}
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
  \item{n}{sample sizes for the first and second groups, specified as c(n1, n2).}
  \item{n.total}{total sample size, which is sum of cell frequencies in 2 x 2 table (f11 + f10 + f01 + f00), or number of rows in a data frame with group variable stacked.}
}

\examples{
# example data for a randomized controlled trial
# subject  group    success
# <int>    <dbl>      <dbl>
#   1        1          1
#   2        0          1
#   3        1          0
#   4        0          1
#   5        1          1
#   ...     ...        ...
#   100      0          0

# prob1 = mean(success | group = 1)
# prob2 = mean(success | group = 0)

# post-hoc exact power
power.exact.fisher(prob1 = 0.60, prob2 = 0.40, n2 = 50)

# we may have 2 x 2 joint probs such as
# -------------------------------------
#             | group (1) | group (0) |
# -------------------------------------
# success (1) |  0.24    |   0.36     |
# -------------------------------------
# success (0) |  0.16    |   0.24     |
# -------------------------------------

# convert joint probs to marginal probs
marginal.probs.2x2(prob11 = 0.24, prob10 = 0.36,
                   prob01 = 0.16, prob00 = 0.24)
}

\references{

Bennett, B. M., & Hsu, P. (1960). On the power function of the exact test for the 2 x 2 contingency table. Biometrika, 47(3/4), 393-398. \doi{10.2307/2333309}

Fisher, R. A. (1935). The logic of inductive inference. Journal of the Royal Statistical Society, 98(1), 39-82. \doi{10.2307/2342435}

}
