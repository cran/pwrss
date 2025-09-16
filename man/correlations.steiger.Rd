\name{power.z.twocors.steiger}
\alias{power.z.twocors.steiger}
\alias{power.z.steiger}

\title{Power Analysis for Dependent Correlations (Steiger's Z-Test)}

\description{
  Calculates power or sample size (only one can be NULL at a time) to test difference between paired correlations (Pearson) using Fisher's Z-transformation.

 Validated via PASS and G*Power.
}

\usage{
power.z.twocors.steiger(rho12, rho13, rho23,
                        rho14 = NULL, rho24 = NULL, rho34 = NULL,
                        n = NULL, power = NULL, alpha = 0.05,
                        alternative = c("two.sided", "one.sided"),
                        pooled = TRUE, common.index = FALSE,
                        ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{rho12}{correlation between variable V1 and V2 (one common index and no common index). Check examples below.}
  \item{rho13}{correlation between variable V1 and V3 (one common index and no common index). Check examples below.}
  \item{rho23}{correlation between variable V2 and V3 (one common index and no common index). Check examples below.}
  \item{rho14}{correlation between variable V1 and V4 (no common index only). Check examples below.}
  \item{rho24}{correlation between variable V2 and V4 (no common index only). Check examples below.}
  \item{rho34}{correlation between variable V3 and V4 (no common index only). Check examples below.}
  \item{n}{integer; sample size.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; the direction or type of the hypothesis test: "two.sided" or "one.sided".}
  \item{pooled}{logical; whether standard error should be pooled. \code{TRUE} by default.}
  \item{common.index}{logical; whether calculations pertain to one common index. \code{TRUE} means calculations involve correlations with a common index (where both correlations share one variable). \code{FALSE} (default) means calculations pertain to correlations with no common index (where all relevant correlations must be explicitly specified). Check examples below.}
  \item{ceiling}{logical; if \code{TRUE} rounds up sample size.}
  \item{verbose}{logical; if \code{FALSE} no output is printed on the console.}
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
  \item{n}{sample size for the first and second groups, in the form of c(n1, n2).}
}

\examples{

# example data for one common index
# compare cor(V1, V2) to cor(V1, V3)

# subject    V1       V2      V3
# <int>    <dbl>    <dbl>    <dbl>
#   1       1.2      2.3      0.8
#   2      -0.0      1.1      0.7
#   3       1.9     -0.4     -2.3
#   4       0.7      1.3      0.4
#   5       2.1     -0.1      0.8
#   ...     ...      ...      ...
#   1000   -0.5      2.7     -1.7

# V1: socio-economic status (common)
# V2: pretest
# V3: post-test

power.z.twocors.steiger(rho12 = 0.35, rho13 = 0.45, rho23 = 0.05,
                        n = 1000, power = NULL, alpha = 0.05,
                        alternative = "two.sided",
                        common.index = TRUE)


# example data for no common index
# compare cor(V1, V2) to cor(V3, V4)

# subject    V1       V2       V3       V4
# <int>    <dbl>    <dbl>    <dbl>    <dbl>
#   1       1.2      2.3      0.8      1.2
#   2      -0.0      1.1      0.7      0.9
#   3       1.9     -0.4     -2.3     -0.1
#   4       0.7      1.3      0.4     -0.3
#   5       2.1     -0.1      0.8      2.7
#   ...     ...      ...      ...      ...
#   1000   -0.5      2.7     -1.7      0.8

# V1: pretest reading
# V2: pretest math
# V3: post-test reading
# V4: post-test math

power.z.twocors.steiger(rho12 = 0.45, rho13 = 0.45, rho23 = 0.50,
                        rho14 = 0.50, rho24 = 0.80, rho34 = 0.55,
                        n = 1000, power = NULL, alpha = 0.05,
                        alternative = "two.sided",
                        common.index = FALSE)
}

\references{
Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87(2), 245-251. \doi{10.1037/0033-2909.87.2.245}
}
