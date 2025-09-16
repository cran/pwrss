\name{power.chisq.gof}
\alias{power.chisq.gof}
\alias{pwrss.chisq.gofit}

\title{Power and Sample Size for Chi-square Goodness-of-Fit or Independence Tests}

\description{
Calculates power or sample size (only one can be NULL at a time) for Chi-square goodness-of-fit or independence tests.

NOTE: The \code{pwrss.chisq.gofit()} function is deprecated. However, it will remain available as a wrapper for the \code{power.chisq.gof()} function.

}

\usage{
power.chisq.gof(w, null.w = 0, df,
                n = NULL, power = NULL, alpha = 0.05,
                ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{

  \item{w}{Cohen's w effect size under alternative. It can be any of Cohen's W, Phi coefficient, or Cramer's V but degrees of freedom should be specified accordingly. Phi coefficient is defined as \code{sqrt(X2/n)} and Cramer's V is defined as \code{sqrt(X2/(n*v))} where \code{v} is \code{min(nrow - 1, ncol - 1)} and X2 is the chi-square statistic.}
  \item{null.w}{Cohen's w effect size under null.}
  \item{df}{integer; degrees of freedom.  Defined as (n.cells - 1) if \code{p1} is a vector, and as (n.rows - 1) * (n.cols - 1) if \code{p1} is a matrix.}
  \item{n}{integer; total sample size.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (Chi-square Test).}
  \item{df}{degrees of freedom.}
  \item{ncp}{non-centrality parameter under alternative.}
  \item{null.ncp}{non-centrality parameter under null.}
  \item{chisq.alpha}{critical value.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n}{total sample size.}
}

\examples{
# ---------------------------------------------------------#
# Example 1: Cohen's W                                     #
# goodness-of-fit test for 1 x k or k x 1 table            #
# How many subjects are needed to claim that               #
# girls choose STEM related majors less than males?       #
# ---------------------------------------------------------#

## Option 1: Use cell probabilities
## from https://www.aauw.org/resources/research/the-stem-gap/
## 28 percent of the  workforce in STEM field is women
prob.vector <- c(0.28, 0.72)
null.prob.vector <- c(0.50, 0.50)
probs.to.w(prob.vector, null.prob.vector)

power.chisq.gof(w = 0.44, df = 1,
                alpha = 0.05, power = 0.80)


# ---------------------------------------------------------#
# Example 2: Phi Coefficient (or Cramer's V or Cohen's W)  #
# test of independence for 2 x 2 contingency tables        #
# How many subjects are needed to claim that               #
# girls are underdiagnosed with ADHD?                      #
# ---------------------------------------------------------#

## from https://time.com/growing-up-with-adhd/
## 5.6 percent of girls and 13.2 percent of boys are diagnosed with ADHD
prob.matrix <- rbind(c(0.056, 0.132),
                     c(0.944, 0.868))
colnames(prob.matrix) <- c("Girl", "Boy")
rownames(prob.matrix) <- c("ADHD", "No ADHD")
prob.matrix

probs.to.w(prob.matrix)

power.chisq.gof(w = 0.1302134, df = 1,
                alpha = 0.05, power = 0.80)


# --------------------------------------------------------#
# Example 3: Cramer's V (or Cohen's W)                    #
# test of independence for j x k contingency tables       #
# How many subjects are needed to detect the relationship #
# between depression severity and gender?                 #
# --------------------------------------------------------#

## from https://doi.org/10.1016/j.jad.2019.11.121
prob.matrix <- cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078),
                     c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101))
rownames(prob.matrix) <- c("Normal", "Mild", "Moderate",
                           "Severe", "Extremely Severe")
colnames(prob.matrix) <- c("Female", "Male")
prob.matrix

probs.to.w(prob.matrix)

power.chisq.gof(w = 0.03022008, df = 4,
                alpha = 0.05, power = 0.80)

}

\references{
Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
