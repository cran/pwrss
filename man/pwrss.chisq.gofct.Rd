\name{Chi-square Goodness-of-Fit or Independence Test (Non-parametric)}
\alias{pwrss.chisq.gofit}


\title{Chi-square Goodness-of-Fit or Independence Test}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test goodness-of-fit or independence.
}

\usage{
pwrss.chisq.gofit(p1 = c(0.50, 0.50),
                  p0 = .chisq.fun(p1)$p0,
                  w = .chisq.fun(p1)$w,
                  df = .chisq.fun(p1)$df,
                  n = NULL, power = NULL,
                  alpha = 0.05, verbose = TRUE)
}

\arguments{
  \item{p1}{a vector or a matrix of cell probabilities under alternative hypothesis}
  \item{p0}{a vector or a matrix of cell probabilities under null hypothesis}
  \item{w}{effect size. By default, it is computed from \code{p1} and \code{p0}; however, when specified by the user it can be any of the Cohen's W, Phi coefficient, Cramer's V, etc. Phi coefficient is defined as \code{sqrt(X2/n)} and Cramer's V is defined as \code{sqrt(X2/(n*v))} where \code{v} is \code{min(nrow - 1, ncol - 1)}}
  \item{df}{degrees of freedom. (Number of cells - 1) for a vector, (Number of rows - 1) * (Number of columns - 1) for a matrix}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{df}{degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{total sample size}
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
prob.mat <- c(0.28, 0.72) # null hypothesis states that c(0.50, 0.50)
pwrss.chisq.gofit(p1 = c(0.28, 0.72),
                  alpha = 0.05, power = 0.80)

## Option 2: Use Cohe's W = 0.44
## df is k - 1 for Cohen's W
pwrss.chisq.gofit(w = 0.44, df = 1,
                  alpha = 0.05, power = 0.80)


# ---------------------------------------------------------#
# Example 2: Phi Coefficient (or Cramer's V or Cohen's W)  #
# test of independence for 2 x 2 contingency tables        #
# How many subjects are needed to claim that               #
# girls are underdiagnosed with ADHD?                      #
# ---------------------------------------------------------#

## Option 1: Use cell probabilities
## from https://time.com/growing-up-with-adhd/
## 5.6 percent of girls and 13.2 percent of boys are diagnosed with ADHD
prob.mat <- rbind(c(0.056, 0.132),
                  c(0.944, 0.868))
colnames(prob.mat) <- c("Girl", "Boy")
rownames(prob.mat) <- c("ADHD", "No ADHD")
prob.mat
pwrss.chisq.gofit(p1 = prob.mat,
                  alpha = 0.05, power = 0.80)

## Option 2: Use Phi coefficient = 0.1302134
## df is 1 for Phi coefficient
pwrss.chisq.gofit(w = 0.1302134, df = 1,
                  alpha = 0.05, power = 0.80)


# --------------------------------------------------------#
# Example 3: Cramer's V (or Cohen's W)                    #
# test of independence for j x k contingency tables       #
# How many subjects are needed to detect the relationship #
# between depression severity and gender?                 #
# --------------------------------------------------------#

## Option 1: Use cell probabilities
## from https://doi.org/10.1016/j.jad.2019.11.121
prob.mat <- cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078),
                  c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101))
rownames(prob.mat) <- c("Normal", "Mild", "Moderate", "Severe", "Extremely Severe")
colnames(prob.mat) <- c("Female", "Male")
prob.mat
pwrss.chisq.gofit(p1 = prob.mat,
                  alpha = 0.05, power = 0.80)

# Option 2: Use Cramer's V = 0.03022008 based on 5 x 2 contingency table
# df is (nrow - 1) * (ncol - 1) for Cramer's V
pwrss.chisq.gofit(w = 0.03022008, df = 4,
                  alpha = 0.05, power = 0.80)

}

\references{
Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
