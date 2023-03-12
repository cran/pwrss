\name{Difference between Two Means (Non-parametric)}
\alias{pwrss.np.2means}

\title{Difference between Two Means (Non-parametric Test for Independent and Paired Samples)}

\description{
  Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two means.
  For standardized mean difference (Cohen's d) set \code{mu1 = d} and use defaults for \code{mu2}, \code{sd1}, and \code{sd2}.
  If pooled standard deviation (psd) is available set \code{sd1 = psd}.

  Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.
}

\usage{
pwrss.np.2means(mu1 = 0.20, mu2 = 0,
                sd1 = ifelse(paired, sqrt(1/(2*(1-paired.r))), 1), sd2 = sd1,
                margin = 0, alpha = 0.05, paired = FALSE, paired.r = 0.50,
                kappa = 1, n2 = NULL, power = NULL,
                alternative = c("not equal", "greater", "less",
                                "non-inferior", "superior", "equivalent"),
                distribution = c("normal", "uniform", "double exponential",
                                 "laplace", "logistic"),
                method = c("guenther", "noether"),
                verbose = TRUE)

}

\arguments{
  \item{mu1}{expected mean in the first group}
  \item{mu2}{expected mean in the second group}
  \item{sd1}{standard deviation in the first group}
  \item{sd2}{standard deviation in the second group}
  \item{kappa}{\code{n1/n2}}
  \item{paired}{if \code{TRUE} paired samples}
  \item{paired.r}{correlation between the two repeated measures (e.g., pretest and posttest)}
  \item{n2}{sample size in the second group (or for the single group in paired samples)}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin: boundry of \code{mu1 - mu2} that is practically insignificant)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{distribution}{parent distribution: "normal", "uniform", "double exponential", "laplace", or "logistic"}
  \item{method}{non-parametric method: "guenther" (default) or "noether"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{df}{degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}


\examples{

# Mann-Whitney U or Wilcoxon rank-sum test
# (a.k.a Wilcoxon-Mann-Whitney test) for independent samples

## difference between group 1 and group 2 is not equal to zero
## estimated difference is Cohen'd = 0.25
pwrss.np.2means(mu1 = 0.25, mu2 = 0, power = 0.80,
                alternative = "not equal")

## difference between group 1 and group 2 is greater than zero
## estimated difference is Cohen'd = 0.25
pwrss.np.2means(mu1 = 0.25, mu2 = 0, power = 0.80,
                alternative = "greater")

## mean of group 1 is practically not smaller than mean of group 2
## estimated difference is Cohen'd = 0.10 and can be as small as -0.05
pwrss.np.2means(mu1 = 0.25, mu2 = 0.15,
                margin = -0.05, power = 0.80,
                alternative = "non-inferior")

## mean of group 1 is practically greater than mean of group 2
## estimated difference is Cohen'd = 0.10 and can be as small as 0.05
pwrss.np.2means(mu1 = 0.25, mu2 = 0.15,
                margin = 0.05, power = 0.80,
                alternative = "superior")

## mean of group 1 is practically same as mean of group 2
## estimated difference is Cohen'd = 0
## and can be as small as -0.05 and as high as 0.05
pwrss.np.2means(mu1 = 0.25, mu2 = 0.25,
                margin = 0.05, power = 0.80,
                alternative = "equivalent")


# Wilcoxon signed-rank test for matched pairs (dependent samples)

## difference between time 1 and time 2 is not equal to zero
## estimated difference between time 1 and time 2 is Cohen'd = -0.25
pwrss.np.2means(mu1 = 0, mu2 = 0.25, power = 0.80,
                paired = TRUE, paired.r = 0.50,
                alternative = "not equal")

## difference between time 1 and time 2 is greater than zero
## estimated difference between time 1 and time 2 is Cohen'd = -0.25
pwrss.np.2means(mu1 = 0, mu2 = 0.25, power = 0.80,
                paired = TRUE, paired.r = 0.50,
                alternative = "greater")

## mean of time 1 is practically not smaller than mean of time 2
## estimated difference is Cohen'd = -0.10 and can be as small as 0.05
pwrss.np.2means(mu1 = 0.15, mu2 = 0.25, margin = 0.05,
                paired = TRUE, paired.r = 0.50, power = 0.80,
                alternative = "non-inferior")

## mean of time 1 is practically greater than mean of time 2
## estimated difference is Cohen'd = -0.10 and can be as small as -0.05
pwrss.np.2means(mu1 = 0.15, mu2 = 0.25, margin = -0.05,
                paired = TRUE, paired.r = 0.50, power = 0.80,
                alternative = "superior")

## mean of time 1 is practically same as mean of time 2
## estimated difference is Cohen'd = 0
## and can be as small as -0.05 and as high as 0.05
pwrss.np.2means(mu1 = 0.25, mu2 = 0.25, margin = 0.05,
                paired = TRUE, paired.r = 0.50, power = 0.80,
                alternative = "equivalent")
}

\references{
Al-Sunduqchi, M. S. (1990). Determining the appropriate sample size for inferences based on the Wilcoxon statistics [Unpublished doctoral dissertation]. University of Wyoming - Laramie

Chow, S. C., Shao, J., Wang, H., and Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Noether, G. E. (1987). Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82(1), 645-647.

Ruscio, J. (2008). A probability-based measure of effect size: Robustness to base rates and other factors. Psychological Methods, 13(1), 19-30.

Ruscio, J., & Mullen, T. (2012). Confidence intervals for the probability of superiority effect size measure and the area under a receiver operating characteristic curve. Multivariate Behavioral Research, 47(2), 201-223.

Zhao, Y.D., Rahardja, D., & Qu, Y. (2008). Sample size calculation for the Wilcoxon-Mann-Whitney test adjusting for ties. Statistics in Medicine, 27(3), 462-468.
}
