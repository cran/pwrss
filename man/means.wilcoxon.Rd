\name{power.np.wilcoxon}
\alias{power.np.wilcoxon}
\alias{power.np.wilcox}
\alias{pwrss.np.2means}
\alias{pwrss.np.2groups}

\title{Power Analysis for Non-parametric Rank-Based Tests (One-Sample, Independent, and Paired Designs)}

\description{
Calculates power or sample size (only one can be NULL at a time) for non-parametric rank-based tests. The following tests and designs are available:

\itemize{
   \item Wilcoxon Signed-Rank Test (One Sample)
   \item Wilcoxon Rank-Sum or Mann-Whitney U Test (Independent Samples)
   \item Wilcoxon Matched-Pairs Signed-Rank Test (Paired Samples)
 }

Use \code{means.to.d()} to convert raw means and standard deviations to Cohen's d, and \code{d.to.cles()} to convert Cohen's d to the probability of superiority. Note that this interpretation is appropriate only when the underlying distribution is approximately normal and the two groups have similar population variances.

Formulas are validated using G*Power and tables in PASS documentation. However, we adopt rounding convention used by G*Power.

Note that R has a partial matching feature which allows you to specify shortened versions of arguments, such as \code{alt} instead of \code{alternative}, or \code{dist} instead of \code{distribution}.

NOTE: \code{pwrss.np.2means()} function is no longer supported. \code{pwrss.np.2groups()} will remain available for some time.

}

\usage{
power.np.wilcoxon(d, null.d = 0, margin = 0,
                  n2 = NULL, n.ratio = 1, power = NULL, alpha = 0.05,
                  alternative = c("two.sided", "one.sided", "two.one.sided"),
                  design = c("independent", "paired", "one.sample"),
                  distribution = c("normal", "uniform", "double.exponential",
                                   "laplace", "logistic"),
                  method = c("guenther", "noether"),
                  ceiling = TRUE, verbose = TRUE, pretty = FALSE)

}

\arguments{
  \item{d}{Cohen's d or Hedges' g.}
  \item{null.d}{Cohen's d or Hedges' g under null, typically 0 (zero).}
  \item{margin}{margin - ignorable \code{d} - \code{null.d} difference.}
  \item{n2}{integer; sample size in the second group (or for the single group in paired samples or one-sample)}
  \item{n.ratio}{\code{n1/n2} ratio (applies to independent samples only)}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{design}{character; "independent" (default), "one.sample", or "paired".}
  \item{alternative}{character; direction or type of the hypothesis test: "two.sided", "one.sided", or "two.one.sided".}
  \item{distribution}{character; parent distribution: "normal", "uniform", "double.exponential", "laplace", or "logistic".}
  \item{method}{character; non-parametric approach: "guenther" (default) or "noether"}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (Z- or T-Test).}
  \item{df}{degrees of freedom (applies when method = 'guenther').}
  \item{ncp}{non-centrality parameter for the alternative (applies when method = 'guenther').}
  \item{null.ncp}{non-centrality parameter for the null (applies when method = 'guenther').}
  \item{t.alpha}{critical value(s) (applies when method = 'guenther').}
  \item{mean}{mean of the alternative (applies when method = 'noether').}
  \item{null.mean}{mean of the null (applies when method = 'noether').}
  \item{sd}{standard deviation of the alternative (applies when method = 'noether').}
  \item{null.sd}{standard deviation of the null (applies when method = 'noether').}
  \item{z.alpha}{critical value(s) (applies when method = 'noether').}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n}{sample size (`n` or `c(n1, n2)` depending on the design.}
}


\examples{
# Mann-Whitney U or Wilcoxon rank-sum test
# (a.k.a Wilcoxon-Mann-Whitney test) for independent samples

## difference between group 1 and group 2 is not equal to zero
## estimated difference is Cohen'd = 0.25
power.np.wilcoxon(d = 0.25,
                power = 0.80)

## difference between group 1 and group 2 is greater than zero
## estimated difference is Cohen'd = 0.25
power.np.wilcoxon(d = 0.25,
                power = 0.80,
                alternative = "one.sided")

## mean of group 1 is practically not smaller than mean of group 2
## estimated difference is Cohen'd = 0.10 and can be as small as -0.05
power.np.wilcoxon(d = 0.10,
                margin = -0.05,
                power = 0.80,
                alternative = "one.sided")

## mean of group 1 is practically greater than mean of group 2
## estimated difference is Cohen'd = 0.10 and can be as small as 0.05
power.np.wilcoxon(d = 0.10,
                margin = 0.05,
                power = 0.80,
                alternative = "one.sided")

## mean of group 1 is practically same as mean of group 2
## estimated difference is Cohen'd = 0
## and can be as small as -0.05 and as high as 0.05
power.np.wilcoxon(d = 0,
                margin = c(-0.05, 0.05),
                power = 0.80,
                alternative = "two.one.sided")


# Wilcoxon signed-rank test for matched pairs (dependent samples)

## difference between time 1 and time 2 is not equal to zero
## estimated difference between time 1 and time 2 is Cohen'd = -0.25
power.np.wilcoxon(d = -0.25,
                power = 0.80,
                design = "paired")

## difference between time 1 and time 2 is greater than zero
## estimated difference between time 1 and time 2 is Cohen'd = -0.25
power.np.wilcoxon(d = -0.25,
                power = 0.80,
                design = "paired",
                alternative = "one.sided")

## mean of time 1 is practically not smaller than mean of time 2
## estimated difference is Cohen'd = -0.10 and can be as small as 0.05
power.np.wilcoxon(d = -0.10,
                margin = 0.05,
                power = 0.80,
                design = "paired",
                alternative = "one.sided")

## mean of time 1 is practically greater than mean of time 2
## estimated difference is Cohen'd = -0.10 and can be as small as -0.05
power.np.wilcoxon(d = -0.10,
                margin = -0.05,
                power = 0.80,
                design = "paired",
                alternative = "one.sided")

## mean of time 1 is practically same as mean of time 2
## estimated difference is Cohen'd = 0
## and can be as small as -0.05 and as high as 0.05
power.np.wilcoxon(d = 0,
                margin = c(-0.05, 0.05),
                power = 0.80,
                design = "paired",
                alternative = "two.one.sided")
}

\references{
Al-Sunduqchi, M. S. (1990). Determining the appropriate sample size for inferences based on the Wilcoxon statistics [Unpublished doctoral dissertation]. University of Wyoming - Laramie

Chow, S. C., Shao, J., Wang, H., and Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Lehmann, E. (1975). Nonparameterics: Statistical methods based on ranks. McGraw-Hill.

Noether, G. E. (1987). Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82(1), 645-647.

Ruscio, J. (2008). A probability-based measure of effect size: Robustness to base rates and other factors. Psychological Methods, 13(1), 19-30.

Ruscio, J., & Mullen, T. (2012). Confidence intervals for the probability of superiority effect size measure and the area under a receiver operating characteristic curve. Multivariate Behavioral Research, 47(2), 201-223.

Zhao, Y.D., Rahardja, D., & Qu, Y. (2008). Sample size calculation for the Wilcoxon-Mann-Whitney test adjusting for ties. Statistics in Medicine, 27(3), 462-468.
}
