\name{power.t.student}
\alias{power.t.student}
\alias{power.t.welch}
\alias{pwrss.t.2means}
\alias{pwrss.z.2means}
\alias{pwrss.t.mean}
\alias{pwrss.z.mean}

\title{Power Analysis for Student's and Welch's T-Tests}

\description{
  Calculates power or sample size (only one can be NULL at a time) for Student's and Welch's T-Tests. Welch's T-Test implementation relies on formulas proposed by Bulus (2024).

  Use \code{means.to.d()} to convert raw means and standard deviations to Cohen's d, and \code{d.to.cles()} to convert Cohen's d to the probability of superiority. Note that this interpretation is appropriate only when the underlying distribution is approximately normal and the two groups have similar population variances.

 In contrast to previous versions, users can now specify whether their claims will be based on raw score mean difference with P-values or standardized mean difference with confidence intervals. While results typically differ by only a few units, these distinctions can be particularly consequential in studies with small sample sizes or high-risk interventions.

  Formulas are validated using Monte Carlo simulations (see Bulus, 2024), G*Power, \url{http://powerandsamplesize.com/}, and tables in PASS documentation. One key difference between PASS and \code{pwrss} lies in how they handle non-inferiority and superiority tests-that is, one-sided tests defined by a negligible effect margin (implemented as of this version). PASS shifts the test statistic so that the null hypothesis assumes a zero effect, treating the negligible margin as part of the alternative hypothesis. As a result, the test statistic is evaluated against a central distribution. In contrast, \code{pwrss} treats the negligible effect as the true null value, and the test statistic is evaluated under a non-central distribution. This leads to slight differences up to third decimal place. To get the same results, reflect the margin in \code{null.d} and specify \code{margin = 0}.

   Equivalence tests are implemented in line with Bulus and Polat (2023), Chow et al. (2018) and Lakens (2017).

  NOTE: The functions \code{pwrss.z.mean()} and \code{pwrss.z.2means()} are no longer supported. The \code{pwrss.t.mean()} and \code{pwrss.t.2means()} functions are deprecated, but they will remain available as wrappers for \code{power.t.student()} or \code{power.t.welch()} during the transition period.
}

\usage{
power.t.student(d, null.d = 0, margin = 0,
                n2 = NULL, n.ratio = 1, power = NULL, alpha = 0.05,
                alternative = c("two.sided", "one.sided", "two.one.sided"),
                design = c("independent", "paired", "one.sample"),
                claim.basis = c("md.pval", "smd.ci"),
                ceiling = TRUE, verbose = TRUE, pretty = FALSE)

power.t.welch(d, null.d = 0, margin = 0,
              var.ratio = 1, n.ratio = 1, n2 = NULL,
              power = NULL, alpha = 0.05,
              alternative = c("two.sided", "one.sided", "two.one.sided"),
              claim.basis = c("md.pval", "smd.ci"),
              ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{d}{Cohen's d or Hedges' g.}
  \item{null.d}{Cohen's d or Hedges' g under null, typically 0(zero).}
  \item{margin}{margin - ignorable \code{d} - \code{null.d} difference.}
  \item{var.ratio}{variance ratio in the form of sd1^2 / sd2^2.}
  \item{n2}{integer; sample size in the second group (or for the single group in paired samples or one-sample).}
  \item{n.ratio}{\code{n1/n2} ratio (applies to independent samples only)}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; the direction or type of the hypothesis test: "one.sided", "two.sided", "two.one.sided". For non-inferiority or superiority tests, add or subtract the margin from the null hypothesis value and use \code{alternative = "one.sided"}.}
  \item{design}{character; "independent", "paired" or "one.sample".}
  \item{claim.basis}{character; "md.pval" when claims are based on raw mean differences and p-values, "smd.ci" when claims are based on standardized mean differences and confidence intervals.}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (T-Test).}
  \item{df}{degrees of freedom.}
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null.}
  \item{t.alpha}{critical value(s).}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n}{sample size (`n` or `c(n1, n2)` depending on the design.}
}


\examples{
#######################
# Independent Samples #
#######################

## difference between group 1 and group 2 is not equal to zero
## targeting minimal difference of Cohen'd = 0.20
## non-parametric
power.np.wilcoxon(d = 0.20,
                  power = 0.80,
                  alternative = "two.sided",
                  design = "independent")

## parametric
power.t.student(d = 0.20,
                power = 0.80,
                alternative = "two.sided",
                design = "independent")

## when sample size ratio and group variances differ
power.t.welch(d = 0.20,
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alternative = "two.sided")


## difference between group 1 and group 2 is greater than zero
## targeting minimal difference of Cohen'd = 0.20
## non-parametric
power.np.wilcoxon(d = 0.20,
                  power = 0.80,
                  alternative = "one.sided",
                  design = "independent")

## parametric
power.t.student(d = 0.20,
                power = 0.80,
                alternative = "one.sided",
                design = "independent")

## when sample size ratio and group variances differ
power.t.welch(d = 0.20,
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alternative = "one.sided")


## mean of group 1 is practically not smaller than mean of group 2
## targeting minimal difference of Cohen'd = 0.20 and can be as small as -0.05
## non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = -0.05,
                  power = 0.80,
                  alternative = "one.sided",
                  design = "independent")

## parametric
power.t.student(d = 0.20,
                margin = -0.05,
                power = 0.80,
                alternative = "one.sided",
                design = "independent")

## when sample size ratio and group variances differ
power.t.welch(d = 0.20,
              margin = -0.05,
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alternative = "one.sided")


## mean of group 1 is practically greater than mean of group 2
## targeting minimal difference of Cohen'd = 0.20 and can be as small as 0.05
## non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = 0.05,
                  power = 0.80,
                  alternative = "one.sided",
                  design = "independent")

## parametric
power.t.student(d = 0.20,
                margin = 0.05,
                power = 0.80,
                alternative = "one.sided",
                design = "independent")

## when sample size ratio and group variances differ
power.t.welch(d = 0.20,
              margin = 0.05,
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alternative = "one.sided")


## mean of group 1 is practically same as mean of group 2
## targeting minimal difference of Cohen'd = 0
## and can be as small as -0.05 or as high as 0.05
## non-parametric
power.np.wilcoxon(d = 0,
                  margin = c(-0.05, 0.05),
                  power = 0.80,
                  alternative = "two.one.sided",
                  design = "independent")

## parametric
power.t.student(d = 0,
                margin = c(-0.05, 0.05),
                power = 0.80,
                alternative = "two.one.sided",
                design = "independent")

## when sample size ratio and group variances differ
power.t.welch(d = 0,
              margin = c(-0.05, 0.05),
              n.ratio = 2,
              var.ratio = 2,
              power = 0.80,
              alternative = "two.one.sided")


##################
# Paired Samples #
##################

## difference between time 1 and time 2 is not equal to zero
## targeting minimal difference of Cohen'd = -0.20
## non-parametric
power.np.wilcoxon(d = -0.20,
                  power = 0.80,
                  alternative = "two.sided",
                  design = "paired")

## parametric
power.t.student(d = -0.20,
                power = 0.80,
                alternative = "two.sided",
                design = "paired")

## difference between time 1 and time 2 is less than zero
## targeting minimal difference of Cohen'd = -0.20
## non-parametric
power.np.wilcoxon(d = -0.20,
                  power = 0.80,
                  alternative = "one.sided",
                  design = "paired")

## parametric
power.t.student(d = -0.20,
                power = 0.80,
                alternative = "one.sided",
                design = "paired")

## mean of time 1 is practically not greater than mean of time 2
## targeting minimal difference of Cohen'd = -0.20 and can be as small as 0.05
## non-parametric
## non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = 0.05,
                  power = 0.80,
                  alternative = "one.sided",
                  design = "paired")

## parametric
power.t.student(d = 0.20,
                margin = 0.05,
                power = 0.80,
                alternative = "one.sided",
                design = "paired")

## mean of time 1 is practically greater than mean of time 2
## targeting minimal difference of Cohen'd = -0.20 and can be as small as -0.05
## non-parametric
power.np.wilcoxon(d = 0.20,
                  margin = -0.05,
                  power = 0.80,
                  alternative = "one.sided",
                  design = "paired")

## parametric
power.t.student(d = 0.20,
                margin = -0.05,
                power = 0.80,
                alternative = "one.sided",
                design = "paired")


## mean of time 1 is practically same as mean of time 2
## targeting minimal difference of Cohen'd = 0
## and can be as small as -0.05 or as high as 0.05
## non-parametric
## non-parametric
power.np.wilcoxon(d = 0,
                  margin = c(-0.05, 0.05),
                  power = 0.80,
                  alternative = "two.one.sided",
                  design = "paired")

## parametric
power.t.student(d = 0,
                margin = c(-0.05, 0.05),
                power = 0.80,
                alternative = "two.one.sided",
                design = "paired")
}

\references{
Bulus, M. (2024). Robust standard errors and confidence intervals for standardized mean difference [Preprint].\doi{10.31219/osf.io/k6mbs}

Bulus, M., & Polat, C. (2023). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi, 24(3), 2207-2328. \doi{10.29299/kefad.1209913}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.

Lakens, D. (2017). Equivalence tests: A practical primer for t tests, correlations, and meta-analyses. Social psychological and personality science, 8(4), 355-362. \doi{10.1177/1948550617697177}
}
