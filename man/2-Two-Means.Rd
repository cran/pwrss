\name{pwrss.t.2means}
\alias{pwrss.t.2means}
\alias{pwrss.z.2means}

\title{Difference between Two Means (t or z Test for Independent or Paired Samples)}

\description{
  Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test difference between two means.
  For standardized mean difference (Cohen's d) set \code{mu1 = d} and use defaults for \code{mu2}, \code{sd1}, and \code{sd2}.
  If pooled standard deviation (psd) is available set \code{sd1 = psd}.

  Formulas are validated using Monte Carlo simulation, G*Power, \url{http://powerandsamplesize.com/}, and tables in PASS documentation.
}

\usage{
pwrss.t.2means(mu1, mu2 = 0, margin = 0,
               sd1 = ifelse(paired, sqrt(1/(2*(1-paired.r))), 1), sd2 = sd1,
               kappa = 1, paired = FALSE, paired.r = 0.50,
               alpha = 0.05, welch.df = FALSE,
               alternative = c("not equal","greater","less",
                               "equivalent","non-inferior","superior"),
               n2 = NULL, power = NULL, verbose = TRUE)

pwrss.z.2means(mu1, mu2 = 0, sd1 = 1, sd2 = sd1, margin = 0,
               kappa = 1, alpha = 0.05,
               alternative = c("not equal", "greater", "less",
                               "equivalent", "non-inferior", "superior"),
               n2 = NULL, power = NULL, verbose = TRUE)

}

\arguments{
  \item{mu1}{expected mean in the first group}
  \item{mu2}{expected mean in the second group}
  \item{sd1}{expected standard deviation in the first group}
  \item{sd2}{expected standard deviation in the second group}
  \item{paired}{if \code{TRUE} paired samples t test}
  \item{paired.r}{correlation between repeated measures for paired samples (e.g., pretest and posttest)}
  \item{n2}{sample size in the second group (or for the single group in paired samples)}
  \item{kappa}{\code{n1/n2} ratio (applies to independent samples only)}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{welch.df}{if \code{TRUE} uses Welch's degrees of freedom adjustment when groups sizes or variances are not equal (applies to independent samples t test only)}
  \item{margin}{non-inferority, superiority, or equivalence margin (margin: boundry of \code{mu1 - mu2} that is practically insignificant)}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent", "non-inferior", or "superior"}
  \item{verbose}{if \code{FALSE} no output is printed on the console}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z or t test)}
  \item{df}{degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}


\examples{
# independent samples t test

## difference between group 1 and group 2 is not equal to zero
## estimated difference is Cohen'd = 0.25
pwrss.t.2means(mu1 = 0.25, mu2 = 0, power = 0.80,
                alternative = "not equal")

## difference between group 1 and group 2 is greater than zero
## estimated difference is Cohen'd = 0.25
pwrss.t.2means(mu1 = 0.25, mu2 = 0, power = 0.80,
                alternative = "greater")

## mean of group 1 is practically not smaller than mean of group 2
## estimated difference is Cohen'd = 0.10 and can be as small as -0.05
pwrss.t.2means(mu1 = 0.25, mu2 = 0.15,
                margin = -0.05, power = 0.80,
                alternative = "non-inferior")

## mean of group 1 is practically greater than mean of group 2
## estimated difference is Cohen'd = 0.10 and can be as small as 0.05
pwrss.t.2means(mu1 = 0.25, mu2 = 0.15,
                margin = 0.05, power = 0.80,
                alternative = "superior")

## mean of group 1 is practically same as mean of group 2
## estimated difference is Cohen'd = 0
## and can be as small as -0.05 and as high as 0.05
pwrss.t.2means(mu1 = 0.25, mu2 = 0.25,
                margin = 0.05, power = 0.80,
                alternative = "equivalent")


#  dependent samples (matched pairs) t test

## difference between time 1 and time 2 is not equal to zero
## estimated difference between time 1 and time 2 is Cohen'd = -0.25
pwrss.t.2means(mu1 = 0, mu2 = 0.25, power = 0.80,
                paired = TRUE, paired.r = 0.50,
                alternative = "not equal")

## difference between time 1 and time 2 is less than zero
## estimated difference between time 1 and time 2 is Cohen'd = -0.25
pwrss.t.2means(mu1 = 0, mu2 = 0.25, power = 0.80,
                paired = TRUE, paired.r = 0.50,
                alternative = "less")

## mean of time 1 is practically not smaller than mean of time 2
## estimated difference is Cohen'd = -0.10 and can be as small as 0.05
pwrss.t.2means(mu1 = 0.15, mu2 = 0.25, margin = 0.05,
                paired = TRUE, paired.r = 0.50, power = 0.80,
                alternative = "non-inferior")

## mean of time 1 is practically greater than mean of time 2
## estimated difference is Cohen'd = -0.10 and can be as small as -0.05
pwrss.t.2means(mu1 = 0.15, mu2 = 0.25, margin = -0.05,
                paired = TRUE, paired.r = 0.50, power = 0.80,
                alternative = "superior")

## mean of time 1 is practically same as mean of time 2
## estimated difference is Cohen'd = 0
## and can be as small as -0.05 and as high as 0.05
pwrss.t.2means(mu1 = 0.25, mu2 = 0.25, margin = 0.05,
                paired = TRUE, paired.r = 0.50, power = 0.80,
                alternative = "equivalent")

}

\references{
Bulus, M., & Polat, C. (in press). pwrss R paketi ile istatistiksel guc analizi [Statistical power analysis with pwrss R package]. Ahi Evran Universitesi Kirsehir Egitim Fakultesi Dergisi. \url{https://osf.io/ua5fc/download/}

Chow, S. C., Shao, J., Wang, H., & Lokhnygina, Y. (2018). Sample size calculations in clinical research (3rd ed.). Taylor & Francis/CRC.

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
