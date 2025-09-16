\name{power.t.regression}
\alias{power.t.regression}
\alias{pwrss.t.regression}
\alias{pwrss.z.regression}
\alias{pwrss.t.reg}
\alias{pwrss.z.reg}
\alias{power.t.reg}

\title{Power Analysis for Linear Regression: Single Coefficient (T-Test)}

\description{
Calculates power or sample size (only one can be NULL at a time) to test a single coefficient in multiple linear regression. The predictor is assumed to be continuous by default. However, one can calculate power or sample size for a binary predictor (such as treatment and control groups in an experimental design) by specifying \code{sd.predictor = sqrt(p*(1-p))} where \code{p} is the proportion of subjects in one of the groups. The sample size in each group would be \code{n*p} and \code{n*(1-p)}. \code{power.t.regression()}\code{pwrss.t.regression()} are the same functions, as well as \code{power.t.reg()} and \code{pwrss.t.reg()}.

Minimal effect and equivalence tests are implemented in line with Hodges and Lehmann (1954), Kim and Robinson (2019), Phillips (1990), and Dupont and Plummer (1998).

Formulas are validated using Monte Carlo simulation, G*Power, tables in PASS documentation, and tables in Bulus (2021).

 NOTE: The \code{pwrss.t.regression()} function and its alias \code{pwrss.z.reg()} are deprecated, but they will remain available as a wrapper for \code{power.t.regression()} during the transition period.

}

\usage{
power.t.regression(beta, null.beta = 0, margin = 0,
                   sd.predictor = 1, sd.outcome = 1,
                   r.squared = (beta * sd.predictor / sd.outcome)^2,
                   k.total = 1, n = NULL, power = NULL, alpha = 0.05,
                   alternative = c("two.sided", "one.sided", "two.one.sided"),
                   ceiling = TRUE, verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{beta}{regression coefficient. One can use standardized regression coefficient, but should keep \code{sd.predictor = 1} and \code{sd.outcome = 1} or leave them out as they are default specifications.}
  \item{null.beta}{regression coefficient under null hypothesis (typically zero). One can use standardized regression coefficient, but should keep \code{sd.predictor = 1} and \code{sd.outcome = 1} or leave them out as they are default specifications.}
  \item{margin}{margin - ignorable \code{beta} - \code{null.beta} difference.}
  \item{sd.predictor}{standard deviation of the predictor. For a binary predictor, \code{sd.predictor = sqrt(p * (1 - p))} where \code{p} is the proportion of subjects in one of the groups.}
  \item{sd.outcome}{standard deviation of the outcome.}
  \item{k.total}{integer; total number of predictors, including the predictor of interest.}
  \item{r.squared}{model R-squared. The default is \code{r.squared = (beta * sd.predictor / sd.outcome)^2} assuming a linear regression with one predictor. Thus, an \code{r.squared} below this value will throw a warning. To consider other covariates in the model provide a value greater than the default \code{r.squared} along with the argument \code{k.total > 1}.}
  \item{n}{integer; sample size.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; the direction or type of the hypothesis test: "two.sided", "one.sided", or "two.one.sided".}
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
  \item{n}{sample size.}
}

\examples{
# continuous predictor x (and 4 covariates)
power.t.regression(beta = 0.20,
            k.total = 5,
            r.squared = 0.30,
            power = 0.80)

# binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
power.t.regression(beta = 0.20,
            sd.predictor = sqrt(p*(1-p)),
            k.total = 5,
            r.squared = 0.30,
            power = 0.80)

# non-inferiority test with binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
power.t.regression(beta = 0.20, # Cohen's d
            margin = -0.05, # non-inferiority margin in Cohen's d unit
            alternative = "one.sided",
            sd.predictor = sqrt(p*(1-p)),
            k.total = 5,
            r.squared = 0.30,
            power = 0.80)

# superiority test with binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
power.t.regression(beta = 0.20, # Cohen's d
            margin = 0.05, # superiority margin in Cohen's d unit
            alternative = "one.sided",
            sd.predictor = sqrt(p*(1-p)),
            k.total = 5,
            r.squared = 0.30,
            power = 0.80)

# equivalence test with binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
power.t.regression(beta = 0, # Cohen's d
            margin = c(-0.05, 0.05), # equivalence bounds in Cohen's d unit
            alternative = "two.one.sided",
            sd.predictor = sqrt(p*(1 - p)),
            k.total = 5,
            r.squared = 0.30,
            power = 0.80)
}

\references{
Bulus, M. (2021). Sample size determination and optimal design of randomized/non-equivalent pretest-post-test control-group designs. Adiyaman University Journal of Educational Sciences, 11(1), 48-69. \doi{10.17984/adyuebd.941434}

Hodges Jr, J. L., & Lehmann, E. L. (1954). Testing the approximate validity of statistical hypotheses. Journal of the Royal Statistical Society Series B: Statistical Methodology, 16(2), 261-268. \doi{10.1111/j.2517-6161.1954.tb00169.x}

Kim, J. H., & Robinson, A. P. (2019). Interval-based hypothesis testing and its applications to economics and finance. Econometrics, 7(2), 21. \doi{10.3390/econometrics7020021}

Phillips, K. F. (1990). Power of the two one-sided tests procedure in bioequivalence. Journal of Pharmacokinetics and Biopharmaceutics, 18(2), 137-144. \doi{10.1007/bf01063556}

Dupont, W. D., and Plummer, W. D. (1998). Power and sample size calculations for studies involving linear regression. Controlled Clinical Trials, 19(6), 589-601. \doi{10.1016/s0197-2456(98)00037-3}
}
