\name{pwrss.z.logreg}
\alias{pwrss.z.logreg}
\alias{pwrss.z.logistic}

\title{Logistic Regression: Single Coefficient (Large Sample Approx. Wald's z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a single coefficient in logistic regression. \code{pwrss.z.logistic()} and \code{pwrss.z.logreg()} are the same functions. The distribution of the predictor variable can be one of the following: \code{c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")} for Demidenko (2007) procedure but only \code{c("normal", "binomial", "bernouilli")} for Hsieh et al. (1998) procedure. The default parameters for these distributions are

\code{distribution = list(dist = "normal", mean = 0, sd = 1)} \cr
\code{distribution = list(dist = "poisson", lambda = 1)} \cr
\code{distribution = list(dist = "uniform", min = 0, max = 1)} \cr
\code{distribution = list(dist = "exponential", rate = 1)} \cr
\code{distribution = list(dist = "binomial", size = 1, prob = 0.50)} \cr
\code{distribution = list(dist = "bernoulli", prob = 0.50)} \cr
\code{distribution = list(dist = "lognormal", meanlog = 0, sdlog = 1)} \cr

Parameters defined in \code{list()} form can be modified, but the names should be kept the same. It is sufficient to use distribution's name for default parameters (e.g. \code{dist = "normal"}).


Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.
}

\usage{
pwrss.z.logreg(p1 = 0.10, p0 = 0.15,
               odds.ratio  = (p1/(1-p1))/(p0/(1-p0)),
               beta0 = log(p0/(1-p0)), beta1 = log(odds.ratio),
               n = NULL, power = NULL, r2.other.x = 0, alpha = 0.05,
               alternative = c("not equal", "less", "greater"),
               method = c("demidenko(vc)", "demidenko", "hsieh"),
               distribution = "normal", verbose = TRUE)


pwrss.z.logistic(p1 = 0.10, p0 = 0.15,
                 odds.ratio  = (p1/(1-p1))/(p0/(1-p0)),
                 beta0 = log(p0/(1-p0)), beta1 = log(odds.ratio),
                 n = NULL, power = NULL, r2.other.x = 0, alpha = 0.05,
                 alternative = c("not equal", "less", "greater"),
                 method = c("demidenko(vc)", "demidenko", "hsieh"),
                 distribution = "normal", verbose = TRUE)
}

\arguments{
  \item{p0}{base probability under null hypothesis (probability that an event occurs without the influence of the predictor X - or when the value of the predictor is zero)}
  \item{p1}{probability under alternative hypothesis (probability that an event occurs when the value of the predictor X is increased by one unit)}
  \item{beta0}{regression coefficient defined as \cr \code{beta0 = log( p0/(1-p0) )}}
  \item{beta1}{regression coefficient for the predictor X defined as \cr \code{beta1 = log( (p1/(1-p1)) / (p0/(1-p0)) )}}
  \item{odds.ratio}{odds ratio defined as \cr \code{odds.ratio = exp(beta1) = (p1/(1-p1)) / (p0/(1-p0))}}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{r2.other.x}{proportion of variance in the predictor X explained by other covariates. Not to be confused with pseudo R-squared}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less"}
  \item{method}{calculation method. \code{"demidenko(vc)"} stands for Demidenko (2007) procedure with variance correction; \code{"demidenko"} stands for Demidenko (2007) procedure without variance correction; \code{"hsieh"} stands for Hsieh et al. (1998) procedure. \code{"demidenko"} and \code{"hsieh"} methods produce similiar results but \code{"demidenko(vc)"} is more precise}
  \item{distribution}{distribution family. Can be one of the \code{c("noramal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")} for Demidenko (2007) procedure but only \code{c("normal", "binomial", "bernouilli")} for Hsieh et al. (1998) procedure}
  \item{verbose}{if \code{FALSE} no output is printed on the console}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z test)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{total sample size}
}

\examples{

# predictor X follows normal distribution

## probability specification
pwrss.z.logreg(p0 = 0.15, p1 = 0.10,
               alpha = 0.05, power = 0.80,
               dist = "normal")

## odds ratio specification
pwrss.z.logreg(p0 = 0.15, odds.ratio = 0.6296,
               alpha = 0.05, power = 0.80,
               dist = "normal")

## regression coefficient specification
pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626,
               alpha = 0.05, power = 0.80,
               dist = "normal")

## change parameters associated with predictor X
dist.x <- list(dist = "normal", mean = 10, sd = 2)
pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626,
               alpha = 0.05, power = 0.80,
               dist = dist.x)


# predictor X follows Bernoulli distribution (such as treatment/control groups)

## probability specification
pwrss.z.logreg(p0 = 0.15, p1 = 0.10,
               alpha = 0.05, power = 0.80,
               dist = "bernoulli")

## odds ratio specification
pwrss.z.logreg(p0 = 0.15, odds.ratio = 0.6296,
               alpha = 0.05, power = 0.80,
               dist = "bernoulli")

## regression coefficient specification
pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626,
               alpha = 0.05, power = 0.80,
               dist = "bernoulli")

## change parameters associated with predictor X
dist.x <- list(dist = "bernoulli", prob = 0.30)
pwrss.z.logreg(p0 = 0.15, beta1 = -0.4626,
               alpha = 0.05, power = 0.80,
               dist = dist.x)
}

\references{
Demidenko, E. (2007). Sample size determination for logistic regression revisited. Statistics in Medicine, 26(18), 3385-3397.

Hsieh, F. Y., Bloch, D. A., & Larsen, M. D. (1998). A simple method of sample size calculation for linear and logistic regression. Statistics in Medicine, 17(4), 1623-1634.
}
