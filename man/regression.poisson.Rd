\name{power.z.poisson}
\alias{power.z.poisreg}
\alias{power.z.poisson}
\alias{pwrss.z.poisreg}
\alias{pwrss.z.poisson}

\title{Power Analysis for Poisson Regression Coefficient (Wald's z Test)}

\description{
Calculates power or sample size (only one can be NULL at a time) to test a single coefficient in poisson regression. \code{power.z.poisson()} and \code{power.z.poisreg()} are the same functions, as well as \code{pwrss.z.poisson()} and \code{pwrss.z.poisreg()}. The distribution of the predictor variable can be one of the following: \code{c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")}. The default parameters for these distributions are

\code{distribution = list(dist = "normal", mean = 0, sd = 1)} \cr
\code{distribution = list(dist = "poisson", lambda = 1)} \cr
\code{distribution = list(dist = "uniform", min = 0, max = 1)} \cr
\code{distribution = list(dist = "exponential", rate = 1)} \cr
\code{distribution = list(dist = "binomial", size = 1, prob = 0.50)} \cr
\code{distribution = list(dist = "bernoulli", prob = 0.50)} \cr
\code{distribution = list(dist = "lognormal", meanlog = 0, sdlog = 1)} \cr

Parameters defined in \code{list()} form can be modified, but the names should be kept the same. It is sufficient to use distribution's name for default parameters (e.g. \code{dist = "normal"}).

Formulas are validated using Monte Carlo simulation, G*Power, and tables in PASS documentation.

NOTE: The \code{pwrss.z.poisson()} and its alias \code{pwrss.z.poisreg()} are deprecated. However, they will remain available as wrappers for the \code{power.z.logistic()} function.
}

\usage{
power.z.poisson(base.rate = NULL, rate.ratio = NULL,
                beta0 = log(base.rate), beta1 = log(rate.ratio),
                n = NULL, power = NULL,
                r.squared.predictor = 0, mean.exposure = 1,
                alpha = 0.05, alternative = c("two.sided", "one.sided"),
                method = c("demidenko(vc)", "demidenko", "signorini"),
                distribution = "normal", ceiling = TRUE,
                verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{base.rate}{the base mean event rate.}
  \item{rate.ratio}{event rate ratio. The relative increase in the mean event rate for one unit increase in the predictor (similar to odds ratio in logistic regression).}
  \item{beta0}{\code{log(base.rate)} or natural logarithm of the base mean event rate.}
  \item{beta1}{\code{log(rate.ratio)} or natural logarithm of the relative increase in the mean event rate for one unit increase in the predictor.}
  \item{mean.exposure}{the mean exposure time (should be > 0). Usually 1}
  \item{n}{integer; sample size.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{r.squared.predictor}{proportion of variance in the predictor accounted for by other covariates. This is not a pseudo R-squared. To compute it, regress the predictor on the covariates and extract the adjusted R-squared from that model.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; direction or type of the hypothesis test: "not equal", "greater", "less".}
  \item{method}{character; calculation method. \code{"demidenko(vc)"} stands for Demidenko (2007) procedure with variance correction; \code{"demidenko"} stands for Demidenko (2007) procedure without variance correction; \code{"signorini"} stands for Signorini (1991) procedure. \code{"demidenko"} and \code{"signorini"} methods produce similar results but \code{"demidenko(vc)"} is more precise.}
  \item{distribution}{character; distribution family. Can be one of the \code{c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")}.}
  \item{ceiling}{logical; whether sample size should be rounded up. \code{TRUE} by default.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (Z-Test).}
  \item{mean}{mean of the alternative distribution.}
  \item{sd}{standard deviation of the alternative distribution.}
  \item{null.mean}{mean of the null distribution.}
  \item{null.sd}{standard deviation of the null distribution.}
  \item{z.alpha}{critical value(s).}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size.}
}

\examples{
# predictor X follows normal distribution

## regression coefficient specification
power.z.poisson(beta0 = 0.50, beta1 = -0.10,
                alpha = 0.05, power = 0.80,
                dist = "normal")

## rate ratio specification
power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = "normal")

## change parameters associated with predictor X
dist.x <- list(dist = "normal", mean = 10, sd = 2)
power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = dist.x)


# predictor X follows Bernoulli distribution (such as treatment/control groups)

## regression coefficient specification
power.z.poisson(beta0 = 0.50, beta1 = -0.10,
                alpha = 0.05, power = 0.80,
                dist = "bernoulli")

## rate ratio specification
power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = "bernoulli")

## change parameters associatied with predictor X
dist.x <- list(dist = "bernoulli", prob = 0.30)
power.z.poisson(base.rate = exp(0.50),
                rate.ratio = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = dist.x)
}

\references{
Demidenko, E. (2007). Sample size determination for logistic regression revisited. Statistics in Medicine, 26(18), 3385-3397. \doi{10.1002/sim.2771}

Signorini, D. F. (1991). Sample size for poisson regression. Biometrika, 78(2), 446-450.
}
