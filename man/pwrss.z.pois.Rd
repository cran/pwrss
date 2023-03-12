\name{Poisson Regression}
\alias{pwrss.z.poisreg}
\alias{pwrss.z.poisson}

\title{Test of a Single Coefficient in Poisson Regression (Large Sample Approx. Wald's z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a single coefficient in poisson regression. \code{pwrss.z.poisson()} and \code{pwrss.z.poisreg()} are the same functions. The distribution of the predictor variable can be one of the following: \code{c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")}. The default parameters for these distributions are

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
pwrss.z.poisreg(exp.beta0 = 1.10, exp.beta1 = 1.16,
                beta0 = log(exp.beta0), beta1 = log(exp.beta1),
                mean.exposure = 1, n = NULL, power = NULL, r2.other.x = 0,
                alpha = 0.05, alternative = c("not equal", "less", "greater"),
                method = c("demidenko(vc)", "demidenko", "signorini"),
                distribution = "normal", verbose = TRUE)


pwrss.z.poisson(exp.beta0 = 1.10, exp.beta1 = 1.16,
                beta0 = log(exp.beta0), beta1 = log(exp.beta1),
                mean.exposure = 1, n = NULL, power = NULL, r2.other.x = 0,
                alpha = 0.05, alternative = c("not equal", "less", "greater"),
                method = c("demidenko(vc)", "demidenko", "signorini"),
                distribution = "normal", verbose = TRUE)
}

\arguments{
  \item{exp.beta0}{the base event rate}
  \item{exp.beta1}{rate ratio: the relative increase of the event rate for one unit increase in the predictor X over the base event rate of exp(beta0) (similiar to odds ratio in logistic regression)}
  \item{beta0}{\code{log(exp.beta0)} or natural logarithm of the base event rate}
  \item{beta1}{\code{log(exp.beta1)} or natural logarithm of the relative increase of the event rate for one unit increase in the predictor X over the base event rate}
  \item{mean.exposure}{the mean exposure time (should be > 0). Usually 1}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{r2.other.x}{proportion of variance in the predictor X explained by other covariates. Not to be confused with explanatory power of covariates in the outcome model (pseudo R-squared)}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less"}
  \item{method}{calculation method. \code{"demidenko(vc)"} stands for Demidenko (2007) procedure with variance correction; \code{"demidenko"} stands for Demidenko (2007) procedure without variance correction; \code{"signorini"} stands for Signorini (1991) procedure. \code{"demidenko"} and \code{"signorini"} methods produce similiar results but \code{"demidenko(vc)"} is more precise}
  \item{distribution}{distribution family. Can be one of the \code{c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")}}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{total sample size}
}

\examples{
# predictor X follows normal distribution

## regression coefficient specification
pwrss.z.poisreg(beta0 = 0.50, beta1 = -0.10,
                alpha = 0.05, power = 0.80,
                dist = "normal")

## rate ratio specification
pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = "normal")

## change parameters associated with predictor X
dist.x <- list(dist = "normal", mean = 10, sd = 2)
pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = dist.x)


# predictor X follows Bernoulli distribution (such as treatment/control groups)

## regression coefficient specification
pwrss.z.poisreg(beta0 = 0.50, beta1 = -0.10,
                alpha = 0.05, power = 0.80,
                dist = "bernoulli")

## rate ratio specification
pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = "bernoulli")

## change parameters associatied with predictor X
dist.x <- list(dist = "bernoulli", prob = 0.30)
pwrss.z.poisreg(exp.beta0 = exp(0.50),
                exp.beta1 = exp(-0.10),
                alpha = 0.05, power = 0.80,
                dist = dist.x)
}

\references{
Demidenko, E. (2007). Sample size determination for logistic regression revisited. Statistics in Medicine, 26(18), 3385-3397.

Hsieh, F. Y., Bloch, D. A., & Larsen, M. D. (1998). A simple method of sample size calculation for linear and logistic regression. Statistics in Medicine, 17(4), 1623-1634.

Signorini, D. F. (1991). Sample size for poisson regression. Biometrika, 78(2), 446-450.
}
