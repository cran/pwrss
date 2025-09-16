\name{power.z.logistic}
\alias{power.z.logistic}
\alias{pwrss.z.logistic}
\alias{power.z.logreg}
\alias{pwrss.z.logreg}

\title{Power Analysis for Logistic Regression Coefficient (Wald's Z-Test)}

\description{
Calculates power or sample size (only one can be NULL at a time) to test a single coefficient in logistic regression. \code{power.z.logistic()} and \code{power.z.logreg()} are the same functions, as well as \code{pwrss.z.logistic()} and \code{pwrss.z.logreg()}.

The distribution of the predictor variable can be one of the following: \code{c("normal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")} for Demidenko (2007) procedure but only \code{c("normal", "binomial", "bernouilli")} for Hsieh et al. (1998) procedure. The default parameters for these distributions are

\code{distribution = list(dist = "normal", mean = 0, sd = 1)} \cr
\code{distribution = list(dist = "poisson", lambda = 1)} \cr
\code{distribution = list(dist = "uniform", min = 0, max = 1)} \cr
\code{distribution = list(dist = "exponential", rate = 1)} \cr
\code{distribution = list(dist = "binomial", size = 1, prob = 0.50)} \cr
\code{distribution = list(dist = "bernoulli", prob = 0.50)} \cr
\code{distribution = list(dist = "lognormal", meanlog = 0, sdlog = 1)} \cr

Parameters defined in \code{list()} form can be modified, but element names should be kept the same. It is sufficient to use distribution's name for default parameters (e.g. \code{dist = "normal"}).

NOTE: The \code{pwrss.z.logistic()} and its alias \code{pwrss.z.logreg()} are deprecated. However, they will remain available as wrappers for the \code{power.z.logistic()} function.

Formulas are validated using G*Power and tables in PASS documentation.
}

\usage{
power.z.logistic(prob = NULL, base.prob = NULL,
                 odds.ratio  = (prob/(1-prob))/(base.prob/(1-base.prob)),
                 beta0 = log(base.prob/(1-base.prob)), beta1 = log(odds.ratio),
                 n = NULL, power = NULL, r.squared.predictor = 0,
                 alpha = 0.05, alternative = c("two.sided", "one.sided"),
                 method = c("demidenko(vc)", "demidenko", "hsieh"),
                 distribution = "normal", ceiling = TRUE,
                 verbose = TRUE, pretty = FALSE)
}

\arguments{
  \item{base.prob}{base probability under null hypothesis (probability that an event occurs without the influence of the predictor - or when the value of the predictor is zero).}
  \item{prob}{probability under alternative hypothesis (probability that an event occurs when the value of the predictor is increased from 0 to 1). Warning: This is base probability + incremental increase.}
  \item{beta0}{regression coefficient defined as \cr \code{beta0 = log(base.prob/(1-base.prob))}}
  \item{beta1}{regression coefficient for the predictor X defined as \cr \code{beta1 = log((prob/(1-prob))/(base.prob/(1-base.prob)))}}
  \item{odds.ratio}{odds ratio defined as \cr \code{odds.ratio = exp(beta1) = (prob/(1-prob))/(base.prob/(1-base.prob))}}
  \item{n}{integer; sample size}
   \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{r.squared.predictor}{proportion of variance in the predictor accounted for by other covariates. This is not a pseudo R-squared. To compute it, regress the predictor on the covariates and extract the adjusted R-squared from that model.}
   \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{alternative}{character; direction or type of the hypothesis test: "not equal", "greater", "less"}
  \item{method}{character; analytic method. \code{"demidenko(vc)"} stands for Demidenko (2007) procedure with variance correction; \code{"demidenko"} stands for Demidenko (2007) procedure without variance correction; \code{"hsieh"} stands for Hsieh et al. (1998) procedure. \code{"demidenko"} and \code{"hsieh"} methods produce similar results but \code{"demidenko(vc)"} is more precise}
  \item{distribution}{character; distribution family. Can be one of the \code{c("noramal", "poisson", "uniform", "exponential", "binomial", "bernouilli", "lognormal")} for Demidenko (2007) procedure but only \code{c("normal", "binomial", "bernouilli")} for Hsieh et al. (1998) procedure.}
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
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n}{sample size.}
}

\examples{

###########################################
# predictor X follows normal distribution #
###########################################

## probability specification
power.z.logistic(base.prob = 0.15, prob = 0.20,
                 alpha = 0.05, power = 0.80,
                 dist = "normal")

## odds ratio specification
power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667,
                 alpha = 0.05, power = 0.80,
                 dist = "normal")

## regression coefficient specification
power.z.logistic(beta0 = -1.734601, beta1 = 0.3483067,
                 alpha = 0.05, power = 0.80,
                 dist = "normal")

## change parameters associated with predictor X
pred.dist <- list(dist = "normal", mean = 10, sd = 2)
power.z.logistic(base.prob = 0.15, beta1 = 0.3483067,
                 alpha = 0.05, power = 0.80,
                 dist = pred.dist)


##############################################
# predictor X follows Bernoulli distribution #
# (such as treatment/control groups)         #
##############################################

## odds ratio specification
power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667,
                 alpha = 0.05, power = 0.80,
                 dist = "bernoulli")

## change parameters associated with predictor X
pred.dist <- list(dist = "bernoulli", prob = 0.30)
power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667,
                 alpha = 0.05, power = 0.80,
                 dist = pred.dist)

####################################
# predictor X is an ordinal factor #
####################################

## generating an ordinal predictor
x.ord <- sample(
  x = c(1, 2, 3, 4), # levels
  size = 1e5, # sample size large enough to get stable estimates
  prob = c(0.25, 0.25, 0.25, 0.25), # category probabilities
  replace = TRUE
)

## dummy coding the ordinal predictor
x.ord <- factor(x.ord, ordered = TRUE)
contrasts(x.ord) <- contr.treatment(4, base = 4)
x.dummy <- model.matrix( ~ x.ord)[,-1]
x.data <- as.data.frame(x.dummy)

## fit linear regression to get multiple r-squared
x.fit <- lm(x.ord1 ~ x.ord2 + x.ord3, data = x.data)

## extract parameters
bern.prob <- mean(x.data$x.ord1)
r.squared.pred <- summary(x.fit)$adj.r.squared

## change parameters associated with predictor X
pred.dist <- list(dist = "bernoulli", prob = bern.prob)
power.z.logistic(base.prob = 0.15, odds.ratio = 1.416667,
               alpha = 0.05, power = 0.80,
               r.squared.pred = r.squared.pred,
               dist = pred.dist)

}

\references{

Demidenko, E. (2007). Sample size determination for logistic regression revisited. Statistics in Medicine, 26(18), 3385-3397. \doi{10.1002/sim.2771}

Hsieh, F. Y., Bloch, D. A., & Larsen, M. D. (1998). A simple method of sample size calculation for linear and logistic regression. Statistics in Medicine, 17(4), 1623-1634.
}
