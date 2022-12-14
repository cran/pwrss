\name{pwrss.z.med}
\alias{pwrss.z.med}

\title{Test of an Indirect Effect in the Mediation Analysis (z Test, Joint Test, and Monte Carlo)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test indirect effects in mediation analysis. One can consider explanatory power of the covariates in the mediation and outcome model via specifying R-squared values accordingly.
}

\usage{
pwrss.z.med(a, b, cp = 0,
            sdx = 1, sdm = 1, sdy = 1,
            r2m.x = a^2 * sdx^2 / sdm^2,
            r2y.mx = (b^2 * sdm^2 + cp^2 * sdx^2) / sdy^2,
            n = NULL, power = NULL, alpha = 0.05,
            alternative = c("not equal", "less", "greater"),
            mc = TRUE, nsims = 1000, ndraws = 1000,
            verbose = TRUE)
}

\arguments{
  \item{a}{expected regression coefficient for X -> M path. One can use standardized regression coefficient, but should use \code{sdx = 1} and \code{sdm = 1} or leave them out as they are default specifications}
  \item{b}{expected regression coefficient for M -> Y path. One can use standardized regression coefficient, but should use \code{sdm = 1} and \code{sdy = 1} or leave them out as they are default specifications}
  \item{cp}{expected regression coefficient for X -> Y path. One can use standardized regression coefficient, but should use \code{sdx = 1} and \code{sdy = 1} or leave them out as they are default specifications}
  \item{sdx}{standard deviation of the predictor (X). For a binary predictor, \code{sdx = sqrt(p*(1-p))} where\code{p} is the proportion of subjects in one of the groups}
  \item{sdm}{standard deviation of the mediator (M)}
  \item{sdy}{standard deviation of the outcome (Y)}
  \item{r2m.x}{expected R-squared value for the mediator model (M ~ X). The default is \code{r2m.x = a^2 * sdx^2 / sdm^2} assuming that X is the only predictor. Thus, an \code{r2m.x} below this value will throw a warning. To consider other covariates in the model provide a value greater than the default.}
  \item{r2y.mx}{expected R-squared value for the outcome model (Y ~ M + X). The default is \code{r2y.mx = (b^2 * sdm^2 + cp^2 * sdx^2) / sdy^2} assuming that M and X are the only predictors. Thus, an \code{r2y.mx} below this value will throw a warning. To consider other covariates in the model provide a value greater than the default.}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent"}
  \item{mc}{logical; TRUE for monte carlo simulation based power}
  \item{nsims}{number of replications, if mc = TRUE}
  \item{ndraws}{number of draws from the distribution of the path coefficients for each replication, if mc = TRUE}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}

\examples{
# with standardized coefficients
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10, n = 200)
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10, power = 0.80)

# with binary predictor X such as treatment/control variable
# in this case standardized coefficients for path a and cp would be Cohen's d values
p <- 0.50 # proportion of subjects in one group
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            sdx = sqrt(p*(1-p)), n = 200)
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            sdx = sqrt(p*(1-p)), power = 0.80)

}

\references{
Add references.
}
