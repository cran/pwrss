\name{Indirect Effect in Mediation Models (z, Joint, Monte Carlo Tests)}
\alias{pwrss.z.med}
\alias{pwrss.z.mediation}

\title{Test of Indirect Effects in Mediation Analysis (z Test, Joint Test, and Monte Carlo Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test indirect effects in mediation analysis. One can consider explanatory power of the covariates in the mediation and outcome model via specifying R-squared values accordingly. \code{pwrss.z.mediation()} and \code{pwrss.z.med()} are the same functions.

Formulas are validated using Monte Carlo simulation.
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
  \item{r2m.x}{expected R-squared value for the mediator model (M ~ X). The default is \code{r2m.x = a^2 * sdx^2 / sdm^2} assuming that X is the only predictor. Thus, an \code{r2m.x} below this value will throw a warning. To consider other covariates in the mediator model provide a value greater than the default}
  \item{r2y.mx}{expected R-squared value for the outcome model (Y ~ M + X). The default is \code{r2y.mx = (b^2 * sdm^2 + cp^2 * sdx^2) / sdy^2} assuming that M and X are the only predictors. Thus, an \code{r2y.mx} below this value will throw a warning. To consider other covariates in the outcome model provide a value greater than the default}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less"}
  \item{mc}{logical; if \code{TRUE}, statistical power based on monte carlo simulation}
  \item{nsims}{number of replications, if \code{mc = TRUE}}
  \item{ndraws}{number of draws from the distribution of the path coefficients for each replication, if \code{mc = TRUE}}
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
# with standardized coefficients

## statistical power
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            alpha = 0.05, n = 200, mc = TRUE)

## minimum required sample size
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            alpha = 0.05, power = 0.80)

## adjust for covariates in the outcome model
pwrss.z.med(a = 0.25, b = 0.25, cp = 0.10,
            r2y.mx = 0.50,
            alpha = 0.05, power = 0.80)

# with binary predictor X such as treatment/control variable
# in this case standardized coefficients for path a and cp would be Cohen's d values

## statistical power
p <- 0.50 # proportion of subjects in one group
pwrss.z.med(a = 0.40, b = 0.25, cp = 0.10,
            sdx = sqrt(p*(1-p)),
            alpha = 0.05, n = 200, mc = TRUE)

## minimum required sample size
pwrss.z.med(a = 0.40, b = 0.25, cp = 0.10,
            sdx = sqrt(p*(1-p)),
            alpha = 0.05, power = 0.80)

## adjust for covariates in outcome model
pwrss.z.med(a = 0.40, b = 0.25, cp = 0.10,
            r2y.mx = 0.50, sdx = sqrt(p*(1-p)),
            alpha = 0.05, power = 0.80)
}

\references{
Aroian, L. A. (1947). The probability function of the product of two normally distributed variables. Annals of Mathematical Statistics, 18(2), 265-271.

Goodman, L. A. (1960). On the exact variance of products. Journal of the American Statistical Association, 55(292), 708-713.

MacKinnon, D. P., & Dwyer, J. H. (1993). Estimating mediated effects in prevention studies. Evaluation Review, 17(2), 144-158.

MacKinnon, D. P., Warsi, G., & Dwyer, J. H. (1995). A simulation study of mediated effect measures. Multivariate Behavioral Research, 30(1), 41-62.

Preacher, K. J., & Hayes, A. F. (2004). SPSS and SAS procedures for estimating indirect effects in simple mediation models. Behavior Research Methods, Instruments, & Computers, 36, 717-731.

Preacher, K. J., & Hayes, A. F. (2008). Asymptotic and resampling strategies for assessing and comparing indirect effects in multiple mediator models. Behavior Research Methods, 40, 879-891.

Sobel, M. E. (1982). Asymptotic intervals for indirect effects in structural equations models. In S. Leinhart (Ed.), Sociological methodology 1982 (pp. 290-312). Jossey-Bass.
}
