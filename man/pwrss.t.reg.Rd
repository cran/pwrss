\name{pwrss.t.reg}
\alias{pwrss.t.reg}
\alias{pwrss.z.reg}

\title{Test of a Single Coefficient in Multiple Linear Regression (t Test or z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a single coefficient in multiple linear regression. The predictor is assumed to be continuous by default. However, one can find statistical power or minimum required sample size for a binary predictor (such as treatment and control groups in experimental designs) by specifying \code{sdx = sqrt(p*(1-p))} where\code{p} is the proportion of subjects in one of the groups. The sample size in each group would be \code{n*p} and \code{n*(1-p)}.
}

\usage{
pwrss.t.reg(beta = 0.25, beta0 = 0,
            sdx = 1, sdy = 1,
            k = 1, r2 = (beta * sdx / sdy)^2,
            alpha = 0.05, n = NULL, power = NULL,
            alternative = c("not equal", "less", "greater"),
            verbose = TRUE)

pwrss.z.reg(beta = 0.25, beta0 = 0,
            sdx = 1, sdy = 1,
            k = 1, r2 = (beta * sdx / sdy)^2,
            alpha = 0.05, n = NULL, power = NULL,
            alternative = c("not equal", "less", "greater"),
            verbose = TRUE)
}

\arguments{
  \item{beta}{expected regression coefficient. One can use standardized regression coefficient, but should use \code{sdx = 1} and \code{sdy = 1} or leave them out as they are default specifications}
  \item{beta0}{regression coefficient under null hypothesis}
  \item{sdx}{standard deviation of the predictor. For a binary predictor, \code{sdx = sqrt(p*(1-p))} where\code{p} is the proportion of subjects in one of the groups}
  \item{sdy}{standard deviation of the outcome}
  \item{k}{(total) number of predictors}
  \item{r2}{expected model R-squared. The default is \code{r2 = (beta * sdx / sdy)^2} assuming a linear regression with one predictor. Thus, an \code{r2} below this value will throw a warning. To consider other covariates in the model provide a value greater than the default along with argument \code{k>1}.}
  \item{n}{sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "equivalent"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{df}{numerator degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{sample size}
}

\examples{
# continuous predictor x
pwrss.t.reg(beta = 0.20,
            k = 5, r2 = 0.30,
            power = 0.80)
pwrss.t.reg(beta = 0.20,
            k = 5, r2 = 0.30,
            n = 140)

# binary predictor x
p <- 0.50 # proportion of subjects in one group
pwrss.t.reg(beta = 0.20, sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            power = 0.80)
pwrss.t.reg(beta = 0.20, sdx = sqrt(p*(1-p)) ,
            k = 5, r2 = 0.30,
            n = 552)
}

\references{
Add references.
}
