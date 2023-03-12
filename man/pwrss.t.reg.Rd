\name{Linear Regression (z or t Test)}
\alias{pwrss.t.reg}
\alias{pwrss.z.reg}
\alias{pwrss.t.regression}
\alias{pwrss.z.regression}

\title{Test of a Single Coefficient in Linear Regression (t Test or z Test)}

\description{
Calculates statistical power or minimum required sample size (only one can be NULL at a time) to test a single coefficient in multiple linear regression. The predictor is assumed to be continuous by default. However, one can find statistical power or minimum required sample size for a binary predictor (such as treatment and control groups in an experimental design) by specifying \code{sdx = sqrt(p*(1-p))} where\code{p} is the proportion of subjects in one of the groups. The sample size in each group would be \code{n*p} and \code{n*(1-p)}. \code{pwrss.t.regression()} and \code{pwrss.t.reg()} are the same functions. \code{pwrss.z.regression()} and \code{pwrss.z.reg()} are the same functions.

Formulas are validated using Monte Carlo simulation, G*Power, tables in PASS documentation, and tables in Bulus (2021).
}

\usage{
pwrss.t.reg(beta1 = 0.25, beta0 = 0, margin = 0,
            sdx = 1, sdy = 1,
            k = 1, r2 = (beta1 * sdx / sdy)^2,
            alpha = 0.05, n = NULL, power = NULL,
            alternative = c("not equal", "less", "greater",
                            "non-inferior", "superior", "equivalent"),
            verbose = TRUE)

pwrss.z.reg(beta1 = 0.25, beta0 = 0, margin = 0,
            sdx = 1, sdy = 1,
            k = 1, r2 = (beta1 * sdx / sdy)^2,
            alpha = 0.05, n = NULL, power = NULL,
            alternative = c("not equal", "less", "greater",
                            "non-inferior", "superior", "equivalent"),
            verbose = TRUE)
}

\arguments{
  \item{beta1}{expected regression coefficient. One can use standardized regression coefficient, but should use \code{sdx = 1} and \code{sdy = 1} or leave them out as they are default specifications}
  \item{beta0}{regression coefficient under null hypothesis (usually zero). Not to be confused with the intercept. One can use standardized regression coefficient, but should use \code{sdx = 1} and \code{sdy = 1} or leave them out as they are default specifications}
  \item{margin}{non-inferiority, superiority, or equivalence margin (margin: boundry of \code{beta1 - beta0} that is practically insignificant). For non-inferiority and superiority tests 'margin' shifts the value of 'beta0' as 'beta0 + margin'. Thus for a positive 'beta0' a negative 'margin' will decrease the value of 'beta0' while a positive 'margin' will increase its value. For non-inferiority and superiority tests 'beta1' is expected be greater than 'beta0 + margin' (when higher values of the outcome is better). For equivalence tests the value of 'beta0' shifts both to the left and right as 'beta0 - margin' and 'beta0 + margin'. For equivalence test 'margin' is always positive and 'beta1' is expected to fall within the interval}
  \item{sdx}{standard deviation of the predictor. For a binary predictor, \code{sdx = sqrt(p*(1-p))} where\code{p} is the proportion of subjects in one of the groups}
  \item{sdy}{standard deviation of the outcome}
  \item{k}{(total) number of predictors}
  \item{r2}{expected model R-squared. The default is \code{r2 = (beta * sdx / sdy)^2} assuming a linear regression with one predictor. Thus, an \code{r2} below this value will throw a warning. To consider other covariates in the model provide a value greater than the default \code{r2} along with the argument \code{k>1}.}
  \item{n}{total sample size}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{alpha}{probability of type I error}
  \item{alternative}{direction or type of the hypothesis test: "not equal", "greater", "less", "non-inferior", "superior", or "equivalent"}
  \item{verbose}{if \code{FALSE} no output is printed on the console. Useful for simulation, plotting, and whatnot}
}

\value{
  \item{parms}{list of parameters used in calculation}
  \item{test}{type of the statistical test (z, t or F?)}
  \item{df}{numerator degrees of freedom}
  \item{ncp}{non-centrality parameter}
  \item{power}{statistical power \eqn{(1-\beta)}}
  \item{n}{total sample size}
}

\examples{
# continuous predictor x (and 4 covariates)
pwrss.t.reg(beta1 = 0.20, alpha = 0.05,
            alternative = "not equal",
            k = 5, r2 = 0.30,
            power = 0.80)

pwrss.t.reg(beta1 = 0.20, alpha = 0.05,
            alternative = "not equal",
            k = 5, r2 = 0.30,
            n = 138)

# binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
pwrss.t.reg(beta1 = 0.20, alpha = 0.05,
            alternative = "not equal",
            sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            power = 0.80)

pwrss.t.reg(beta1 = 0.20, alpha = 0.05,
            alternative = "not equal",
            sdx = sqrt(p*(1-p)) ,
            k = 5, r2 = 0.30,
            n = 550)

# non-inferiority test with binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
pwrss.t.reg(beta1 = 0.20, beta0 = 0.10, margin = -0.05,
            alpha = 0.05, alternative = "non-inferior",
            sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            power = 0.80)

pwrss.t.reg(beta1 = 0.20, beta0 = 0.10, margin = -0.05,
            alpha = 0.05, alternative = "non-inferior",
            sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            n = 770)

# superiority test with binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
pwrss.t.reg(beta1 = 0.20, beta0 = 0.10, margin = 0.01,
            alpha = 0.05, alternative = "superior",
            sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            power = 0.80)

pwrss.t.reg(beta1 = 0.20, beta0 = 0.10, margin = 0.01,
            alpha = 0.05, alternative = "superior",
            sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            n = 2138)

# equivalence test with binary predictor x (and 4 covariates)
p <- 0.50 # proportion of subjects in one group
pwrss.t.reg(beta1 = 0.20, beta0 = 0.20, margin = 0.05,
            alpha = 0.05, alternative = "equivalent",
            sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            power = 0.80)

pwrss.t.reg(beta1 = 0.20, beta0 = 0.20, margin = 0.05,
            alpha = 0.05, alternative = "equivalent",
            sdx = sqrt(p*(1-p)),
            k = 5, r2 = 0.30,
            n = 9592)
}

\references{
Bulus, M. (2021). Sample size determination and optimal design of randomized/non-equivalent pretest-posttest control-group designs. Adiyaman Univesity Journal of Educational Sciences, 11(1), 48-69.

Phillips, K. F. (1990). Power of the two one-Sided tests procedure in bioequivalence. Journal of Pharmacokinetics and Biopharmaceutics, 18(2), 137-144.

Dupont, W. D., and Plummer, W. D. (1998). Power and sample size calculations for studies involving linear regression. Controlled Clinical Trials, 19(6), 589-601.
}
