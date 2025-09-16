\name{pwrss.t.contrasts}
\alias{pwrss.t.contrasts}
\alias{pwrss.t.contrast}
\alias{power.t.contrasts}
\alias{power.t.contrast}

\title{Power Analysis for One-, Two-, Three-Way ANCOVA Contrasts and Multiple Comparisons (T-Tests)}

\description{
 Calculates power or sample size for one-, two-, three-Way ANCOVA contrasts and multiple comparisons. The \code{pwrss.t.contrast()} function allows a single contrast. For multiple contrast testing (multiple comparisons) use the \code{pwrss.t.contrasts()} function. The \code{pwrss.t.contrasts()} function also allows adjustment to alpha due to multiple testing. All arguments are the same as with the earlier function except that it accepts an object returned from the \code{pwrss.f.ancova.shieh()} function for convenience. Beware that, in this case, all other arguments are ignored except \code{alpha} and \code{adjust.alpha}.

 The \code{pwrss.t.contrasts()} function returns a data frame with as many rows as number of contrasts and eight columns with names 'contrast', 'comparison', 'psi', 'd', 'ncp', 'df', 'n.total', and 'power'. 'psi' is the contrast estimate. 'd'
 is the standardized contrast estimate. Remaining parameters are explained elsewhere.

 Note that R has a partial matching feature which allows you to specify shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of \code{k.covariates}.

  Formulas are validated using examples and tables in Shieh (2017).
}

\usage{

power.t.contrast(mu.vector,
                 sd.vector,
                 contrast.vector,
                 n.vector = NULL, p.vector = NULL,
                 r.squared = 0, k.covariates = 1,
                 power = NULL, alpha = 0.05,
                 tukey.kramer = FALSE,
                 ceiling = TRUE, verbose = TRUE, pretty = FALSE)

power.t.contrasts(x = NULL,
                  mu.vector = NULL,
                  sd.vector = NULL,
                  n.vector = NULL, p.vector = NULL,
                  r.squared = 0, k.covariates = 1,
                  contrast.matrix = NULL,
                  power = NULL, alpha = 0.05,
                  adjust.alpha = c("none", "tukey", "bonferroni",
                                   "holm", "hochberg", "hommel",
                                   "BH", "BY", "fdr"),
                  ceiling = TRUE, verbose = TRUE, pretty = FALSE)

}

\arguments{
  \item{x}{object; an object returned from \code{pwrss.f.ancova.shieh()} function.}
  \item{mu.vector}{vector; adjusted means (or estimated marginal means) for each level of a factor. Ignored when 'x' is specified.}
  \item{sd.vector}{vector; unadjusted standard deviations for each level of a factor. Ignored when 'x' is specified.}
  \item{n.vector}{vector; sample sizes for each level of a factor. Ignored when 'x' is specified.}
  \item{p.vector}{vector; proportion of total sample size in each level of a factor. These proportions should sum to one. Ignored when 'x' is specified.}
  \item{r.squared}{explanatory power of covariates (R-squared) in the ANCOVA model. Ignored when 'x' is specified.}
  \item{k.covariates}{Number of covariates in the ANCOVA model. Ignored when 'x' is specified.}
  \item{contrast.vector}{vector; a single contrast in the form of a vector with as many elements as number of levels or groups (or cells in factorial designs). Ignored when 'x' is specified.}
 \item{contrast.matrix}{vector or matrix; contrasts should not be confused with the model (design) matrix. Rows of contrast matrix indicate independent vector of contrasts summing to zero. The default contrast matrix is constructed using deviation coding scheme (a.k.a. effect coding). Columns in the contrast matrix indicate number of levels or groups (or cells in factorial designs). Ignored when 'x' is specified.}
  \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}. Ignored when 'x' is specified.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}. Note that this should be specified even if 'x' is specified. The 'alpha' value within the 'x' object pertains to the omnibus test, NOT the test of contrasts.}
  \item{tukey.kramer}{logical; \code{TRUE} by default. If \code{FALSE} no adjustment will be made to control Type 1 error.}
  \item{adjust.alpha}{character; one of the methods in c("none", "tukey", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr") to control Type 1 error. See \code{?stats::p.adjust}.}
  \item{ceiling}{logical;  \code{TRUE} by default. If \code{FALSE} sample size in each cell is NOT rounded up.}
  \item{verbose}{logical; \code{TRUE} by default. If \code{FALSE} no output is printed on the console.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (T-Test).}
  \item{psi}{contrast-weighted mean difference.}
  \item{d}{contrast-weighted standardized mean difference.}
  \item{df}{degrees of freedom.}
  \item{t.alpha}{critical values.}
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n.total}{total sample size.}
}


\examples{

  ###################################################################
  #######################  planned contrasts  #######################
  ###################################################################

  #########################
  ## dummy coding scheme ##
  #########################

  contrast.object <- factorial.contrasts(factor.levels = 3, # one factor w/ 3 levels
                                         coding = "treatment") # use dummy coding scheme

  # get contrast matrix from the contrast object
  contrast.matrix <- contrast.object$contrast.matrix

  # calculate sample size given design characteristics
  ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30,  0.20), # marginal means
                                        sd.vector = c(1, 1, 1), # unadjusted standard deviations
                                        p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
                                        contrast.matrix = contrast.matrix,
                                        r.squared = 0.50,
                                        k.covariates = 1,
                                        alpha = 0.05, power = 0.80)

  # power of planned contrasts, adjusted for alpha level
  power.t.contrasts(ancova.design, adjust.alpha = "fdr")

  ###########################
  ## Helmert coding scheme ##
  ###########################

  contrast.object <- factorial.contrasts(factor.levels = 3, # one factor w/ 4 levels
                                         coding = "helmert") # use helmert coding scheme

  # get contrast matrix from the contrast object
  contrast.matrix <- contrast.object$contrast.matrix

  # calculate sample size given design characteristics
  ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30,  0.20), # marginal means
                                        sd.vector = c(1, 1, 1), # unadjusted standard deviations
                                        p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
                                        contrast.matrix = contrast.matrix,
                                        r.squared = 0.50,
                                        k.covariates = 1,
                                        alpha = 0.05, power = 0.80)

  # power of planned contrasts
  power.t.contrasts(ancova.design)

  ##############################
  ## polynomial coding scheme ##
  ##############################

  contrast.object <- factorial.contrasts(factor.levels = 3, # one factor w/ 4 levels
                                         coding = "poly") # use polynomial coding scheme

  # get contrast matrix from the contrast object
  contrast.matrix <- contrast.object$contrast.matrix

  # calculate sample size given design characteristics
  ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                                        sd.vector = c(1, 1, 1), # unadjusted standard deviations
                                        p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
                                        contrast.matrix = contrast.matrix,
                                        r.squared = 0.50,
                                        k.covariates = 1,
                                        alpha = 0.05, power = 0.80)

  # power of the planned contrasts
  power.t.contrasts(ancova.design)

  ######################
  ## custom contrasts ##
  ######################

  # custom contrasts
  contrast.matrix <- rbind(
    cbind(A1 = 1, A2 = -0.50, A3 = -0.50),
    cbind(A1 = 0.50, A2 = 0.50, A3 = -1)
  )
  # labels are not required for custom contrasts,
  # but they make it easier to understand power.t.contrasts() output

  # calculate sample size given design characteristics
  ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                                        sd.vector = c(1, 1, 1), # unadjusted standard deviations
                                        p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
                                        contrast.matrix = contrast.matrix,
                                        r.squared = 0.50,
                                        k.covariates = 1,
                                        alpha = 0.05, power = 0.80)

  # power of the planned contrasts
  power.t.contrasts(ancova.design)

}

\references{
  Shieh, G. (2017). Power and sample size calculations for contrast analysis in ANCOVA. Multivariate Behavioral Research, 52(1), 1-11. \doi{10.1080/00273171.2016.1219841}
}
