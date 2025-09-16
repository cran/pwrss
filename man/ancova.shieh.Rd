\name{power.f.ancova.shieh}
\alias{power.f.ancova.shieh}

\title{Power Analysis for One-, Two-, Three-Way ANCOVA Using Means, Standard Deviations, and (Optionally) Contrasts (F test)}

\description{
  Calculates power or sample size for one-, two-, three-way ANCOVA. For factorial designs, use the argument \code{factor.levels} but note that unique combination of levels (cells in this case) should follow a specific order for the test of interaction. The order of marginal means and standard deviations  is printed as a warning message.

  Note that R has a partial matching feature which allows you to specify shortened versions of arguments, such as \code{mu} or \code{mu.vec} instead of \code{mu.vector}, or such as \code{k} or \code{k.cov} instead of \code{k.covariates}.

  Formulas are validated using examples and tables in Shieh (2020).
}

\usage{
power.f.ancova.shieh(mu.vector, sd.vector,
                     n.vector = NULL, p.vector = NULL,
                     factor.levels = length(mu.vector),
                     r.squared = 0, k.covariates = 1,
                     contrast.matrix = NULL,
                     power = NULL, alpha = 0.05,
                     ceiling = TRUE, verbose = TRUE,
                     pretty = FALSE)
}

\arguments{
  \item{mu.vector}{vector; adjusted means (or estimated marginal means) for each level of a factor.}
  \item{sd.vector}{vector; unadjusted standard deviations for each level of a factor. If a pooled standard deviation is provided, repeat its value to match the number of group means. A warning will be issued if group standard deviations differ substantially beyond what is expected due to sampling error.}
  \item{n.vector}{vector; sample sizes for each level of a factor.}
  \item{p.vector}{vector; proportion of total sample size in each level of a factor. These proportions should sum to one.}
  \item{factor.levels}{integer; number of levels or groups in each factor. For example, for two factors each having two levels or groups use e.g. c(2, 2), for three factors each having two levels or groups use e.g. c(2, 2, 2).}
  \item{r.squared}{explanatory power of covariates (R-squared) in the ANCOVA model.}
  \item{k.covariates}{integer; number of covariates in the ANCOVA model.}
  \item{contrast.matrix}{vector or matrix; contrasts should not be confused with the model (design) matrix. Rows of contrast matrix indicate independent vector of contrasts summing to zero. The default contrast matrix is constructed using deviation coding scheme (a.k.a. effect coding). Columns in the contrast matrix indicate number of levels or groups (or cells in factorial designs).}
   \item{power}{statistical power, defined as the probability of correctly rejecting a false null hypothesis, denoted as \eqn{1 - \beta}.}
  \item{alpha}{type 1 error rate, defined as the probability of incorrectly rejecting a true null hypothesis, denoted as \eqn{\alpha}.}
  \item{ceiling}{logical;  \code{TRUE} by default. If \code{FALSE} sample size in each cell is NOT rounded up.}
  \item{verbose}{logical; \code{TRUE} by default. If \code{FALSE} no output is printed on the console.}
  \item{pretty}{logical; whether the output should show Unicode characters (if encoding allows for it). \code{FALSE} by default.}
}

\value{
  \item{parms}{list of parameters used in calculation.}
  \item{test}{type of the statistical test (F-Test)}
  \item{eta.squared}{(partial) eta-squared.}
  \item{f}{Cohen's f statistic.}
  \item{df1}{numerator degrees of freedom.}
  \item{df2}{denominator degrees of freedom.}
  \item{ncp}{non-centrality parameter for the alternative.}
  \item{null.ncp}{non-centrality parameter for the null.}
  \item{power}{statistical power \eqn{(1-\beta)}.}
  \item{n.total}{total sample size.}
}

\examples{
###################################################################
##########################  main effect  ##########################
###################################################################

# power for one-way ANCOVA (two levels or groups)
power.f.ancova.shieh(mu.vector = c(0.20, 0), # marginal means
                     sd.vector = c(1, 1), # unadjusted standard deviations
                     n.vector = c(150, 150), # sample sizes
                     r.squared = 0.50, # proportion of variance explained by covariates
                     k.covariates = 1, # number of covariates
                     alpha = 0.05)


# sample size for one-way ANCOVA (two levels or groups)
power.f.ancova.shieh(mu.vector = c(0.20, 0), # marginal means
                     sd.vector = c(1, 1), # unadjusted standard deviations
                     p.vector = c(0.50, 0.50), # allocation, should sum to 1
                     r.squared = 0.50,
                     k.covariates = 1,
                     alpha = 0.05,
                     power = 0.80)

###################################################################
#######################  interaction effect  ######################
###################################################################

# sample size for two-way ANCOVA (2 x 2)
power.f.ancova.shieh(mu.vector = c(0.20, 0.25, 0.15, 0.05), # marginal means
                     sd.vector = c(1, 1, 1, 1), # unadjusted standard deviations
                     p.vector = c(0.25, 0.25, 0.25, 0.25), # allocation, should sum to 1
                     factor.levels = c(2, 2), # 2 by 2 factorial design
                     r.squared = 0.50,
                     k.covariates = 1,
                     alpha = 0.05,
                     power = 0.80)
# Elements of `mu.vector`, `sd.vector`, `n.vector` or `p.vector` should follow this specific order:
#  A1:B1  A1:B2  A2:B1  A2:B2

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
ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                                      sd.vector = c(1, 1, 1), # unadjusted standard deviations
                                      p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
                                      contrast.matrix = contrast.matrix,
                                      r.squared = 0.50,
                                      k.covariates = 1,
                                      alpha = 0.05,
                                      power = 0.80)

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
ancova.design <- power.f.ancova.shieh(mu.vector = c(0.15, 0.30, 0.20), # marginal means
                                      sd.vector = c(1, 1, 1), # unadjusted standard deviations
                                      p.vector = c(1/3, 1/3, 1/3), # allocation, should sum to 1
                                      contrast.matrix = contrast.matrix,
                                      r.squared = 0.50,
                                      k.covariates = 1,
                                      alpha = 0.05,
                                      power = 0.80)

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
                                      alpha = 0.05,
                                      power = 0.80)

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
                                      alpha = 0.05,
                                      power = 0.80)

# power of the planned contrasts
power.t.contrasts(ancova.design)
}

\references{
Shieh, G. (2020). Power analysis and sample size planning in ANCOVA designs. Psychometrika, 85(1), 101-120. \doi{10.1007/s11336-019-09692-3}
}
