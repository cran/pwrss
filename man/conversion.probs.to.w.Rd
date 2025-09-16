\name{probs.to.w}
\alias{probs.to.w}


\title{Conversion from Probabilities to Cohen's w}

\description{
  Helper function to convert (multinomial or product-multinomial) probabilities to Cohen's w.
}

\usage{
 probs.to.w(prob.matrix,
           null.prob.matrix = NULL,
           verbose = TRUE)
}

\arguments{
\item{prob.matrix}{a vector or matrix of cell probabilities under alternative hypothesis}
\item{null.prob.matrix}{a vector or matrix of cell probabilities under null hypothesis. Calculated automatically when \code{prob.matrix} is specified. The default can be overwritten by the user via providing a vector of the same size or matrix of the same dimensions as \code{prob.matrix}}
\item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
}

\value{
    \item{w}{Cohen's w effect size. It can be any of Cohen's W, Phi coefficient, Cramer's V. Phi coefficient is defined as \code{sqrt(X2/n)} and Cramer's V is defined as \code{sqrt(X2/(n*v))} where \code{v} is \code{min(nrow - 1, ncol - 1)} and X2 is the chi-square statistic.}
    \item{df}{degrees of freedom.}
}

\examples{
  # ---------------------------------------------------------#
  # Example 1: Cohen's W                                     #
  # goodness-of-fit test for 1 x k or k x 1 table            #
  # How many subjects are needed to claim that               #
  # girls choose STEM related majors less than males?       #
  # ---------------------------------------------------------#

  ## from https://www.aauw.org/resources/research/the-stem-gap/
  ## 28 percent of the  workforce in STEM field is women
  prob.vector <- c(0.28, 0.72)
  null.prob.vector <- c(0.50, 0.50)
  probs.to.w(prob.vector, null.prob.vector)

  power.chisq.gof(w = 0.44, df = 1,
                  alpha = 0.05, power = 0.80)


  # ---------------------------------------------------------#
  # Example 2: Phi Coefficient (or Cramer's V or Cohen's W)  #
  # test of independence for 2 x 2 contingency tables        #
  # How many subjects are needed to claim that               #
  # girls are underdiagnosed with ADHD?                      #
  # ---------------------------------------------------------#

  ## from https://time.com/growing-up-with-adhd/
  ## 5.6 percent of girls and 13.2 percent of boys are diagnosed with ADHD
  prob.matrix <- rbind(c(0.056, 0.132),
                       c(0.944, 0.868))
  colnames(prob.matrix) <- c("Girl", "Boy")
  rownames(prob.matrix) <- c("ADHD", "No ADHD")
  prob.matrix

  probs.to.w(prob.matrix)

  power.chisq.gof(w = 0.1302134, df = 1,
                  alpha = 0.05, power = 0.80)


  # --------------------------------------------------------#
  # Example 3: Cramer's V (or Cohen's W)                    #
  # test of independence for j x k contingency tables       #
  # How many subjects are needed to detect the relationship #
  # between depression severity and gender?                 #
  # --------------------------------------------------------#

  ## from https://doi.org/10.1016/j.jad.2019.11.121
  prob.matrix <- cbind(c(0.6759, 0.1559, 0.1281, 0.0323, 0.0078),
                       c(0.6771, 0.1519, 0.1368, 0.0241, 0.0101))
  rownames(prob.matrix) <- c("Normal", "Mild", "Moderate",
                             "Severe", "Extremely Severe")
  colnames(prob.matrix) <- c("Female", "Male")
  prob.matrix

  probs.to.w(prob.matrix)

  power.chisq.gof(w = 0.03022008, df = 4,
                  alpha = 0.05, power = 0.80)

}

\references{
  Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
