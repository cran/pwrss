\name{factorial.contrasts}
\alias{factorial.contrasts}

\title{Factorial Contrasts}

\description{

Helper function to construct the default contrast coefficients for various coding schemes.

Note that R has a partial matching feature which allows you to specify shortened versions of arguments, such as \code{coding} instead of \code{coding.scheme}.

Validated using \code{lm()} and \code{aov()} functions.

}

\usage{

factorial.contrasts(factor.levels = c(3, 2),
                    coding.scheme = rep("deviation", length(factor.levels)),
                    base = factor.levels,
                    intercept = FALSE,
                    verbose = TRUE)
}

\arguments{
 \item{factor.levels}{Integer. Number of levels or groups in each factor. For example, for two factors each having two levels or groups use e.g. c(2, 2), for three factors each having two levels or groups use e.g. c(2, 2, 2)}
  \item{coding.scheme}{Character vector. Coding scheme for each factor. "sum" for deviation or effect coding, "treatment" for dummy coding, "helmert" for Helmert type of coding, and "poly" for polynomial coding. Each factor can have their own coding scheme. If a single character value is provided, it will be copied to other factors}
  \item{base}{Integer vector. Specifies which group is considered the baseline group. Ignored for coding schemes other than "treatment"}
  \item{intercept}{Logical. \code{FALSE} by default. If \code{TRUE} contrast matrix includes the intercept}
  \item{verbose}{Logical. \code{TRUE} by default. If \code{FALSE} no output is printed on the console}
}

\value{
  \item{factor.levels}{Number of levels (or groups) in each factor}
  \item{factor.data}{Unique combination of factor levels}
  \item{model.matrix}{Model (design) matrix based on unique combination of factor levels}
  \item{contrast.matrix}{Contrast matrix}
}


\examples{

###################################################
############### dummy coding scheme ###############
####################################################

# one factor w/ 3 levels
contrast.object <- factorial.contrasts(factor.levels = 3,
                                       coding = "treatment")
# model matrix
contrast.object$model.matrix

# contrast matrix
contrast.object$contrast.matrix

#######################################################
###############  deviation coding scheme ##############
#######################################################

# especially useful for factorial designs
# two factors w/ 2 and 3 levels, respectively
contrast.object <- factorial.contrasts(factor.levels = c(2, 3),
                                       coding = "sum")

# model matrix
contrast.object$model.matrix

# contrast matrix
contrast.object$contrast.matrix


######################################################
###############  Helmert coding scheme ###############
######################################################

# one factor w/ 3 levels
contrast.object <- factorial.contrasts(factor.levels = 3,
                                       coding = "helmert")

# model matrix
contrast.object$model.matrix

# contrast matrix
contrast.object$contrast.matrix

#########################################################
###############  polynomial coding scheme ###############
#########################################################

# one factor w/ 3 levels
contrast.object <- factorial.contrasts(factor.levels = 3,
                                       coding = "poly")

# model matrix
contrast.object$model.matrix

# contrast matrix
contrast.object$contrast.matrix

}

