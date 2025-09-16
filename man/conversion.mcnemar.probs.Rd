\name{joint.probs.2x2}
\alias{joint.probs.2x2}
\alias{marginal.probs.2x2}

\title{Conversion Between Joint and Marginal Probabilities for the McNemar Test}

\description{

  Helper function to converts joint probabilities to marginal probabilities (and vice versa) for the McNemar test applied to paired binary data.

}

\usage{

joint.probs.2x2(prob1, prob2, rho = 0.50, verbose = TRUE)

marginal.probs.2x2(prob11, prob10, prob01, prob00, verbose = TRUE)

}

\arguments{

  \item{prob1}{(marginal) probability of success in case group (or after).}
  \item{prob2}{(marginal) probability of success in matched-control group (or before).}
  \item{rho}{the correlation between case and matched-control, or after and before (phi coefficient).}
  \item{prob11}{(joint) probability of success in both groups. 'prob11' and 'prob00' are known as concordant probs.}
  \item{prob10}{(joint) probability of success in case (or after) but failure in matched control (or before). 'prob10' and 'prob01' are known as discordant probs.}
  \item{prob01}{(joint) probability of failure in case (or after) but success in matched control (or before). prob10' and 'prob01' are known as discordant probs.}
  \item{prob00}{(joint) probability of failure in both groups. 'prob11' and 'prob00' are known as concordant probs.}
  \item{verbose}{if \code{FALSE} no output is printed on the console.}

}

\value{

  \item{parms}{list of parameters used in calculation.}
  \item{prob1}{(marginal) probability of success in case group (or after).}
  \item{prob2}{(marginal) probability of success in matched-control group (or before).}
  \item{rho}{the correlation between case and matched-control, or after and before (phi coefficient).}
  \item{prob11}{(joint) probability of success in both groups. 'prob11' and 'prob00' are known as concordant probs.}
  \item{prob10}{(joint) probability of success in case (or after) but failure in matched control (or before). 'prob10' and 'prob01' are known as discordant probs.}
  \item{prob01}{(joint) probability of failure in case (or after) but success in matched control (or before). prob10' and 'prob01' are known as discordant probs.}
  \item{prob00}{(joint) probability of failure in both groups. 'prob11' and 'prob00' are known as concordant probs.}

}


\examples{

# example data for a matched case-control design
# subject  case     control
# <int>    <dbl>    <dbl>
#   1        1        1
#   2        0        1
#   3        1        0
#   4        0        1
#   5        1        1
#   ...     ...      ...
#   100      0        0

# example summary stats
# prob1 = mean(case) which is 0.55
# prob2 = mean(control) which is 0.45
# rho = cor(case, control) which is 0.4141414


# example data for a before-after design
# subject  before   after
# <int>    <dbl>    <dbl>
#   1        1        1
#   2        0        1
#   3        1        0
#   4        0        1
#   5        1        1
#   ...     ...      ...
#   100      0        0

# example summary stats
# prob1 = mean(after) which is 0.55
# prob2 = mean(before) which is 0.45
# rho = cor(after, before) which is 0.4141414

# convert to a 2 x 2 frequency table
freqs <- matrix(c(30, 10, 20, 40), nrow = 2, ncol = 2)
colnames(freqs) <- c("control_1", "control_0")
rownames(freqs) <- c("case_1", "case_0")
freqs

# convert to a 2 x 2 proportion table
props <- freqs / sum(freqs)
props

# discordant pairs (0 and 1, or 1 and 0) in 'props' matrix
# are the sample estimates of prob01 and prob10


# we may not have 2 x 2 joint probs
# convert marginal probs to joint probs using summary stats
jp <- joint.probs.2x2(prob1 = 0.55, # mean of case (or after)
                          prob2 = 0.45, # mean of matched control (or before)
                          # correlation b/w matched case-control / before-after
                          rho = 0.4141414)

# required sample size for exact test
# assuming prob01 and prob10 are population parameters
power.exact.mcnemar(prob01 = jp$prob01,
                    prob10 = jp$prob10,
                    power = 0.80, alpha = 0.05,
                    method = "exact")

# convert joint probs to marginal probs and calc phi coefficient (rho)
# these values can be used in other procedures
marginal.probs.2x2(prob11 = 0.35, # mean of case (or after)
                    prob10 = 0.20, # mean of matched control (or before)
                    prob01 = 0.10,
                    prob00 = 0.35)
}

\references{

Zhang, S., Cao, J., and Ahn, C. (2017). Inference and sample size calculation for clinical trials with incomplete observations of paired binary outcomes. Statistics in Medicine, 36(4), 581-591. \doi{10.1002/sim.7168}

}
