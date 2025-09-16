\name{d.to.cles}
\alias{d.to.cles}
\alias{cles.to.d}

\title{Conversion from Cohen's d to Common Language Effect Size}

\description{

Helper function to convert Cohen's d to common language effect size (or vice versa). The result is the probability of superiority for independent samples. It can be interpreted as the probability that a randomly selected observation from Group 1 exceeds a randomly selected observation from Group 2. The rationale is the same for paired-samples and one-sample designs, but the interpretation differs: For paired samples, it can be interpreted as the probability that the difference score (i.e., the score under Condition 1 minus the score under Condition 2) is greater than zero for a randomly selected individual. For a one-sample design, it can be interpreted as the probability that a randomly selected observation is greater than the reference value (e.g., 0).

}

\usage{

  d.to.cles(d, design = c("independent", "paired", "one.sample"), verbose = TRUE)

  cles.to.d(cles, design = c("independent", "paired", "one.sample"), verbose = TRUE)

}

\arguments{
  \item{d}{Cohen's d}
  \item{design}{character; one of the "independent", "paired", or "one.sample". The default is "independent".}
  \item{cles}{common language effect size.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
}

\value{
  \item{d}{Cohen's d}
  \item{cles}{common language effect size.}
}


\examples{

d.to.cles(0.20) # small
d.to.cles(0.50) # medium
d.to.cles(0.80) # large

cles.to.d(0.5562315)
cles.to.d(0.6381632)
cles.to.d(0.7141962)

}

