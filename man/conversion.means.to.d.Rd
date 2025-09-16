\name{means.to.d}
\alias{means.to.d}

\title{Conversion from Means and Standard Deviations to Cohen's d}

\description{

  Helper function to convert means and standard deviations to Cohen's d.

}

\usage{
means.to.d(mu1, mu2 = 0,
           sd1 = 1, sd2 = 1,
           n2, n.ratio = 1,
           paired = FALSE,
           rho.paired = 0.50,
           verbose = TRUE)
}

\arguments{

  \item{mu1}{mean of the first group.}
  \item{mu2}{mean of the second group.}
  \item{sd1}{standard deviation of the first group.}
  \item{sd2}{standard deviation of the second group.}
  \item{n.ratio}{\code{n1/n2} ratio (applies to independent samples only).}
  \item{paired}{if \code{TRUE} paired samples}
  \item{rho.paired}{correlation between repeated measures for paired samples (e.g., pretest and post-test).}
  \item{n2}{integer; sample size in the second group (or for the single group in paired samples).}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
}

\value{

  \item{d}{Cohen's d}

}


\examples{

# means and standard deviations from independent samples
means.to.d(mu1 = 20, mu2 = 17.5,
           sd1 = 5, sd2 = 15,
           n2 = 30, n.ratio = 1)

# means and standard deviations from paired samples
means.to.d(mu1 = 20, mu2 = 17.5,
           sd1 = 5, sd2 = 15,
           n2 = 30, n.ratio = 1,
           paired = TRUE,
           rho.paired = 0.50)
}

