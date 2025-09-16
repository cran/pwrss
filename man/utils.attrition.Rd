\name{inflate.sample}
\alias{inflate.sample}

\title{Inflate Sample Size for Attrition}

\description{

  Helper function to inflate sample size for attrition.

}

\usage{
  inflate.sample(n, rate = 0.05,
               ceiling = TRUE,
               verbose = TRUE)
}

\arguments{
  \item{n}{sample size.}
  \item{rate}{attrition rate.}
  \item{ceiling}{rounds-up the inflated sample size.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
}

\value{inflated sample size.}

\examples{

inflate.sample(n = 100, rate = 0.05)

}

