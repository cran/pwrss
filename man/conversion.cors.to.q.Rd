\name{cors.to.q}
\alias{cors.to.q}
\alias{q.to.cors}
\alias{cor.to.z}

\title{Conversion from Correlation Difference to Cohen's q}

\description{

  Helper function to convert correlation difference to Cohen's q (and vice versa). \code{cor.to.z()} function applies Fisher's z transformation.

}

\usage{
  cor.to.z(rho, verbose = TRUE)

  cors.to.q(rho1, rho2, verbose = TRUE)

  q.to.cors(q, rho1 = NULL, rho2 = NULL, verbose = TRUE)

}

\arguments{
  \item{rho}{correlation.}
  \item{rho1}{first correlation.}
  \item{rho2}{second correlation.}
  \item{q}{Cohen's q effect size.}
  \item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
}

\value{
  \item{rho1}{first correlation.}
  \item{rho2}{second correlation.}
  \item{delta}{correlation difference: rho1 - rho2.}
  \item{q}{Cohen's q effect size.}
}


\examples{

q.to.cors(q = 0.10, rho1 = 0.50)
q.to.cors(q = 0.30, rho1 = 0.50)
q.to.cors(q = 0.50, rho1 = 0.50)

cors.to.q(rho2 = 0.5712027, rho1 = 0.50)
cors.to.q(rho2 = 0.6907068, rho1 = 0.50)
cors.to.q(rho2 = 0.7815365, rho1 = 0.50)

}

