\name{f.to.etasq}
\alias{f.to.etasq}
\alias{etasq.to.f}

\title{Conversion from Cohen's f to Eta-squared}

\description{
Helper function to convert between Cohen's f and Eta-squared.
}

\usage{
f.to.etasq(f, verbose = TRUE)

etasq.to.f(eta.squared, verbose = TRUE)
}

\arguments{
\item{f}{Cohen's f.}
\item{eta.squared}{(Partial) Eta-squared.}
\item{verbose}{logical; whether the output should be printed on the console. \code{TRUE} by default.}
}

\value{
\item{f}{Cohen's f.}
\item{eta.squared}{(Partial) Eta-squared.}
}

\examples{
f.to.etasq(f = 0.10) # small
f.to.etasq(f = 0.25) # medium
f.to.etasq(f = 0.40) # large
}

\references{
  Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Lawrence Erlbaum Associates.
}
