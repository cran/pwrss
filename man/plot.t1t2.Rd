\name{plot}
\alias{plot.pwrss}

\title{Type I and Type II Error Plots for t, z, and F Tests}

\description{
Plots Type I (alpha) and Type II (beta) errors for t, z, and F tests.
}

\usage{
\method{plot}{pwrss}(x, ...)
}

\arguments{
  \item{x}{an object of the type "pwrss" returned from one of the \code{pwrss} functions}
  \item{...}{for S3 generic/method consistency}
}

\value{no return value at the moment}

\examples{
design <- pwrss.f.ancova(n.levels = c(3,3),
               n = 50, eta2 = 0.10)
plot(design)
}
