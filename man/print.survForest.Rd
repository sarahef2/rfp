% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.survForest.r
\name{print.survForest}
\alias{print.survForest}
\title{print.survForest}
\usage{
\method{print}{survForest}(x, ...)
}
\arguments{
\item{x}{A fitted survForest object}

\item{...}{...}
}
\description{
Print a survForest object
}
\examples{
x = matrix(rnorm(1000), ncol = 10)
y = exp(rowMeans(x))
c = rbinom(100, 1, 0.9)
fit = survForest(x, y, c, ntrees = 10)
fit
}
