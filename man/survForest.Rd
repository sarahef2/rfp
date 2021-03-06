% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survForest.r
\name{survForest}
\alias{survForest}
\title{survForest}
\usage{
survForest(x, y, censor, ntrees = 500, mtry = max(1, as.integer(ncol(x)/3)),
  split.gen = "random", split.rule = "logrank", nsplit = 1,
  nmin = max(1, as.integer(log(nrow(x)))), alpha = 0, replacement = FALSE,
  resample.prob = 0.632, subject.weight = NULL, variable.weight = NULL,
  importance = TRUE, nimpute = 1, use.cores = 0, verbose = FALSE)
}
\arguments{
\item{x}{A matrix or data.frame for features}

\item{y}{Response variable, a numeric/factor vector or a Surv object}

\item{censor}{The censoring indicator if survival model is used}

\item{ntrees}{Number of trees default is \code{ntrees = 500}}

\item{mtry}{Number of variables used at each internal node}

\item{split.gen}{How the cutting points are generated}

\item{split.rule}{How to compare the splits}

\item{nsplit}{Number of random cutting points to compare for each variable at an internal node}

\item{nmin}{Minimum number of observations reqired in an internal node to perform a split. Set this to twice of the desired terminal node size.}

\item{alpha}{Minimum number of observations required for each child node as a portion of the parent node. Must be within \code{(0, 0.5]}.}

\item{replacement}{Whether the in-bag samples are sampled with replacement}

\item{resample.prob}{Proportion of in-bag samples}

\item{subject.weight}{Subject weights}

\item{variable.weight}{Variable weights when randomly sample \code{mtry} to select the splitting rule}

\item{importance}{Should importance measures be calculated}

\item{use.cores}{Number of cores for parallel computing}

\item{verbose}{Printing additional information}
}
\value{
A \code{survForest} object;
}
\description{
Main function to fit survival forests
}
\examples{
x = matrix(rnorm(1000), ncol = 10)
y = exp(rowMeans(x))
c = rbinom(100, 1, 0.9)
fit = survForest(x, y, c, ntrees = 10)
}
