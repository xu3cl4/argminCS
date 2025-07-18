% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_hypothesis_test.R
\name{argmax.HT}
\alias{argmax.HT}
\title{A wrapper to perform argmax hypothesis test.}
\usage{
argmax.HT(data, r = NULL, method = "softmin.LOO", ...)
}
\arguments{
\item{data}{(1) A n by p matrix of raw samples (for GTA), or
(2) A n by (p-1) difference matrix (for SML, HML, NS, MT). Each row is a sample.}

\item{r}{The dimension of interest for testing; defaults to NULL. Required for GTA.}

\item{method}{A string indicating the method to use. Defaults to 'softmin.LOO'.
See **Details** for supported methods and abbreviations.}

\item{...}{Additional arguments passed to \code{\link{argmin.HT.LOO}},
\code{\link{argmin.HT.MT}}, \code{\link{argmin.HT.nonsplit}}, or \code{\link{argmin.HT.gupta}}.}
}
\value{
A character string: 'Accept' or 'Reject', indicating whether the dimension could be an argmax, and relevant statistics.
}
\description{
This function performs a hypothesis test to evaluate whether a given dimension may be the argmax.
It internally negates the data and reuses the implementation from \code{\link{argmin.HT}}.
}
\details{
The supported methods include:\tabular{ll}{
  \code{softmin.LOO (SML)} \tab Leave-one-out algorithm using exponential weighting. \cr
  \tab \cr
  \code{argmin.LOO (HML)} \tab A variant of SML that uses hard argmin instead of exponential weighting. Not recommended. \cr
  \tab \cr
  \code{nonsplit (NS)} \tab Variant of SML without data splitting. Requires a fixed lambda value. Not recommended. \cr
  \tab \cr
  \code{Bonferroni (MT)} \tab Multiple testing using Bonferroni correction. \cr
  \tab \cr
  \code{Gupta (GTA)} \tab The method from \insertRef{gupta.1965}{argminCS}. \cr
}
}
\examples{
set.seed(108)
n <- 200
p <- 20
mu <- (1:p)/p
cov <- diag(p)
data <- MASS::mvrnorm(n, mu, cov)

## Define the dimension of interest
r <- 4

## Construct difference matrix for dimension r
difference.matrix.r <- matrix(rep(data[, r], p - 1), ncol = p - 1, byrow = FALSE) - data[, -r]

## softmin.LOO (SML)
argmax.HT(difference.matrix.r)

## use seed
argmax.HT(difference.matrix.r, seed=19)

## With known true difference
true.mean.diff <- mu[r] - mu[-r]
argmax.HT(difference.matrix.r, true.mean = true.mean.diff)

## Without scaling
argmax.HT(difference.matrix.r, scale.input = FALSE)

## With a user-specified lambda
argmax.HT(difference.matrix.r, lambda = sqrt(n) / 2.5)

## Add a seed for reproducibility
argmax.HT(difference.matrix.r, seed = 17)

## argmin.LOO (HML)
argmax.HT(difference.matrix.r, method = "HML")

## nonsplit method
argmax.HT(difference.matrix.r, method = "NS", lambda = sqrt(n)/2.5)

## Bonferroni method (choose t test for normal data)
argmax.HT(difference.matrix.r, method = "MT", test = "t")

## Gupta method (pass full data matrix)
critical.val <- get.quantile.gupta.selection(p = length(mu))
argmax.HT(data, r, method = "GTA", critical.val = critical.val)

}
\references{
\insertRef{cck.many.moments}{argminCS}

\insertRef{gupta.1965}{argminCS}

\insertRef{futschik.1995}{argminCS}
}
