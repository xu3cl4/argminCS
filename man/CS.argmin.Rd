% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/construct_confidence_set.R
\name{CS.argmin}
\alias{CS.argmin}
\title{Construct a discrete confidence set for argmin.}
\usage{
CS.argmin(data, method = "softmin.LOO", alpha = 0.05, ...)
}
\arguments{
\item{data}{A n by p data matrix; each of its row is a p-dimensional sample.}

\item{method}{A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
For the list of supported methods and their abbreviations, see Details.}

\item{alpha}{The significance level; defaults to 0.05. The function produces a (1-alpha) confidence set.}

\item{...}{Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}, \link{argmin.HT.fold}, \link{is.lambda.feasible.fold}, \link{argmin.HT.GU}, \link{argmin.HT.SN}, \link{argmin.HT.bootstrap}, \link{argmin.HT.MT}, \link{argmin.HT.gupta}.
A correct argument name needs to be specified if it is used.}
}
\value{
A vector of indices (0-based) representing the (1 - alpha) confidence set.
}
\description{
This is a wrapper to construct a discrete confidence set for argmin.. Multiple methods are supported.
}
\details{
The supported methods include:\tabular{ll}{
   \code{softmin.LOO (SML)} \tab LOO (leave-one-out) algorithm, using the exponential weightings. \cr
   \tab \cr
   \code{argmin.LOO (HML)} \tab A variant of SML, but it uses (hard) argmin rather than exponential weighting.
   The method is not recommended. \cr
   \tab \cr
   \code{nonsplit (NS)} \tab  A variant of SML, but no splitting is involved.
   One needs to pass a fixed lambda value as a required additional argument.\cr
   \tab \cr
   \code{fold (FD)} \tab A n fold version of SML; defaults to two folds. \cr
   \tab \cr
   \code{GU} \tab The method in \insertCite{dey.2024}{argminCS}. \cr
   \tab \cr
   \code{CCK.self.normalization (SN)} \tab Modified from the self-normalization method in \insertCite{cck.many.moments}{argminCS}. \cr
   \tab \cr
   \code{CCK.bootstrap (CB)} \tab Modified from the bootstrap method in \insertCite{cck.many.moments}{argminCS}. See also \insertCite{lei.cvc}{argminCS}. \cr
   \tab \cr
   \code{Bonferroni (MT)} \tab Multiple testing with Bonferroni's correction. \cr
   \tab \cr
   \code{Gupta (GTA)} \tab The method in \insertCite{gupta.1965}{argminCS}. \cr
   \tab \cr
   \code{Futschik (FCHK)} \tab The two-step method in \insertCite{futschik.1995}{argminCS} \cr
}
If computation is a concern, use 'MT'. Otherwise, 'SML' is recommended.
}
\examples{
r <- 4
n <- 200
mu <- (1:20)/20
cov <- diag(length(mu))
set.seed(108)
data <- MASS::mvrnorm(n, mu, cov)
sample.mean <- colMeans(data)

## softmin.LOO
CS.argmin(data)

## argmin.LOO
CS.argmin(data, method='HML')

## nonsplit
CS.argmin(data, method='NS', lambda=sqrt(n)/2.5)

### fold
## defaults to 2 fold
CS.argmin(data, method='FD')
## 5 fold
CS.argmin(data, method='FD', n.fold=5)

## GU
CS.argmin(data, method='GU')

## self-normalization
CS.argmin(data, method='SN')

## bootstrap
CS.argmin(data, method='CB')

## Bonferroni (choose t test because of normal data)
CS.argmin(data, method='MT', test='t')

## Gupta
CS.argmin(data, method='GTA')

## Futschik two-step method
# default alpha.1, alpha.2
CS.argmin(data, method='FCHK')

alpha.1 <- 0.0005
alpha.2 <- 1 - (0.95/(1 - alpha.1))
CS.argmin(data, method='FCHK', alpha.1=0.0005, alpha.2=alpha.2)

}
\references{
{
  \insertRef{cck.many.moments}{argminCS}

  \insertRef{lei.cvc}{argminCS}

  \insertRef{dey.2024}{argminCS}

 \insertRef{gupta.1965}{argminCS}

 \insertRef{futschik.1995}{argminCS}
}
}
