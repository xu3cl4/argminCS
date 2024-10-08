% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_hypothesis_test.R
\name{argmin.HT}
\alias{argmin.HT}
\title{A wrapper to perform argmin hypothesis test.}
\usage{
argmin.HT(data, r, method = "softmin.LOO", ...)
}
\arguments{
\item{data}{A n by p data matrix; each of its row is a p-dimensional sample.}

\item{r}{The dimension of interest for hypothesis test.}

\item{method}{A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
For the list of supported methods and their abbreviations, see Details.}

\item{...}{Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}, \link{argmin.HT.GU}, \link{argmin.HT.SN}, \link{argmin.HT.bootstrap}, \link{argmin.HT.MT}, \link{argmin.HT.gupta}.
A correct argument name needs to be specified if it is used.}
}
\value{
'Accept' or 'Reject'. A string indicating whether the given dimension could be an argmin (Accept) or not (Reject).
}
\description{
This is a wrapper to perform hypothesis test to see if a given dimension may be an argmin. Multiple methods are supported.
}
\details{
The supported methods include:\tabular{ll}{
   \code{softmin.LOO (SML)} \tab LOO (leave-one-out) algorithm, using the exponential weightings. \cr
   \tab \cr
   \code{argmin.LOO (HML)} \tab A variant of SML, but it uses (hard) argmin rather than exponential weighting.
   The method is not recommended because its type 1 error is not controlled. \cr
   \tab \cr
   \code{nonsplit (NS)} \tab  A variant of SML, but no splitting is involved.
   One needs to pass a fixed lambda value as a required additional argument.\cr
   \tab \cr
   \code{fold (FD)} \tab A n fold version of SML. \cr
   \tab \cr
   \code{GU} \tab The method in \insertCite{dey.2024}{argminCS}. \cr
   \tab \cr
   \code{CCK.self.normalization (SN)} \tab Modified from the self-normalization method in \insertCite{cck.many.moments}{argminCS}. \cr
   \tab \cr
   \code{CCK.bootstrap (CB)} \tab Modified from the bootstrap method in \insertCite{cck.many.moments}{argminCS}. See also \insertCite{lei.cvc}{argminCS}. \cr
   \tab \cr
   \code{Bonferroni (MT)} \tab Multiple testing with Bonferroni's correction. \cr
   \tab \cr
   \code{Gupta (GTA)} \tab Multiple testing with Bonferroni's correction. \cr
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
argmin.HT(data, r)
# provide centered test statistic
argmin.HT(data, r, true.mean=mu)

## argmin.LOO
argmin.HT(data, r, method='HML')

## nonsplit
argmin.HT(data, r, method='NS', lambda=sqrt(n)/2.5)

### fold
## defaults to 2 fold
argmin.HT(data, r, method='FD')
argmin.HT(data, r, method='FD', min.algor='argmin')
## 5 fold
argmin.HT(data, r, method='FD', n.fold=5)

### GU
## calculate omega first
omega <- omega.bootstrap(data, alpha=0.05)
omega

## let the function to perform data-splitting
## and other necessary hypothesis testing preparation
argmin.HT(data, r, method='GU', omega=omega)
# one can also let the function to tune omega automatically
argmin.HT(data, r, method='GU')

## split the data by ourselves
set.seed(32)
indices.training <- sample(n, n/2, replace=FALSE)
data.training <- data[indices.training,]
data.testing <- data[-indices.training,]
estimated.minimum.mean <- mean(data.testing[,which.min(colMeans(data.training))])
argmin.HT(data, 1, method='GU', omega=omega,
estimated.minimum.mean=estimated.minimum.mean, mean.r=mean(data.testing[,r]))

## self-normalization
argmin.HT(data, r, method='SN')

## bootstrap
argmin.HT(data, r, method='CB')

## Bonferroni (choose t test because of normal data)
argmin.HT(data, r, method='MT', test='t')

## Gupta (choose t test because of normal data)
critical.val <- get.quantile.gupta.selection(p=length(mu))
argmin.HT(data, r, method='GTA', critical.val=critical.val)
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
