#' Construct a discrete confidence set for argmin.
#'
#' This is a wrapper to construct a discrete confidence set for argmin.. Multiple methods are supported.
#'
#' @details The supported methods include:\tabular{ll}{
#'    \code{softmin.LOO (SML)} \tab LOO (leave-one-out) algorithm, using the exponential weightings. \cr
#'    \tab \cr
#'    \code{argmin.LOO (HML)} \tab A variant of SML, but it uses (hard) argmin rather than exponential weighting.
#'    The method is not recommended. \cr
#'    \tab \cr
#'    \code{nonsplit (NS)} \tab  A variant of SML, but no splitting is involved.
#'    One needs to pass a fixed lambda value as a required additional argument.\cr
#'    \tab \cr
#'    \code{fold (FD)} \tab A n fold version of SML.
#'    One needs to pass a fixed lambda value as a required additional argument.\cr
#'    \tab \cr
#'    \code{GU} \tab The method in \insertCite{dey.2024}{argminCS}. \cr
#'    \tab \cr
#'    \code{CCK.self.normalization (SN)} \tab Modified from the self-normalization method in \insertCite{cck.many.moments}{argminCS}. \cr
#'    \tab \cr
#'    \code{CCK.bootstrap (CB)} \tab Modified from the bootstrap method in \insertCite{cck.many.moments}{argminCS}. See also \insertCite{lei.cvc}{argminCS}. \cr
#'    \tab \cr
#'    \code{Bonferroni (MT)} \tab Multiple testing with Bonferroni's correction. \cr
#' }
#' If computation is a concern, use 'MT'. Otherwise, 'SML' is recommended.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param method A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
#' For the list of supported methods and their abbreviations, see Details.
#' @param alpha The significance level; defaults to 0.05. The function produces a (1-alpha) confidence set.
#' @param ... Additional arguments to \link{argmin.HT.LOO}, \link{argmin.HT.nonsplit}, \link{argmin.HT.fold}, \link{argmin.HT.MT}, \link{argmin.HT.SN}, \link{argmin.HT.bootstrap}.
#' A correct argument name needs to be specified if it is used.
#'
#' @return A vector of indices (0-based) representing the (1 - alpha) confidence set.
#' @export
#'
#' @examples
#' r <- 4
#' n <- 200
#' mu <- (1:20)/20
#' cov <- diag(length(mu))
#' set.seed(108)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#'
#' ## softmin.LOO
#' CS.argmin(data)
#'
#' ## argmin.LOO
#' CS.argmin(data, method='HML')
#'
#' ## nonsplit
#' CS.argmin(data, method='NS', lambda=sqrt(n)/2.5)
#'
#' ### fold
#' ## defaults to 2 fold
#' CS.argmin(data, method='FD', lambda=sqrt(n)/2.5)
#' ## 5 fold
#' CS.argmin(data, method='FD', lambda=sqrt(n)/2.5, n.fold=5)
#'
#' ## self-normalization
#' CS.argmin(data, method='SN')
#'
#' ## bootstrap
#' CS.argmin(data, method='CB')
#'
#' ## Bonferroni (choose t test because of normal data)
#' CS.argmin(data, method='MT', test='t')
#'
#' @importFrom Rdpack reprompt
#' @references{
#'   \insertRef{cck.many.moments}{argminCS}
#'
#'   \insertRef{lei.cvc}{argminCS}
#'
#'   \insertRef{dey.2024}{argminCS}
#' }
#'
CS.argmin <- function(data, method='softmin.LOO', alpha=0.05, ...){

  p <- ncol(data)

  if (method == 'softmin.LOO' | method == 'SML'){
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.LOO(data, r, sample.mean=sample.mean, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'argmin.LOO' | method == 'HML') {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.LOO(
      data, r, min.algor=getMin.argmin.LOO, sample.mean=sample.mean, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'nonsplit' | method == 'NS') {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.nonsplit(data, r, sample.mean=sample.mean, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'fold' | method == 'FD') {
    res <- sapply(1:p, function(r) {argmin.HT.fold(data, r, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'GU') {

  } else if (method == 'CCK.self.normalization' | method == 'SN'){
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.SN(data, r, sample.mean=sample.mean, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'CCK.bootstrap' | method == 'CB') {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.bootstrap(data, r, sample.mean=sample.mean, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'Bonferroni' | method == 'MT') {
    sample.mean <- colMeans(data) # np
    min.indices <- which(sample.mean == min(sample.mean))
    r.min <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
    r.min.sec <- find.sub.argmin(sample.mean, r.min)
    res <- sapply(1:p, function(r) {argmin.HT.MT(data, r, r.min=r.min, r.min.sec=r.min.sec, ...)$ans})
    return (which(res == 'Accept'))

  } else {
    stop("'method' should be one of 'softmin.LOO' (SML), 'argmin.LOO' (HML),
         'nonsplit' (NS), 'fold' (FD), 'GU', 'CCK.bootstrap' (CB), 'Bonferroni' (MT)")
  }
}

# CS.GU <- function(data, alpha, omega=NULL){
#
#   n <- nrow(data)
#   p <- ncol(data)
#
#   ## split the data
#   withr::with_seed(p*n*data[1,1], {
#     idx.tr <- sample(1:n, n/2, replace=F)
#   })
#   data.tr <- data[idx.tr,]
#   data.tt <- data[-idx.tr,]
#
#   # get the best model over training set
#   sample.mean.tr <- colMeans(data.tr)
#   min.indices <- which(sample.mean.tr == min(sample.mean.tr))
#   idx.min.tr <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
#
#   # evaluate the best model from the training set over the testing set
#   sample.mean.tt <- colMeans(data.tt)
#   risk.idx.min.tr <- sample.mean.tt[idx.min.tr]
#
#   # compute the learning rate
#   if (is.null(omega)){
#     omega <- omega.bootstrap(data, alpha)
#   }
#
#   idx <- 1:p
#   res <- rep(NA, p)
#
#   res.s <- idx[sapply(1:length(idx), function(j) argmin.HT.GU(
#     sample.mean.tt[j], alpha, risk.idx.min.tr, omega, idx=j))]
#   res.s <- na.omit(res.s)
#   res[res.s] <- res.s
#   return (res)
# }


