#' Perform argmin hypothesis test
#'
#' This is a wrapper to perform hypothesis test to see if a given dimension may be an argmin. Multiple methods are supported.
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
#'    \code{fold (FD)} \tab A n fold version of SML. \cr
#'    \tab \cr
#'    \code{GU} \tab The method in \insertCite{dey.2024}{argminCS}. \cr
#'    \tab \cr
#'    \code{CCK.self.normalization (SN)} \tab Modified from the self-normalization method in \insertCite{cck.many.moments}{argminCS}. \cr
#'    \tab \cr
#'    \code{CCK.bootstrap (CB)} \tab Modified from the bootstrap method in \insertCite{cck.many.moments}{argminCS}. See also \insertCite{lei.cvc}{argminCS}. \cr
#'    \tab \cr
#'    \code{Bonferroni (MT)} \tab Multiple testing with Bonferroni's correction. \cr
#' }
#' If computation is a concern, use 'SN' or 'MT'. Otherwise, 'SML' is recommended.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param method A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
#' For the list of supported methods and their abbreviations, see Details.
#' @param ... Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible}.
#' A correct argument name needs to be specified if it is used.
#'
#' @return 'Accept' or 'Reject'. A string indicating whether the given dimension could be an argmin (Accept) or not (Reject).
#' @export
#'
#' @examples
#'
#' @references{
#'   \insertRef{cck.many.moments}{argminCS}
#'
#'   \insertRef{lei.cvc}{argminCS}
#'
#'   \insertRef{dey.2024}{argminCS}
#' }
#'
argmin.HT <- function(data, r, method='softmin.LOO', ...){
  if (method == 'softmin.LOO' | method == 'SML'){
    return (argmin.HT.LOO(data, r, ...)$ans)

  } else if (method == 'argmin.LOO' | method == 'HML') {
    return (argmin.HT.LOO(data, r, min.algor=getMin.argmin.LOO, ...)$ans)

  } else if (method == 'nonsplit' | method == 'NS') {
    return (argmin.HT.nonsplit(data, r, ...)$ans)

  } else if (method == 'fold' | method == 'FD') {
    return (argmin.HT.fold(data, r, ...)$ans)

  } else if (method == 'GU') {

  } else if (method == 'CCK.self.normalization' | method == 'SN') {

  } else if (method == 'CCK.bootstrap' | method == 'CB') {

  } else if (method == 'Bonferroni' | method == 'MT') {
    return (argmin.HT.MT(data, r, ...)$ans)

  } else {
    stop("'method' should be one of 'softmin.LOO' (SML), 'argmin.LOO' (HML),
         'nonsplit' (NS), 'fold' (FD), 'GU', 'CCK.bootstrap' (CB), 'Bonferroni' (MT)")
  }
}

#' Perform argmin hypothesis test, using the LOO (leave-one-out) algorithm.
#'
#' Test if a dimension may be argmin, using the LOO (leave-one-out) algorithm in Zhang et al 2024.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(data) to speed up computation.
#' @param min.algor The algorithm to find the minimum excluding the r-th dimension; defaults to \link{getMin.softmin.LOO}. The other option is \link{getMin.argmin.LOO}.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin); defaults to NULL.
#' If lambda=NULL (recommended), the function would determine a lambda value in a data-driven way.
#' @inheritParams lambda.adaptive
#' @param enlarge A boolean value indicating if the data-driven lambda should be determined via an iterative enlarging algorithm; defaults to TRUE.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param ... Additional arguments to \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible}.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat.scale} \tab The scaled test statistic \cr
#'    \tab \cr
#'    \code{std} \tab The standard deviation estimate. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
argmin.HT.LOO <- function(data, r, sample.mean=NULL, min.algor=getMin.softmin.LOO,
                          lambda=NULL, const=2.5, enlarge=TRUE, alpha=0.05, ...){

  n <- nrow(data)
  p <- ncol(data)
  val.critical <- qnorm(1-alpha, 0, 1)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  if (is.null(lambda)){
    lambda <- lambda.adaptive(data, r, sample.mean=sample.mean, const=const)
    if (enlarge) {
      lambda <- lambda.adaptive.enlarge(
        lambda, data, r, sample.mean=sample.mean, ...)
    }
  }

  Qs <- NULL
  if (p == 2){
    Qs <- data[,-r]
  } else {
    Qs <- sapply(1:n, function(i) return (min.algor(i, r, data, lambda, sample.mean=sample.mean))) # 4np
  }
  diffs <- data[,r] - Qs # n
  sigma <- sd(diffs) # 2n
  test.stat <- sqrt(n)*(sample.mean[r] - mean(Qs))
  test.stat.scale <- test.stat/sigma
  ans <- ifelse(test.stat.scale <= val.critical, 'Accept', 'Reject')
  return (list(test.stat.scale=test.stat.scale, std=sigma, ans=ans))
}


#' Perform argmin hypothesis test without any splitting.
#'
#' Test if a dimension may be argmin withtout any splitting.
#'
#' @details
#' This method is not recommended, given its poor performance when p is small.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(data) to speed up computation.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat.scale} \tab The scaled test statistic \cr
#'    \tab \cr
#'    \code{std} \tab The standard deviation estimate. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
argmin.HT.nonsplit <- function(data, r, lambda, sample.mean=NULL, alpha=0.05){
  n <- nrow(data)
  val.critical <- qnorm(1-alpha, 0, 1)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  weights <- LDATS::softmax(-lambda*sample.mean[-r])
  Qs <- data[,-r] %*% weights
  diffs <- data[,r] - Qs

  sigma <- sd(diffs) # 2n
  test.stat <- sqrt(n)*(sample.mean[r] - mean(Qs))
  test.stat.scale <- test.stat/sigma

  ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')
  return (list(test.stat.scale=test.stat.scale, std=sigma, ans=ans))
}


#' Perform argmin hypothesis test by splitting into folds.
#'
#' Test if a dimension may be argmin by splitting into folds.
#'
#' @details
#' This method does not support any data-driven way to tune lambda for now. Try to use argmin.HT.LOO instead.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param n.fold The number of folds.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat.scale} \tab The scaled test statistic \cr
#'    \tab \cr
#'    \code{std} \tab The standard deviation estimate. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
argmin.HT.fold <- function(data, r, lambda, alpha=0.05, n.fold=2){
  n <- nrow(data)
  p <- ncol(data)
  val.critical <- qnorm(1-alpha, 0, 1)

  if (n.fold == n){
    return (argmin.HT.LOO(data, r, sample.mean=colMeans(data), lambda=lambda))
  }

  if (p == 2){
    diffs <- data[,j] - data[,-r]
    var.smp <- sum((diffs - mean(diffs))^2)/(n-1)
    test.stat <- sqrt(n)*mean(diffs)
    test.stat.scale <- test.stat/sqrt(var.smp)
    ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')
    return (list(test.stat.scale=test.stat.scale, std=sigma, ans=ans))
  } else{
    ### create folds
    set.seed(13*n.fold)
    flds <- caret::createFolds(1:n, k=n.fold, list=T, returnTrain=F)

    diffs <- lapply(1:n.fold, function(fold) {
      mu.out.fold <- colMeans(data[setdiff(1:n, flds[[fold]]),])
      in.fold <- data[flds[[fold]],]
      weights <- LDATS::softmax(-lambda*mu.out.fold[-r])
      Qs <- in.fold[,-r] %*% weights
      diffs.fold <- in.fold[,j] - Qs
      return (diffs.fold)
    })
    diffs <- do.call(c, diffs)
    var.smp <- sum((diffs - mean(diffs))^2)/(n-1)

    test.stat <- sqrt(n)*mean(diffs)
    test.stat.scale <- test.stat/sqrt(var.smp)

    ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')
    return (list(test.stat.scale=test.stat.scale, std=sigma, ans=ans))
  }
}

## TO DO
omega.bootstrap <- function(X, alpha, B=100){
  # data: n by p
  # alpha: a desired significance level
  n <- nrow(X)
  p <- nrow(X)

  risks <- colMeans(X)
  idx.theta.hat <- which.min(risks)

  omegas <- 1:100
  coverages <- lapply(1:B, function(i){
    set.seed(13*i)
    smp <- sample(1:n, n, replace=T)
    X.boot <- X[smp,]

    # split the data
    set.seed(i)
    idx.tr <- sample(1:n, n/2, replace=F)
    X.tr <- X.boot[idx.tr,]
    X.tt <- X.boot[-idx.tr,]

    # get the best model over training set
    risks.tr <- colMeans(X.tr)
    idx.min.tr <- which.min(risks.tr)

    # evaluate the best model from the training set over the testing set
    risk.idx.min.tr <- mean(X.tt[idx.min.tr,])

    # give any idx (it's simply a dummy value)
    risk.theta.hat <- mean(X.tt[idx.theta.hat,])
    res <- sapply(omegas, function(omega)
      argmin.HT.GU(risk.theta.hat, alpha, risk.idx.min.tr, omega, idx=1))
    res <- ifelse(is.na(res), 0, 1)
    return (res)
  })

  coverages <- as.matrix(do.call('rbind', coverages))
  coverages <- colSums(coverages)
  return (omegas[which.min(abs(coverages - (1 - alpha)*B))])
}

argmin.HT.GU <- function(risk.theta, alpha, risk.best.idx.tr, omega, idx=1){

  test.stat <- exp(-omega*(risk.best.idx.tr - risk.theta))
  # reject the test stat if test >= 1/alpha; otherwise, keep it
  return (if (test.stat < 1/alpha) idx else NA)
}

#' Perform argmin hypothesis test, using multiple testing with Bonferroni's correction.
#'
#' Test if a dimension may be argmin, using multiple testing with Bonferroni's correction.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param r.min The sample argmin of the data; defaults to NULL. It can be calculated via which.min(colMeans(data)).
#' @param test The test to perform: 't' or 'z' test; defaults to 'z'.
#' If data are believed to be normally distributed, use 't'; otherwise 'z'.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{p.val} \tab p value without Bonferroni's correction. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
argmin.HT.MT <- function(data, r, r.min=NULL, test='z', alpha=0.05){

  val.critical <- alpha/(p-1)

  if (is.null(r.min)){
    mean <- colMeans(Xs) # np
    r.min <- which.min(mean) #p
  }

  p.val <- NULL
  if (test == 't'){
    # t test
    p.val <- stats::t.test(Xs[,j]-Xs[,r.min], alternative='greater')$p.value
  } else{
    # z.test
    diffs <- Xs[,j] - Xs[,r.min]
    p.val <- BSDA::z.test(diffs, sigma.x=sd(diffs), alternative='greater')$p.value
  }

  ans <- ifelse(p.val > critical, 'Accept', 'Reject')
  return (list(p.val=p.val, ans=ans))
}
