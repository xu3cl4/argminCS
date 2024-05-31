argmin.HT <- function(Xs, mean.smp, j, lambda, alpha, y.algor){
  # data: n by p
  # j : the index to test
  # lambda: parameter for exponential weighting
  n <- nrow(Xs)
  val.critical <- qnorm(1-alpha, 0, 1)

  ys <- sapply(1:n, function(i) return (y.algor(i, j, Xs, mean.smp, lambda)))
  diffs <- Xs[,j] - ys
  sigma <- sd(diffs)
  test.stat <- sqrt(n)*(mean.smp[j] - mean(ys))
  test.stat.scale <- test.stat/sigma
  return (if (test.stat.scale < val.critical) j else NA)
}

argmin.HT.adaptive <- function(Xs, mean.smp, j, alpha, const, y.algor, enlarge=F){
  # data: n by p
  # j : the index to test
  # lambda: parameter for exponential weighting
  n <- nrow(Xs)
  val.critical <- qnorm(1-alpha, 0, 1)

  lambda <- lambda.adaptive(Xs, mean.smp, j, const) # n + p
  if (enlarge) {
    lambda <- lambda.adaptive.enlarge(lambda, Xs, mean.smp, j, const=0.05)
  }
  ys <- sapply(1:n, function(i) return (y.algor(i, j, Xs, mean.smp, lambda))) # 4np
  diffs <- Xs[,j] - ys # n
  sigma <- sd(diffs) # 2n
  test.stat <- sqrt(n)*(mean.smp[j] - mean(ys))
  test.stat.scale <- test.stat/sigma
  # if (j == 39){
  #   print(glue('{j}: {list(round(diffs, 2))} {sigma} {test.stat} {test.stat.scale}'))
  # }
  return (if (test.stat.scale < val.critical) j else NA)
}

## non-split argmin algorithm
argmin.HT.nonsplit <- function(Xs, mean.smp, j, lambda, alpha, y.algor){
  # data: n by p
  # j : the index to test
  # lambda: parameter for exponential weighting
  n <- nrow(Xs)
  val.critical <- qnorm(1-alpha, 0, 1)

  weights <- softmax(-lambda*mean.smp[-j])
  ys <- Xs[,-j] %*% weights
  diffs <- Xs[,j] - ys

  var.smp <- sum((diffs - mean(diffs))^2)/(n-1)

  test.stat <- sqrt(n)*mean(diffs)
  test.stat.scale <- test.stat/sqrt(var.smp)

  return (if (test.stat.scale < val.critical) j else NA)
}

## split argmin algorithm
argmin.HT.fold <- function(Xs, j, lambda, alpha, n.fold=2){
  # data: n by p
  # j : the index to test
  # lambda: parameter for exponential weighting
  n <- nrow(Xs)
  p <- ncol(Xs)
  val.critical <- qnorm(1-alpha, 0, 1)

  if (p == 2){
    diffs <- Xs[,j] - Xs[,-j]
    var.smp <- sum((diffs - mean(diffs))^2)/(n-1)
    test.stat <- sqrt(n)*mean(diffs)
    test.stat.scale <- test.stat/sqrt(var.smp)
    return (if (test.stat.scale < val.critical) j else NA)
  } else{
    ### create folds
    set.seed(13*n.fold)
    flds <- createFolds(1:n, k=n.fold, list=T, returnTrain=F)

    diffs <- lapply(1:n.fold, function(fold) {
      mu.out.fold <- colMeans(Xs[setdiff(1:n, flds[[fold]]),])
      in.fold <- Xs[flds[[fold]],]
      weights <- softmax(-lambda*mu.out.fold[-j])
      ys <- in.fold[,-j] %*% weights
      diffs.fold <- in.fold[,j] - ys
      return (diffs.fold)
    })
    diffs <- do.call(c, diffs)
    var.smp <- sum((diffs - mean(diffs))^2)/(n-1)

    test.stat <- sqrt(n)*mean(diffs)
    test.stat.scale <- test.stat/sqrt(var.smp)

    return (if (test.stat.scale < val.critical) j else NA)
  }
}

omega.bootstrap <- function(X, alpha, B=100){
  # Xs: n by p
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
  ## note that we may not need this subroutine
  ## we may try to simply calculate the test statistic all at once

  test.stat <- exp(-omega*(risk.best.idx.tr - risk.theta))
  # reject the test stat if test >= 1/alpha; otherwise, keep it
  return (if (test.stat < 1/alpha) idx else NA)
}
