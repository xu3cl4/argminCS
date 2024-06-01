CS.argmin <- function(data, method, alpha=0.05){

}

CS <- function(Xs, lambda, alpha, y.algor=getY.softmax){
  # Xs: n by p
  # lambda: parameter for exponential weighting
  # alpha: a desired significance level
  p <- ncol(Xs)
  idx <- 1:p
  mean.smp <- colMeans(Xs)

  res <- rep(NA, p)
  if (length(idx) == 1) {
    res[idx] <- idx
    return (res)
  }
  res.s <- idx[sapply(1:length(idx), function(j) argmin.HT(Xs, mean.smp, j, lambda, alpha, y.algor))]
  res.s <- na.omit(res.s)
  res[res.s] <- res.s
  return (res)
}

CS.adaptive <- function(Xs, alpha, const, y.algor=getY.softmax, enlarge=F){
  # Xs: n by p
  # lambda: parameter for exponential weighting
  # alpha: a desired significance level
  # 8 np^2
  p <- ncol(Xs)
  idx <- 1:p
  mean.smp <- colMeans(Xs)
  if (pre.screening){
    idx <- screening(Xs, mean.smp)
    Xs <- Xs[,idx]
    mean.smp <- mean.smp[idx]
  }

  res <- rep(NA, p)
  if (length(idx) == 1){
    res[idx] <- idx
    return (res)
  }
  # p
  res.s <- idx[sapply(1:length(idx), function(j)
    argmin.HT.adaptive(Xs, mean.smp, j, alpha, const, y.algor, enlarge))]
  res.s <- na.omit(res.s)
  res[res.s] <- res.s
  return (res)
}

CS.GU <- function(X, alpha, omega=NULL){

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(omega)){
    omega <- omega.bootstrap(X, alpha)
  }

  ## split the data
  set.seed(p*n)
  idx.tr <- sample(1:n, n/2, replace=F)
  X.tr <- X[idx.tr,]
  X.tt <- X[-idx.tr,]

  # get the best model over training set
  risks.tr <- colMeans(X.tr)
  idx.min.tr <- which.min(risks.tr)

  # evaluate the best model from the training set over the testing set
  risk.idx.min.tr <- mean(X.tt[idx.min.tr,])

  # compute the learning rate
  if (is.null(omega)){
    omega <- omega.bootstrap(X, alpha)
  }

  idx <- 1:p
  res <- rep(NA, p)
  smp.mean.tt <- colMeans(X.tt)
  res.s <- idx[sapply(1:length(idx), function(j) argmin.HT.GU(
    smp.mean.tt[j], alpha, risk.idx.min.tr, omega, idx=j))]
  res.s <- na.omit(res.s)
  res[res.s] <- res.s
  return (res)
}

CS.non.split <- function(Xs, lambda, alpha, y.algor=getY.softmax, pre.screening=F, n.fold=2){
  # Xs: n by p
  # lambda: parameter for exponential weighting
  # alpha: a desired significance level
  p <- ncol(Xs)
  idx <- 1:p
  mean.smp <- colMeans(Xs)
  if (pre.screening){
    idx <- screening(Xs, mean.smp)
    Xs <- Xs[,idx]
    mean.smp <- mean.smp[idx]
  }

  res <- rep(NA, p)
  if (length(idx) == 1) {
    res[idx] <- idx
    return (res)
  }
  res.s <- idx[sapply(1:length(idx), function(j) argmin.HT.nonsplit(Xs, mean.smp, j, lambda, alpha, y.algor))]
  res.s <- na.omit(res.s)
  res[res.s] <- res.s
  return (res)
}

CS.split <- function(Xs, lambda, alpha, pre.screening=F, n.fold=2){
  # Xs: n by p
  # lambda: parameter for expoential weighting
  # alpha: a desired significance level
  p <- ncol(Xs)
  idx <- 1:p
  mean.smp <- NULL
  if (pre.screening){
    mean.smp <- colMeans(Xs)
    idx <- screening(Xs, mean.smp)
    Xs <- Xs[,idx]
    mean.smp <- mean.smp[idx]
  }

  res <- rep(NA, p)
  if (length(idx) == 1) {
    res[idx] <- idx
    return (res)
  }
  res.s <- idx[sapply(1:length(idx), function(j) argmin.HT.fold(Xs, j, lambda, alpha, n.fold=n.fold))]
  res.s <- na.omit(res.s)
  res[res.s] <- res.s
  return (res)
}

CS.naive <- function(Xs, alpha){
  # Xs: n by p
  # alpha: the desired significance level
  p <- ncol(Xs)
  n <- nrow(Xs)
  val.critical <- qnorm(1-alpha, 0, 1)

  mean <- colMeans(Xs) # np
  j.star <- which.min(mean) #p
  # np
  res <- sapply(1:p, function(j){
    if (j == j.star){
      return (j.star)
    }
    else{
      test.stat <- sqrt(n)*(mean[j] - mean[j.star])
      diffs <- Xs[,j] - Xs[,j.star] #n
      sigma <- sd(diffs) # n
      test.stat.scale <- test.stat/sigma
    }
    return (if (test.stat.scale < val.critical) j else NA)
  })
  return (res)
}

CS.MT <- function(Xs, alpha, test='t'){
  # Xs:    n by p
  # alpha: the desired significance level to be Bonferroni corrected
  # test:  the specified test for hypothesis testing
  #        two options: (1) 't': t test (2) 'z': z test
  p <- ncol(Xs)
  n <- nrow(Xs)
  val.critical <- alpha/(p-1)

  mean <- colMeans(Xs) # np
  j.star <- which.min(mean) #p

  res <- NULL
  if (test == 't'){
    # t test
    res <- sapply(1:p, function(j){
      if (j == j.star){
        return (j.star)
      } else{
        p.val <- t.test(Xs[,j]-Xs[,j.star], alternative='greater')$p.value
        return (if (p.val > val.critical) j else NA)
      }
    })
  } else{
    # z test
    res <- sapply(1:5, function(j){
      if (j == j.star){
        return (j.star)
      } else{
        #print(sd(Xs[,j]))
        #print(sd(Xs[,j.star]))
        p.val <- z.test(Xs[,j], sigma.x=sd(Xs[,j]),
                        Xs[,j.star], sigma.y=sd(Xs[,j.star]),
                        alternative='greater')$p.value
        return (if (p.val > val.critical) j else NA)
      }
    })
    res <- c(res, rep(NA, p-s-1))
  }
  return (res)
}

CS.CCK.self.normalized <- function(Xs, alpha){
  ## approximately 4*p^2*n
  n <- nrow(Xs)
  p <- ncol(Xs)

  norm.quantile <- qnorm(1 - alpha/(p-1))
  val.critical <- ifelse(norm.quantile^2 >= n,Inf,norm.quantile/sqrt(1-norm.quantile^2/n))

  mean.smp <- colMeans(Xs)
  res <- sapply(1:p, function(j){
    mean.diffs <- mean.smp[j] - mean.smp[-j] # p
    diffs <- matrix(rep(Xs[,j], p-1), nrow=n, byrow=F) - Xs[,-j] # np
    sd.diffs <- apply(diffs, 2, sd) # 3np
    ratios <- mean.diffs/sd.diffs
    test.stat <- sqrt(n)*max(ratios)
    return (if (test.stat <= val.critical) j else NA)
  })

  return (res)
}

CS.CCK.bootstrap <- function(data, alpha){

}
