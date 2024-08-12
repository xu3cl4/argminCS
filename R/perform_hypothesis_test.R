#' A wrapper to perform argmin hypothesis test.
#'
#' This is a wrapper to perform hypothesis test to see if a given dimension may be an argmin. Multiple methods are supported.
#'
#' @details The supported methods include:\tabular{ll}{
#'    \code{softmin.LOO (SML)} \tab LOO (leave-one-out) algorithm, using the exponential weightings. \cr
#'    \tab \cr
#'    \code{argmin.LOO (HML)} \tab A variant of SML, but it uses (hard) argmin rather than exponential weighting.
#'    The method is not recommended because its type 1 error is not controlled. \cr
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
#'    \tab \cr
#'    \code{Gupta (GTA)} \tab Multiple testing with Bonferroni's correction. \cr
#' }
#' If computation is a concern, use 'MT'. Otherwise, 'SML' is recommended.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param method A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
#' For the list of supported methods and their abbreviations, see Details.
#' @param ... Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}, \link{argmin.HT.GU}, \link{argmin.HT.SN}, \link{argmin.HT.bootstrap}, \link{argmin.HT.MT}, \link{argmin.HT.gupta}.
#' A correct argument name needs to be specified if it is used.
#'
#' @return 'Accept' or 'Reject'. A string indicating whether the given dimension could be an argmin (Accept) or not (Reject).
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
#' argmin.HT(data, r)
#' # provide centered test statistic
#' argmin.HT(data, r, true.mean=mu)
#'
#' ## argmin.LOO
#' argmin.HT(data, r, method='HML')
#'
#' ## nonsplit
#' argmin.HT(data, r, method='NS', lambda=sqrt(n)/2.5)
#'
#' ### fold
#' ## defaults to 2 fold
#' argmin.HT(data, r, method='FD')
#' argmin.HT(data, r, method='FD', min.algor='argmin')
#' ## 5 fold
#' argmin.HT(data, r, method='FD', n.fold=5)
#'
#' ### GU
#' ## calculate omega first
#' omega <- omega.bootstrap(data, alpha=0.05)
#' omega
#'
#' ## let the function to perform data-splitting
#' ## and other necessary hypothesis testing preparation
#' argmin.HT(data, r, method='GU', omega=omega)
#' # one can also let the function to tune omega automatically
#' argmin.HT(data, r, method='GU')
#'
#' ## split the data by ourselves
#' set.seed(32)
#' indices.training <- sample(n, n/2, replace=FALSE)
#' data.training <- data[indices.training,]
#' data.testing <- data[-indices.training,]
#' estimated.minimum.mean <- mean(data.testing[,which.min(colMeans(data.training))])
#' argmin.HT(data, 1, method='GU', omega=omega,
#' estimated.minimum.mean=estimated.minimum.mean, mean.r=mean(data.testing[,r]))
#'
#' ## self-normalization
#' argmin.HT(data, r, method='SN')
#'
#' ## bootstrap
#' argmin.HT(data, r, method='CB')
#'
#' ## Bonferroni (choose t test because of normal data)
#' argmin.HT(data, r, method='MT', test='t')
#'
#' ## Gupta (choose t test because of normal data)
#' critical.val <- get.quantile.gupta.selection(p=length(mu))
#' argmin.HT(data, r, method='GTA', critical.val=critical.val)
#' @importFrom Rdpack reprompt
#' @references{
#'   \insertRef{cck.many.moments}{argminCS}
#'
#'   \insertRef{lei.cvc}{argminCS}
#'
#'   \insertRef{dey.2024}{argminCS}
#'
#'   \insertRef{gupta.1965}{argminCS}
#'
#'   \insertRef{futschik.1995}{argminCS}
#' }
argmin.HT <- function(data, r, method='softmin.LOO', ...){
  if (method == 'softmin.LOO' | method == 'SML'){
    return (argmin.HT.LOO(data, r, ...))

  } else if (method == 'argmin.LOO' | method == 'HML') {
    return (argmin.HT.LOO(data, r, min.algor=getMin.argmin.LOO, ...))

  } else if (method == 'nonsplit' | method == 'NS') {
    return (argmin.HT.nonsplit(data, r, ...))

  } else if (method == 'fold' | method == 'FD') {
    return (argmin.HT.fold(data, r, ...))

  } else if (method == 'GU') {
    return (argmin.HT.GU(data, r, ...))

  } else if (method == 'CCK.self.normalization' | method == 'SN'){
    return (argmin.HT.SN(data, r, ...))

  } else if (method == 'CCK.bootstrap' | method == 'CB') {
    return (argmin.HT.bootstrap(data, r, ...))

  } else if (method == 'Bonferroni' | method == 'MT') {
    return (argmin.HT.MT(data, r, ...))

  } else if (method == 'Gupta' | method == 'GTA' | method=='gupta') {
    return (argmin.HT.gupta(data, r, ...))

  } else {
    stop("'method' should be one of 'softmin.LOO' (SML), 'argmin.LOO' (HML),
         'nonsplit' (NS), 'fold' (FD), 'GU', 'CCK.bootstrap' (CB), 'Bonferroni' (MT)")
  }
}

#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, using the LOO (leave-one-out) algorithm in Zhang et al 2024.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, compute colMeans(data) and pass it to sample.mean to speed up computation.
#' @param min.algor The algorithm to find the minimum excluding the r-th dimension; defaults to \link{getMin.softmin.LOO}. The other option is \link{getMin.argmin.LOO}.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin); defaults to NULL.
#' If lambda=NULL (recommended), the function would determine a lambda value in a data-driven way.
#' @inheritParams lambda.adaptive
#' @param enlarge A boolean value indicating if the data-driven lambda should be determined via an iterative enlarging algorithm; defaults to TRUE.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param true.mean The population mean landscape; defaults to NULL. If a vector were provided, the centered test statistic would be outputted.
#' It is only useful for a simulation purpose.
#' @param ... Additional arguments to \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat.scale} \tab The scaled test statistic \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being greater than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{std} \tab The standard deviation estimate. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#'    \tab \cr
#'    \code{lambda} \tab The lambda used in the hypothesis testing. \cr
#'    \tab \cr
#'    \code{test.stat.centered} \tab (Optional) The centered test statistic. Outputted only when true.mean is not NULL. \cr
#' }
argmin.HT.LOO <- function(data, r, sample.mean=NULL, min.algor=getMin.softmin.LOO,
                          lambda=NULL, const=2.5, enlarge=TRUE, alpha=0.05, true.mean=NULL, ...){

  n <- nrow(data)
  p <- ncol(data)
  val.critical <- stats::qnorm(1-alpha, 0, 1)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  if (is.null(lambda) & identical(deparse(min.algor), deparse(getMin.softmin.LOO))){
    lambda <- lambda.adaptive(data, r, sample.mean=sample.mean, const=const)
    if (enlarge) {
      lambda <- lambda.adaptive.enlarge(
        lambda, data, r, algorithm='LOO', sample.mean=sample.mean, ...)
    }
  }

  Qs <- NULL
  if (p == 2){
    Qs <- data[,-r]
    ## centered test statistic
    Qs.true.mean <- NULL
    if (!is.null(true.mean)){
      Qs.true.mean <- true.mean[1] - true.mean[2]
    }
  } else {
    #Qs <- sapply(1:n, function(i) return (min.algor(i, r, data, lambda, sample.mean=sample.mean))) # 4np
    ## centered test statistic
    res <- lapply(1:n, function(i) return (min.algor(i, r, data, lambda, sample.mean=sample.mean, true.mean=true.mean))) # 4np
    res <- do.call(rbind, res)
    Qs <- unlist(res[,1])
    Qs.true.mean <- unlist(res[,2])
  }
  diffs <- data[,r] - Qs # n
  sigma <- stats::sd(diffs) # 2n
  test.stat <- sqrt(n)*(sample.mean[r] - mean(Qs))
  test.stat.scale <- test.stat/sigma
  ans <- ifelse(test.stat.scale <= val.critical, 'Accept', 'Reject')

  if (is.null(true.mean)) {
    return (list(test.stat.scale=test.stat.scale, critical.value=val.critical, std=sigma, ans=ans,
                 lambda=lambda))
  } else {
    ## provide test.stat.centered
    random.center <- (sqrt(n)*(true.mean[r] - mean(Qs.true.mean)))/sigma
    test.stat.centered <- test.stat.scale - random.center
    return (list(test.stat.scale=test.stat.scale, critical.value=val.critical, std=sigma, ans=ans,
                 lambda=lambda, test.stat.centered=test.stat.centered))
  }
}


#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin without any splitting.
#'
#' @details
#' This method is not recommended, given its poor performance when p is small.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, compute colMeans(data) and pass it to sample.mean to speed up computation.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat.scale} \tab The scaled test statistic \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being greater than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{std} \tab The standard deviation estimate. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
argmin.HT.nonsplit <- function(data, r, lambda, sample.mean=NULL, alpha=0.05){
  n <- nrow(data)
  val.critical <- stats::qnorm(1-alpha, 0, 1)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  weights <- LDATS::softmax(-lambda*sample.mean[-r])
  Qs <- data[,-r] %*% weights
  diffs <- data[,r] - Qs

  sigma <- stats::sd(diffs) # 2n
  test.stat <- sqrt(n)*(sample.mean[r] - mean(Qs))
  test.stat.scale <- test.stat/sigma

  ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')
  return (list(test.stat.scale=test.stat.scale, critical.value=val.critical, std=sigma, ans=ans))
}


#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin by splitting into folds.
#'
#' @details
#' This method does not support any data-driven way to tune lambda for now. Try to use argmin.HT.LOO instead.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param n.fold The number of folds.
#' @param flds A list of row position integers corresponding to folds; defaults to NULL.
#' @param min.algor The algorithm to compute the minimum: 'softmin' or 'argmin'; defaults to 'softmin'
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin); defaults to NULL.
#' If lambda=NULL (recommended), the function would determine a lambda value in a data-driven way.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' It is used to tune lambda in a data-driven way.
#' @inheritParams lambda.adaptive.fold
#' @param enlarge A boolean value indicating if the data-driven lambda should be determined via an iterative enlarging algorithm; defaults to TRUE.
#' @param true.mean The population mean landscape; defaults to NULL. If a vector were provided, the centered test statistic would be outputted.
#' It is only useful for a simulation purpose.
#' @param ... Additional arguments to \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.fold}.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat.scale} \tab The scaled test statistic \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being greater than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{std} \tab The standard deviation estimate. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#'    \tab \cr
#'    \code{lambda} \tab The lambda used in the hypothesis testing. \cr
#'    \tab \cr
#'    \code{test.stat.centered} \tab (Optional) The centered test statistic. Outputted only when true.mean is not NULL. \cr
#' }
argmin.HT.fold <- function(data, r, alpha=0.05, n.fold=2, flds=NULL, sample.mean=NULL,
                           min.algor='softmin', lambda=NULL, const=2.5, enlarge=TRUE, true.mean=NULL, ...){
  n <- nrow(data)
  p <- ncol(data)
  val.critical <- stats::qnorm(1-alpha, 0, 1)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  if (n.fold == n){
    # use the LOO method
    return (argmin.HT.LOO(data, r, sample.mean=sample.mean, lambda=lambda, true.mean=true.mean))
  }

  test.stat.centered <- NULL
  if (p == 2){
    diffs <- data[,r] - data[,-r]
    sigma <- stats::sd(diffs)
    test.stat <- sqrt(n)*mean(diffs)
    test.stat.scale <- test.stat/sigma
    ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')
    ## ouptput centered test statistic
    if (!is.null(true.mean)){
      test.stat.centered <- test.stat.scale - sqrt(n)*(true.mean[r] - true.mean[-r])/sigma
    }
    return (list(test.stat.scale=test.stat.scale,
                 std=sigma, ans=ans, test.stat.centered=test.stat.centered))
  } else{
    ### create folds if not given
    if (is.null(flds)){
      seed <- ceiling(abs(13*r*data[1,r]*n.fold)) + r
      if (seed >  2^31 - 1){
        seed <- seed %% (2^31 - 1)
      }
      withr::with_seed(seed, {
        flds <- caret::createFolds(1:n, k=n.fold, list=TRUE, returnTrain=FALSE)
      })
    }

    if (is.null(lambda)){
      lambda <- lambda.adaptive.fold(data, r, flds, const=const)

      if (enlarge) {
        lambda <- lambda.adaptive.enlarge(
          lambda, data, r, algorithm='fold', sample.mean=sample.mean, flds=flds, ...)
      }
    }

    diffs <- lapply(1:n.fold, function(fold) {
      mu.out.fold <- colMeans(data[setdiff(1:n, flds[[fold]]),])
      in.fold <- data[flds[[fold]],]
      if (min.algor == 'softmin'){
        # try softmin
        weights <- LDATS::softmax(-lambda*mu.out.fold[-r])
        # testing
        Qs <- in.fold[,-r] %*% weights
        diffs.fold <- in.fold[,r] - Qs

        ## output centered test statistic
        Q.true.mean <- sum(true.mean[-r]*weights)
        diffs.fold.centered <- diffs.fold - rep(true.mean[r] - Q.true.mean, nrow(in.fold))

      } else if (min.algor == 'argmin'){
        # try argmin
        idx.min <- find.sub.argmin(mu.out.fold, r)
        diffs.fold <- in.fold[,r] - in.fold[,idx.min]

        ## output centered test statistic
        diffs.fold.centered <- diffs.fold - rep(true.mean[r] - true.mean[idx.min], nrow(in.folds))
      } else {
        # error
        stop("'min.algor' should be either 'softmin' or 'argmin'")
      }
      #return (diffs.fold)
      return (data.frame(diffs.fold=diffs.fold, diffs.fold.centered=diffs.fold.centered))
    })
    #diffs <- do.call(c, diffs)
    res <- do.call(rbind,diffs)
    diffs <- res[,1]
    sigma <- stats::sd(diffs)

    test.stat <- sqrt(n)*mean(diffs)
    test.stat.scale <- test.stat/sigma

    ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')

    test.stat.centered <- NULL
    if (!is.null(true.mean)) {
      diffs.centered <- res[,2]
      test.stat.centered <- test.stat.scale - sqrt(n)*mean(diffs.centered)/sigma
    }
    return (list(test.stat.scale=test.stat.scale, critical.value=val.critical, std=sigma, ans=ans,
                 lambda=lambda, test.stat.centered=test.stat.centered))
  }
}

#' tune the omega used for the universal inference method by \insertCite{dey.2024}{argminCS}.
#'
#' tune the omega used for the universal inference method by \insertCite{dey.2024}{argminCS} by bootstrapping.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param B The number of bootstrap iterations; defaults to 100.
#' @param omegas A vector including the omegas to perform bootstrapping over; defaults to 1:100.
#'
#' @return An omega value for hypothesis testing, using the universal inference algorithm.
#' @export
#'
#' @examples
#' n <- 200
#' mu <- (1:20)/20
#' cov <- diag(length(mu))
#' set.seed(38)
#' data <- MASS::mvrnorm(n, mu, cov)
#' omega <- omega.bootstrap(data)
#' omega
#'
#' @references{
#'   \insertRef{dey.2024}{argminCS}
#' }
omega.bootstrap <- function(data, alpha=0.05, B=100, omegas=1:100){
  # data: n by p
  # alpha: a desired significance level

  n <- nrow(data)
  p <- ncol(data)

  sample.mean <- colMeans(data)
  min.indices <- which(sample.mean == min(sample.mean))
  seed.argmin.all <- ceiling(23*data[n,p]*p)
  if (seed.argmin.all > (2^31-1)){
    seed.argmin.all <- seed.argmin.all %%  (2^31-1)
  }
  withr::with_seed(seed.argmin.all, {
    idx.min <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
  })

  coverages <- lapply(1:B, function(i){

    seed <- ceiling(i*seed.argmin.all + i)
    if (seed > 2^31 - 1){
      seed <- seed %% (2^31 - 1)
    }
    withr::with_seed(seed, {
      indices <- sample(1:n, n, replace=T)
    })
    data.boot <- data[indices,]

    # split the data
    seed.splitting <- seed*11+i
    if (seed.splitting > 2^31 - 1){
      seed.splitting <- seed.splitting %% (2^31 - 1)
    }
    withr::with_seed(seed.splitting, {
      indices.boot.training <- sample(1:n, floor(n/2), replace=FALSE)
    })
    data.boot.training <- data.boot[indices.boot.training,]
    data.boot.testing <- data.boot[-indices.boot.training,]

    # get the best dimension over the bootstrapped training set
    sample.mean.boot.training <- colMeans(data.boot.training)
    min.indices.boot <- which(sample.mean.boot.training == min(sample.mean.boot.training))
    seed.argmin <- 3*seed+i
    if (seed.argmin >  2^31-1){
      seed.argmin <- seed.argmin %%  2^31-1
    }
    withr::with_seed(seed.argmin, {
      idx.min.boot.training <- ifelse((length(min.indices.boot) > 1),
                                      sample(c(min.indices.boot), 1),
                                      min.indices.boot[1])
    })

    estimated.minimum.mean.boot <- mean(data.boot.testing[,idx.min.boot.training])

    # give any idx (it's simply a dummy value)
    res <- sapply(omegas, function(omega){
      argmin.HT.GU(data, 1, omega, estimated.minimum.mean=estimated.minimum.mean.boot,
                   mean.r=mean(data.boot.testing[,idx.min]), alpha=alpha)$ans} )

    # reward the omegas that fail to reject idx.min
    res <- ifelse(res=='Reject', 0, 1)
    return (res)
  })

  coverages <- as.matrix(do.call('rbind', coverages))
  coverages <- colSums(coverages)

  distances.from.nominals <- abs(coverages - (1 - alpha)*B)
  omega.indices.candidates <- which(distances.from.nominals == min(distances.from.nominals))
  seed.final <- floor(seed.argmin.all*data[n,p])
  if (seed.final > 2^31 - 1){
    seed.final <- seed.final %% (2^31 - 1)
  }
  withr::with_seed(seed.final, {
    omega.index <- ifelse((length(omega.indices.candidates) > 1),
                          sample(c(omega.indices.candidates), 1),
                          omega.indices.candidates[1])
  })
  return (omegas[omega.index])
}

#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, employing the universal inference framework by \insertCite{dey.2024}{argminCS}.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param omega The tuning parameter in the statistic; defaults to NULL. Recommend to compute it using \link{omega.bootstrap}.
#' @param estimated.minimum.mean Upon a splitting over data, the mean of r'-th dimension over the testing set, where r' is the argmin over the training set; defaults to NULL.
#' Highly recommend to provide a value to improve the computation speed.
#' @param mean.r Upon a splitting over data, the mean of the r-th dimension over the testing set; defaults to NULL.
#' Highly recommend to provide a value to improve the computation speed.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param ... Additional arguments to \link{omega.bootstrap}.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.statistic} \tab The off-line test statistic in \insertCite{dey.2024}{argminCS}. \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being greater than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
#'
#' @references{
#'   \insertRef{dey.2024}{argminCS}
#' }
argmin.HT.GU <- function(data, r, omega=NULL, estimated.minimum.mean=NULL, mean.r=NULL, alpha=0.05, ...){

  n <- nrow(data)
  p <- ncol(data)
  n.testing <- ceiling(n/2)

  if (is.null(omega)){
    omega <- omega.bootstrap(data, alpha=alpha, ...)
  }

  # hypothesis testing preparation
  if (is.null(estimated.minimum.mean) | is.null(mean.r)){
    seed <- ceiling(data[n,p]*n*p + p)
    if (seed > 2^31 - 1){
      seed <- seed %% 2^31 -1
    }
    withr::with_seed(seed, {
      indices.training <- sample(1:n, floor(n/2), replace=FALSE)
      })

    data.training <- data[indices.training,]
    data.testing  <- data[-indices.training,]

    # get the best dimension over training set
    sample.mean.training <- colMeans(data.training)
    min.indices <- which(sample.mean.training == min(sample.mean.training))
    seed.argmin <- 3*seed
    if (seed.argmin >  (2^31-1)){
      seed.argmin <- seed.argmin %% (2^31-1)
    }
    withr::with_seed(seed.argmin, {
      idx.min.training <- ifelse((length(min.indices) > 1),
                                 sample(c(min.indices), 1),
                                 min.indices[1])
    })

    estimated.minimum.mean <- mean(data.testing[,idx.min.training])
    mean.r <- mean(data.testing[,r])
  }

  # perform hypothesis testing
  test.statistic <- exp(-omega*n.testing*(estimated.minimum.mean - mean.r))
  critical.value <- 1/alpha
  ans <- ifelse(test.statistic >= critical.value, 'Reject', 'Accept')
  return (list(test.statistic=test.statistic, critical.value=critical.value, ans=ans))
}

#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, using multiple testing with Bonferroni's correction.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param r.min The sample argmin of the sample mean of data; defaults to NULL.
#' @param r.min.sec The index of the second smallest dimension of the sample mean of the data; defaults to NULL.
#' @param test The test to perform: 't' or 'z' test; defaults to 'z'.
#' If data are believed to be normally distributed, use 't'; otherwise 'z'.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{p.val} \tab p value without Bonferroni's correction. \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being less than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
argmin.HT.MT <- function(data, r, test='z', r.min=NULL, r.min.sec=NULL, alpha=0.05){

  p <- ncol(data)
  val.critical <- alpha/(p-1)

  if (is.null(r.min) | is.null(r.min.sec)){
    sample.mean <- colMeans(data) # np
    min.indices <- which(sample.mean == min(sample.mean))
    seed <- data[1,r]*r*37
    if (seed >  2^31 - 1){
      seed <- seed %%  2^31 - 1
    }
    withr::with_seed(seed, {
      r.min <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
    })
    r.min.sec <- find.sub.argmin(sample.mean, r.min)
  }

  p.val <- NULL
  if (test == 't'){
    # t test
    if (r == r.min){
      p.val <- stats::t.test(data[,r]-data[,r.min.sec], alternative='greater')$p.value
    } else {
      p.val <- stats::t.test(data[,r]-data[,r.min], alternative='greater')$p.value
    }
  } else{
    # z.test
    if (r == r.min){
      diffs <- data[,r] - data[,r.min.sec]
      p.val <- BSDA::z.test(diffs, sigma.x=stats::sd(diffs), alternative='greater')$p.value
    } else {
      diffs <- data[,r] - data[,r.min]
      p.val <- BSDA::z.test(diffs, sigma.x=stats::sd(diffs), alternative='greater')$p.value
    }
  }

  ans <- ifelse(p.val > val.critical, 'Accept', 'Reject')
  return (list(p.val=p.val, critical.value=val.critical, ans=ans))
}

#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, adapting the one-step self-normalization method in \insertCite{cck.many.moments}{argminCS}.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, compute colMeans(data) and pass it to sample.mean to speed up computation.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat} \tab The test statistic \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being greater than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
#'
#' @importFrom Rdpack reprompt
#' @references{
#'  \insertRef{cck.many.moments}{argminCS}
#' }
argmin.HT.SN <- function(data, r, sample.mean=NULL, alpha=0.05){

  n <- nrow(data)
  p <- ncol(data)

  norm.quantile <- stats::qnorm(1 - alpha/(p-1))
  val.critical <- ifelse(norm.quantile^2 >= n,Inf,norm.quantile/sqrt(1-norm.quantile^2/n))

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  diffs <- matrix(rep(data[,r], p-1), nrow=n, byrow=F) - data[,-r] # np
  sd.diffs <- apply(diffs, 2, stats::sd) # 3np
  mean.diffs <- NULL
  if (is.null(sample.mean)){
    # even without this step, the algorithm still requires np operations
    mean.diffs <- colMeans(diffs)
  } else {
    mean.diffs <- sample.mean[r] - sample.mean[-r] # p
  }
  ratios <- mean.diffs/sd.diffs
  test.stat <- sqrt(n)*max(ratios)
  ans <- ifelse(test.stat <= val.critical, 'Accept', 'Reject')
  return (list(test.stat=test.stat, critical.value=val.critical, ans=ans))
}

#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, adapting the one-step bootstrap method in \insertCite{cck.many.moments}{argminCS}.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, compute colMeans(data) and pass it to sample.mean to speed up computation.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param B The number of bootstrap samples; defaults to 200.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{p.val} \tab p value of the test \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
#'
#' @importFrom Rdpack reprompt
#' @references{
#'  \insertRef{cck.many.moments}{argminCS}
#'
#'  \insertRef{lei.cvc}{argminCS}
#' }
argmin.HT.bootstrap <- function(data, r, sample.mean=NULL, alpha=0.05, B=200){

  n <- nrow(data)
  p <- ncol(data)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  diffs <- matrix(rep(data[,r], p-1), nrow=n, byrow=F) - data[,-r] # np
  sd.diffs <- apply(diffs, 2, stats::sd) # 3np
  mean.diffs <- NULL
  if (is.null(sample.mean)){
    # even without this step, the algorithm still requires np operations
    mean.diffs <- colMeans(diffs)
  } else {
    mean.diffs <- sample.mean[r] - sample.mean[-r] # p
  }
  ratios <- mean.diffs/sd.diffs
  test.stat <- sqrt(n)*max(ratios)

  # n by p matrix
  diffs.centered <- diffs - matrix(rep(mean.diffs, n), nrow=n, byrow=T)
  test.stat.MBs <- sapply(1:B,
                          function(i){
                            seed <- ceiling(abs(53*i*r*data[1,r]*sample.mean[r]))+i
                            if (seed >  2^31 - 1){
                              seed <- seed %%  2^31 - 1
                            }
                            withr::with_seed(seed, {
                              Gaussian.vec <- stats::rnorm(n, 0, 1)
                             })
                            test.stat.MB <- sqrt(n)*max(colMeans(diffs.centered*Gaussian.vec)/sd.diffs)
                            return (test.stat.MB)
                            })

  p.val <- mean(test.stat.MBs > test.stat)
  ans <- ifelse(p.val > alpha, 'Accept', 'Reject')
  return (list(p.val=p.val, test.stat=test.stat, ans=ans))
}

#' @title Generate the quantile used for the selection procedure in \insertCite{gupta.1965}{argminCS}.
#'
#' @description Generate the quantile used for the selection procedure in \insertCite{gupta.1965}{argminCS} by Monte Carlo estimation.
#'
#' @note The quantile is pre-calculated for some common configurations of (p, alpha)
#'
#' @param p The number of dimensions in your data matrix.
#' @param alpha The level of the upper quantile; defaults to 0.05 (95\% percentile).
#' @param N The number of Monte Carlo repetitions; defaults to 100000.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{critica.val} \tab The 1 - alpha upper quantile. \cr
#' }
#' @export
#'
#' @examples
#' get.quantile.gupta.selection(p=10)
#'
#' get.quantile.gupta.selection(p=100)
#'
#' @importFrom Rdpack reprompt
#' @references{
#'  \insertRef{gupta.1965}{argminCS}
#'
#'  \insertRef{futschik.1995}{argminCS}
#' }
get.quantile.gupta.selection <- function(p, alpha=0.05, N=100000){
  quantile <- quantiles.gupta[quantiles.gupta$p == p, as.character(alpha)]

  if (p == 1){
    return (stats::qnorm(1-alpha))

  }
  else if (p == 2){
    return (stats::qnorm(1-alpha, 0, sqrt(2)))

  } else {
    if (length(quantile) == 1){
      return (quantile)
    } else{
      mult.normals <- MASS::mvrnorm(n=N, mu=rep(0, p), Sigma=diag(p))
      diffs.with.max <- apply(mult.normals, 1, function(row) {max(row[-p]) - row[p]})
      return (stats::quantile(diffs.with.max, 1-alpha))
    }
  }
}

#' @title Perform argmin hypothesis test.
#'
#' @description Test if a dimension may be argmin, using the method in \insertCite{gupta.1965}{argminCS}.
#'
#' @note This method requires independence among the dimensions.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, compute colMeans(data) and pass it to sample.mean to speed up computation.
#' @param stds The equal dimension-wise ( population) standard deviations of the n samples in data; defaults a vector of 1's.
#' @param critical.val The quantile for the hypothesis test; defaults to NULL. It can be calculated via \link{get.quantile.gupta.selection}.
#' If your experiment involves hypothesis testing over more than one dimension, pass a quantile to speed up computation.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param ... Additional argument to \link{get.quantile.gupta.selection}.
#' A correct argument name needs to be specified if it is used.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat} \tab The test statistic \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being greater than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
#'
#' @importFrom Rdpack reprompt
#' @references{
#'  \insertRef{gupta.1965}{argminCS}
#'
#'  \insertRef{futschik.1995}{argminCS}
#' }
argmin.HT.gupta <- function(data, r, sample.mean=NULL, stds=NULL, critical.val=NULL, alpha=0.05, ...){

  # note that the implementation can, in turn, asks for scaled.sample.mean, its min and second min
  # to speed up the computation, in case that there are many dimensions involving.
  n <- nrow(data)
  p <- ncol(data)

  if (is.null(critical.val)){
    critical.val <- get.quantile.gupta.selection(p=p, alpha=alpha, ...)
  }

  scaled.sample.mean <- sample.mean
  if (is.null(scaled.sample.mean)){
    scaled.sample.mean <- colMeans(data)
  }
  if (is.null(stds)){
    stds <- rep(1, p)
  }
  scaled.sample.mean <- scaled.sample.mean/stds

  # Gupta's paper sets up the selection rule to find argmax
  # we thus have to modify the 'data' by -1 if intended to follow their procedure closely
  # note that the above sample mean is calculated without a multiplication of -1
  # the list of standard deviations would not be affected by a multiplication of -1
  # Gupta's test stat: sqrt(n)*(max(-scaled.sample.mean[-r]) - (-scaled.sample.mean[r]))
  # it is equal to the following
  test.stat <- sqrt(n)*(scaled.sample.mean[r] - min(scaled.sample.mean[-r]))
  ans <- ifelse(test.stat <= critical.val, 'Accept', 'Reject')
  return (list(test.stat=test.stat, critical.val=critical.val, ans=ans))
}
