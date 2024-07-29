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
#'    \code{fold (FD)} \tab A n fold version of SML; defaults to two folds. \cr
#'    \tab \cr
#'    \code{GU} \tab The method in \insertCite{dey.2024}{argminCS}. \cr
#'    \tab \cr
#'    \code{CCK.self.normalization (SN)} \tab Modified from the self-normalization method in \insertCite{cck.many.moments}{argminCS}. \cr
#'    \tab \cr
#'    \code{CCK.bootstrap (CB)} \tab Modified from the bootstrap method in \insertCite{cck.many.moments}{argminCS}. See also \insertCite{lei.cvc}{argminCS}. \cr
#'    \tab \cr
#'    \code{Bonferroni (MT)} \tab Multiple testing with Bonferroni's correction. \cr
#'    \tab \cr
#'    \code{Gupta (GTA)} \tab The method in \insertCite{gupta.1965}{argminCS}. \cr
#'    \tab \cr
#'    \code{Futschik (FCHK)} \tab The two-step method in \insertCite{futschik.1995}{argminCS} \cr
#' }
#' If computation is a concern, use 'MT'. Otherwise, 'SML' is recommended.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param method A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
#' For the list of supported methods and their abbreviations, see Details.
#' @param alpha The significance level; defaults to 0.05. The function produces a (1-alpha) confidence set.
#' @param ... Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}, \link{argmin.HT.fold}, \link{is.lambda.feasible.fold}, \link{argmin.HT.GU}, \link{argmin.HT.SN}, \link{argmin.HT.bootstrap}, \link{argmin.HT.MT}, \link{argmin.HT.gupta}.
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
#' CS.argmin(data, method='FD')
#' ## 5 fold
#' CS.argmin(data, method='FD', n.fold=5)
#'
#' ## GU
#' CS.argmin(data, method='GU')
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
#' ## Gupta
#' CS.argmin(data, method='GTA')
#'
#' ## Futschik two-step method
#' # default alpha.1, alpha.2
#' CS.argmin(data, method='FCHK')
#'
#' alpha.1 <- 0.0005
#' alpha.2 <- 1 - (0.95/(1 - alpha.1))
#' CS.argmin(data, method='FCHK', alpha.1=0.0005, alpha.2=alpha.2)
#'
#' @importFrom Rdpack reprompt
#' @references{
#'   \insertRef{cck.many.moments}{argminCS}
#'
#'   \insertRef{lei.cvc}{argminCS}
#'
#'   \insertRef{dey.2024}{argminCS}
#'
#'  \insertRef{gupta.1965}{argminCS}
#'
#'  \insertRef{futschik.1995}{argminCS}
#' }
#'
CS.argmin <- function(data, method='softmin.LOO', alpha=0.05, ...){

  n <- nrow(data)
  p <- ncol(data)

  if (method == 'softmin.LOO' | method == 'SML'){
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.LOO(data, r, sample.mean=sample.mean, alpha=alpha, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'argmin.LOO' | method == 'HML') {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.LOO(
      data, r, min.algor=getMin.argmin.LOO, sample.mean=sample.mean, alpha=alpha, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'nonsplit' | method == 'NS') {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.nonsplit(data, r, sample.mean=sample.mean, alpha=alpha, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'fold' | method == 'FD') {
    additional.arguments <- list(...)
    if (methods::hasArg(n.fold)){
      n.fold <- additional.arguments$n.fold
    } else{
      n.fold <- 2 # defaults
    }
    seed <- ceiling(abs(13*data[n,p]*n.fold)) + p
    if (seed >  2^31 - 1){
      seed <- seed %%  2^31 - 1
    }
    withr::with_seed(seed, {
      flds <- caret::createFolds(1:n, k=n.fold, list=TRUE, returnTrain=FALSE)
    })
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.fold(data, r, alpha=alpha, flds=flds, sample.mean=sample.mean,...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'GU') {
    omega.provided <- methods::hasArg(omega)
    if (!omega.provided){
      omega <- omega.bootstrap(data)
    }
    seed <- floor(data[n,p]*p + p)
    withr::with_seed(seed, {
      indices.training <- sample(n, floor(n/2), replace=F)
    })
    data.training <- data[indices.training,]
    sample.mean.training <- colMeans(data.training)
    min.indices <- which(sample.mean.training == min(sample.mean.training))
    seed.argmin <- 3*seed
    if (seed.argmin >  (2^31-1)){
      seed.argmin <- seed.argmin %% (2^31-1)
    }
    withr::with_seed(seed.argmin, {
      idx.min.training <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
    })

    data.testing <- data[-indices.training,]
    sample.mean.testing <- colMeans(data.testing)

    res <- NULL
    if (omega.provided){
      res <- sapply(1:p, function(r) {
      argmin.HT.GU(data, r,
                   estimated.minimum.mean=sample.mean.testing[idx.min.training],
                   mean.r=sample.mean.testing[r], alpha=alpha, ...)$ans})
    } else {
      res <- sapply(1:p, function(r) {
        argmin.HT.GU(data, r, omega = omega,
                     estimated.minimum.mean=sample.mean.testing[idx.min.training],
                     mean.r=sample.mean.testing[r], alpha=alpha, ...)$ans})
    }
    return (which(res == 'Accept'))
  } else if (method == 'CCK.self.normalization' | method == 'SN'){
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.SN(data, r, sample.mean=sample.mean, alpha=alpha, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'CCK.bootstrap' | method == 'CB') {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {argmin.HT.bootstrap(data, r, sample.mean=sample.mean, alpha=alpha, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'Bonferroni' | method == 'MT') {
    sample.mean <- colMeans(data) # np
    min.indices <- which(sample.mean == min(sample.mean))
    seed <- ceiling(107*sample.mean[1]*data[n,p])
    if (seed >  2^31-1){
      seed <- seed %%  2^31-1
    }
    withr::with_seed(seed, {
      r.min <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
    })
    r.min.sec <- find.sub.argmin(sample.mean, r.min)
    res <- sapply(1:p, function(r) {argmin.HT.MT(data, r, r.min=r.min, r.min.sec=r.min.sec, alpha=alpha, ...)$ans})
    return (which(res == 'Accept'))

  } else if (method == 'Gupta' | method == 'gupta' | method == 'GTA') {
    p <- ncol(data)
    critical.val <- get.quantile.gupta.selection(p=p, alpha=alpha)
    sample.mean <- colMeans(data)
    stds <- NULL
    if (methods::hasArg(std)){
      additional.arguments <- list(...)
      std <- additional.arguments$std
      stds <- rep(std, p)
    ## testing
    # if (methods::hasArg(stds)){
    #   additional.arguments <- list(...)
    #   stds <- additional.arguments$stds
    #   res <- sapply(1:p, function(r) {argmin.HT.gupta(
    #     data, r, critical.val=critical.val, sample.mean=sample.mean, alpha=alpha, ...)$ans})
    } else {
      stds <- rep(1, p)
      #stds <- apply(data, 2, stats::sd)
    }
    res <- sapply(1:p, function(r) {argmin.HT.gupta(
      data, r, critical.val=critical.val, sample.mean=sample.mean, stds=stds, alpha=alpha, ...)$ans})

    return (which(res == 'Accept'))

  } else if (method == 'futschik' | method == 'Futschik' | method == 'FCHK'){

    additional.arguments <- list(...)
    if (methods::hasArg(alpha.1) & methods::hasArg(alpha.2)){
      alpha.1 <- additional.arguments$alpha.1
      alpha.2 <- additional.arguments$alpha.2
      if (alpha.1 < 0 | alpha.1 > 1 | alpha.2 < 0 | alpha.2 > 1){
        stop('invalid configurations of (alpha.1, alpha.2)')
      }
    } else{
      alpha.1 <- alpha/10
      alpha.2 <- 1 - (1 - alpha)/(1 - alpha.1)
    }

    # # comment out to show the violation of validity
    # if (methods::hasArg(std)){
    #   std <- additional.arguments$std
    #   stds <- rep(std, p)
    # } else {
    #   stds <- rep(1, p)
    # }
    stds <- apply(data, 2, stats::sd)

    sample.mean <- colMeans(data)
    # step 1
    p.1 <- ncol(data)
    critical.val.1 <- get.quantile.gupta.selection(p=p.1, alpha=alpha.1)
    res.1 <- sapply(1:p, function(r) {argmin.HT.gupta(
      data, r, critical.val=critical.val.1, sample.mean=sample.mean, stds=stds, alpha=alpha.1, ...)$ans})
    res.1 <- unname(res.1)

    # step 2
    p.2 <- sum(res.1 == 'Accept') - 1
    if (p.2 == 0){
      return (which(res.1 == 'Accept'))
    } else {
      critical.val.2 <- get.quantile.gupta.selection(p=p.2, alpha=alpha.2)
      res.2 <- sapply(1:p, function(r) {argmin.HT.gupta(
        data, r, critical.val=critical.val.2, sample.mean=sample.mean, stds=stds, alpha=alpha.2, ...)$ans})
      res.2 <- unname(res.2)

      return (which(res.1 == 'Accept' & res.2 == 'Accept'))
    }

  } else {
    stop("'method' should be one of 'softmin.LOO' (SML), 'argmin.LOO' (HML),
         'nonsplit' (NS), 'fold' (FD), 'GU', 'CCK.bootstrap' (CB), 'Bonferroni' (MT), Gupta (GTA), 'Futschik (FCHK)")
  }
}

