#' Construct a discrete confidence set for argmin.
#'
#' This is a wrapper to construct a discrete confidence set for argmin. Multiple methods are supported.
#'
#' @importFrom Rdpack reprompt
#' @details The supported methods include:\tabular{ll}{
#'   \code{softmin.LOO (SML)} \tab Leave-one-out algorithm using exponential weighting. \cr
#'   \tab \cr
#'   \code{argmin.LOO (HML)} \tab A variant of SML that uses hard argmin instead of exponential weighting. Not recommended. \cr
#'   \tab \cr
#'   \code{nonsplit (NS)} \tab A variant of SML without data splitting. Requires a fixed lambda value as an additional argument. Not recommended\cr
#'   \tab \cr
#'   \code{Bonferroni (MT)} \tab Multiple testing using Bonferroni correction. \cr
#'   \tab \cr
#'   \code{Gupta (GTA)} \tab The method proposed by \insertCite{gupta.1965;textual}{argminCS}.
#'   Requires independence and the same population standard deviation for all dimensions. \cr
#'   \tab \cr
#'   \code{Futschik (FCHK)} \tab A two-step method from \insertCite{futschik.1995;textual}{argminCS}.
#'   Requires independence and the same population standard deviation for all dimensions. \cr
#' }
#'
#' @param data A n by p data matrix; each row is a p-dimensional sample.
#' @param method A string indicating the method used to construct the confidence set. Defaults to 'softmin.LOO'.
#' Can be abbreviated (e.g., 'SML' for 'softmin.LOO'). See **Details** for available methods and abbreviations.
#' @param alpha The significance level; defaults to 0.05. The function produces a \eqn{1 - \alpha} confidence set.
#' @param ... Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}, \link{argmin.HT.MT}, \link{argmin.HT.gupta}.
#' A correct argument name needs to be specified if it is used.
#'
#' @return A vector of indices (1-based) representing the (1 - alpha) confidence set.
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
#' ## use seed
#' CS.argmin(data, seed=13)
#'
#' ## argmin.LOO
#' CS.argmin(data, method='HML')
#'
#' ## nonsplit
#' CS.argmin(data, method='NS', lambda=sqrt(n)/2.5)
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
#'  \insertRef{gupta.1965}{argminCS}
#'
#'  \insertRef{futschik.1995}{argminCS}
#' }
#'
CS.argmin <- function(data, method='softmin.LOO', alpha=0.05, ...){

  n <- nrow(data)
  p <- ncol(data)

  if (method %in% c('softmin.LOO', 'SML', 'sml')){
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {
      difference.matrix <- get.difference.matrix(data, r)
      sample.mean.r <- get.sample.mean.r(sample.mean, r)
      return (argmin.HT.LOO(difference.matrix, sample.mean=sample.mean.r, alpha=alpha, ...)$ans)
      })
    return (which(res == 'Accept'))

  } else if (method %in% c('argmin.LOO', 'HML', 'hml')) {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {
      difference.matrix <- get.difference.matrix(data, r)
      sample.mean.r <- get.sample.mean.r(sample.mean, r)
      return (argmin.HT.LOO(difference.matrix, sample.mean=sample.mean.r, min.algor='argmin', alpha=alpha, ...)$ans)
    })
    return (which(res == 'Accept'))

  } else if (method %in% c('nonsplit', 'NS', 'ns')) {
    sample.mean <- colMeans(data)
    res <- sapply(1:p, function(r) {
      difference.matrix <- get.difference.matrix(data, r)
      sample.mean.r <- get.sample.mean.r(sample.mean, r)
      return (argmin.HT.nonsplit(difference.matrix, sample.mean=sample.mean.r, alpha=alpha, ...)$ans)
      })
    return (which(res == 'Accept'))

  } else if (method %in% c('Bonferroni', 'bonferroni', 'MT', 'mt')) {
    sample.mean <- colMeans(data) # np
    res <- sapply(1:p, function(r) {
      difference.matrix <- get.difference.matrix(data, r)
      sample.mean.r <- get.sample.mean.r(sample.mean, r)
      return (argmin.HT.MT(difference.matrix, sample.mean=sample.mean.r, alpha=alpha, ...)$ans)
      })
    return (which(res == 'Accept'))

  } else if (method %in% c('Gupta', 'gupta', 'GTA', 'gta')) {
    p <- ncol(data)
    critical.val <- get.quantile.gupta.selection(p=p, alpha=alpha)
    sample.mean <- colMeans(data)
    stds <- NULL
    if (methods::hasArg(std)){
      additional.arguments <- list(...)
      std <- additional.arguments$std
      stds <- rep(std, p)

    } else {
      stds <- rep(1, p)
      # the method does not support the use of sample standard deviations
      # stds <- apply(data, 2, stats::sd)
    }
    res <- sapply(1:p, function(r) {argmin.HT.gupta(
      data, r, critical.val=critical.val, sample.mean=sample.mean, stds=stds, alpha=alpha, ...)$ans})

    return (which(res == 'Accept'))

  } else if (method %in% c('futschik', 'FCHK', 'fchk')) {

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
    if (methods::hasArg(std)){
      std <- additional.arguments$std
      stds <- rep(std, p)
    } else {
      stds <- rep(1, p)
      ## the method does not support the use of sample standard deviations
      # stds <- apply(data, 2, stats::sd)
    }

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
       'nonsplit' (NS), 'Bonferroni' (MT), 'Gupta' (GTA), 'Futschik' (FCHK)")
  }
}

#' Construct a discrete confidence set for argmax.
#'
#' This is a wrapper to construct a confidence set for the argmax by negating the input and reusing \code{\link{CS.argmin}}.
#'
#' @importFrom Rdpack reprompt
#'
#' @details The supported methods include:\tabular{ll}{
#'   \code{softmin.LOO (SML)} \tab Leave-one-out algorithm using exponential weighting. \cr
#'   \code{argmin.LOO (HML)} \tab Variant of SML that uses hard argmin instead of soft weighting. Not recommended. \cr
#'   \code{nonsplit (NS)} \tab Variant of SML without data splitting. Requires a fixed lambda value. Not recommended. \cr
#'   \code{Bonferroni (MT)} \tab Multiple testing using Bonferroni correction. \cr
#'   \code{Gupta (GTA)} \tab The method of \insertRef{gupta.1965}{argminCS}. \cr
#'   \code{Futschik (FCHK)} \tab A two-step method from \insertRef{futschik.1995}{argminCS}. \cr
#' }
#'
#' @param data An \eqn{n \times p} matrix; each row is a p-dimensional sample.
#' @param method A string indicating the method to use; defaults to 'softmin.LOO'.
#' Can be abbreviated (e.g., 'SML' for 'softmin.LOO'). See Details for full list.
#' @param alpha Significance level. The function returns a \eqn{1 - \alpha} confidence set.
#' @param ... Additional arguments passed to corresponding testing functions.
#'
#' @return A vector of indices (1-based) representing the confidence set for the argmax.
#' @export
#'
#' @examples
#' set.seed(108)
#' n <- 200
#' p <- 20
#' mu <- (1:p)/p
#' cov <- diag(p)
#' data <- MASS::mvrnorm(n, mu, cov)
#'
#' ## softmin.LOO (SML)
#' CS.argmax(data)
#'
#' ## argmin.LOO (HML)
#' CS.argmax(data, method = "HML")
#'
#' ## nonsplit (NS) - requires lambda
#' CS.argmax(data, method = "NS", lambda = sqrt(n)/2.5)
#'
#' ## Bonferroni (MT) - t test default
#' CS.argmax(data, method = "MT", test = "t")
#'
#' ## Gupta (GTA)
#' CS.argmax(data, method = "GTA")
#'
#' ## Futschik (FCHK) with default alpha.1 and alpha.2
#' CS.argmax(data, method = "FCHK")
#'
#' ## Futschik (FCHK) with user-specified alpha.1 and alpha.2
#' alpha.1 <- 0.001
#' alpha.2 <- 1 - (0.95 / (1 - alpha.1))
#' CS.argmax(data, method = "FCHK", alpha.1 = alpha.1, alpha.2 = alpha.2)
#'
#' @references
#' \insertRef{gupta.1965}{argminCS}
#'
#' \insertRef{futschik.1995}{argminCS}
#'
#' \insertRef{cck.many.moments}{argminCS}
CS.argmax <- function(data, method = "softmin.LOO", alpha = 0.05, ...) {
  negated.data <- -data
  return(CS.argmin(negated.data, method = method, alpha = alpha, ...))
}

#' Construct a difference matrix for argmin hypothesis testing
#'
#' Given a data matrix and a reference column index, construct the difference matrix
#' used in hypothesis testing procedures. Each column represents the difference
#' between the reference dimension and one of the remaining dimensions.
#'
#' @param data A \code{n} by \code{p} data matrix; each row is a \code{p}-dimensional sample.
#' @param r An integer between 1 and \code{p}, indicating the reference column (dimension).
#'
#' @return A \code{n} by \code{(p-1)} matrix where each row is the difference between the \code{r}-th column and the remaining columns.
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(50), nrow = 10)
#' diff.mat <- get.difference.matrix(data, r = 2)
#'
#' @export
#'
#' @keywords internal
get.difference.matrix <- function(data, r) {
  difference.matrix <- matrix(rep(data[, r], ncol(data) - 1), ncol = ncol(data) - 1, byrow = FALSE) - data[, -r]
  return (difference.matrix)
}

#' Compute sample mean differences for hypothesis testing
#'
#' Computes the vector of differences between the sample mean at the reference index
#' and the remaining dimensions.
#'
#' @param sample.mean A vector of length \code{p} containing the sample means of each dimension.
#' @param r An integer between 1 and \code{p}, indicating the reference dimension.
#'
#' @return A vector of length \code{p - 1} giving the differences: sample.mean[r] - sample.mean[-r].
#'
#' @examples
#' sample.mean <- 1:5
#' get.sample.mean.r(sample.mean, r = 3)
#'
#' @export
#'
#' @keywords internal
get.sample.mean.r <- function(sample.mean, r) {
  return (sample.mean[r] - sample.mean[-r])
}



