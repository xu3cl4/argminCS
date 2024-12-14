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
#'    \code{Bonferroni (MT)} \tab Multiple testing with Bonferroni's correction. \cr
#'    \tab \cr
#'    \code{Gupta (GTA)} \tab The method in \insertRef{gupta.1965}{argminCS}. \cr
#' }
#'
#' @param data (1) A n by p data matrix for (GTA); each of its row is a p-dimensional sample, or
#' (2) A n by (p-1) difference matrix for (SML, HML, NS, MT); each of its row is a (p-1)-dimensional sample differences
#' @param r The dimension of interest for hypothesis test.
#' @param method A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
#' For the list of supported methods and their abbreviations, see Details.
#' @param ... Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}, \link{argmin.HT.MT}, \link{argmin.HT.gupta}.
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
#' difference.matrix.r <- matrix(rep(data[,r], p-1), ncol=p-1, byrow=FALSE) - data[,-r]
#' argmin.HT(difference.matrix, r)
#' # provide centered test statistic (to simulate asymptotic normality)
#' true.mean.difference.r <- mu[r] - mu[-r]
#' argmin.HT(difference.matrix, r, true.mean=true.mean.difference.r)
#'
#' ## argmin.LOO
#' argmin.HT(difference.matrix, r, method='HML')
#'
#' ## nonsplit
#' argmin.HT(difference.matrix, r, method='NS', lambda=sqrt(n)/2.5)
#'
#' ## Bonferroni (choose t test because of normal data)
#' argmin.HT(difference.matrix, r, method='MT', test='t')
#'
#' ## Gupta
#' critical.val <- get.quantile.gupta.selection(p=length(mu))
#' argmin.HT(data, r, method='GTA', critical.val=critical.val)
#' @importFrom Rdpack reprompt
#' @references{
#'   \insertRef{cck.many.moments}{argminCS}
#'
#'   \insertRef{gupta.1965}{argminCS}
#'
#'   \insertRef{futschik.1995}{argminCS}
#' }
argmin.HT <- function(data, r, method='softmin.LOO', ...){
  if (method == 'softmin.LOO' | method == 'SML'){
    return (argmin.HT.LOO(difference.matrix, ...))

  } else if (method == 'argmin.LOO' | method == 'HML') {
    return (argmin.HT.LOO(difference.matrix, min.algor=getMin.argmin.LOO, ...))

  } else if (method == 'nonsplit' | method == 'NS') {
    return (argmin.HT.nonsplit(difference.matrix, ...))

  } else if (method == 'Bonferroni' | method == 'MT') {
    return (argmin.HT.MT(difference.matrix, ...))

  } else if (method == 'Gupta' | method == 'GTA' | method=='gupta') {
    return (argmin.HT.gupta(data, r, ...))

  } else {
    stop("'method' should be one of 'softmin.LOO' (SML), 'argmin.LOO' (HML),
         'nonsplit' (NS), 'Bonferroni' (MT), 'Gupta' (GTA)")
  }
}

#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, using the LOO (leave-one-out) algorithm in Zhang et al 2024.
#'
#' @param difference.matrix A n by (p-1) difference data matrix (reference dimension - the rest);
#' each of its row is a (p-1)-dimensional vector of differences.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(difference.matrix).
#' If your experiment involves hypothesis testing over more than one dimension, compute colMeans(data) and pass it to sample.mean to speed up computation.
#' @param min.algor The algorithm to find the minimum excluding the r-th dimension; defaults to 'softmin'. The other option is 'argmin'.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin); defaults to NULL.
#' If lambda=NULL (recommended), the function would determine a lambda value in a data-driven way.
#' @inheritParams lambda.adaptive.enlarge
#' @param enlarge A boolean value indicating if the data-driven lambda should be determined via an iterative enlarging algorithm; defaults to TRUE.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param true.mean.difference The population mean of the differences; defaults to NULL.
#' If a vector were provided, the centered test statistic would be outputted.
#' It is only useful for a simulation purpose.
#' @param output.weights A boolean variable specifying whether the exponential weights should be outputted; defaults to FALSE.
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
#'    \code{lambda.capped} \tab Boolean variable indicating the data-driven lambda has reached the large threshold n^5 \cr
#'    \tab \cr
#'    \code{residual.slepian} \tab The final approximate first order stability term for the data-driven lambda.
#'    \tab \cr
#'    \code{variance.bound} \tab The final variance bound for the data-driven lambda.
#'    \tab \cr
#'    \code{test.stat.centered} \tab (Optional) The centered test statistic. Outputted only when true.mean is not NULL. \cr
#'    \tab \cr
#'    \code{exponential.weights} \tab (Optional) A (n.fold by p) matrix storing the exponential weightings in the test statistic. \cr
#' }
argmin.HT.LOO <- function(difference.matrix, sample.mean=NULL, min.algor='softmin',
                          lambda=NULL, const=2.5, enlarge=TRUE, alpha=0.05, true.mean.difference=NULL, output.weights=FALSE, ...){

  ## scale the difference.matrix (pre-processing step)
  sd.difference.matrix <- apply(difference.matrix, 2, stats::sd)
  non.identical.columns <- which(!(sd.difference.matrix == 0 & difference.matrix[1,] == 0))
  scaled.difference.matrix <- sweep(difference.matrix[,non.identical.columns],
                                    2, sd.difference.matrix[non.identical.columns], FUN='/')

  ## sample mean
  sample.mean <- colMeans(scaled.difference.matrix)

  if (!is.null(true.mean.difference)){
    scaled.true.mean.difference <- true.mean.difference/sd.difference.matrix
    scaled.true.mean.difference <- scaled.true.mean.difference[non.identical.columns]
  }

  ## parameters
  n <- nrow(scaled.difference.matrix)
  p.minus.1 <- ncol(scaled.difference.matrix)
  p <- p.minus.1 + 1
  val.critical <- stats::qnorm(1-alpha, 0, 1)

  ## determine lambda if needed
  capped <- FALSE
  if (is.null(lambda) & min.algor=='softmin'){
    lambda <- lambda.adaptive.LOO(scaled.difference.matrix, sample.mean=sample.mean, const=const)
    if (enlarge) {
      res <- lambda.adaptive.enlarge(
        lambda, scaled.difference.matrix, algorithm='LOO', sample.mean=sample.mean, ...)
      lambda <- res$lambda
      capped <- res$capped # A boolean variable indicating if the resulting lambda reaches the capped value
      residual.slepian <- res$residual.slepian
      variance.bound <- res$variance.bound
    }
  }

  if (output.weights){
    # n.fold by (p-1)
    exponential.weights <- matrix(0, n, p - 1)
  }

  test.stat.centered <- NULL
  if (p == 2){
    ## the algorithm reduces to the pairwise t-test
    sigma <- 1
    test.stat <- sqrt(n)*mean(scaled.difference.matrix[,1]) ## already scaled
    ans <- ifelse(test.stat < val.critical, 'Accept', 'Reject')

    ## ouptput centered test statistic
    if (!is.null(true.mean.difference)){
      test.stat.centered <- test.stat - sqrt(n)*true.mean.difference/sd.difference.matrix[1]
      ## in this case, true.mean.difference is a scalar
    }
  } else {
    diffs.weighted <- rep(NA, n)
    diffs.weighted.centered <- rep(NA, n)
    for (i in 1:n){
      sample.mean.noi <- (sample.mean*n - scaled.difference.matrix[i,])/(n-1)
      if (min.algor == 'softmin'){
        weights <- LDATS::softmax(lambda*sample.mean.noi)
        if (output.weights){
          exponential.weights[i,] <- weights
        }
        diffs.weighted[i] <- sum(scaled.difference.matrix[i,]*weights)

        ## output centered test statistic
        if (!is.null(true.mean.difference)){
          diffs.weighted.centered[i] <- diffs.weighted[i] - sum(scaled.true.mean.difference*weights)
        }
      } else if (min.algor == 'argmin') {
        min.indices <- which(sample.mean.noi == max(sample.mean.noi))
        seed.argmin.i <- ceiling(abs(23*sample.mean.noi[p.minus.1]*p) + i) %% (2^31-1)
        withr::with_seed(seed.argmin.i, {
          idx.min <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
        })
        diffs.weighted[i] <- scaled.difference.matrix[i,idx.min]

        ## output centered test statistic
        if (!is.null(true.mean.difference)){
          diffs.weighted.centered[i] <- diffs.weighted[i] - scaled.true.mean.difference[idx.min]
        }
      } else {
        # error
        stop("'min.algor' should be either 'softmin' or 'argmin'")
      }
    }

    sigma <- stats::sd(diffs.weighted)
    test.stat <- sqrt(n)*(mean(diffs.weighted))
    test.stat.scale <- test.stat/sigma
    ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')

    if (!is.null(true.mean.difference)) {
      test.stat.centered <- sqrt(n)*mean(diffs.weighted.centered)/sigma
    }
  }

  if(output.weights){
    return (list(test.stat.scale=test.stat.scale, critical.value=val.critical, std=sigma, ans=ans,
                 lambda=lambda, lambda.capped=capped,
                 residual.slepian=residual.slepian,
                 variance.bound=variance.bound,
                 test.stat.centered=test.stat.centered,
                 exponential.weights=exponential.weights))
  } else {
    return (list(test.stat.scale=test.stat.scale, critical.value=val.critical, std=sigma, ans=ans,
                 lambda=lambda, lambda.capped=capped,
                 residual.slepian=residual.slepian,
                 variance.bound=variance.bound,
                 test.stat.centered=test.stat.centered))
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
    seed <- (data[1,r]*r*37) %% (2^31 - 1)
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
