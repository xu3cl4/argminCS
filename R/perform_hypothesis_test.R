#' A wrapper to perform argmin hypothesis test.
#'
#' This is a wrapper to perform hypothesis test to see if a given dimension may be an argmin. Multiple methods are supported.
#'
#' @importFrom Rdpack reprompt
#' @details The supported methods include:\tabular{ll}{
#'    \code{softmin.LOO (SML)} \tab LOO (leave-one-out) algorithm, using the exponential weightings. \cr
#'    \tab \cr
#'    \code{argmin.LOO (HML)} \tab A variant of SML, but it uses (hard) argmin rather than exponential weighting.
#'    The method is not recommended because its type 1 error is not controlled. \cr
#'    \tab \cr
#'    \code{nonsplit (NS)} \tab  A variant of SML, but no splitting is involved.
#'    One needs to pass a fixed lambda value as a required additional argument.
#'    The method is not recommended because its type 1 error is not controlled. \cr
#'    \tab \cr
#'    \code{Bonferroni (MT)} \tab Multiple testing with Bonferroni's correction. \cr
#'    \tab \cr
#'    \code{Gupta (GTA)} \tab The method in \insertRef{gupta.1965}{argminCS}. \cr
#' }
#'
#' @param data (1) A n by p data matrix for (GTA); each of its row is a p-dimensional sample, or
#' (2) A n by (p-1) difference matrix for (SML, HML, NS, MT); each of its row is a (p-1)-dimensional sample differences
#' @param r The dimension of interest for hypothesis test; defaults to NULL. (Only needed for GTA)
#' @param method A string indicating the method for hypothesis test; defaults to 'softmin.LOO'. Passing an abbreviation is allowed.
#' For the list of supported methods and their abbreviations, see Details.
#' @param ... Additional arguments to \link{argmin.HT.LOO}, \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}, \link{argmin.HT.MT}, \link{argmin.HT.gupta}.
#' A correct argument name needs to be specified if it is used.
#'
#' @return 'Accept' or 'Reject'. A string indicating whether the given dimension could be an argmin (Accept) or not (Reject), and relevant statistics.
#' @export
#'
#' @examples
#' r <- 4
#' n <- 200
#' p <- 20
#' mu <- (1:p)/p
#' cov <- diag(length(mu))
#' set.seed(108)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#'
#' ## softmin.LOO
#' difference.matrix.r <- matrix(rep(data[,r], p-1), ncol=p-1, byrow=FALSE) - data[,-r]
#' argmin.HT(difference.matrix.r)
#'
#' ## use seed
#' argmin.HT(difference.matrix.r, seed=19)
#'
#' # provide centered test statistic (to simulate asymptotic normality)
#' true.mean.difference.r <- mu[r] - mu[-r]
#' argmin.HT(difference.matrix.r, true.mean=true.mean.difference.r)
#'
#' # keep the data unstandardized
#' argmin.HT(difference.matrix.r, scale.input=FALSE)
#'
#' # use an user-specified lambda
#' argmin.HT(difference.matrix.r, lambda=sqrt(n)/2.5)
#'
#' # add a seed
#' argmin.HT(difference.matrix.r, seed=19)
#'
#' ## argmin.LOO/hard min
#' argmin.HT(difference.matrix.r, method='HML')
#'
#' ## nonsplit
#' argmin.HT(difference.matrix.r, method='NS', lambda=sqrt(n)/2.5)
#'
#' ## Bonferroni (choose t test because of normal data)
#' argmin.HT(difference.matrix.r, method='MT', test='t')
#' ## z test
#' argmin.HT(difference.matrix.r, method='MT', test='z')
#'
#' ## Gupta
#' critical.val <- get.quantile.gupta.selection(p=length(mu))
#' argmin.HT(data, r, method='GTA', critical.val=critical.val)
#'
#' @importFrom Rdpack reprompt
#' @references{
#'   \insertRef{zhang2024winners}{argminCS}
#'
#'   \insertRef{cck.many.moments}{argminCS}
#'
#'   \insertRef{gupta.1965}{argminCS}
#'
#'   \insertRef{futschik.1995}{argminCS}
#'}
#'
argmin.HT <- function(data, r = NULL, method = 'softmin.LOO', ...) {
  method <- tolower(method)  # Case-insensitive matching
  method <- match.arg(method,
                      choices = c('softmin.loo', 'sml',
                                  'argmin.loo', 'hml',
                                  'nonsplit', 'ns',
                                  'bonferroni', 'mt',
                                  'gupta', 'gta', 'gupta'))

  # Precompute difference.matrix if needed
  methods.needing.differences <- c('softmin.loo', 'sml', 'argmin.loo', 'hml', 'nonsplit', 'ns', 'bonferroni', 'mt')
  if (method %in% methods.needing.differences) {
    if (is.null(r)) {
      difference.matrix <- data
    } else {
      p <- ncol(data)
      difference.matrix <- matrix(rep(data[, r], p - 1), ncol = p - 1, byrow = FALSE) - data[, -r]
    }
  }

  # Dispatch based on method
  switch(method,
         'softmin.loo' =, 'sml' = argmin.HT.LOO(difference.matrix, ...),
         'argmin.loo' =, 'hml' = argmin.HT.LOO(difference.matrix, min.algor = 'argmin', ...),
         'nonsplit' =, 'ns' = argmin.HT.nonsplit(difference.matrix, ...),
         'bonferroni' =, 'mt' = argmin.HT.MT(difference.matrix, ...),
         'gupta' =, 'gta' = argmin.HT.gupta(data, r, ...),
         stop("'method' should be one of: 'softmin.LOO' (SML), 'argmin.LOO' (HML), 'nonsplit' (NS), 'Bonferroni' (MT), or 'Gupta' (GTA)")
  )
}

#' A wrapper to perform argmax hypothesis test.
#'
#' This function performs a hypothesis test to evaluate whether a given dimension may be the argmax.
#' It internally negates the data and reuses the implementation from \code{\link{argmin.HT}}.
#'
#' @importFrom Rdpack reprompt
#'
#' @details The supported methods include:\tabular{ll}{
#'   \code{softmin.LOO (SML)} \tab Leave-one-out algorithm using exponential weighting. \cr
#'   \tab \cr
#'   \code{argmin.LOO (HML)} \tab A variant of SML that uses hard argmin instead of exponential weighting. Not recommended. \cr
#'   \tab \cr
#'   \code{nonsplit (NS)} \tab Variant of SML without data splitting. Requires a fixed lambda value. Not recommended. \cr
#'   \tab \cr
#'   \code{Bonferroni (MT)} \tab Multiple testing using Bonferroni correction. \cr
#'   \tab \cr
#'   \code{Gupta (GTA)} \tab The method from \insertRef{gupta.1965}{argminCS}. \cr
#' }
#'
#' @param data (1) A n by p matrix of raw samples (for GTA), or
#' (2) A n by (p-1) difference matrix (for SML, HML, NS, MT). Each row is a sample.
#' @param r The dimension of interest for testing; defaults to NULL. Required for GTA.
#' @param method A string indicating the method to use. Defaults to 'softmin.LOO'.
#' See **Details** for supported methods and abbreviations.
#' @param ... Additional arguments passed to \code{\link{argmin.HT.LOO}},
#' \code{\link{argmin.HT.MT}}, \code{\link{argmin.HT.nonsplit}}, or \code{\link{argmin.HT.gupta}}.
#'
#' @return A character string: 'Accept' or 'Reject', indicating whether the dimension could be an argmax, and relevant statistics.
#'
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
#' ## Define the dimension of interest
#' r <- 4
#'
#' ## Construct difference matrix for dimension r
#' difference.matrix.r <- matrix(rep(data[, r], p - 1), ncol = p - 1, byrow = FALSE) - data[, -r]
#'
#' ## softmin.LOO (SML)
#' argmax.HT(difference.matrix.r)
#'
#' ## use seed
#' argmax.HT(difference.matrix.r, seed=19)
#'
#' ## With known true difference
#' true.mean.diff <- mu[r] - mu[-r]
#' argmax.HT(difference.matrix.r, true.mean = true.mean.diff)
#'
#' ## Without scaling
#' argmax.HT(difference.matrix.r, scale.input = FALSE)
#'
#' ## With a user-specified lambda
#' argmax.HT(difference.matrix.r, lambda = sqrt(n) / 2.5)
#'
#' ## Add a seed for reproducibility
#' argmax.HT(difference.matrix.r, seed = 17)
#'
#' ## argmin.LOO (HML)
#' argmax.HT(difference.matrix.r, method = "HML")
#'
#' ## nonsplit method
#' argmax.HT(difference.matrix.r, method = "NS", lambda = sqrt(n)/2.5)
#'
#' ## Bonferroni method (choose t test for normal data)
#' argmax.HT(difference.matrix.r, method = "MT", test = "t")
#'
#' ## Gupta method (pass full data matrix)
#' critical.val <- get.quantile.gupta.selection(p = length(mu))
#' argmax.HT(data, r, method = "GTA", critical.val = critical.val)
#'
#' @references
#' \insertRef{cck.many.moments}{argminCS}
#'
#' \insertRef{gupta.1965}{argminCS}
#'
#' \insertRef{futschik.1995}{argminCS}
argmax.HT <- function(data, r = NULL, method = "softmin.LOO", ...) {
  args <- list(...)

  # Negate data
  negated.data <- -data

  # If true.mean is provided, negate it
  if ("true.mean" %in% names(args)) {
    args$true.mean <- -args$true.mean
  }

  # Call argmin.HT with modified data and possibly modified true.mean
  do.call(argmin.HT, c(list(data = negated.data, r = r, method = method), args))
}


#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, using the LOO (leave-one-out) algorithm in Zhang et al 2024.
#'
#' @param difference.matrix A n by (p-1) difference data matrix (reference dimension - the rest);
#' each of its row is a (p-1)-dimensional vector of differences.
#' @param sample.mean The sample mean of differences; defaults to NULL. It can be calculated via colMeans(difference.matrix).
#' @param min.algor The algorithm to compute the test statistic by weighting across dimensions; 'softmin' uses exponential weighting,
#' while 'argmin' picks the largest mean coordinate directly. Defaults to 'softmin'.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin); defaults to NULL.
#' If lambda=NULL (recommended), the function would determine a lambda value in a data-driven way.
#' @param const The scaling constant for initial data-driven lambda
#' @param enlarge A boolean value indicating if the data-driven lambda should be determined via an iterative enlarging algorithm; defaults to TRUE.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param true.mean.difference The population mean of the differences. (Optional); used to compute a centered test statistic for simulation or diagnostic purposes.
#' @param output.weights A boolean variable specifying whether the exponential weights should be outputted; defaults to FALSE.
#' @param scale.input A boolean variable specifying whether the input difference matrix should be standardized. Defaults to TRUE
#' @param seed (Optional) If provided, used to seed the random sampling (for reproducibility).
#' @param ... Additional arguments to \link{lambda.adaptive.enlarge}, \link{is.lambda.feasible.LOO}.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{test.stat.scale} \tab The scaled test statistic \cr
#'    \tab \cr
#'    \code{critical.value} \tab The critical value for the hypothesis test. Being greater than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{std} \tab The standard deviation estimate. \cr
#'    \tab \cr
#'    \code{ans} \tab A character string: either 'Reject' or 'Accept', depending on the test outcome. \cr
#'    \tab \cr
#'    \code{lambda} \tab The lambda used in the hypothesis testing. \cr
#'    \tab \cr
#'    \code{lambda.capped} \tab Boolean variable indicating the data-driven lambda has reached the large threshold n^5 \cr
#'    \tab \cr
#'    \code{residual.slepian} \tab The final approximate first order stability term for the data-driven lambda. \cr
#'    \tab \cr
#'    \code{variance.bound} \tab The final variance bound for the data-driven lambda. \cr
#'    \tab \cr
#'    \code{test.stat.centered} \tab (Optional) The centered test statistic, computed only if \code{true.mean.difference} is provided. \cr
#'    \tab \cr
#'    \code{exponential.weights} \tab (Optional) A (n by p-1) matrix storing the exponential weightings in the test statistic. \cr
#' }
argmin.HT.LOO <- function(difference.matrix, sample.mean=NULL, min.algor='softmin',
                          lambda=NULL, const=2.5, enlarge=TRUE, alpha=0.05,
                          true.mean.difference=NULL, output.weights=FALSE, scale.input=TRUE, seed=NULL, ...){

  ## Pre-processing: scale the difference.matrix or not
  if (scale.input) {
    sd.difference.matrix <- apply(difference.matrix, 2, stats::sd)

    # Exclude columns with zero standard deviation
    non.identical.columns <- which(sd.difference.matrix != 0)

    # Scale the matrix while preserving matrix structure in case of single-column extraction
    scaled.difference.matrix <- sweep(difference.matrix[, non.identical.columns, drop=FALSE],
                                      2, sd.difference.matrix[non.identical.columns], FUN='/')

    # Adjust sample.mean and true.mean.difference accordingly
    if (is.null(sample.mean)) {
      sample.mean <- colMeans(scaled.difference.matrix)
    } else {
      sample.mean <- sample.mean[non.identical.columns]/sd.difference.matrix[non.identical.columns]
    }

    if (!is.null(true.mean.difference)) {
      scaled.true.mean.difference <- true.mean.difference[non.identical.columns]/sd.difference.matrix[non.identical.columns]
    }
  } else {
    # Use the raw matrix
    scaled.difference.matrix <- difference.matrix

    if (is.null(sample.mean)) {
      sample.mean <- colMeans(difference.matrix)
    }

    if (!is.null(true.mean.difference)) {
      scaled.true.mean.difference <- true.mean.difference
    }
  }

  ## parameters
  n <- nrow(scaled.difference.matrix)
  p.minus.1 <- ncol(scaled.difference.matrix)
  p <- p.minus.1 + 1
  val.critical <- stats::qnorm(1-alpha, 0, 1)

  ## some initializations
  residual.slepian <- NULL
  variance.bound <- NULL
  capped <- FALSE
  test.stat.centered <- NULL

  ## determine lambda if needed
  if (is.null(lambda) & min.algor=='softmin'){
    lambda <- lambda.adaptive.LOO(scaled.difference.matrix, sample.mean=sample.mean, const=const, seed=seed)
    if (enlarge) {
      res <- lambda.adaptive.enlarge(
        lambda, scaled.difference.matrix, sample.mean=sample.mean, seed=seed, ...)
      lambda <- res$lambda
      capped <- res$capped # A boolean variable indicating if the resulting lambda reaches the capped value
      residual.slepian <- res$residual.slepian
      variance.bound <- res$variance.bound
    }
  }

  if (output.weights){
    # n by (p-1)
    exponential.weights <- matrix(0, n, p - 1)
  }

  if (p == 2){
    ## the algorithm reduces to the pairwise t-test
    sigma <- stats::sd(scaled.difference.matrix[,1])
    test.stat <- sqrt(n)*mean(scaled.difference.matrix[,1]) ## already scaled
    test.stat.scale <- test.stat/sigma
    ans <- ifelse(test.stat < val.critical, 'Accept', 'Reject')

    ## output centered test statistic
    if (!is.null(true.mean.difference)){
      test.stat.centered <- test.stat - sqrt(n)*true.mean.difference/sd.difference.matrix[1]
      ## in this case, true.mean.difference is a scalar
    }

    if (output.weights){
      exponential.weights <- matrix(1, n, 1)
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
        if (!is.null(seed)) {
          # seed.argmin.i <- ceiling(abs(seed*sample.mean.noi[p.minus.1]*p) + i) %% (2^31-1)
          seed.argmin.i <- (seed*p + i) %% (2^31-1)
          withr::with_seed(seed.argmin.i, {
            idx.min <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
          })
        } else {
          idx.min <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
        }
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

  out <- list(test.stat.scale = test.stat.scale,
              critical.value = val.critical,
              std = sigma,
              ans = ans,
              lambda = lambda,
              lambda.capped = capped,
              residual.slepian = residual.slepian,
              variance.bound = variance.bound)

  if (output.weights) {
    out$exponential.weights <- exponential.weights
  }

  if (!is.null(true.mean.difference)) {
    out$test.stat.centered <- test.stat.centered
  }
  return (out)
}


#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin without any splitting.
#'
#' @details
#' This method is not recommended, given its poor performance when p is small.
#'
#' @param difference.matrix A n by (p-1) difference data matrix (reference dimension - the rest);
#' each of its row is a (p-1)-dimensional vector of differences.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param sample.mean The sample mean of differences; defaults to NULL. It can be calculated via colMeans(difference.matrix).
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#' @param scale.input A boolean variable specifying whether the input difference matrix should be standardized defaults to TRUE
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
argmin.HT.nonsplit <- function(difference.matrix, lambda, sample.mean=NULL, alpha=0.05, scale.input=TRUE){

  n <- nrow(difference.matrix)
  p <- ncol(difference.matrix)

  if (scale.input) {
    # Scale the difference.matrix (pre-processing step)
    sd.difference.matrix <- apply(difference.matrix, 2, stats::sd)
    # Filter out columns with zero sd and zero in first row (non-identical)
    non.identical.columns <- which(!(sd.difference.matrix == 0 & difference.matrix[1,] == 0))
    scaled.difference.matrix <- sweep(difference.matrix[,non.identical.columns, drop=FALSE],
                                      2, sd.difference.matrix[non.identical.columns], FUN='/')

    # Adjust sample.mean accordingly
    if (is.null(sample.mean)){
      sample.mean <- colMeans(scaled.difference.matrix)
    } else {
      sample.mean <- sample.mean[non.identical.columns]/sd.difference.matrix[non.identical.columns]
    }
  } else {
    # Use the raw matrix without scaling
    scaled.difference.matrix <- difference.matrix
    if (is.null(sample.mean)) {
      sample.mean <- colMeans(difference.matrix)
    }
  }

  val.critical <- stats::qnorm(1-alpha, 0, 1)

  weights <- LDATS::softmax(lambda*sample.mean)
  diffs <- scaled.difference.matrix %*% weights

  sigma <- stats::sd(diffs)
  test.stat <- sqrt(n)*(mean(diffs))
  test.stat.scale <- test.stat/sigma

  ans <- ifelse(test.stat.scale < val.critical, 'Accept', 'Reject')
  return (list(test.stat.scale=test.stat.scale, critical.value=val.critical, std=sigma, ans=ans))
}

#' Perform argmin hypothesis test.
#'
#' Test if a dimension may be argmin, using multiple testing with Bonferroni's correction.
#'
#' @param difference.matrix A n by (p-1) difference data matrix (reference dimension - the rest);
#' each of its row is a (p-1)-dimensional vector of differences.
#' @param sample.mean The sample mean of differences; defaults to NULL. It can be calculated via colMeans(difference.matrix).
#' @param test The test to perform: 't' or 'z'; defaults to 'z'.
#' If the data are assumed normally distributed, use 't'; otherwise 'z'.
#' @param alpha The significance level of the hypothesis test; defaults to 0.05.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{p.val} \tab p value without Bonferroni's correction. \cr
#'    \tab \cr
#'.   \code{critical.value} \tab The critical value for the hypothesis test. Being less than it leads to a rejection. \cr
#'    \tab \cr
#'    \code{ans} \tab 'Reject' or 'Accept' \cr
#' }
argmin.HT.MT <- function(difference.matrix, sample.mean=NULL, test='z', alpha=0.05){

  test <- match.arg(test, choices = c('t', 'z'))

  p <- ncol(difference.matrix) + 1
  val.critical <- alpha/(p-1)

  sd.difference.matrix <- apply(difference.matrix, 2, stats::sd)
  non.identical.columns <- which(!(sd.difference.matrix == 0 & difference.matrix[1,] == 0))

  if (is.null(sample.mean)){
    mean.difference <- colMeans(difference.matrix[,non.identical.columns,drop = FALSE])
  } else {
    mean.difference <- sample.mean[non.identical.columns]
  }

  scaled.mean.differences <- mean.difference/sd.difference.matrix[non.identical.columns]
  argmax <- which.max(scaled.mean.differences)

  difference.matrix.non.identical <- difference.matrix[,non.identical.columns,drop = FALSE]
  if (test == 't') {
    res <- stats::t.test(difference.matrix.non.identical[,argmax], alternative='greater')
  } else {
    ## z test

    # ## hard-coded version
    # sample.size <- nrow(difference.matrix)
    # wald.test.stat <-
    #   mean(difference.matrix.non.identical[,argmax])/(stats::sd(difference.matrix.non.identical[,argmax])/sqrt(sample.size))
    # p.value <- 1 - pnorm(wald.test.stat)
    # res <- list(p.value=p.value, test.stat=wald.test.stat)

    res <- BSDA::z.test(difference.matrix.non.identical[,argmax],
                        sigma.x=stats::sd(difference.matrix.non.identical[,argmax]),
                          alternative='greater')
  }
  p.val <- res$p.value
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
#' @references
#' \insertRef{gupta.1965}{argminCS}
#'
#' \insertRef{futschik.1995}{argminCS}
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

#' @title Perform argmin hypothesis test using Gupta's method.
#'
#' @importFrom Rdpack reprompt
#' @description Test whether a dimension is the argmin, using the method in \insertCite{gupta.1965}{argminCS}.
#'
#' @note This method requires independence among the dimensions.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If performing multiple tests across dimensions, pre-computing \code{sample.mean} and \code{critical.val}
#' can significantly reduce computation time.
#' @param stds A vector of the same (population) standard deviations for all dimensions; defaults to a vector of 1's.
#' These are used to standardize the sample means.
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
#' @references
#' \insertRef{gupta.1965}{argminCS}
#'
#' \insertRef{futschik.1995}{argminCS}
#'
argmin.HT.gupta <- function(data, r, sample.mean=NULL, stds=NULL, critical.val=NULL, alpha=0.05, ...){

  # note that the implementation can, in turn, asks for scaled.sample.mean, its min and second min
  # to speed up the computation, in case that there are many dimensions involving.
  n <- nrow(data)
  p <- ncol(data)

  if (r < 1 || r > p) stop("Parameter 'r' must be an integer between 1 and the number of columns in 'data'.")

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
