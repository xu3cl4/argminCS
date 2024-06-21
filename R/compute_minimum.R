#' Compute the minimum, using softmin.
#'
#' Compute the minimum, using softmin, under the leave-one-out (LOO) scheme;
#' more specifically, the quantity Q in Zhang et al 2024.
#'
#' @param i An index of the sample to be left out.
#' @param r The dimension of interest for hypothesis test; excluded from the calculation of softmin.
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(data) to speed up computation.
#' @param print.weights A boolean specifying if the exponential weightings are printed to the console; defaults to False.
#'
#' @return The minimum calculated from softmin under the leave-one-out scheme with the r-th dimension excluded (The quantity Q in Zhang et al).
#' @export
#'
#' @importFrom MASS mvrnorm
#' @examples
#' n <- 100
#' mu <- (1:10)/10
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' ## let the function calculate the sample.mean
#' getMin.softmin.LOO(7, 1, data, sqrt(n))
#'
#' ## calculate the sample.mean outside the function.
#' ## This may foster LOO implementation computation-wise
#' getMin.softmin.LOO(7, 1, data, sqrt(n), sample.mean=colMeans(data))
#'
#  ## require the function to print out the exponential weights to the console
#' getMin.softmin.LOO(7, 1, data, sqrt(n), sample.mean=colMeans(data), print.weights=TRUE)

getMin.softmin.LOO <- function(i, r, data, lambda, sample.mean=NULL, print.weights=FALSE){
  n <- nrow(data)
  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }
  mu.hat.noi <- (sample.mean*n - data[i,])/(n-1) #3p
  weights <- LDATS::softmax(-lambda*mu.hat.noi[-r]) #4p
  if (print.weights){
    print(glue::glue('weights are {list(weights)}'))
  }
  return (sum(weights*data[i,-r])) #2p
}

#' Compute the minimum, using argmin.
#'
#' Compute the minimum, using (hard)argmin, under the leave-one-out (LOO) scheme;
#' more specifically, the quantity Q in Zhang et al 2024.
#'
#' @param i An index of the sample to be left out.
#' @param r The dimension of interest for hypothesis test; excluded from the calculation of softmin.
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param lambda NULL; in-place to ease the implementation.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(data) to speed up computation.
#' @param ties.method A string indicating the method to tackle with ties: 'average' (or simply 'a'), 'random' (or simple 'r); defaults to average.
#' @param seed An integer seed in case that 'random' is chosen to tackle with ties. If no value is given, the seed would be set to \eqn{i*r + 11}.
#'
#' @return The minimum calculated from argmin under the leave-one-out scheme with the r-th dimension excluded (The quantity Q in Zhang et al).
#' @export
#'
#' @importFrom MASS mvrnorm
#' @examples
#' n <- 100
#' mu <- (1:10)/10
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' ## let the function calculate the sample.mean
#' getMin.argmin.LOO(7, 1, data, sqrt(n))
#'
#' ## calculate the sample.mean outside the function.
#' ## This may foster LOO implementation computation-wise
#' getMin.argmin.LOO(7, 1, data, sqrt(n), sample.mean=colMeans(data))

getMin.argmin.LOO <- function(i, r, data, lambda=NULL, sample.mean=NULL, ties.method='average', seed=NULL){
  n <- nrow(data)
  p <- ncol(data)
  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }
  mu.hat.noi <- (sample.mean*n - data[i,])/(n-1)
  min.indices <- c(1:(p-1))[which(mu.hat.noi[-r] == min(mu.hat.noi[-r]))]

  if (ties.method == 'random' | ties.method == 'r'){
    #randomly select an argmin
    if (is.null(seed)) {
      seed <- ceiling(abs(i*r*11*data[1,1]*n*p))
    }
    if (seed > 2^32){
      seed <- seed %% 2^32
    }
    withr::with_seed(seed, {
      r.i.hat <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
    })
    return (data[i,-r][r.i.hat])
  } else if (ties.method == 'average' | ties.method == 'a'){
    # average over all argmins
    return (sum((data[i,-r])[min.indices])/length(min.indices))
  } else {
    stop("'ties.method' should be one of 'average', 'random'")
  }
}
