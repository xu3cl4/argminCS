#' Check the feasibility of a tuning parameter \eqn{\lambda}
#'
#' Check the feasibility of a tuning parameter \eqn{\lambda} by examining
#' whether its resulting \eqn{\nabla_i K_j} is less than a threshold value,
#' i.e., the first order stability is likely achieved.
#' For further details, we refer to the paper Zhang et al 2024.
#'
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(data) to speed up computation.
#' @param threshold A threshold value to examine if the first order stability is likely achieved; defaults to 0.05. As its value gets smaller, the first order stability tends to increase while power might decrease.
#' @param n.pairs The number of \eqn{(i,j)} pairs for estimation; defaults to 100.
#' @param seed An integer-valued seed for subsampling. If no value is given, the seed would be set, using the value of other arguments.
#'
#' @return A boolean value indicating if the given \eqn{\lambda} likely gives the first order stability.
#' @export
#'
#' @importFrom MASS mvrnorm
#' @examples
#' n <- 300
#' mu <- (1:10)/10
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#'
#' ## compute a data-driven lambda, and check its feasibility
#' lambda <- lambda.adaptive(data, 1, sample.mean=sample.mean)
#' is.lambda.feasible(lambda, data, 1, sample.mean=sample.mean)
#'
#' ## want to ensure a greater stability
#' is.lambda.feasible(lambda, data, 1, sample.mean=sample.mean, threshold=0.01)
#'
#' ## smaller n.pairs to speed up computation
#' is.lambda.feasible(lambda, data, 1, sample.mean=sample.mean, n.pairs=50)
is.lambda.feasible <- function(lambda, data, r, sample.mean=NULL, threshold=0.05, n.pairs=100, seed=NULL){

  n <- nrow(data)
  p <- ncol(data)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  ## use at most 100 pairs of samples for estimation
  n.pairs <- 100
  if (n < 200){
    if (n %% 2 == 0){
      n.pairs <- as.integer(n/2)
    } else {
      n.pairs <- as.integer((n - 1)/2)
    }
  }

  ## subsample from the given sample
  #if (is.null(seed)){
  #  seed <- ceiling(abs(data[1,1]*sample.mean[1]*lambda))
  #}
  #withr::with_seed(seed, {
  sample.indices <- sample(n, 2*n.pairs)
  #})
  index.pairs <- cbind(sample.indices[1:n.pairs], sample.indices[(n.pairs+1):(2*n.pairs)])

  differences.by.perturbing.one <- sapply(1:n.pairs, function(i){
    # i is the perturbed index
    index.first <- index.pairs[i,1]
    index.second <- index.pairs[i,2]

    # get an index different from index.first and index.second
    j <- max(index.first, index.second)+1
    if (j > n){
      index.candidates <- setdiff(1:3, c(index.first, index.second))
      if (length(index.candidates) > 1){
        #withr::with_seed(ceiling(i*seed), {
        j <- sample(index.candidates, 1)
        #})
      } else {
        j <- index.candidates[1]
      }
    }

    mu.no.j.i <- (sample.mean*n - data[j,] - data[index.second,])/(n-2)
    mu.no.j.i.perturbed <- (sample.mean*n - data[j,] - data[index.first,])/(n-2)

    weights.1 <- LDATS::softmax(-lambda*mu.no.j.i[-r])
    weights.2 <- LDATS::softmax(-lambda*mu.no.j.i.perturbed[-r])

    difference <- sum((weights.1 - weights.2)*(data[j, -r] - sample.mean[-r]))
    return (difference)
  })

  difference.by.perturbing.one.squared <- mean(differences.by.perturbing.one^2)

  # estimate variance by leaving one out
  Qs <- sapply(index.pairs[,1], function(x) return (getMin.softmin.LOO(x, r, data, lambda, sample.mean)))
  diffs <- data[index.pairs[,1],r] - Qs
  variance <- stats::var(diffs)

  scaled.difference.by.perturbing.one.squared <- difference.by.perturbing.one.squared/variance
  return (ifelse(n*scaled.difference.by.perturbing.one.squared < threshold, T, F))
}

#' Iteratively enlarge a tuning parameter \eqn{\lambda} in a data-driven way.
#'
#' Iteratively enlarge a tuning parameter \eqn{\lambda} to enhance the power of hypothesis testing.
#' The iterative algorithm ends when an enlarged \eqn{\lambda} unlikely yields the first order stability.
#'
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(data) to speed up computation.
#' @param mult.factor In each iteration, \eqn{\lambda} would be multiplied by mult.factor to yield an enlarged \eqn{\lambda}; defaults to 2.
#' @param verbose A boolean value indicating if the number of iterations should be printed to console; defaults to FALSE.
#' @param ... Additional arguments to \link{is.lambda.feasible}.
#'
#' @return An (potentially) enlarged \eqn{\lambda}.
#' @export
#'
#' @seealso \link{is.lambda.feasible}
#'
#' @importFrom MASS mvrnorm
#' @examples
#' n <- 300
#' mu <- (1:10)/5
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#'
#' ## compute a data-driven lambda, and check its feasibility
#' lambda <- lambda.adaptive(data, 1, sample.mean=sample.mean)
#' lambda.adaptive.enlarge(lambda, data, 1, sample.mean=sample.mean)
#'
#' ## want to ensure a greater stability
#' lambda.adaptive.enlarge(lambda, data, 1, sample.mean=sample.mean, threshold=0.001)
#'
#' ## print out the number of iterations
#' lambda.adaptive.enlarge(lambda, data, 1, sample.mean=sample.mean, verbose=TRUE)
lambda.adaptive.enlarge <- function(lambda, data, r, sample.mean=NULL, mult.factor=2, verbose=FALSE, ...){

  n <- nrow(data)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  lambda.curr <- lambda
  lambda.next <- mult.factor*lambda
  feasible <- is.lambda.feasible(lambda.next, data, r, sample.mean=sample.mean, ...)
  count <- 1
  while (feasible & lambda.next < n){
    lambda.curr <- lambda.next
    lambda.next <- mult.factor*lambda.next
    feasible <- is.lambda.feasible(lambda.next, data, r, sample.mean=sample.mean, ...)
    count <- count + 1
  }
  if (verbose){
    print(glue::glue('before: {lambda}, after: {lambda.curr}, iteration: {count}'))
  }
  return (lambda.curr)
}

#' Generate a data-driven \eqn{\lambda}.
#'
#' Generate a data-driven \eqn{\lambda} motivated by the derivation of the first order stability.
#' For its precise definition, we refer to the paper Zhang et al 2024.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' @param const A scaling constant for the data driven \eqn{\lambda}; defaults to 2.5. As its value gets larger, the first order stability tends to increase while power might decrease.
#'
#' @return A data-driven \eqn{\lambda}.
#' @export
#'
#' @importFrom MASS mvrnorm
#' @examples
#' n <- 300
#' mu <- (1:10)/10
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#'
#' ## compute a data-driven lambda, and check its feasibility
#' lambda.adaptive(data, 1, sample.mean=sample.mean)
#'
#' ## want to ensure a greater stability
#' lambda.adaptive(data, 1, sample.mean=sample.mean, const=3)
lambda.adaptive <- function(data, r, sample.mean=NULL, const=2.5){
  n <- nrow(data)
  p <- ncol(data)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }
  Xs.min <- sapply(1:n, function(i){
    mu.hat.noi <- (sample.mean*n - data[i,])/(n-1)

    min.indices <- which(mu.hat.noi[-r] == min(mu.hat.noi[-r]))

    #seed <- ceiling(abs(i*r*data[1,1]*sample.mean[1]))
    #withr::with_seed(seed, {
    min.idx <- ifelse(length(min.indices) > 1,
                        sample(c(min.indices), 1), min.indices[1])
    #})

    X.min <- data[i, -r][min.idx]
    return (X.min)
  })
  return (sqrt(n)/(const*stats::sd(Xs.min)))
}
