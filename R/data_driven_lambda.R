#' Iteratively enlarge a tuning parameter \eqn{\lambda} in a data-driven way.
#'
#' Iteratively enlarge a tuning parameter \eqn{\lambda} to enhance the power of hypothesis testing.
#' The iterative algorithm ends when an enlarged \eqn{\lambda} unlikely yields the first order stability.
#'
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param algorithm 'LOO' or 'fold'.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(data) to speed up computation.
#' @param flds A list of row position integers corresponding to folds. It is for the fixed fold algorithm.
#' @param mult.factor In each iteration, \eqn{\lambda} would be multiplied by mult.factor to yield an enlarged \eqn{\lambda}; defaults to 2.
#' @param verbose A boolean value indicating if the number of iterations should be printed to console; defaults to FALSE.
#' @param ... Additional arguments to \link{is.lambda.feasible.LOO}.
#'
#' @return An (potentially) enlarged \eqn{\lambda}.
#' @export
#'
#' @seealso \link{is.lambda.feasible.LOO}
#'
#' @importFrom MASS mvrnorm
#' @importFrom caret createFolds
#' @examples
#' n <- 300
#' mu <- (1:10)/5
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#'
#' ## compute a data-driven lambda, and check its feasibility (LOO algorithm)
#' lambda <- lambda.adaptive(data, 1, sample.mean=sample.mean)
#' lambda.adaptive.enlarge(lambda, data, 1, algorithm='LOO',
#'         sample.mean=sample.mean)
#'
#' ## want to ensure a greater stability
#' lambda.adaptive.enlarge(lambda, data, 1, algorithm='LOO',
#'         sample.mean=sample.mean, threshold=0.001)
#'
#' ## print out the number of iterations
#' lambda.adaptive.enlarge(lambda, data, 1, algorithm='LOO',
#'         sample.mean=sample.mean, verbose=TRUE)
#'
#' ## compute a data-driven lambda, and check its feasibility (fold algorithm)
#' flds <- caret::createFolds(1:n, k=2, list=TRUE, returnTrain=FALSE)
#' lambda <- lambda.adaptive(data, 1, sample.mean=sample.mean)
#' lambda.adaptive.enlarge(lambda, data, 1, algorithm='fold',
#'           sample.mean=sample.mean, flds=flds)
#'
#' ## want to ensure a greater stability
#' lambda.adaptive.enlarge(lambda, data, 1, algorithm='fold',
#'           sample.mean=sample.mean, flds=flds, threshold=0.001)
#'
#' ## print out the number of iterations
#' lambda.adaptive.enlarge(lambda, data, 1, algorithm='fold',
#'           sample.mean=sample.mean, flds=flds, verbose=TRUE)
lambda.adaptive.enlarge <- function(lambda, data, r, algorithm, sample.mean=NULL, flds=NULL,
                                    mult.factor=2, verbose=FALSE, ...){

  n <- nrow(data)

  if (is.null(sample.mean) & algorithm=='LOO'){
    sample.mean <- colMeans(data)
  }

  if (is.null(flds) & algorithm=='fold'){
    stop("lambda.adaptive.enlarge: needs to provide 'flds' when tuning the parameter for fixed fold algorithm")
  }

  lambda.curr <- lambda
  lambda.next <- mult.factor*lambda
  if (algorithm == 'LOO'){
    feasible <- is.lambda.feasible.LOO(lambda.next, data, r, sample.mean=sample.mean, ...)
    threshold <- n
    # threshold <- n^2
  } else if (algorithm == 'fold'){
    feasible <- is.lambda.feasible.fold(lambda.next, data, r, sample.mean=sample.mean, flds=flds, ...)
    ## experiments over the threshold
    n.fold <- length(flds)
    if (n.fold == 2){
      threshold <- n^2
    } else if (n.fold== 5){
      threshold <- n^(3/2)
    } else{
      stop("lambda.adaptive.enlarge: only supports LOO, 2-fold and 5-fold algorithms for now")
    }
  } else {
    stop("'algorithm' should be either 'LOO' or 'fold'")
  }
  count <- 1
  while (feasible & lambda.next < threshold){
    lambda.curr <- lambda.next
    lambda.next <- mult.factor*lambda.next
    if (algorithm == 'LOO'){
      feasible <- is.lambda.feasible.LOO(lambda.next, data, r, sample.mean=sample.mean, ...)
    } else {
      feasible <- is.lambda.feasible.fold(lambda.next, data, r, sample.mean=sample.mean, flds=flds, ...)
    }
    count <- count + 1
  }
  if (verbose){
    print(glue::glue('before: {lambda}, after: {lambda.curr}, iteration: {count}'))
  }
  return (lambda.curr)
  # return (lambda.next)
}

#' Generate a data-driven \eqn{\lambda} for LOO algorithm.
#'
#' Generate a data-driven \eqn{\lambda} for LOO algorithm motivated by the derivation of the first order stability.
#' For its precise definition, we refer to the paper Zhang et al 2024.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param sample.mean The sample mean of the n samples in data; defaults to NULL. It can be calculated via colMeans(data).
#' @param const A scaling constant for the data driven \eqn{\lambda}; defaults to 2.5. As its value gets larger, the first order stability tends to increase while power might decrease.
#'
#' @return A data-driven \eqn{\lambda} for LOO algorithm.
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

    seed <- ceiling(abs(19*i*r*data[r,r]*sample.mean[r])+i) %% (2^31 - 1)
    withr::with_seed(seed, {
      min.idx <- ifelse(length(min.indices) > 1,
                        sample(c(min.indices), 1), min.indices[1])
    })

    X.min <- data[i,-r][min.idx]
    return (X.min)
  })
  lambda.suggested <- sqrt(n)/(const*stats::sd(Xs.min))
  return (lambda.suggested)
}

#' Generate a data-driven \eqn{\lambda} for fixed fold algorithm.
#'
#' Generate a data-driven \eqn{\lambda} for fixed fold algorithm motivated by the derivation of the first order stability.
#'
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param flds A list of row position integers corresponding to folds.
#' @param const A scaling constant for the data driven \eqn{\lambda}; defaults to 2.5. As its value gets larger, the first order stability tends to increase while power might decrease.
#'
#' @return A data-driven \eqn{\lambda} for fixed fold algorithm.
#' @export
#'
#' @importFrom caret createFolds
#' @importFrom MASS mvrnorm
#' @examples
#' n <- 300
#' mu <- (1:10)/10
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' flds <- caret::createFolds(1:n, k=2, list=TRUE, returnTrain=FALSE)
#'
#' ## compute a data-driven lambda, and check its feasibility
#' lambda.adaptive.fold(data, 1, flds)
#'
#' ## want to ensure a greater stability
#' lambda.adaptive.fold(data, 1, flds, const=3)
lambda.adaptive.fold <- function(data, r, flds, const=2.5){
  n <- nrow(data)
  n.fold <- length(flds)

  Xs.min <- lapply(1:n.fold, function(fold){
    mu.out.fold <- colMeans(data[setdiff(1:n, flds[[fold]]),])
    min.idx <- find.sub.argmin(mu.out.fold, r)
    Xs.min.fold <- data[flds[[fold]], min.idx]
    return (Xs.min.fold)
  })
  Xs.min <- do.call(c, Xs.min)

  return (sqrt(n)/(const*stats::sd(Xs.min)))
}
