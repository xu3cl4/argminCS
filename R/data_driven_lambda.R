#' Iteratively enlarge a tuning parameter \eqn{\lambda} in a data-driven way.
#'
#' Iteratively enlarge a tuning parameter \eqn{\lambda} to enhance the power of hypothesis testing.
#' The iterative algorithm ends when an enlarged \eqn{\lambda} unlikely yields the first order stability.
#'
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param scaled.difference.matrix A n by (p-1) difference scaled.difference.matrix matrix after column-wise scaling (reference dimension - the rest);
#' each of its row is a (p-1)-dimensional vector of differences.
#' @param sample.mean The sample mean of the n samples in scaled.difference.matrix; defaults to NULL. It can be calculated via colMeans(scaled.difference.matrix).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(scaled.difference.matrix) to speed up computation.
#' @param mult.factor In each iteration, \eqn{\lambda} would be multiplied by mult.factor to yield an enlarged \eqn{\lambda}; defaults to 2.
#' @param verbose A boolean value indicating if the number of iterations should be printed to console; defaults to FALSE.
#' @param seed (Optional) If provided, used to seed for tie-breaking (for reproducibility).
#' @param ... Additional arguments to \link{is.lambda.feasible.LOO}.
#'
#' @return A list containing:
#' \tabular{ll}{
#'   \code{lambda} \tab The final (enlarged) lambda that is still feasible. \cr
#'   \tab \cr
#'   \code{capped} \tab Logical, \code{TRUE} if the enlargement was capped due to reaching the threshold. \cr
#'   \tab \cr
#'   \code{residual.slepian} \tab Residual value from the feasibility check at the final lambda. \cr
#'   \tab \cr
#'   \code{variance.bound} \tab Variance bound used in the final feasibility check. \cr
#' }
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' r <- 4
#' n <- 200
#' mu <- (1:20)/20
#' cov <- diag(length(mu))
#' set.seed(108)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#' diff.mat <- get.difference.matrix(data, r)
#' sample.mean.r <- get.sample.mean.r(sample.mean, r)
#' lambda <- lambda.adaptive.LOO(diff.mat, sample.mean=sample.mean.r)
#'
#' # Run the enlargement algorithm
#' res <- lambda.adaptive.enlarge(lambda, diff.mat, sample.mean=sample.mean.r, seed=3)
#' res
#' # with a seed
#' res <- lambda.adaptive.enlarge(lambda, diff.mat, sample.mean=sample.mean.r, seed=3)
#' res
#'
#' @export
lambda.adaptive.enlarge <- function(lambda, scaled.difference.matrix, sample.mean=NULL,
                                    mult.factor=2, verbose=FALSE, seed=NULL, ...){

  if (lambda < 0) {
    stop('The weighting parameter lambda needs to be non-negative.')
  }

  n <- nrow(scaled.difference.matrix)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(scaled.difference.matrix)
  }

  lambda.curr <- lambda
  lambda.next <- mult.factor*lambda
  threshold <- n^5 # any large value

  # keep track of the feasibility criterion bounds
  residual.slepian <- 0
  variance.bound <- 0

  # check feasibility of the next lambda
  res <- is.lambda.feasible.LOO(lambda.next, scaled.difference.matrix, sample.mean=sample.mean, seed=seed, ...)
  feasible <- res$feasible
  residual.slepian.next <- res$residual.slepian
  variance.bound.next <- res$variance.bound

  count <- 1
  while (feasible & lambda.next <= threshold){
    ## updating rule
    lambda.curr <- lambda.next
    residual.slepian <- residual.slepian.next
    variance.bound <- variance.bound.next
    lambda.next <- mult.factor*lambda.next

    ## check feasibility
    res <- is.lambda.feasible.LOO(lambda.next, scaled.difference.matrix, sample.mean=sample.mean, seed=seed, ...)
    feasible <- res$feasible
    residual.slepian.next <- res$residual.slepian
    variance.bound.next <- res$variance.bound
    count <- count + 1
  }
  if (residual.slepian == 0 & variance.bound == 0){
    # just to store how it fails to perform the first iteration
    residual.slepian <- residual.slepian.next
    variance.bound <- variance.bound.next
  }

  if (verbose){
    print(glue::glue('before: {lambda}, after: {lambda.curr}, iteration: {count}'))
  }

  capped <- ifelse(lambda.next <= threshold, FALSE, TRUE)
  return (list(lambda=lambda.curr, capped=capped,
               residual.slepian=residual.slepian,
               variance.bound=variance.bound))
}

#' Generate a scaled.difference.matrix-driven \eqn{\lambda} for LOO algorithm.
#'
#' Generate a scaled.difference.matrix-driven \eqn{\lambda} for LOO algorithm motivated by the derivation of the first order stability.
#' For its precise definition, we refer to the paper Zhang et al 2024.
#'
#' @param scaled.difference.matrix A n by (p-1) difference scaled.difference.matrix matrix after column-wise scaling (reference dimension - the rest);
#' each of its row is a (p-1)-dimensional vector of differences.
#' @param sample.mean The sample mean of the n samples in scaled.difference.matrix; defaults to NULL. It can be calculated via colMeans(scaled.difference.matrix).
#' @param const A scaling constant for the scaled.difference.matrix driven \eqn{\lambda}; defaults to 2.5. As its value gets larger, the first order stability tends to increase while power might decrease.
#' @param seed (Optional) If provided, used to seed for tie-breaking (for reproducibility).
#'
#' @return A scaled.difference.matrix-driven \eqn{\lambda} for LOO algorithm.
#'
#' @examples
#' # Simulate data
#' set.seed(123)
#' r <- 4
#' n <- 200
#' mu <- (1:20)/20
#' cov <- diag(length(mu))
#' set.seed(108)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#' diff.mat <- get.difference.matrix(data, r)
#' sample.mean.r <- get.sample.mean.r(sample.mean, r)
#' lambda <- lambda.adaptive.LOO(diff.mat, sample.mean=sample.mean.r)
#'
#' @export
#' @importFrom MASS mvrnorm
lambda.adaptive.LOO <- function(scaled.difference.matrix, sample.mean=NULL, const=2.5, seed=NULL){
  n <- nrow(scaled.difference.matrix)
  p.minus.1 <- ncol(scaled.difference.matrix)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(scaled.difference.matrix)
  }
  Xs.min <- sapply(1:n, function(i){
    mu.hat.noi <- (sample.mean*n - scaled.difference.matrix[i,])/(n-1)

    min.indices <- which(mu.hat.noi == max(mu.hat.noi))
    if (!is.null(seed)) {
      seed <- ceiling(abs(seed*scaled.difference.matrix[n,p.minus.1] + i*seed)) %% (2^31 - 1)
      withr::with_seed(seed, {
        min.idx <- ifelse(length(min.indices) > 1,
                          sample(c(min.indices), 1), min.indices[1])
      })
    } else {
      min.idx <- ifelse(length(min.indices) > 1,
                        sample(c(min.indices), 1), min.indices[1])
    }

    X.min <- scaled.difference.matrix[i,min.idx]
    return (X.min)
  })

  return (sqrt(n)/(const*stats::sd(Xs.min)))
}
