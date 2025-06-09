#' Check the feasibility of a tuning parameter \eqn{\lambda} for LOO algorithm.
#'
#' Check the feasibility of a tuning parameter \eqn{\lambda} for LOO algorithm by examining
#' whether its resulting \eqn{\nabla_i K_j} is less than a threshold value,
#' i.e., the first order stability is likely achieved.
#' For further details, we refer to the paper Zhang et al 2024.
#'
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param scaled.difference.matrix A n by (p-1) difference scaled.difference.matrix matrix after column-wise scaling (reference dimension - the rest);
#' each of its row is a (p-1)-dimensional vector of differences.
#' @param sample.mean The sample mean of the n samples in scaled.difference.matrix; defaults to NULL. It can be calculated via colMeans(scaled.difference.matrix).
#' If your experiment involves hypothesis testing over more than one dimension, pass sample.mean=colMeans(scaled.difference.matrix) to speed up computation.
#' @param threshold.1 A threshold value to examine if the first order stability is likely achieved; defaults to 0.05. As its value gets smaller, the first order stability tends to increase while power might decrease.
#' @param n.pairs The number of \eqn{(i,j)} pairs for estimation; defaults to 100.
#' @param seed (Optional) An integer-valued seed for subsampling.
#'
#' @return A boolean value indicating if the given \eqn{\lambda} likely gives the first order stability.
#' @export
#'
#' @importFrom MASS mvrnorm
is.lambda.feasible.LOO <- function(lambda, scaled.difference.matrix, sample.mean=NULL,
                                   threshold.1=0.05, n.pairs=100, seed=NULL){

  n <- nrow(scaled.difference.matrix)
  p.minus.1 <- ncol(scaled.difference.matrix)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(scaled.difference.matrix)
  }

  ## estimate the expectation of nabla i K 1 squares
  ## using leave-two-out estimates
  ## use at most 100 pairs of samples for estimation
  n.pairs <- min(n.pairs, n)

  # sub-sample from the given sample
  if (!is.null(seed)){
    seed <- ceiling(abs(seed*scaled.difference.matrix[n,p.minus.1]*lambda + n.pairs*seed)) %% (2^31 - 1)
    withr::with_seed(seed, {
    # we simply shuffle and sub-sample the indices 1:n to get triplets (i,j,k)
    indices <- sample(n, n.pairs, replace=FALSE)
    })
  } else {
    indices <- sample(n, n.pairs, replace=FALSE)
  }

  differences.by.perturbing.one <- rep(NA, n.pairs)
  diffs.weighted <- rep(NA, n.pairs)
  # diffs.weighted.true.mean <- rep(NA, n.pairs)
  for (repeat_idx in (1:n.pairs)) {
    index.j <- indices[repeat_idx]
    index.i <- indices[ifelse((repeat_idx + 1) > n.pairs, (repeat_idx + 1) %% n.pairs, (repeat_idx + 1))]
    index.k <- indices[ifelse((repeat_idx + 2) > n.pairs, (repeat_idx + 2) %% n.pairs, (repeat_idx + 2))]

    mu.no.j.k <- (sample.mean*n - scaled.difference.matrix[index.j,] - scaled.difference.matrix[index.k,])/(n-2)
    mu.no.j.i <- (sample.mean*n - scaled.difference.matrix[index.j,] - scaled.difference.matrix[index.i,])/(n-2)

    weights.1 <- LDATS::softmax(lambda*mu.no.j.k)
    weights.2 <- LDATS::softmax(lambda*mu.no.j.i)

    difference <- sum((weights.1 - weights.2)*(scaled.difference.matrix[index.j,] - sample.mean))
    differences.by.perturbing.one[repeat_idx] <- difference

    mu.no.j <- (sample.mean*n - scaled.difference.matrix[index.j,])/(n-1)
    weights <- LDATS::softmax(lambda*mu.no.j)
    diffs.weighted[repeat_idx] <- sum(weights*scaled.difference.matrix[index.j,])
    # diffs.weighted.true.mean[repeat_idx] <- sum(weights*sample.mean) # use sample mean to estimate true mean
  }
  difference.by.perturbing.one.squared <- mean(differences.by.perturbing.one^2)
  residual.slepian <- n*difference.by.perturbing.one.squared

  # estimate variance by leaving one out
  variance <- stats::var(diffs.weighted)
  upper.bound.1 <- threshold.1*variance

  return (list(feasible=ifelse(residual.slepian < upper.bound.1, TRUE, FALSE),
               residual.slepian=residual.slepian, variance.bound=upper.bound.1))
}
