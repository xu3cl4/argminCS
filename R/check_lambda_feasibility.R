#' Check the feasibility of a tuning parameter \eqn{\lambda} for LOO algorithm.
#'
#' Check the feasibility of a tuning parameter \eqn{\lambda} for LOO algorithm by examining
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
#' is.lambda.feasible.LOO(lambda, data, 1, sample.mean=sample.mean)
#'
#' ## want to ensure a greater stability
#' is.lambda.feasible.LOO(lambda, data, 1, sample.mean=sample.mean, threshold=0.01)
#'
#' ## smaller n.pairs to speed up computation
#' is.lambda.feasible.LOO(lambda, data, 1, sample.mean=sample.mean, n.pairs=50)
is.lambda.feasible.LOO <- function(lambda, data, r, sample.mean=NULL, threshold=0.05, n.pairs=100, seed=NULL){

  n <- nrow(data)
  p <- ncol(data)

  if (is.null(sample.mean)){
    sample.mean <- colMeans(data)
  }

  ## estimate the expectation of nabla i K 1 squares
  ## using leave-two-out estimates
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
  if (is.null(seed)){
    seed <- ceiling(abs(17*data[1,r]*sample.mean[r]*lambda*r + r))
  }
  if (seed >  2^31 - 1){
    seed <- seed %%  2^31 - 1
  }
  withr::with_seed(seed, {
    sample.indices <- sample(n, 2*n.pairs)
  })
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
        new.seed <- 37*i*seed + i
        if (new.seed >  2^31 - 1){
          new.seed <- new.seed %%  2^31 - 1
        }
        withr::with_seed(new.seed, {
          j <- sample(index.candidates, 1)
        })
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


#' Check the feasibility of a tuning parameter \eqn{\lambda} for LOO algorithm.
#'
#' Check the feasibility of a tuning parameter \eqn{\lambda} for LOO algorithm by examining
#' whether its resulting \eqn{\nabla_i K_j} is less than a threshold value,
#' i.e., the first order stability is likely achieved.
#' For further details, we refer to the paper Zhang et al 2024.
#'
#' @param lambda The real-valued tuning parameter for exponential weightings (the calculation of softmin).
#' @param data A n by p data matrix; each of its row is a p-dimensional sample.
#' @param r The dimension of interest for hypothesis test.
#' @param flds A list of row position integers corresponding to folds.
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
#' @importFrom caret createFolds
#' @examples
#' n <- 300
#' mu <- (1:10)/10
#' cov <- diag(length(mu))
#' set.seed(31)
#' data <- MASS::mvrnorm(n, mu, cov)
#' sample.mean <- colMeans(data)
#' flds <- caret::createFolds(1:n, k=2, list=TRUE, returnTrain=FALSE)
#'
#' ## compute a data-driven lambda, and check its feasibility
#' lambda <- lambda.adaptive.fold(data, 1, flds=flds)
#' is.lambda.feasible.fold(lambda, data, 1, sample.mean=sample.mean, flds=flds)
#'
#' ## want to ensure a greater stability
#' is.lambda.feasible.fold(lambda, data, 1, sample.mean=sample.mean, flds=flds, threshold=0.01)
#'
#' ## smaller n.pairs to speed up computation
#' is.lambda.feasible.fold(lambda, data, 1, sample.mean=sample.mean, flds=flds, n.pairs=50)
is.lambda.feasible.fold <- function(lambda, data, r, flds, sample.mean=NULL, threshold=0.05, n.pairs=100, seed=NULL){

  n <- nrow(data)
  p <- ncol(data)
  n.fold <- length(flds)

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
  n.pairs.per.fold <- floor(n.pairs/n.fold)

  if (is.null(seed)){
    seed <- ceiling(abs(17*data[1,r]*lambda*r + r))
  }

  differences.by.perturbing.one <- c()
  diffs <- c()
  for (fold in 1:n.fold){
    out.fold.sample <- data[setdiff(1:n, flds[[fold]]),]
    in.fold.sample <- data[flds[[fold]],]
    n.out.fold <- nrow(out.fold.sample)
    n.in.fold <- nrow(in.fold.sample)
    mu.out.fold <- colMeans(out.fold.sample)
    n.pairs <- n.pairs.per.fold

    # subsample pairs of indices from the out.fold.sample for leave two out estimates
    # subsample indices from in fold sample
    seed.fold <- n.fold*(seed) + n.fold
    if (seed.fold >  2^31 - 1){
      seed.fold <- seed.fold %% (2^31 - 1)
    }
    if (2*n.pairs > n.out.fold){
      n.pairs <- n.out.fold %/% 2
      # in case that createFolds does not split the data evenly (which may be often the case)
    }
    withr::with_seed(seed.fold, {
      sample.indices <- sample(n.out.fold , 2*n.pairs)
      in.fold.indices <- sample(n.in.fold, n.pairs)
    })
    index.pairs <- cbind(sample.indices[1:n.pairs],
                         sample.indices[(n.pairs+1):(2*n.pairs)])

    differences.by.perturbing.one.fold <- sapply(1:n.pairs, function(i){
      # i is the perturbed index
      index.first <- index.pairs[i,1]
      index.second <- index.pairs[i,2]

      # get an index in fold
      index.in.fold <- in.fold.indices[i]

      mu.i <- (mu.out.fold*n.out.fold - out.fold.sample[index.second,])/(n.out.fold-1)
      mu.i.perturbed <- (mu.out.fold*n.out.fold - out.fold.sample[index.first,])/(n.out.fold-1)

      weights.1 <- LDATS::softmax(-lambda*mu.i[-r])
      weights.2 <- LDATS::softmax(-lambda*mu.i.perturbed[-r])

      difference <- sum((weights.1 - weights.2)*(in.fold.sample[index.in.fold, -r] - sample.mean[-r]))
      return (difference)
    })
    differences.by.perturbing.one <- c(differences.by.perturbing.one, differences.by.perturbing.one.fold)

    # estimate differences for variance
    weights <- LDATS::softmax(-lambda*mu.out.fold[-r])
    Qs <- in.fold.sample[,-r] %*% weights
    diffs.fold <- in.fold.sample[,r] - Qs
    diffs <- c(diffs, diffs.fold)
  }

  difference.by.perturbing.one.squared <- mean(differences.by.perturbing.one^2)
  variance <- stats::var(diffs)

  scaled.difference.by.perturbing.one.squared <- difference.by.perturbing.one.squared/variance
  return (ifelse(n*scaled.difference.by.perturbing.one.squared < threshold, T, F))
}
