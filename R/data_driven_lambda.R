is.lambda.feasible <- function(lambda, data, mean.smp, r, const=0.1){
  #### the subroutine is used to check whether the given lambda satisfies the first
  #### order stability, using the leave-2-out estimates

  # lambda   : a proposed lambda to be checked if it satisfies the first order stability
  # data     : n by p
  # mean.smp : sample mean of data (p by 1), calculated by colMeans(data)
  # r        : the dimension to be tested

  # output : a boolean value indicating if the given lambda satisfies the first order stability

  n <- nrow(data)
  p <- ncol(data)
  val.critical <- const

  ## use at most 100 pairs of samples for estimation
  n.pairs <- 100
  if (n < 200){
    if (n %% 2 ==0){
      n.pairs <- as.integer(n/2)
    } else {
      n.pairs <- as.integer((n - 1)/2)
    }
  }

  ## subsample from the given sample
  set.seed(ceiling(const*lambda))
  sample.indices <- sample(n, 2*n.pairs)
  index.pairs <- cbind(sample.indices[1:n.pairs], sample.indices[(n.pairs+1):(2*n.pairs)])

  differences.by.perturbing.one <- sapply(1:n.pairs, function(i){

    index.first <- index.pairs[i,1]
    index.second <- index.pairs[i,2]

    # get an index different from index.first and index.second
    j <- max(index.first, index.second)+1
    if (j > n){
      index.candidates <- setdiff(1:3, c(index.first, index.second))
      set.seed(ceiling(i*lambda))
      j <- sample(index.candidates, 1)
    }

    mu.no.j.i <- (mean.smp*n - data[j,] - data[index.second,])/(n-2)
    mu.no.j.i.perturbed <- (mean.smp*n - data[j,] - data[index.first,])/(n-2)

    weights.1 <- softmax(-lambda*mu.no.j.i[-r])
    weights.2 <- softmax(-lambda*mu.no.j.i.perturbed[-r])

    difference <- sum((weights.1 - weights.2)*(data[j, -r] - mean.smp[-r]))
    return (difference)
  })

  difference.by.perturbing.one.squared <- mean(differences.by.perturbing.one^2)

  # estimate variance by leaving one out
  Qs <- sapply(index.pairs[,1], function(x) return (getY.softmax(x, r, data, mean.smp, lambda)))
  diffs <- data[index.pairs[,1],r] - Qs
  variance <- var(diffs)

  scaled.difference.by.perturbing.one.squared <- difference.by.perturbing.one.squared/variance
  return (ifelse(n*scaled.difference.by.perturbing.one.squared < val.critical, T, F))
}

lambda.adaptive.enlarge <- function(lambda, data, mean.smp, r, const=0.1){
  # lambda   : initial lambda
  # data     : n by p
  # mean.smp : sample mean of data (p by 1), calculated by colMeans(data)
  # r        : the dimension to be tested
  n <- nrow(data)

  lambda.curr <- lambda
  lambda.next <- 2*lambda
  feasible <- is.lambda.feasible(lambda.next, data, mean.smp, r, const)
  count <- 1
  while (feasible & lambda.next < n){
    lambda.curr <- lambda.next
    lambda.next <- 2*lambda.next
    feasible <- is.lambda.feasible(lambda.next, data, mean.smp, r, const)
    #print(glue('{lambda.next} {feasible}'))
    count <- count + 1
  }
  print(glue('before: {lambda}, after: {lambda.curr}, iteration: {count}'))
  return (lambda.curr)
}

lambda.adaptive <- function(data, mean.smp, r, const=1){
  # r is the index of the 'argmin'
  # data: n by p
  # const: scaling parameter for lambda
  n <- nrow(data)
  p <- ncol(data)
  Xs.min <- sapply(1:n, function(i){
    mu.hat.noi <- (mean.smp*n - data[i,])/(n-1)

    # find argmin among all dimensions except r
    # s <- 0
    # val.min <- 1e4 # any sufficiently large value
    # for (t in 1:p){
    #   if (t != r & mu.hat.noi[t] < val.min){
    #     s <- t
    #     val.min <- mu.hat.noi[t]
    #   }
    # }
    # return (data[i,s])
    min.indices <- which(mu.hat.noi[-r] == min(mu.hat.noi[-r]))

    set.seed(i)
    min.idx <- ifelse(length(min.indices) > 1,
                      sample(c(min.indices), 1), min.indices[1])

    X.min <- data[i, -r][min.idx]
    return (X.min)
  })
  return (sqrt(n)/(const*sd(Xs.min)))
}
