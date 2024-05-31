HT.performance <- function(j.star, HT.res){
  # j.star is the true argmin
  # HT.res n.sim by p
  # a data matrix. Each of its row is a confidence set from a simulation
  n.sim <- nrow(HT.res)
  coverage <- sum(HT.res[,j.star])/n.sim
  lengths <- rowSums(HT.res)
  overcover.rate <- length(which((lengths > 1) & (HT.res[,j.star]==1)))/n.sim
  length.avg <- mean(lengths)
  return (list(coverage=coverage, overcover.rate=overcover.rate, length.avg=length.avg))
}
