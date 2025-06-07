#' Get the index of the smallest dimension apart from an index
#'
#' @param nums A vector of numbers
#' @param idx An index to be excluded
#' @param seed (Optional) If provided, used to seed the random sampling (for reproducibility).
#'
#' @return The index of the second smallest dimension (as an integer).
#' @export
#'
#' @examples
#' nums <- c(1,3,2)
#' find.sub.argmin(nums,1)
#' ## return 3
#'
#' nums <- c(1,1,2)
#' find.sub.argmin(nums,1)
#' ## return 2
find.sub.argmin <- function(nums, idx, seed=NULL){

  min.val <- min(nums[-idx])
  min.indices <- setdiff(which(nums == min.val), c(idx))
  if (!is.null(seed)) {
    withr::with_seed(seed, {
      if (length(min.indices) > 1) {
        min.idx.sec <- sample(min.indices, 1)
      } else {
        min.idx.sec <- min.indices[1]
      }
    })
  } else {
    if (length(min.indices) > 1) {
      min.idx.sec <- sample(min.indices, 1)
    } else {
      min.idx.sec <- min.indices[1]
    }
  }
  return (min.idx.sec)
}
