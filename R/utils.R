#' Get the index of the smallest dimension apart from an index
#'
#' @param nums A vector of numbers
#' @param idx An index to be excluded
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
find.sub.argmin <- function(nums, idx){
  p <- length(nums)

  min.idx.sec <- -1
  min.val <- 1e308
  for (s in 1:p){
    if (s != idx){
      if (nums[s] < min.val){
        min.idx.sec <- s
        min.val <- nums[s]
      }
    }
  }
  return (min.idx.sec)
}
