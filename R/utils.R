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

  min.val <- min(nums[-idx])
  min.indices <- setdiff(which(nums == min.val), c(idx))
  # withr::with_seed(ceiling(abs(idx*7+nums[idx]) + idx), {
    min.idx.sec <- ifelse((length(min.indices) > 1), sample(c(min.indices), 1), min.indices[1])
  # })
  return (min.idx.sec)
}
