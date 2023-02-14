#' @title
#' Scale to
#'
#' @description
#' Scale a column of numbers between a minimum and maximum value.
#'
#' @param values_in Vector of values to be scaled.
#' @param min_range Lower bound of range to scale, defaults to -1.
#' @param max_range Upper bound of range to scale, defaults to 1.
#' @param min_val Minimum value to arbitrarily tie the range to. If missing, the lowest number in <values_in> will be used.
#' @param max_val Maximum value to arbitrarily tie the range to. If missing, the highest number in <values_in> will be used.
#'
#' @export

scale_to        <- function(values_in,min_range=-1,max_range=1,min_val,max_val){
  if(missing(min_val)){
    min_val <- min(values_in)
  }
  if(missing(max_val)){
    max_val <- max(values_in)
  }
  (max_range - min_range) * (values_in - min_val) / (max_val - min_val) + min_range
}
