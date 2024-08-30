#' @title
#' Thin column
#'
#' @description
#' Given a vector of column values, convert any value that is equal to the value
#' preceding it to empty space. Always returns a character vector, generally
#' useful for tables to display.
#'
#' For instance, the column A: "Apples","Apples","Eggs","Eggs","Eggs" will
#' become "Apples","","Eggs","","".
#'
#' @param vals_in Vector of values to interpret.
#'
#' @export

thin_column <- function(vals_in){
  vals_sft  <- vals_in[-1]
  idx <- which(vals_sft == vals_in[-length(vals_in)]) + 1
  vals_out  <- vals_in
  vals_out[idx] <- ""
  return(vals_out)
}
