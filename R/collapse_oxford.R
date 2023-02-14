#' @title Oxford collapse
#' 
#' @description Collapse a list of things into a string with the correct Oxford
#' comma (Things, more things, and final things.)
#' 
#' @param text_list List of character values to collapse.
#' 
#' @import stringr
#' 
#' @export

oxford_collapse <- function(text_list){
  paste0(
    paste(text_list[-length(text_list)],collapse=", "),
    ", and ",text_list[length(text_list)])
}
