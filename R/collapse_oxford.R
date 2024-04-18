#' @title Oxford collapse
#'
#' @description Collapse a list of things into a string with the correct Oxford
#' comma (Things, more things, and final things.) If the list is of length 1 or
#' 2, text will reflect this without commas.
#'
#' "eggs"
#' "eggs and bacon"
#' "eggs, bacon, and jury duty."
#'
#' @param text_list List of character values to collapse.
#'
#' @import stringr
#'
#' @export

oxford_collapse <- function(text_list){
  if(length(text_list) == 1){
    txt <- text_list[1]
  }else if(length(text_list) == 2){
    txt <- paste(text_list,collapse=" and ")
  }else{
    txt <- paste0(
      paste(text_list[-length(text_list)],collapse=", "),
      ", and ",text_list[length(text_list)])
  }
  return(txt)
}
