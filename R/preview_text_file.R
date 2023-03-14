#' @title Preview text file.
#'
#' @description Given a text file name (such as a script or README), print out
#' the top <lines_out> lines followed by a final <suffix_line> using cat().
#' Useful for previewing text documents within markdowns.
#'
#' @param file_name Name of text file to preview.
#' @param lines_out Number of lines to preview; defaults to 12.
#' @param suffix_line Line to append to the tail of the output. Defaults to "...".
#'
#' @import magrittr
#'
#' @export

preview_text_file   <- function(file_name,lines_out=12,suffix_line="..."){
  scan(file_name,what=character(),sep = "\n",quiet = TRUE) %>%
    head(n=lines_out) %>%
    c("...") %>%
    paste(collapse="\n") %>%
    cat
}
