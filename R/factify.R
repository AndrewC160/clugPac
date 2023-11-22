#' @title
#' Factify
#'
#' @description
#' Convert tibble columns to factors if they meet specific requirements. Factors
#' already present in the table are ignored, as are logical columns.
#'
#' For character columns, factors must:
#' 1) Have fewer than nrow() or <max_unique_values> unique values in the table.
#' 2) Have fewer than <max_char> characters.
#'
#' For numeric columns, factors must:
#' 1) Be integers.
#' 2) Have fewer than nrow() or <max_unique_values> unique values in the table.
#' 3) Have fewer than <max_digits> digits.
#'
#' @param table_in Input table.
#' @param max_unique_vals Maximum unique values allowed in larger tables. Defaults to 1,000.
#' @param max_char Maximum number of characters allowed (all values in the column) in character columns. Defaults to 15.
#' @param max_digits Maximum number of digits allowed (all values in the column) in integer columns. Defaults to 3.
#'
#' @import dplyr
#' @import magrittr
#'
#' @export

factify <- function(table_in,max_unique_vals=1000,max_char=15,max_digits=3){
  m_unq <- min(nrow(table_in),max_unique_vals)
  table_in %>%
    mutate_if(function(x) !is.factor(x) & !is.logical(x) & length(unique(x)) <= m_unq & !is.numeric(x) & max(nchar(x)) <= max_char,factor) %>%
    mutate_if(function(x) !is.factor(x) & !is.logical(x) & length(unique(x)) <= m_unq &  is.integer(x) & max(nchar(as.character(x))) <= max_digits,factor) %>%
    return()
}
