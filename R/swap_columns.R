#' @title Swap Columns
#'
#' @description
#' Given a table with at least two named columns, swap the values in columnA
#' (<colA>) with those of column B (<colB>) and vice versa, then return the
#' table. Swapping can be made conditional either by providing a function to
#' mapply which evaluates to TRUE/FALSE for each column value (for instance,
#' condition=function(x,y) x>y will swap values if column A is greater than
#' columnB). For more complex conditions or conditions which depend on other
#' columns in the table, a boolean column can be generated beforehand and its
#' name specified with a string.
#'
#' @param tb_in Tibble to modify.
#' @param colA Column name (as a character vector) of first column.
#' @param colB Column name of second column.
#' @param condition Function of colA and colB values to determine if they should be swapped, or the name of a third column with boolean values. Defaults to NULL, i.e. swap all values.
#'
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import magrittr
#'
#' @examples
#'library(magrittr)
#'library(tibble)
#'library(dplyr)
#'library(tidyr)
#'
#'# Dummy table.
#'tb <- tibble(idx=c(1:29)) %>%
#'  mutate(x=sample(1:1000,size=n()),
#'         y=sample(1:1000,size=n()))
#'
#'swap_columns(tb,colA = "x",colB = "y")
#'
#'# Swap columns if x < y
#'tb %>%
#'  swap_columns("x","y",function(x,y) x<y) %>%
#'  mutate(XgtY = x>y)
#'
#'# Swap columns if y is even.
#'tb %>%
#'  swap_columns("x","y",function(x,y) y %% 2 == 0)
#'
#'# Swap columns based on a third column.
#'tb %>%
#'  mutate(cond = sample(c(TRUE,FALSE),size=n(),replace=T)) %>%
#'  swap_columns("x","y","cond")
#'
#' @export

swap_columns  <- function(tb_in,colA,colB,condition=NULL){
  valsA1 <- valsA2 <- tb_in[,colA] %>% unlist
  valsB1 <- valsB2 <- tb_in[,colB] %>% unlist
  valsC <- NULL

  if(is.null(condition)){
    valsC <- rep(TRUE,length(valsA1))
  }else if(is.character(condition)){
    valsC <- tb_in[,condition] %>% unlist
  }else if(is.function(condition)){
    valsC <- mapply(condition,valsA1,valsB1)
  }else{
    stop("Condition argument must be a column name, a function to mapply across colA and colB, or NULL.")
  }

  valsA2[valsC] <- valsB1[valsC]
  valsB2[valsC] <- valsA1[valsC]

  tb_in %>%
    mutate(!!as.name(colA) := valsA2,
           !!as.name(colB) := valsB2) %>%
    return()
}



