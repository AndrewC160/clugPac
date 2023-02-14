#' @title
#' Truncate tibble.
#'
#' @description
#' Given a table input, print out <n_row> and <n_col> rows and columns,
#' respectively, with an added row/column of ellipses denoting how many have
#' been truncated. Useful for previewing large tables in an HTML.
#'
#' @param table_in Table to be truncated.
#' @param n_row Number of rows to include.
#' @param n_col Number of columns to include.
#' @param row_suffix Suffix to describe truncated rows. Defaults to "x...", i.e. "23x...".
#' @param col_suffix Suffix to describe truncated columns in heading. Defaults to "y...", i.e. "23y...".
#' @param sig_figs Number of significant figures to round numbers to. Defaults to 2.
#'
#' @import tibble
#' @import dplyr
#' @import magrittr
#'
#' @export

trunc_tibble    <- function(table_in, n_row = 10,n_col=10,row_suffix = "x...",col_suffix="x...",sig_figs=2){
  select  <- dplyr::select
  filter  <- dplyr::filter
  rename  <- dplyr::rename

  dm  <- dim(table_in)
  n_row <- min(c(n_row,dm[1]))
  n_col <- min(c(n_col,dm[2]))

  tb_out<- table_in %>%
    filter(row_number() <= n_row) %>%
    select(c(1:n_col)) %>%
    mutate_if(is.double,round,digits=sig_figs) %>%
    mutate_if(is.factor,as.character)

  if(dm[2] != ncol(tb_out)){
    col_nm <- paste0("(",prettyNum(dm[2],big.mark=","),col_suffix,")")
    tb_out <- tb_out %>%
      cbind(trunc_col = rep("...",nrow(.))) %>%
      rename(!!as.name(col_nm):=trunc_col)
  }

  if(dm[1] != nrow(tb_out)){
    trunc_row <- c(paste0("(",prettyNum(dm[1],big.mark=","),row_suffix,")"),
                   rep("...",ncol(tb_out)-1))
    tb_out <- rbind(tb_out,trunc_row)
  }
  return(tb_out)
}
