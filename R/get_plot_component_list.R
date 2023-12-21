#' @title Get plot component list
#' 
#' @description
#' Wrapper for cowplot::get_plot_component() which always returns a list and can
#' omit zero grobs (extremely reasonable yet somehow not standard behavior). 
#' 
#' @param gtable_in GTable object normally fed to get_plot_component.
#' @param pattern String pattern provided to get_plot_component used to select grobs ('guide-box-bottom', for instance).
#' @param return_all Should all matching grobs be returned (in a list)? Defaults to TRUE.
#' @param drop_zero_grobs Should zero (empty) grobs be ignored? Defaults to TRUE.
#' 
#' @import cowplot
#' 
#' @export

get_plot_component_list <- function(gtable_in,pattern,return_all=TRUE,drop_zero_grobs=TRUE){
  g_tb_out <- get_plot_component(gtable_in,pattern,return_all=return_all)
  if(inherits(g_tb_out,"grob")){
    g_tb_out<- list(g_tb_out)
  }
  if(drop_zero_grobs){
    g_tb_out[sapply(g_tb_out,inherits,"zeroGrob")] <- NULL
  }
  return(g_tb_out)
}
