#' @title
#' Resize facets
#'
#' @description
#' Given a GGplot with horizontal facets (typically barplots), build that plot
#' using ggplot_build in order to resize facets to fit data points. Has not
#' been fully implemented yet, so only works in cases of barplots or simple
#' scales, and does not work with vertical facets. Work in progress. Note that
#' this function should be a LAST step as things like themes stop working.
#'
#' @param plot_in GGplot object.
#' @param silent Boolean, should warning messages be supressed? Defaults to FALSE.
#'
#' @import ggplot2
#' @import ggplotify
#'
#' @export


resize_facets   <- function(plot_in,silent=FALSE){
  p_blt   <- ggplot_build(plot_in)
  p_grd   <- ggplot_gtable(p_blt)
  fct_szs <- p_blt$data[[1]] %>%
    as_tibble() %>%
    group_by(PANEL) %>%
    summarize(sz = n(),.groups="drop") %>%
    vectify(sz,PANEL)
  fct_szs <- fct_szs / max(fct_szs)
  fct_ps  <- p_blt$layout$facet_params
  fct_arr <- case_when(is.null(fct_ps$ncol) & is.null(fct_ps$nrow) ~ "none",
                       !is.null(fct_ps$ncol) & !is.null(fct_ps$nrow) ~ "too many",
                       is.null(fct_ps$nrow) ~ "vertical",
                       TRUE ~ "horizontal")

  if(fct_arr == "none"){
    return(plot_in)
  }else if(fct_arr == "too many"){
    if(!silent){message("resize_facets() only works with one-dimensional facets, i.e. one column or one row.")}
    return(plot_in)
  }else if(fct_arr == "vertical"){
    grbs  <- which(as.character(p_grd$heights) == "1null")
  }else{
    grbs  <- which(as.character(p_grd$widths) == "1null")
  }
  if(length(fct_szs) != length(grbs)){
    if(!silent){message("Incongruent facet dimensions...did you align vertically / horizontally and fail to correctly specify?")}
    return(plot_in)
  }
  if(fct_arr == "vertical"){
    p_grd$heights[grbs] <- p_grd$heights[grbs] * fct_szs
  }else{
    p_grd$widths[grbs]  <- p_grd$widths[grbs] * fct_szs
  }
  return(ggplotify::as.ggplot(p_grd))
}
