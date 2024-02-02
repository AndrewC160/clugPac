#' @title
#' Resize facets 2D
#'
#' @description
#' Given a GGplot with horizontal and/or vertical discrete facets, resize facets
#' to fit data points based on X and Y ranges in each. Facets can be resized
#' vertically (<vertical>=TRUE) and/or horizontally (<horizontal>=TRUE), with
#' both being the default. Dimensions are calculated based on the x and y max/
#' min values in each row/column of panels, so all row/column panels will have
#' the same respective dimensions (i.e., all panels in a row will have the same
#' height, all panels in a column will have the same width). Plots are built and
#' converted back to ggplot objects before they are returned, so don't expect
#' most ggplot behaviors to be applicable after applying this function (in other
#' words, plan to apply themes/scales/etc. prior to applying this function).
#'
#' @param plot_in GGplot object.
#' @param vertical Boolean; should panels be scaled vertically? Defaults to TRUE.
#' @param horizontal Boolean; should panels be scaled horizontally? Defaults to TRUE.
#' @param min_x_frac Minimum fraction to scale panels to horizontally; defaults to 0.1, in which case a panel will never scale down to less than 10% the size of the largest panel.
#' @param min_y_frac Minimum fraction to scale panels to vertically; defaults to 0.1, in which case a panel will never scale down to less than 10% the size of the tallest panel.
#'
#' @import ggplot2
#' @import ggplotify
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import stringr
#'
#' @export

resize_facets_2d  <- function(plot_in,vertical=TRUE,horizontal=TRUE,min_x_frac = 0.1,min_y_frac=0.1){
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  arrange <- dplyr::arrange

  p_blt   <- ggplot_build(plot_in)
  p_grd   <- ggplot_gtable(p_blt)

  p_szs   <- as_tibble(p_blt$data[[1]]) %>%
    group_by(PANEL) %>%
    summarize(x_min = min(xmin),
              x_max = max(xmax),
              y_min = min(ymin),
              y_max = max(ymax))

  tb_layout<- as_tibble(p_grd$layout) %>%
    filter(grepl("panel",name)) %>%
    mutate(x_idx = as.integer(as.factor(l)),
           y_idx = as.integer(as.factor(t))) %>%
    arrange(y_idx,x_idx) %>%
    mutate(PANEL=as.factor(row_number())) %>%
    inner_join(p_szs,by="PANEL")

  if(vertical){
    height_fracs <- tb_layout %>%
      group_by(y_idx) %>%
      summarize(ymin = min(y_min),
                ymax = max(y_max),
                .groups="drop") %>%
      mutate(range = ymax - ymin,
             frac = range / max(range),
             frac = ifelse(frac < min_x_frac,min_x_frac,frac)) %>%
      vectify(frac,y_idx)

    grbs_v  <- which(as.character(p_grd$heights) == "1null")
    p_grd$heights[grbs_v] <- p_grd$heights[grbs_v] * height_fracs
  }
  if(horizontal){
    width_fracs <- tb_layout %>%
      group_by(x_idx) %>%
      summarize(xmin = min(x_min),
                xmax = max(x_max),
                .groups="drop") %>%
      mutate(range = xmax - xmin,
             frac = range / max(range),
             frac = ifelse(frac < min_y_frac,min_y_frac,frac)) %>%
      vectify(frac,x_idx)

    grbs_h  <- which(as.character(p_grd$widths) == "1null")
    p_grd$widths[grbs_h]  <- p_grd$widths[grbs_h] * width_fracs
  }
  return(ggplotify::as.ggplot(p_grd))
}


# resize_facets   <- function(plot_in,silent=FALSE){
#   p_blt   <- ggplot_build(plot_in)
#   p_grd   <- ggplot_gtable(p_blt)
#   fct_szs <- p_blt$data[[1]] %>%
#     as_tibble() %>%
#     group_by(PANEL) %>%
#     summarize(sz = n(),.groups="drop") %>%
#     vectify(sz,PANEL)
#   fct_szs <- fct_szs / max(fct_szs)
#   fct_ps  <- p_blt$layout$facet_params
#   fct_arr <- case_when(is.null(fct_ps$ncol) & is.null(fct_ps$nrow) ~ "none",
#                        !is.null(fct_ps$ncol) & !is.null(fct_ps$nrow) ~ "too many",
#                        is.null(fct_ps$nrow) ~ "vertical",
#                        TRUE ~ "horizontal")
#
#   if(fct_arr == "none"){
#     return(plot_in)
#   }else if(fct_arr == "too many"){
#     if(!silent){message("resize_facets() only works with one-dimensional facets, i.e. one column or one row.")}
#     return(plot_in)
#   }else if(fct_arr == "vertical"){
#     grbs  <- which(as.character(p_grd$heights) == "1null")
#   }else{
#     grbs  <- which(as.character(p_grd$widths) == "1null")
#   }
#   if(length(fct_szs) != length(grbs)){
#     if(!silent){message("Incongruent facet dimensions...did you align vertically / horizontally and fail to correctly specify?")}
#     return(plot_in)
#   }
#   if(fct_arr == "vertical"){
#     p_grd$heights[grbs] <- p_grd$heights[grbs] * fct_szs
#   }else{
#     p_grd$widths[grbs]  <- p_grd$widths[grbs] * fct_szs
#   }
#   return(ggplotify::as.ggplot(p_grd))
# }

