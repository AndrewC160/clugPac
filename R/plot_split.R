#' @title Plot split.
#'
#' @description
#' Split plot into basic top/right/bottom/left/main panels. Useful when
#' plot_grid fails and when more specific controls are needed for plot component
#' heights/widths. Currently useful for broken cowplot::plot_grid() problem, but
#' very much a work in progress.
#'
#' @param plot_in GGPlot object to split.
#'
#' @import grid
#' @import cowplot
#' @import magrittr
#' @import gridExtra
#'
#' @export

plot_split  <- function(plot_in){
  g_tb <- ggplot_gtable(ggplot_build(plot_in)) %>% suppressWarnings

  #Top panel.
  p_top     <- get_plot_component_list(g_tb,"-t",return_all=TRUE)
  len_ax <- length(p_top)
  panel_top <- grid.arrange(grobs=p_top,layout_matrix=matrix(1:len_ax,ncol=1))
  if(length(panel_top) == 0) panel_top <- NULL

  #Left panel.
  p_left_ttl<- get_plot_component_list(g_tb,"ylab-l")
  p_left    <- get_plot_component_list(g_tb,"axis-l",return_all=TRUE)
  if(!is.null(p_left_ttl)){
    len_axs   <- length(p_left)
    title_pos <- ceiling(len_axs/2)
    title_rows<- rep(NA,length=len_axs)
    title_rows[title_pos] <- len_axs+1
    mtx   <- cbind(title_rows,c(1:len_axs))
    p_left  <- c(p_left,p_left_ttl)
    panel_left<-grid.arrange(grobs=p_left,layout_matrix=mtx)
  }
  if(length(panel_left) == 0) panel_left <- NULL

  #Right panel.
  p_axis  <- get_plot_component_list(g_tb,"axis-r",return_all=TRUE)
  p_strip <- get_plot_component_list(g_tb,"strip-r",return_all=TRUE)
  p_guide <- get_plot_component_list(g_tb,"guide-box-right")
  p_right <- c(p_axis,p_strip)
  if(length(p_axis) == length(p_strip)){
    mtx   <- matrix(1:length(p_right),ncol=2)
  }else if(length(p_strip) > length(p_axis)){
    mtx   <- matrix(1:length(p_right),ncol=1)
  }else{
    mtx   <- NULL
  }
  if(!is.null(p_guide)){
    len_axs   <- nrow(mtx)
    grob_num  <- length(mtx)
    guide_pos <- ceiling(len_axs/2)
    guide_rows<- rep(NA,length=len_axs) %>% suppressWarnings
    guide_rows[guide_pos] <- grob_num + 1
    mtx   <- cbind(mtx,guide_rows)
    p_right   <- c(p_right,p_guide)
  }
  if(!any(is.na(mtx))){
    panel_right <- grid.arrange(grobs=p_right,layout_matrix=mtx)
  }else{
    panel_right <- p_right
  }
  if(length(panel_right) == 0) panel_right <- NULL

  #Bottom panel.
  ## Need to fix bottom panel to allow for multiple axes (for instance when facet_wrap() is in the x-direction).
  p_axis  <- get_plot_component_list(g_tb,"axis-b",return_all=TRUE)
  p_ttl   <- get_plot_component_list(g_tb,"xlab-b",return_all=TRUE)
  p_guide <- get_plot_component_list(g_tb,"guide-box-bottom")
  p_bottom<- c(p_axis,p_ttl)
  if(length(p_guide) > 0){
    mtx   <- matrix(1:3,ncol=1)
    p_bottom<- c(p_bottom,p_guide)
  }else{
    mtx   <- matrix(1:2,ncol=1)
  }
  panel_bottom<- grid.arrange(grobs=p_bottom,layout_matrix=mtx)

  if(length(panel_bottom) == 0) panel_bottom <- NULL

  #Main panel.
  panel_main  <- get_panel(g_tb,return_all=TRUE) %>%
    grid.arrange(grobs=.,layout_matrix = matrix(1:length(.),ncol = 1))

  if(length(panel_main) == 0) panel_main <- NULL

  return(list(main=panel_main,
              top=panel_top,
              right=panel_right,
              bottom=panel_bottom,
              left=panel_left))
}
