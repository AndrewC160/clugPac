#' @title
#' Genome plot
#'
#' @description
#' Create a base plot for a genomic region. Produces a list containing GGplot
#' <plot> as well as a list of features for use with subsequent plot additions,
#' including the original grange <gr_win>,the x range, and a scale for
#' the X-axis.
#'
#' @param gr_window Single contiguous genomic range to display.
#' @param x_name Name to display on the x-axis. Defaults to NULL, in which case output of grange_desc() is used.
#' @param y_name Name to display on the y-axis. Defaults to "Pileup."
#'
#' @import GenomicRanges
#' @import IRanges
#' @import tidyr
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import BiocGenerics
#' @import ggplot2
#'
#' @export

genome_plot <- function(gr_window,x_name=NULL,y_name="Pileup"){
  if(length(gr_window) > 1) warning("Only the first <gr_window> region will be used.")
  gr_win<- gr_window[1]
  tb_win<- as_tibble(gr_win)

  x_rng <- c(BiocGenerics::start(gr_win),
             BiocGenerics::end(gr_win))
  y_rng <- c(0,1)
  if(is.null(x_name)) x_name <- grange_desc(gr_win)
  if(width(gr_window) > 1E5){
    lab_func <- prettyBP
  }else{
    lab_func <- comma
  }
  x_scale<- scale_x_continuous(name=x_name,limits=x_rng,expand=c(0,0),labels=lab_func,oob=scales::oob_keep)
  p <- ggplot(tb_win,mapping=aes(xmin=start,xmax=end,ymin=0,ymax=1)) +
    x_scale +
    scale_fill_identity() +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_linetype_identity() +
    scale_linewidth_identity() +
    theme(plot.background = element_rect(fill="white",color=NA),
          panel.background = element_blank(),
          panel.grid.major.x = element_line(linewidth=0.5,color="gray65",linetype="dotted"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.y = unit(0,"lines"),
          strip.background = element_blank(),
          strip.text.y.right = element_text(angle=0),
          legend.position = "none")
  return(list(plot=p,gr_win=gr_win,x_rng=x_rng,y_rng=y_rng,x_scale=x_scale))
}
