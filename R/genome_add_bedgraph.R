#' @title
#' Genome: add bedgraph
#'
#' @description
#' Add a layer of bedGraph signal to a genome_plot() list and return list along
#' with the new plot and additional tables for future layers. Function will also
#' use the clugPac::downsample_bdg_table() function to downsample output tables
#' and plots to keep approximately 2,000 bins on the X-axis. This number can be
#' adjusted using the <downsample_bins> argument, and to disable downsampling
#' all together use <downsample_bins) = NULL.
#'
#' @param plot_base_in Output stemming from <genome_plot()>.
#' @param bdg_files List of BedGraph files to be included. If named, each name will correspond to a single facet.
#' @param y_name Name to display on the y-axis. Defaults to "Pileup."
#' @param bdg_fills Fill value(s) to use for each bedGraph. Should be either length 1 for all facets or the same length as the number of bedGraph files. Defaults to "black".
#' @param bdg_alphas Alpha value(s) to use for each bedGraph. Defaults to 1.
#' @param bdg_multiplier If a multiplier value should be used to adjust bedGraph scores (such as with spike normalization), include the adjustment values here. Defaults to 1.
#' @param facet_scales Argument for ggplot2::facet_wrap(facet_scales). Defaults to "free_y".
#'
#' @import GenomicRanges
#' @import IRanges
#' @import tidyr
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#'
#' @export

genome_add_bedgraph <- function(plot_base_in,bdg_files,y_name=NULL,bdg_fills="black",bdg_alphas=1,bdg_multiplier=1,facet_scales="free_y",downsample_bins=2000,coarse_downsample_lines=1E5){
  if(is.null(names(bdg_files))) names(bdg_files)  <- gsub(".bdg.gz|.bdg|.bedGraph","",basename(bdg_files))
  if(is.null(y_name)){
    y_name  <- plot_base_in$y_name
    if(is.null(plot_base_in$y_name)){
      y_name<- "Score"
    }
  }
  tb_fls <-
    tibble(name=names(bdg_files),
           fill=bdg_fills,
           alpha=bdg_alphas,
           bdg=bdg_files,
           weight=bdg_multiplier)

  nm_vec  <- vectify(tb_fls,name,bdg)
  tb_bdg  <- vectify(tb_fls,bdg,name) %>%
    scan_bdg(gr_regions=plot_base_in$gr_win,
             name_vec = nm_vec,
             coarse_downsample_lines = coarse_downsample_lines) %>%
    group_by(name) %>%
    downsample_bdg_table(bin_num  = downsample_bins) %>%
    ungroup %>%
    left_join(tb_fls,by="name") %>%
    mutate(score = score * weight)

  tb_scales <- tb_bdg %>%
    group_by(name) %>%
    summarize(ymin=min(score),
              ymax=max(score),
              .groups="drop")

  y_scale <- scale_y_continuous(name=y_name,expand=expansion(mult=c(0,0.1)),labels=comma)

  p <- plot_base_in$plot +
    scale_y_continuous(name=y_name,expand=expansion(mult=c(0,0.1)),labels=comma) +
    facet_wrap(.~name,ncol=1,strip.position = "right",scales=facet_scales) +
    geom_rect(data=tb_bdg,mapping=aes(xmin=start,xmax=end,ymin=0,ymax=score,
                                      fill=fill,alpha=alpha)) +
    theme(strip.background = element_blank(),
          strip.text.y.right = element_text(angle=0))

  plot_base_out <- plot_base_in
  plot_base_out$plot <- p
  plot_base_out$y_scale <- y_scale
  plot_base_out$y_name  <- y_name
  plot_base_out$scale_y_limits <- tb_scales
  return(plot_base_out)
}

